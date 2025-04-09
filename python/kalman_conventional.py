import numpy as np
from scipy.linalg import qr
from covariance_matrix import CovarianceMatrix  # Assuming this is implemented elsewhere
from kalman_base import KalmanBase  # Assuming KalmanBase is implemented as discussed earlier

class KalmanConventional(KalmanBase):
    """
    KalmanConventional: A conventional Kalman filter and RTS smoother.
    Supports both linear and extended filtering/smoothing based on the type of evolution and observation operators.
    """
    
    def __init__(self, options=None):
        if options is None:
            options = {}
        super().__init__(options)
    
    def rollback_step_with_index(self, i):
        self.current = {}
        step_data = self.steps.data[i]
        self.current["dimension"] = step_data["dimension"]
        self.current["step"] = step_data["step"]
        self.current["predictedState"] = step_data["predictedState"]
        self.current["predictedCovariance"] = step_data["predictedCovariance"]
        self.current["F"] = step_data["F"]
        self.current["H"] = step_data["H"]
    
    def evolve(self, n_i, H_i, F_i, c_i, K_i):
        self.current = {"dimension": n_i}
        
        if len(self.steps.data) == 0:
            self.current["step"] = 0
            return
        
        prev_step = self.steps.last()
        self.current["step"] = prev_step["step"] + 1
        
        if H_i is None:
            H_i = np.eye(F_i.shape[0])
        
        self.current["F"] = F_i
        self.current["H"] = H_i
        
        if callable(F_i):
            predicted_state, jacobian = F_i(prev_step["assimilatedState"])
        else:
            predicted_state = F_i @ prev_step["assimilatedState"]
            jacobian = F_i
        
        H_inv = np.linalg.inv(H_i)
        self.current["predictedState"] = H_inv @ (predicted_state + c_i)
        #print( H_inv.shape )
        #print( jacobian.shape )
        #print( prev_step["assimilatedCovariance"].shape )
        #print( K_i.explicit().shape )
        self.current["predictedCovariance"] = H_inv @ jacobian @ prev_step["assimilatedCovariance"] @ jacobian.T @ H_inv.T + H_inv @ K_i.explicit() @ H_inv.T
        
        self.current["estimatedState"] = self.current["predictedState"]
        self.current["estimatedCovariance"] = CovarianceMatrix(self.current["predictedCovariance"], 'C')
    
    def observe(self, G_i=None, o_i=None, C_i=None):		
        if self.current["step"] == 0:
            W_i_G_i = C_i.weigh(G_i)
            W_i_o_i = C_i.weigh(o_i)
            Q, R = qr(W_i_G_i, mode='economic')
            self.current["assimilatedState"] = np.linalg.solve(R, Q.T @ W_i_o_i)
            self.current["assimilatedCovariance"] = np.linalg.inv(R.T @ R)
        else:
            if o_i is None:
                self.current["assimilatedState"] = self.current["predictedState"]
                self.current["assimilatedCovariance"] = self.current["predictedCovariance"]
            else:
                if callable(G_i):
                    predicted_observations, jacobian = G_i(self.current["predictedState"])
                else:
                    predicted_observations = G_i @ self.current["predictedState"]
                    jacobian = G_i
                
                S = jacobian @ self.current["predictedCovariance"] @ jacobian.T + C_i.explicit()
                gain = self.current["predictedCovariance"] @ jacobian.T @ np.linalg.inv(S)
                innovation = o_i - predicted_observations
                
                self.current["assimilatedState"] = self.current["predictedState"] + gain @ innovation
                self.current["assimilatedCovariance"] = (
                    self.current["predictedCovariance"] - gain @ jacobian @ self.current["predictedCovariance"]
                )
        
        self.current["estimatedState"] = self.current["assimilatedState"]
        self.current["estimatedCovariance"] = CovarianceMatrix(self.current["assimilatedCovariance"], 'C')
        
        self.steps.append(self.current.copy())
    
    def smooth(self):
        l = len(self.steps.data)
        self.steps.data[l-1]["smoothedState"] = self.steps.data[l-1]["assimilatedState"]
        self.steps.data[l-1]["smoothedCovariance"] = self.steps.data[l-1]["assimilatedCovariance"]
        
        for i in range(l-2, -1, -1):
            next_step = self.steps.data[i+1]
            curr_step = self.steps.data[i]
            
            F_next = next_step["F"]
            if callable(F_next):
                _, next_evolution_matrix = F_next(curr_step["assimilatedState"])
            else:
                next_evolution_matrix = F_next
            
            H_next = next_step["H"]
            next_evolution_matrix = np.linalg.inv(H_next) @ next_evolution_matrix
            
            backward_innovation = (
                curr_step["assimilatedCovariance"] @ next_evolution_matrix.T
                @ np.linalg.inv(next_step["predictedCovariance"])
            )
            
            curr_step["smoothedState"] = (
                curr_step["assimilatedState"] + backward_innovation @ (
                    next_step["smoothedState"] - next_step["predictedState"]
                )
            )
            curr_step["smoothedCovariance"] = (
                curr_step["assimilatedCovariance"] + backward_innovation @ (
                    next_step["smoothedCovariance"] - next_step["predictedCovariance"]
                ) @ backward_innovation.T
            )
            
            curr_step["estimatedState"] = curr_step["smoothedState"]
            curr_step["estimatedCovariance"] = CovarianceMatrix(curr_step["smoothedCovariance"], 'C')
