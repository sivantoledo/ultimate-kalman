import numpy as np
from scipy.linalg import qr
from covariance_matrix import CovarianceMatrix
from kalman_base import KalmanBase
from flexible_array import FlexibleArray

class KalmanUltimate(KalmanBase):
    """
    KalmanUltimate: An implementation of the Paige-Saunders Kalman filter and smoother.
    
    Inherits from KalmanBase and implements evolve, observe, and smooth methods.
    
    Attributes:
        options (dict): Configuration options, including 'covarianceEstimates' 
                        ('PaigeSaunders' or 'SelInv').
    """

    def __init__(self, options=None):
        """
        Initialize the KalmanUltimate object.
        
        Args:
            options (dict, optional): Configuration options. Defaults to empty dict.
        """
        if options is None:
            options = {}
        super().__init__(options)
        self.steps = FlexibleArray()  # Initialize with FlexibleArray

    def rollback_step_with_index(self, i):
        """
        Reconstruct the current state to just after evolve() at step i.
        
        Args:
            i (int): Index of the step to roll back to.
        """
        self.current = {}
        step_data = self.steps.data[i]
        self.current["dimension"] = step_data["dimension"]
        self.current["step"] = step_data["step"]
        #self.current["Rbar"] = step_data["Rbar"]
        #self.current["ybar"] = step_data["ybar"]
        self.current["Rbar"] = step_data.get("Rbar", None)  # Default to None if missing
        self.current["ybar"] = step_data.get("ybar", None)  # Default to None if missing
    
    def evolve(self, n_i, H_i, F_i, c_i, K_i):
        """
        Evolve the state using the given linear recurrence.
        
        Args:
            n_i (int): Dimension of the current state u_i.
            H_i (np.ndarray): Matrix relating u_i to the evolution equation.
            F_i (np.ndarray): Evolution matrix from u_{i-1} to u_i.
            c_i (np.ndarray): Constant term in the evolution equation.
            K_i (CovarianceMatrix): Covariance matrix of the error term epsilon.
        """
        self.current = {"dimension": n_i}

        if self.steps.size() == 0:  # First step
            self.current["step"] = 0
            return

        ptr_imo = self.steps.last_index()  # Pointer to previous step
        self.current["step"] = self.steps.data[ptr_imo]["step"] + 1
        n_imo = self.steps.data[ptr_imo]["dimension"]  # Dimension of u_{i-1}

        l_i = F_i.shape[0]  # Row dimension of F_i
        if H_i.shape[0] != l_i:  # Handle empty or mismatched H_i
            if l_i == n_i:
                # print("generating H_i=I")
                H_i = np.eye(l_i)
            else:
                H_i = np.hstack([np.eye(l_i), np.zeros((l_i, n_i - l_i))])

        V_i_F_i = -K_i.weigh(F_i)
        V_i_c_i = K_i.weigh(c_i)
        V_i_H_i = K_i.weigh(H_i)
        
        #print(V_i_F_i)
        #print(V_i_c_i)
        #print(V_i_H_i)

        # Construct matrices based on previous step's Rdiag
        if "Rdiag" in self.steps.data[ptr_imo] and self.steps.data[ptr_imo]["Rdiag"] is not None:
            #print("has Rdiag in previous step")
            z_i = self.steps.data[ptr_imo]["Rdiag"].shape[0]
            A = np.vstack([self.steps.data[ptr_imo]["Rdiag"], V_i_F_i])
            B = np.vstack([np.zeros((z_i, n_i)), V_i_H_i])
            y = np.vstack([self.steps.data[ptr_imo]["y"], V_i_c_i])
        else:
            #print("no Rdiag in previous step")
            A = V_i_F_i
            B = V_i_H_i
            y = V_i_c_i
            
        #print("A")
        #print(A)
        #print("B")
        #print(B)
        #print("y")
        #print(y)

        # QR decomposition
        Q, R = qr(A, mode='full')
        #print("Q")
        #print(Q)
        #print("R")
        #print(R)
        B = Q.T @ B
        y = Q.T @ y

        # Update previous step
        self.steps.data[ptr_imo]["Rdiag"] = R[:min(A.shape[0], n_imo), :]
        self.steps.data[ptr_imo]["Rsupdiag"] = B[:min(B.shape[0], n_imo), :]
        self.steps.data[ptr_imo]["y"] = y[:min(y.shape[0], n_imo), :]

        # Store leftovers in current Rbar and ybar
        #print("going to update current")
        #print(B.shape[0])
        #print(n_imo)
        if B.shape[0] > n_imo:
            self.current["Rbar"] = B[n_imo:, :]
            self.current["ybar"] = y[n_imo:, :]
        #print(self.current)

    def observe(self, G_i=None, o_i=None, C_i=None):
        """
        Provide observations of the current state.
        
        Args:
            G_i (np.ndarray, optional): Observation matrix.
            o_i (np.ndarray, optional): Observation vector.
            C_i (CovarianceMatrix, optional): Covariance matrix of observation error.
        """
        n_i = self.current["dimension"]

        if o_i is None:  # No observations
            if "Rbar" in self.current and self.current["Rbar"] is not None:
                A = self.current["Rbar"]
                y = self.current["ybar"]
            else:
                A = np.array([])
                y = np.array([])
        else:  # With observations
            W_i_G_i = C_i.weigh(G_i)
            W_i_o_i = C_i.weigh(o_i)

            if "Rbar" in self.current and self.current["Rbar"] is not None:
                if W_i_G_i.size > 0 and W_i_G_i.shape[0] > 0:  # Check if non-empty and has rows
                    A = np.vstack([self.current["Rbar"], W_i_G_i])
                    y = np.vstack([self.current["ybar"], W_i_o_i])
                else:
                    A = self.current["Rbar"]
                    y = self.current["ybar"]
            else:
                A = W_i_G_i
                y = W_i_o_i

        # Process A and y with QR decomposition if A is not empty
        if A.size > 0:
            if A.shape[0] >= A.shape[1]:
                Q, R = qr(A, mode='economic')
                y = Q.T @ y
                n = min(n_i, A.shape[0])
                self.current["Rdiag"] = R[:n, :]
                self.current["y"] = y[:n, :]
            else:
                self.current["Rdiag"] = A
                self.current["y"] = y

        # Compute estimates if Rdiag is full rank
        if "Rdiag" in self.current and self.current["Rdiag"].shape[0] == n_i:
            self.current["estimatedState"] = np.linalg.solve(self.current["Rdiag"], self.current["y"])
            self.current["estimatedCovariance"] = CovarianceMatrix(self.current["Rdiag"], 'W')

        #print("appending current")
        #print(self.current)
        self.steps.append(self.current.copy())

    def smooth(self):
        """
        Compute smooth estimates of all stored states.
        """
        l = self.steps.size()

        # Backward pass to compute smoothed states
        v = None
        for i in range(l - 1, -1, -1):
            step = self.steps.data[i]
            if i == l - 1:
                v = step["y"]
            else:
                v = step["y"] - step["Rsupdiag"] @ v
            v = np.linalg.solve(step["Rdiag"], v)
            step["estimatedState"] = v

        # Compute covariance estimates based on options
        if "covarianceEstimates" not in self.options or self.options["covarianceEstimates"] == "PaigeSaunders":
            R = None
            for i in range(l - 1, -1, -1):
                step = self.steps.data[i]
                if i == l - 1:
                    R = step["Rdiag"]
                else:
                    n_ipo = R.shape[0]
                    n_i = step["Rdiag"].shape[0]
                    Q, _ = qr(np.vstack([step["Rsupdiag"], R]), mode='full')
                    S = Q.T @ np.vstack([step["Rdiag"], np.zeros((n_ipo, step["Rdiag"].shape[1]))])
                    R = S[n_ipo:n_ipo + n_i, :n_i]
                    step["estimatedCovariance"] = CovarianceMatrix(R, 'W')

        elif self.options.get("covarianceEstimates") == "SelInv":
            C = None
            for i in range(l - 1, -1, -1):
                step = self.steps.data[i]
                if i == l - 1:
                    C = step["estimatedCovariance"].explicit()
                else:
                    n_i = step["Rdiag"].shape[0]
                    Rdiag = step["Rdiag"]
                    Rsupdiag = step["Rsupdiag"]
                    C = np.linalg.solve(Rdiag, (np.eye(n_i) + Rsupdiag @ C @ Rsupdiag.T) @ np.linalg.inv(Rdiag.T))
                    step["estimatedCovariance"] = CovarianceMatrix(C, 'C')