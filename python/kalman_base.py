import numpy as np
from covariance_matrix import CovarianceMatrix
from abc import ABC, abstractmethod
from collections import OrderedDict
from flexible_array import FlexibleArray  # Importing the previously created FlexibleArray

class KalmanBase(ABC):
    """
    KalmanBase: Abstract base class for Kalman filters and smoothers.
    
    Attributes:
        steps (FlexibleArray): Stores time step structures.
        current (dict): Represents the current time step.
        options (dict): Stores configuration options.
    """

    def __init__(self, options):
        """
        Initialize the KalmanBase object.
        
        Args:
            options (dict): Options for the Kalman filter.
        """
        self.steps = FlexibleArray()  # Using FlexibleArray to store steps
        self.current = None  # Current state, moved to steps at the end of observe
        self.options = options

    @abstractmethod
    def evolve(self, n_i, H_i, F_i, c_i, K_i):
        """Evolve the state using a linear matrix equation."""
        pass

    @abstractmethod
    def observe(self, G_i, o_i, C_i):
        """Provide observations of the current state."""
        pass

    @abstractmethod
    def smooth(self):
        """Compute smooth estimates of all stored states."""
        pass

    @abstractmethod
    def rollback_step_with_index(self, i):
        """Reconstructs the object to the state just after evolve() at step i."""
        pass

    def index_of_step(self, s):
        """
        Get the index of step s in the FlexibleArray.

        Args:
            s (int): Step number.

        Returns:
            int: Index in steps.
        
        Raises:
            ValueError: If step is not found or state is invalid.
        """
        #if self.steps.size() == 0:
        #    raise ValueError("Invalid internal state: No steps stored.")

        #earliest = self.steps.first()["step"]
        #i = s - earliest
        #print('xxx')
        #print(f">>> {i} {s} {earliest} {self.steps.size()}")
        #print( self.steps.data.keys() )

        #if i < 0 or i >= self.steps.size():
        #    raise ValueError("Invalid internal state: Step out of range.")

        #return i
        # trivial in this implementation
        return s

    def earliest(self):
        """
        Returns the index of the earliest step that has not been forgotten.

        Returns:
            int: Index of the earliest step, or -1 if no steps exist.
        """
        return self.steps.first()["step"] if self.steps.size() > 0 else -1

    def latest(self):
        """
        Returns the index of the latest step that was observed.

        Returns:
            int: Index of the latest step, or -1 if no steps exist.
        """
        return self.steps.last()["step"] if self.steps.size() > 0 else -1

    def forget(self, s=None):
        """
        Forget the oldest steps to save memory.

        Args:
            s (int, optional): Forget all steps <= s. Defaults to keeping only the last step.
        """
        if self.steps.size() == 0:
            return

        if s is None:
            s=self.latest() - 1
        
        while self.earliest() != -1 and self.earliest() <= s:
            #print(f"rollback calling drop_last s={s} earliest={self.earliest()}")
            self.steps.drop_first()

    def rollback(self, s):
        """
        Rollback the object to just after the evolution of step s (but before observation of step s).

        Args:
            s (int): Step to roll back to.
        """

        ptr_s = self.index_of_step(s)
        self.rollback_step_with_index(ptr_s)

        while self.latest() != -1 and self.latest() >= s:
            #print(f"rollback calling drop_last s={s} latest={self.latest()}")
            self.steps.drop_last();
			
        #self.steps.drop_last(self.steps.size() - ptr_s)
        #print('rollback!')
        #print( self.current )

    def estimate(self, s=None):
        """
        Obtain the most recent estimate of the state of a step.

        Args:
            s (int, optional): Step index. Defaults to latest step.

        Returns:
            tuple: (estimate, covariance matrix).
        """
        if self.steps.size() == 0:
            raise ValueError("No steps available for estimation.")

        if s is None:
            s = self.latest()
            
        ptr_s = self.index_of_step(s)
        #print(f"estimate step={s} ptr_s={ptr_s} earliest={self.earliest()} latest={self.latest()}")
        #print(self.steps.data.keys())
        step_data = self.steps.data[ptr_s]

        if "estimatedState" in step_data and step_data["estimatedState"] is not None:
            return step_data["estimatedState"], step_data["estimatedCovariance"]
        else:
            #import numpy as np
            #dim = step_data["dimension"]
            #return np.full((dim, 1), np.nan), np.full((dim, dim), np.nan)
            dim = step_data["dimension"]
            estimate = np.full((dim, 1), np.nan)
            # Return a CovarianceMatrix object instead of a raw NumPy array
            cov = CovarianceMatrix(np.eye(dim) * np.nan, 'C')
            #print("Warning: state is currently underdetermined")
            return estimate, cov            

    def gather(self):
        """
        Return all estimates concatenated into a single column vector.

        Returns:
            np.ndarray: Concatenated estimates.
        """
        import numpy as np

        N = sum(step["dimension"] for step in self.steps.data.values())
        u = np.zeros((N, 1))

        j = 0
        for step in self.steps.data.values():
            dim = step["dimension"]
            u[j:j + dim] = self.estimate(step["step"])[0]  # Extract state estimate
            j += dim

        return u
