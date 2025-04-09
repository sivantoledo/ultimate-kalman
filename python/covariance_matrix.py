import numpy as np
from scipy.linalg import cholesky, inv

class CovarianceMatrix:
    def __init__(self, Z, matrix_type):
        """
        Initialize a CovarianceMatrix object.

        :param Z: Matrix representing the covariance in various forms.
        :param matrix_type: String specifying the type of representation ('C', 'I', 'W', 'w', 'M', 'U').
        """
        self.type = matrix_type
        self.invW = None
        self.W = None
        self.w = None
        self.C = None  # For lazy factorization

        if matrix_type == 'C':
            self.type = 'C'
            self.C = Z
        elif matrix_type == 'I':
            self.type = 'W'
            self.W = cholesky(Z, lower=True).T
        elif matrix_type == 'M':
            self.type = 'M' if not np.allclose(Z, np.triu(Z)) else 'U'
            self.invW = Z
        elif matrix_type == 'W':
            self.type = 'W'
            self.W = Z
        elif matrix_type == 'w':
            self.type = 'w'
            self.w = Z
        else:
            raise ValueError("Illegal covariance-matrix representation")

    def weigh(self, A):
        """
        Multiply a given matrix or column vector by W.

        :param A: Input matrix or vector.
        :return: Weighted matrix or vector.
        """
        if self.type in {'U', 'M', 'F'}:
            return np.linalg.solve(self.invW, A)
        elif self.type == 'W':
            return self.W @ A
        elif self.type == 'w':
            #print( type(self.w))
            #print( self.w )
            #print (self.w.ravel() )
            #print( self.w.shape )
            #print( self.w.ravel().shape )
         
            #print( type(A))
            #print( A )
            #print( A.shape )
            # ravel() forces self.w to be a 1D array (dimensions (1,1), not (1,) when it is a scalar)
            return np.diag(self.w.ravel()) @ A
        elif self.type == 'C':
            self.type = 'F'  # Lower Cholesky factorized
            self.invW = cholesky(self.C, lower=True)
            return np.linalg.solve(self.invW, A)
        else:
            raise ValueError("Illegal covariance type")

    def rep(self):
        """
        Return the matrix representation of the covariance.

        :return: Tuple (Matrix, Type)
        """
        if self.type in {'U', 'M'}:
            return self.invW, self.type
        elif self.type == 'W':
            return self.W, self.type
        elif self.type == 'w':
            return self.w, self.type
        elif self.type == 'C':
            return self.C, self.type
        elif self.type == 'F':
            return self.invW, self.type
        else:
            raise ValueError("Illegal covariance type")

    def explicit(self):
        """
        Return the explicit dense covariance matrix C.

        :return: Explicit covariance matrix.
        """
        if self.type in {'U', 'M'}:
            return self.invW.T @ self.invW
        elif self.type == 'W':
            return inv(self.W.T @ self.W)
        elif self.type == 'w':
            #print( self.w.flatten().shape )
            #print( self.w )
            #print( self.w.flatten() ** -2 )
            #print( np.diag( self.w.flatten() ** -2 ))
            return np.diag(self.w.flatten() ** -2)
        elif self.type in {'C', 'F'}:
            return self.C
        else:
            raise ValueError("Illegal covariance type")
