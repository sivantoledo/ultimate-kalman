package sivantoledo.kalman;

import java.util.Arrays;

import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.DiagonalCovarianceMatrix.Representation;

public class RealCovarianceMatrix implements CovarianceMatrix {
  
  public static enum Representation {
    COVARIANCE_MATRIX, // covariance matrix C
    FACTOR,            // ??? 
    INVERSE_FACTOR     // W such that W'*W = inv(C)     
  }

  private RealMatrix W;
  
  public int dimension() { return W.getColumnDimension(); }
  
  public RealCovarianceMatrix copy() {
    return new RealCovarianceMatrix(W, Representation.INVERSE_FACTOR);
  }
  
  public RealCovarianceMatrix(RealMatrix v, Representation rep) {
    switch (rep) {
    case COVARIANCE_MATRIX:
      //Matrix.print(v, "going to chol...");
      //System.out.printf("diff %.2e - %.2e = %.2e\n",v.getEntry(0, 1),v.getEntry(1, 0),v.getEntry(0, 1)-v.getEntry(1, 0));
      try {
        CholeskyDecomposition chol = new CholeskyDecomposition(v);
        W = MatrixUtils.inverse(chol.getL());
      } catch (NonPositiveDefiniteMatrixException npdme) {
        EigenDecomposition eig = new EigenDecomposition(v);
        double[] eigenvalues = eig.getRealEigenvalues();
        System.out.printf("covariance matrix is not symmetric positive definite\n");
        Matrix.print(v, "cov");
        for (int i=0; i<eigenvalues.length; i++) {
          System.out.printf("eigenvalue %.2e + i * %.2e\n", eig.getRealEigenvalue(i),eig.getImagEigenvalue(i));
        }
        double mev = Arrays.stream(eigenvalues).min().getAsDouble();
        RealMatrix p = v.copy().add(MatrixUtils.createRealIdentityMatrix(eigenvalues.length).scalarMultiply(2*Math.abs(mev)));
        Matrix.print(v, "perturbed_cov");
        // Try again
        try {
          CholeskyDecomposition chol = new CholeskyDecomposition(p);
          W = MatrixUtils.inverse(chol.getL());
          System.out.printf("Perturbed to positive definiteness\n");
        } catch (NonPositiveDefiniteMatrixException still) {
          System.out.printf("Perturbed, but still not positive definite\n");
        }
      }
      break;
    case FACTOR:
      W = MatrixUtils.inverse(v);
      break;
    case INVERSE_FACTOR:
      W = v.copy();
      break;      
    }
  }

  @Override
  public RealVector weigh(RealVector v) {
    return W.operate(v);
  }

  @Override
  public RealMatrix weigh(RealMatrix A) {
    return W.multiply(A);
  }

  @Override
  public String toString() { return "W="+Matrix.toString(W.getData(),"%.3e"); };

}
