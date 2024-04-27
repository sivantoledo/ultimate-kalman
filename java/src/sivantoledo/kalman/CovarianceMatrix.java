package sivantoledo.kalman;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * Representation of a covariance matrix C.
 * 
 * @author Sivan Toledo
 */

public interface CovarianceMatrix {
  /**
   * Returns the dimension of this square matrix.
   * 
   * @return the dimension of C
   */
  public int dimension();

  /**
   * Multiply the argument by a matrix W such that W'*W = inv(C)
   * 
   * (using Matlab notation, where W' is the transpose of W).
   * 
   * @param v a vector to multiply by W
   * @return W*v
   */
  public RealVector weigh(RealVector v);

  /**
   * Multiply the argument by a matrix W such that W'*W = inv(C)
   * 
   * (using Matlab notation, where W' is the transpose of W).
   * 
   * @param A a matrix to multiply by W
   * @return W*A
   */
  public RealMatrix weigh(RealMatrix A);
  
  /**
   * Returns a copy of C.
   * 
   * @return a copy of C
   */
  public CovarianceMatrix copy();
  
  /**
   * Returns an explicit representation of C.
   * 
   * @return an explicit representation of C.
   */
  public default RealMatrix get() {
    RealMatrix W = weigh(MatrixUtils.createRealIdentityMatrix(dimension()));
    //RealMatrix invW = MatrixUtils.inverse(W);
    //RealMatrix C = invW.multiply(invW.transpose());
    RealMatrix C = MatrixUtils.inverse(W.transpose().multiply(W));
    return C;
  }
}