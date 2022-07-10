package sivantoledo.kalman;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public interface CovarianceMatrix {
  public RealVector weigh(RealVector v);
  public RealMatrix weigh(RealMatrix A);
  public CovarianceMatrix copy();
  public int dimension();
  
  public default RealMatrix get() {
    RealMatrix W = weigh(MatrixUtils.createRealIdentityMatrix(dimension()));
    //RealMatrix invW = MatrixUtils.inverse(W);
    //RealMatrix C = invW.multiply(invW.transpose());
    RealMatrix C = MatrixUtils.inverse(W.transpose().multiply(W));
    return C;
  }
}