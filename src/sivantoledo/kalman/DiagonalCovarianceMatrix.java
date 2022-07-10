package sivantoledo.kalman;

import java.util.Arrays;

import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class DiagonalCovarianceMatrix implements CovarianceMatrix {
  
  public static enum Representation {
    COVARIANCE_MATRIX,
    DIAGONAL_VARIANCES,
    DIAGONAL_STANDARD_DEVIATIONS,
    DIAGONAL_INVERSE_STANDARD_DEVIATIONS,
  }

  private RealVector     weights;
  private DiagonalMatrix W;
  
  public int dimension() { 
    if (W != null) return W.getColumnDimension();
    return weights.getDimension();
  }

  
  public DiagonalCovarianceMatrix copy() {
    return new DiagonalCovarianceMatrix(weights, Representation.DIAGONAL_INVERSE_STANDARD_DEVIATIONS);
  }
  
  public DiagonalCovarianceMatrix(RealVector variances) {
    this(variances, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES);
  }

  public DiagonalCovarianceMatrix(double[] v, DiagonalCovarianceMatrix.Representation rep) {
    this(MatrixUtils.createRealVector(v),rep);
  }

  public DiagonalCovarianceMatrix(RealVector v, DiagonalCovarianceMatrix.Representation rep) {
    switch (rep) {
    case DIAGONAL_VARIANCES:
      weights = v.map( (v_i) -> (1.0/Math.sqrt(v_i)) );
      break;
    case DIAGONAL_STANDARD_DEVIATIONS:
      weights = v.map( (v_i) -> (1.0/v_i) );
      break;
    case DIAGONAL_INVERSE_STANDARD_DEVIATIONS:
      weights = v.copy();
      break;      
    }
    W = new DiagonalMatrix(weights.toArray());
  }

  @Override
  public RealVector weigh(RealVector v) {
    return weights.ebeMultiply(v);
  }

  @Override
  public RealMatrix weigh(RealMatrix A) {
    return W.multiply(A);
  }
  
  @Override
  public String toString() { return String.format("DiagonalCovarianceMatrix(weights=%s)",Arrays.toString(weights.toArray())); };

}
