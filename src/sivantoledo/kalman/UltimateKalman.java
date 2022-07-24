package sivantoledo.kalman;

import java.util.ArrayList;

//import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
//import org.apache.commons.math3.exception.MathArithmeticException;
//import org.apache.commons.math3.linear.MatrixUtils;
//import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * 
 * An Ultimate Kalman filter/smoother.
 * 
 * @author Sivan Toledo
 *
 */
public class UltimateKalman {
  
  private static class Step {
    public long       step;
    public int        dimension;
    //public StepState  stepState;
    
    public RealMatrix Rdiag;    // diagonal block of the R factor of WA
    public RealMatrix Rsupdiag; // the next block in the same block row of WA as Rdiag
    public RealVector y;      // the transformed right-hand side for the same blcok row
    
    public RealMatrix Rbar;
    public RealVector ybar;

    public RealVector state;
    public RealMatrix covariance; // the W factor of the covariance matrix
  }
  
  private Step current = null;
  private ArrayList<Step> steps = new ArrayList<Step>();
  private long first = -1; // step number of the first element in the array list

  public UltimateKalman() {
  }
  
  
  /**
   * Index of the first (oldest) step that is still retained (has not been forgotten).
   * 
   * @return index of the first (oldest) step that has not been forgotten
   */
  public long earliest() {
    return first;
  }
  
  /**
   * Index of the latest (newest) step that has not been rolled back.
   * 
   * @return of the latest (newest) step that has not been rolled back
   */
  public long latest() {
    return first + steps.size() - 1;
  }
  
  /**
   * Creates a new step with the given state dimension and provides
   * the evolution equations.
   * 
   * The matrix equation is H_i*u_i = F_i*u_{i-1} + c_i + eps_i
   * 
   * where eps_i represents Gaussian random errors with zero mean and covariance K_i. 
   * 
   * The number of constraints can be less than dimension of the current step.
   *
   * If there are unconstrained elements in the state vector (no evolution
   * from previous state) then they must appear last in the vector.
   * 
   * @param H
   * @param F
   * @param b
   * @param C
   */
  public void evolve(int n_i, RealMatrix H_i, RealMatrix F_i, RealVector c_i, CovarianceMatrix K_i) {
    current = new Step();
    current.dimension = n_i;
    
    if (steps.size() == 0) {
      current.step = 0;
      first = 0;
     return;
    }
    
    Step imo = steps.get( steps.size()-1 ); // previous step
    assert(imo != null);
    current.step = imo.step + 1;
    
    assert(H_i != null);
    assert(F_i != null);
    assert(c_i != null);
    assert(K_i != null);

    RealMatrix V_i_H_i = K_i.weigh(H_i);
    RealMatrix V_i_F_i = K_i.weigh(F_i).scalarMultiply(-1.0);
    RealVector V_i_c_i = K_i.weigh(c_i);
    
    RealMatrix z = (imo.Rdiag==null) ? null : zeros(imo.Rdiag.getRowDimension(),n_i);
    RealMatrix A = vconcat( imo.Rdiag, V_i_F_i );
    RealMatrix B = vconcat( z,         V_i_H_i );
    RealVector y = vconcat( imo.y,     V_i_c_i );
    
    QRDecomposition qr = new QRDecomposition(A);
    
    B = qr.getQT().multiply(B);
    y = qr.getQT().operate(y);
    
    imo.Rdiag     = qr.getR().getSubMatrix(0, Integer.min(imo.dimension,A.getRowDimension())-1,   
                                           0, A.getColumnDimension()-1);
    imo.Rsupdiag = B.getSubMatrix(0, Integer.min(imo.dimension,B.getRowDimension())-1,   
                                  0, B.getColumnDimension()-1);
    imo.y         = y.getSubVector(0, Integer.min(imo.dimension,y.getDimension()));
    
    if (B.getRowDimension() > imo.dimension) {
      current.Rbar = B.getSubMatrix(imo.dimension, B.getRowDimension()-1, 0, B.getColumnDimension()-1);
      current.ybar = y.getSubVector(imo.dimension, y.getDimension() - imo.dimension);
    }
  }

  /**
   * A super simplified version for the first step (or call one of the other overloaded
   * versions; the arguments are ignored).
   * 
   * @param F
   * @param b
   * @param C
   */
  public void evolve(int n_i) {
    evolve(n_i,null,null,null,null);
  }

  /**
   * A simplified version with H=I.
   * 
   * @param F
   * @param b
   * @param C
   */
  public void evolve(int dimension, RealMatrix F_i, RealVector c_i, CovarianceMatrix K_i) {
    int n_i = dimension;
    int l_i = F_i.getRowDimension();
    RealMatrix H_i = MatrixUtils.createRealMatrix(l_i, n_i);
    for (int i=0; i<Integer.min(l_i,n_i); i++) H_i.setEntry(i, i, 1.0);
    evolve(dimension,H_i,F_i,c_i,K_i);   
  }

  /**
   * Adds observation equations to the current step,
   * 
   *   o_i = G_i*u_i + delta_i
   *   
   * where delta_i is a Gaussian random vector with zero mean and covariance C_i.
   * 
   * @param G_i
   * @param o_i
   * @param C_i
   */
  public void observe(RealMatrix G_i, RealVector o_i, CovarianceMatrix C_i) {
    RealMatrix A = null;
    RealVector y = null;
    
    if (o_i == null) {  // no observations
      if (current.Rbar != null) { 
        A = current.Rbar.copy();
        y = current.ybar.copy();
      }
    } else {
    
      RealMatrix W_i_G_i = C_i.weigh(G_i);
      RealVector W_i_o_i = C_i.weigh(o_i);
    
      A = vconcat( current.Rbar, W_i_G_i );
      y = vconcat( current.ybar, W_i_o_i );
    }
    
    if (A != null) {
      if (A.getRowDimension() >= A.getColumnDimension()) {
        QRDecomposition qr = new QRDecomposition(A);    
        y = qr.getQT().operate(y);
        
        int n = Integer.min(current.dimension, A.getRowDimension());
        
        current.Rdiag = qr.getR().getSubMatrix(0, n-1,   
                                               0, A.getColumnDimension()-1);
        current.y     = y.getSubVector(0, n);       
      } else {
        current.Rdiag = A;
        current.y     = y;       
      }
    }
    
    if (current.Rdiag != null && current.Rdiag.getRowDimension() == current.dimension) {
      current.state = current.y.copy();
      MatrixUtils.solveUpperTriangularSystem(current.Rdiag, current.state);
      current.covariance = current.Rdiag.copy();
    }

    steps.add(current);
    current = null;
  }

  /**
   * A simplified version for steps with no observations at all.
   */
  public void observe() {
    observe(null,null,null);
  }
    
  /**
   * The most up to date estimate of step si
   * 
   * @param si the index of a step
   * @return the most up to date estimate of the latest step
   */

  public RealVector estimate(long si) {
    //System.out.printf("estimate si = %d\n", si);
    if (si == -1) si = latest();
    long index = si - first;
    //System.out.printf("estimate si = %d first = %d index = %d size = %d\n", si,first,index,steps.size());
    if (index < 0 || index >= steps.size()) return null;
    Step step = steps.get( (int) (si - first) );
    
    if (step.state != null) return step.state;
    return MatrixUtils.createRealVector(new double[ step.dimension ]).mapMultiply(Double.NaN);
  }

  /**
   * The most up to date estimate of the latest step
   * 
   * @return the most up to date estimate of the latest step
   */

  public RealVector estimate() {
    return estimate(-1);
  }
    
  /**
   * The covariance matrix of the most up to date estimate of the latest step.  
   *
   * @return covariance matrix
   */
  
  public CovarianceMatrix covariance(long si) {
    if (si == -1) si = latest();
    long index = si - first;
    if (index < 0 || index >= steps.size()) return null;
    Step step = steps.get( (int) (si - first) );

    if (step.covariance==null) return null;
    return new RealCovarianceMatrix(step.covariance,RealCovarianceMatrix.Representation.INVERSE_FACTOR);
  }

  /**
   * The covariance matrix of the most up to date estimate of the latest step.  
   *
   * @return covariance matrix
   */
  
  public RealMatrix covarianceInverseFactor(long si) {
    if (si == -1) si = latest();
    long index = si - first;
    if (index < 0 || index >= steps.size()) return null;
    Step step = steps.get( (int) (si - first) );

    if (step.covariance==null) return null;
    return step.covariance.copy();
  }

  /**
   * The covariance matrix of the most up to date estimate of the latest step.  
   *
   * @return covariance matrix
   */
  
  public CovarianceMatrix covariance() {
    return covariance(-1);
  }

  /**
   * Computes the state vectors of all the steps that have not been dropped.
   */

  public void smooth() {
    
    RealVector v = null;
    for (int index = steps.size()-1; index >= 0; index--) {
      Step step = steps.get(index);
      if (v == null) {     // last step
        v = step.y.copy();
      } else {
        v = step.y.subtract( step.Rsupdiag.operate(v) );
      }
      
      MatrixUtils.solveUpperTriangularSystem(step.Rdiag, v);
      step.state = v;
    }
    
    RealMatrix R = null;
    for (int index = steps.size()-1; index >= 0; index--) {
      Step step = steps.get(index);
      if (R == null) {     // last step
        R = step.Rdiag.copy();
      } else {
        int n_ipo = R.getRowDimension();
        RealMatrix A = vconcat( step.Rsupdiag, R );
        RealMatrix B = vconcat( step.Rdiag, MatrixUtils.createRealMatrix(n_ipo, step.Rdiag.getColumnDimension()) );
        
        QRDecomposition qr = new QRDecomposition(A);
        R = qr.getQT().multiply(B);
        R = R.getSubMatrix(n_ipo, R.getRowDimension()-1, 0, step.dimension-1);
        step.covariance = R;
      }      
    }

  }

  /**
   * Drops from memory all but the current step.
   */
  
  public void forget() {
    forget(-1);
  }
  
  /**
   * Drops from memory all but the current step.
   */
  
  public void forget(long si) {
    if (si == -1) si = latest()-1;
    if (si < earliest() || si > latest()-1) return; // nothing to do 
    while (steps.get(0).step <= si) {
      steps.remove(0);
      first++;
    }
  }
  
  public void rollback() {
    rollback(-1);
  }
  
  /**
   * Drops from memory all but the current step.
   */
  
  public void rollback(long si) {
    if (si == -1) si = latest(); // the default is to roll back to to last step
    if (si < earliest() || si > latest()) return; // nothing to do 
    while (steps.get( steps.size()-1 ).step > si) {
      steps.remove( steps.size()-1 );
    }    
    assert( steps.get( steps.size()-1 ).step == si );
    current = steps.get( steps.size()-1 );
    steps.remove( steps.size()-1 );
    
    // now remove all the fields that were added in observe
    current.Rdiag      = null;
    current.Rsupdiag   = null;
    current.state      = null;
    current.covariance = null;
    current.y          = null;
  }
  
  public double[] perftest(RealMatrix H, RealMatrix F, RealVector c, CovarianceMatrix K,
      RealMatrix G, RealVector o, CovarianceMatrix C,
      long count, long decimation) {

    int n = G.getColumnDimension();
    int m = G.getRowDimension();
    
    double[] t = new double[ (int) Math.floor(count/decimation) ];
    
    int j = 0;
    long start = System.nanoTime();
    
    for (int i=0; i<count; i++) {
      evolve(n,H,F,c,K);
      observe(G,o,C);
      estimate();
      forget();
      if (i % decimation == (decimation-1)) {
        // reporting is in seconds, not nanoseconds
        t[ j++ ] = 1e-9 * (double) (System.nanoTime() - start) / (double) decimation;
        start = System.nanoTime();
      }
    }
    
    return t;
  }

  
  /* 
   * utilities 
   */

  private static RealMatrix vconcat(RealMatrix T, RealMatrix B) {
    //int rows, cols;
    
    if (T== null && B == null) return null;
    if (T == null) return B.copy();
    if (B == null) return T.copy();
    
    assert(T.getColumnDimension() == B.getColumnDimension());
    
    RealMatrix C = MatrixUtils.createRealMatrix(T.getRowDimension()+B.getRowDimension(), T.getColumnDimension());
    C.setSubMatrix(T.getData(), 0,                   0);
    C.setSubMatrix(B.getData(), T.getRowDimension(), 0);
    
    return C;
  }
  
  private static RealVector vconcat(RealVector T, RealVector B) {
    if (T== null && B == null) return null;
    if (T == null) return B.copy();
    if (B == null) return T.copy();
    
    RealVector C = MatrixUtils.createRealVector(new double[ T.getDimension() + B.getDimension() ]);
    C.setSubVector(0, T);
    C.setSubVector(T.getDimension(), B);
    
    return C;
  }

  private static RealMatrix zeros(int rows, int cols) {
    RealMatrix C = MatrixUtils.createRealMatrix(rows,cols);
    return C;
  }
  
}