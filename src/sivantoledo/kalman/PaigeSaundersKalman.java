package sivantoledo.kalman;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/*
 * Kalman filtering and smoothing using the Paige-Saunders algorithm.
 */

public class PaigeSaundersKalman {
  
  public static class Step {
    public int        evolutionCount;    // number of evolution equations predicting this state
    public int        observationCount;  // number of observation equations constraining this state
    public int        stateDimension;    // number of state variables
    //public RealMatrix WG;
    //public RealVector Wbo;
    //public RealMatrix WF;
    //public RealMatrix WI;
    //public RealVector Wbe;
    //public RealMatrix Rtilde;
    public RealMatrix Rdiag;
    public RealMatrix Rsupdiag;
    public RealVector QTb;
    
    public RealVector       estimatedState;
    public CovarianceMatrix estimatedCovariance;
    //public RealMatrix estimatedCovariance;
    
    public Step copy() {
      Step copy = new Step();
      copy.evolutionCount = this.evolutionCount;
      copy.observationCount = this.observationCount;
      copy.stateDimension   = this.stateDimension;
      copy.Rdiag            = this.Rdiag.copy();
      copy.Rsupdiag         = Rsupdiag==null ? null : this.Rsupdiag.copy();
      copy.QTb              = this.QTb.copy();
      copy.estimatedState   = estimatedState==null ? null : this.estimatedState.copy();
      copy.estimatedCovariance = estimatedCovariance==null ? null : this.estimatedCovariance.copy();
      return copy;
    }
  }
  
  //private int k = 0;
  private LinkedList<Step> steps = new LinkedList<>();
  private Step current = null; 
  
  public PaigeSaundersKalman copy() {
    //System.out.printf("PSK copy %d %s %s\n", steps.size(),current,steps);
    PaigeSaundersKalman copy = new PaigeSaundersKalman();
    copy.current = current.copy();
    for (Step s: this.steps) copy.steps.add(s.copy());
    return copy;
  }
  
  /*
   * User must call advance before every step (including the first)
   * 
   * The dimension is the dimension of the state vector for that state.
   * The dimension can vary among states.
   * 
   * If there are unconstrained elements in the state vector (no evolution
   * from previous state) then they must appear last in the vector.
   */
  public void advance(int dimension) {
    current = new Step();
    current.stateDimension = dimension;
  }
  
  /*
   * Evolution equations.
   * 
   * The number of constraints can be less than current.stateDimension.
   * If it is less, then some current variables cannot be predicted from
   * the previous state. We assume that they are always the last ones.
   */
  public void evolve(RealMatrix F, RealVector b, CovarianceMatrix C) {
    int evolutionCount = F.getRowDimension();
    RealMatrix H;
    if (evolutionCount == current.stateDimension) {
      H = MatrixUtils.createRealIdentityMatrix(evolutionCount);
    } else {
      H = MatrixUtils.createRealMatrix(evolutionCount,current.stateDimension);
      H.setSubMatrix(MatrixUtils.createRealIdentityMatrix(evolutionCount).getData(), 0, 0);
    }
    
    evolve(H,F,b,C);
  }
  
  public void evolve(RealMatrix H, RealMatrix F, RealVector b, CovarianceMatrix C) {
    current.evolutionCount = F.getRowDimension();
    
    RealMatrix WF  = C.weigh(F).scalarMultiply(-1.0);
    RealVector Wbe = C.weigh(b);
    RealMatrix WI  = C.weigh(H);
    /*
    if (current.evolutionCount == current.stateDimension) {
      //System.out.printf("square evolution\n");
      WI  = C.weigh(MatrixUtils.createRealIdentityMatrix(F.getRowDimension()));
      //WI  = MatrixUtils.createRealIdentityMatrix(F.getRowDimension());
    } else {
      //System.out.printf("incomplete evolution\n");
      WI  = MatrixUtils.createRealMatrix(current.evolutionCount,current.stateDimension);     
      WI.setSubMatrix(C.weigh(MatrixUtils.createRealIdentityMatrix(current.evolutionCount)).getData(), 0, 0);
    }
    */
    
    /*
     * At this point, Rdiag from last step is Rtilde from the book.
     * It can be square, if there were enough observations so far
     * to fully characterize the state, but it can also be rectangular
     * with fewer rows than columns.
     */
    Step last = steps.getLast();
    //RealMatrix prevRdiag     = last.Rdiag;
    //RealVector prevQTb       = last.QTb;
    //int        prevDimension = last.stateDimension;
    
    //int d = b.getDimension();                // dimension of state vector
    //int l = prevRdiag.getRowDimension(); // normally l=d but might be smaller for initial state
    
    //System.out.printf("rdiag %d %d WF dims %d %d\n", last.Rdiag.getRowDimension(), last.Rdiag.getColumnDimension(), WF.getRowDimension(), WF.getColumnDimension());
    /*
     * Construct block column
     */
    RealMatrix K = MatrixUtils.createRealMatrix(last.Rdiag.getRowDimension()+WF.getRowDimension(), last.Rdiag.getColumnDimension());
    K.setSubMatrix(last.Rdiag.getData(), 0,                               0);
    K.setSubMatrix(WF.getData()        , last.Rdiag.getRowDimension(), 0);
    
    int lastRdiagRowDim = last.Rdiag.getRowDimension();
    
    /*
     * Factor the block column.
     * The R factor becomes Rdiag of the last step.
     * The Q factor is applied to to the next block column and to the RHS.
     */
    QRDecomposition qr = new QRDecomposition(K);
    last.Rdiag    = qr.getR().getSubMatrix(0, last.stateDimension-1,                        0, last.stateDimension-1);
    //RealMatrix Q  = qr.getQ().getSubMatrix(0, last.stateDimension+current.evolutionCount-1, 0, last.stateDimension-1); // thin Q factor
    
    RealMatrix Z  = MatrixUtils.createRealMatrix(lastRdiagRowDim+WF.getRowDimension(), current.stateDimension);
    Z.setSubMatrix(WI.getData(), lastRdiagRowDim, 0);
    
    //System.out.printf(" Q %d %d Z %d %d\n", qr.getQ().getRowDimension(),qr.getQ().getColumnDimension(),Z.getRowDimension(),Z.getColumnDimension());
    //RealMatrix Y  = Q.transpose().multiply(Z);
    RealMatrix Y  = qr.getQ().transpose().multiply(Z);
    last.Rsupdiag = Y.getSubMatrix(0,                   last.stateDimension-1,                        0, current.stateDimension-1);
    //System.out.printf("Y is %d by %d Z %d %d\n", Y.getRowDimension(), Y.getColumnDimension(), Z.getRowDimension(), Z.getColumnDimension());
    current.Rdiag = Y.getSubMatrix(last.stateDimension, Y.getRowDimension()-1, 0, current.stateDimension-1);
    
    RealVector rhs = last.QTb.append(Wbe);
    rhs.append(Wbe);
    RealVector v   = qr.getQ().transpose().operate(rhs);
    //System.out.printf("rhs length %d %d %d %d\n", rhs.getDimension(),v.getDimension(),Wbe.getDimension(),last.QTb.getDimension());
    last.QTb       = v.getSubVector(0,                   last.stateDimension); // start and length
    //System.out.printf(">> %d %d\n", last.stateDimension, v.getDimension()-1 );
    current.QTb    = v.getSubVector(last.stateDimension, v.getDimension() - last.stateDimension); // start and length
    
    //System.out.printf("PSK e rhs %d last qtb %d current qtb %d\n",v.getDimension(),last.QTb.getDimension(),current.QTb.getDimension());

    //Matrix.print(current.Rdiag, "R3");
  }


  /*
   * At the end of each step, the user must call observe, even if there 
   * are no actual observations.
   */

  public void observe() {
    current.observationCount = 0; // no observations in this step
    QRDecomposition qr = new QRDecomposition(current.Rdiag);
    //RealMatrix Q  = qr.getQ().getSubMatrix(0, current.observationCount-1, 0, current.stateDimension-1); // thin Q factor
    current.Rdiag = qr.getR().getSubMatrix(0, current.stateDimension-1,   0, current.stateDimension-1);
    current.QTb   = qr.getQ().transpose().operate(current.QTb); 
    
    //RealMatrix RdiagInverse = MatrixUtils.inverse(current.Rdiag);
    //current.estimatedCovariance = RdiagInverse.multiply(RdiagInverse.transpose());
    current.estimatedCovariance = new RealCovarianceMatrix(current.Rdiag,RealCovarianceMatrix.Representation.INVERSE_FACTOR);

    //Matrix.print(current.Rdiag, "R1");
    steps.addLast(current);
  }

  public void observe(RealMatrix G, RealVector b, CovarianceMatrix C) {
    current.observationCount = b.getDimension();
    
    RealMatrix WG  = C.weigh(G);
    RealVector Wbo = C.weigh(b);
    
    if (steps.size()==0 || current.Rdiag==null) {               // normally the first step
      //if (current.observationCount <= current.stateDimension) { // not enough observations to warrant a QR facotrization
      //  current.Rdiag = WG;
      //  current.QTb   = Wbo;
      //} else {
        QRDecomposition qr = new QRDecomposition(WG);
        //RealMatrix Q  = qr.getQ().getSubMatrix(0, current.observationCount-1, 0, current.stateDimension-1); // thin Q factor
        current.Rdiag = qr.getR().getSubMatrix(0, current.stateDimension-1,   0, current.stateDimension-1);
        current.QTb   = qr.getQ().transpose().operate(Wbo); 
      //}
    } else {
      // stack Rdiag (from previous evolve) and WG, QTb (from previous evolve) and Wbo
      RealMatrix K = MatrixUtils.createRealMatrix(current.Rdiag.getRowDimension()+WG.getRowDimension(), current.stateDimension);
      K.setSubMatrix(current.Rdiag.getData(), 0,                               0);
      K.setSubMatrix(WG.getData()           , current.Rdiag.getRowDimension(), 0);
      
      QRDecomposition qr = new QRDecomposition(K);
      //RealMatrix Q  = qr.getQ().getSubMatrix(0,current.Rdiag.getRowDimension()+WG.getRowDimension()-1,0,current.stateDimension-1); // thin Q factor
      current.Rdiag = qr.getR().getSubMatrix(0,current.stateDimension-1,0,current.stateDimension-1);

      RealVector k = current.QTb.append(Wbo); // current RHS  
      current.QTb   = qr.getQ().getSubMatrix(0, K.getRowDimension()-1, 0, K.getColumnDimension()-1).transpose().operate(k);       // transform; equivalent to Q.transpose().operate(k);
    }
    
    //System.out.printf("PSK o current qtb %d\n",current.QTb.getDimension());

    ///RealMatrix RdiagInverse = MatrixUtils.inverse(current.Rdiag);
    //current.estimatedCovariance = RdiagInverse.multiply(RdiagInverse.transpose());
    current.estimatedCovariance = new RealCovarianceMatrix(current.Rdiag,RealCovarianceMatrix.Representation.INVERSE_FACTOR);
    
    //Matrix.print(current.Rdiag, "R2");
    steps.addLast(current);
  }

  
  public void smooth() {
    Iterator<Step> i = steps.descendingIterator();
    
    RealVector v = null;
    while (i.hasNext()) {
      Step s = i.next();
      //Matrix.print(s.Rdiag, "Rdiag");
      
      if (v == null || s.Rsupdiag==null) v = s.QTb.copy();
      else                               v = s.QTb.subtract( s.Rsupdiag.operate(v) );
      
      MatrixUtils.solveUpperTriangularSystem(s.Rdiag, v);
      
      s.estimatedState = v;
      
      //Matrix.print(v.toArray(), "state");      
    }
    //return null;
  }
  
  public CovarianceMatrix covariance() {
    if (steps.size()>1) {
      System.err.printf("covariances for Kalman smoothing not implemented yet (%d steps retained)\n",steps.size());
      System.exit(1);
    }
    return steps.getLast().estimatedCovariance;
  }
  
  public void drop() {
    while (steps.size()>1) steps.removeFirst();
  }
  
  public void filter() throws DimensionMismatchException {
    
    //System.out.printf("PSK QTb %d Rdiag %d %d\n", steps.getLast().QTb.getDimension(), steps.getLast().Rdiag.getRowDimension(), steps.getLast().Rdiag.getColumnDimension());
    if (steps.getLast().Rdiag.getRowDimension() < steps.getLast().Rdiag.getColumnDimension()) throw new DimensionMismatchException(steps.getLast().Rdiag.getRowDimension(),steps.getLast().Rdiag.getColumnDimension());
    
    //Matrix.print(steps.getLast().Rdiag, "Rdiag");
    
    RealVector solution = steps.getLast().QTb.copy(); 
    try {
      MatrixUtils.solveUpperTriangularSystem(steps.getLast().Rdiag, solution);
    } catch (MathArithmeticException mae) {
      System.out.printf("MathArithmeticException: %s\n", mae.getMessage());
      Matrix.print(steps.getLast().Rdiag,"Rdiag (exception)");
      throw mae;
    }

    steps.getLast().estimatedState = solution; // keep the filter's estimate, not sure if useful
    //return solution;
  }
  
  public RealVector state() {
    return steps.getLast().estimatedState;
  }
    
  public static void main(String[] args) {
    
    double dt = 0.1;
    double g = 9.80665;
    double verticalAccel = -g;

    int d = 4;

    RealMatrix initialObservationMatrix = MatrixUtils.createRealIdentityMatrix(d);
    RealVector initialStateExpectation  = MatrixUtils.createRealVector(new double[] {0, 0, 20, 20 });
    RealVector initialObsStdDevs        = MatrixUtils.createRealVector(new double[] { 1e-6, 1e-6, 1e-6, 1e-6 });
    CovarianceMatrix initialObsVariances      = new DiagonalCovarianceMatrix( initialObsStdDevs.map( (s) -> (s*s) ) );
        
    RealMatrix constantAccelertion = MatrixUtils.createRealMatrix(new double[][] {
        { 1, 0, dt,  0 }, 
        { 0, 1,  0, dt },
        { 0, 0,  1,  0 },
        { 0, 0,  0,  1 } });
    
    RealVector transitionStdDevs  =  MatrixUtils.createRealVector(new double[] { 1e-6, 1e-6,  0.1,  0.1 } );
    CovarianceMatrix transitionVariances      = new DiagonalCovarianceMatrix( transitionStdDevs.map( (s) -> (s*s) ) );

    RealVector transformedControl =  MatrixUtils.createRealVector(new double[] { 0, 0, 0, dt*verticalAccel });

    //transformedControl = zeros(6,1);
      
    RealMatrix positionObservation = MatrixUtils.createRealMatrix(new double[][] {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 } });
    RealVector positionObsStdDevs  = MatrixUtils.createRealVector(new double[] { 0.1, 0.1 });
    CovarianceMatrix positionObsVariances = new DiagonalCovarianceMatrix( positionObsStdDevs.map( (s) -> (s*s) ) );
    
    PaigeSaundersKalman kalman = new PaigeSaundersKalman();

    kalman.advance(d);
    kalman.observe(initialObservationMatrix,initialStateExpectation,initialObsVariances);
    
    //kalman.drop();
    kalman.filter();
    kalman.drop();
    //System.out.println( Arrays.toString( kalman.filter().toArray() ));
    
    for (int i=0; i<50; i++) {
      //System.out.println(i);

      kalman.advance(d);      
      kalman.evolve(constantAccelertion, transformedControl, transitionVariances);
      kalman.observe();

      kalman.drop();
      kalman.filter();
      //kalman.covariance();

      System.out.println( Arrays.toString( kalman.state().toArray() ));
      //kalman.drop();
      
      //Matrix.print(kalman.state().toArray(), "rt_state");
      //Matrix.print(kalman.covariance(), "rt_cov");
     
    }
    
    kalman.smooth();
    
    for (Step s: kalman.steps) Matrix.print(s.estimatedState.toArray(), "state");

    //for (Step s: kalman.steps) Matrix.print(s.estimatedCovariance, "cov");

    
  }

}
