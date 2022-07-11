package sivantoledo.kalman;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathArithmeticException;
//import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/*
 * Kalman filtering and smoothing using the Paige-Saunders algorithm.
 */

public class PaigeSaundersKalman implements UltimateKalman {
  
  public static enum StepState {
    ADVANCED,
    EVOLVED,
    OBSERVED,
    SMOOTHED
  }
  
  public static class Step {
    public final int        index;
    public final double     timestamp;         // just a label
    public StepState  stepState;
    
    public int        evolutionCount;    // number of evolution equations predicting this state
    public int        observationCount;  // number of observation equations constraining this state
    public int        stateDimension;    // number of state variables
    //public RealMatrix WG;
    //public RealVector Wbo;
    //public RealMatrix WF;
    //public RealMatrix WI;
    //public RealVector Wbe;
    //public RealMatrix Rtilde;
    public RealMatrix Rdiag;    // diagonal block of the R factor of WA
    public RealMatrix Rsupdiag; // the next block in the same block row of WA as Rdiag
    public RealVector QTb;      // the transformed right-hand side for the same blcok row
    
    public RealVector       estimatedState;
    public CovarianceMatrix estimatedCovariance;
    
    public Step(int index, double timestamp, int dimension) {
      this.index          = index;
      this.timestamp      = timestamp;
      this.stateDimension = dimension;
      /*
       * evolve() cannot be called at index==0, so we consider step 0 as evolved already.
       */
      if (index == 0) this.stepState      = StepState.EVOLVED;
      else            this.stepState      = StepState.ADVANCED;
    }
    
    public Step copy() {
      Step copy = new Step(this.index, this.timestamp, this.stateDimension);
      //copy.stateDimension   = this.stateDimension;
      copy.evolutionCount   = this.evolutionCount;
      copy.observationCount = this.observationCount;
      copy.Rdiag            = this.Rdiag.copy();
      copy.Rsupdiag         = Rsupdiag==null ? null : this.Rsupdiag.copy();
      copy.QTb              = this.QTb.copy();
      copy.estimatedState   = estimatedState==null ? null : this.estimatedState.copy();
      copy.estimatedCovariance = estimatedCovariance==null ? null : this.estimatedCovariance.copy();
      return copy;
    }
  }
  
  
  
  //private int k = 0;
  private boolean          started = false;
  private LinkedList<Step> steps   = new LinkedList<>();
  private Step             current = null; 
  
  private Step step(long i) {
    long indexInSteps = i - firstIndex();
    if (indexInSteps < 0 || indexInSteps >= steps.size()) throw new IllegalArgumentException("step "+i+" is not in memory");
    return steps.get((int) indexInSteps);
  }
  
  @Override
  public long currentIndex() {
    if (!started) return -1;
    return current.index;
  }

  @Override
  public long firstIndex() {
    if (started==false || steps==null || steps.size()<1) throw new IllegalArgumentException("no steps have been created yet");
    return steps.getFirst().index;
  }

  @Override
  public double timestamp(long i) {
    return step(i).timestamp;
  }

  
  
  // looks like it must be ccalled on an object that has already started.
  @Override
  public UltimateKalman copy() {
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
  @Override
  public void advance(int dimension, double timestamp) {
    int index;
    if (started) {
      index = steps.getLast().index + 1;
    } else {
      index = 0;
      started = true;
    }
    current = new Step(index, timestamp, dimension);
    //current.stateDimension = dimension;
  }
  
  /*
   * Evolution equations.
   * 
   * The number of constraints can be less than current.stateDimension.
   * If it is less, then some current variables cannot be predicted from
   * the previous state. We assume that they are always the last ones.
   */
  @Override
  public void evolve(RealMatrix F, RealVector b, CovarianceMatrix C) {
    
    if (current.stepState != StepState.ADVANCED) throw new IllegalStateException("evolve(...) must be called after advance(...)");
    
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
  
  @Override
  public void evolve(RealMatrix H, RealMatrix F, RealVector b, CovarianceMatrix C) {
    
    if (current.stepState != StepState.ADVANCED) throw new IllegalStateException("evolve(...) must be called after advance(...)");

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
    K.setSubMatrix(last.Rdiag.getData(), 0,                            0);
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
    rhs.append(Wbe);                        // XXX not sure why we append again; seems like a bug...
    RealVector v   = qr.getQ().transpose().operate(rhs);
    //System.out.printf("rhs length %d %d %d %d\n", rhs.getDimension(),v.getDimension(),Wbe.getDimension(),last.QTb.getDimension());
    last.QTb       = v.getSubVector(0,                   last.stateDimension); // start and length
    //System.out.printf(">> %d %d\n", last.stateDimension, v.getDimension()-1 );
    current.QTb    = v.getSubVector(last.stateDimension, v.getDimension() - last.stateDimension); // start and length
    
    //System.out.printf("PSK e rhs %d last qtb %d current qtb %d\n",v.getDimension(),last.QTb.getDimension(),current.QTb.getDimension());

    //Matrix.print(current.Rdiag, "R3");
    
    current.stepState = StepState.EVOLVED;
  }


  /*
   * At the end of each step, the user must call observe, even if there 
   * are no actual observations.
   */

  @Override
  public void observe() {
    
    if (current.stepState != StepState.EVOLVED) throw new IllegalStateException("observe(...) must be called after evolve(...)");

    current.observationCount = 0; // no observations in this step
    QRDecomposition qr = new QRDecomposition(current.Rdiag);
    //RealMatrix Q  = qr.getQ().getSubMatrix(0, current.observationCount-1, 0, current.stateDimension-1); // thin Q factor
    current.Rdiag = qr.getR().getSubMatrix(0, current.stateDimension-1,   0, current.stateDimension-1);
    current.QTb   = qr.getQ().transpose().operate(current.QTb); 
    
    //RealMatrix RdiagInverse = MatrixUtils.inverse(current.Rdiag);
    //current.estimatedCovariance = RdiagInverse.multiply(RdiagInverse.transpose());
    current.estimatedCovariance = new RealCovarianceMatrix(current.Rdiag,RealCovarianceMatrix.Representation.INVERSE_FACTOR);

    //Matrix.print(current.Rdiag, "R1");
    
    current.stepState = StepState.OBSERVED;

    steps.addLast(current);
  }

  @Override
  public void observe(RealMatrix G, RealVector b, CovarianceMatrix C) {
    
    //System.out.printf("obs %d %s\n", current.index, current.stepState);

    if (current.stepState != StepState.EVOLVED) throw new IllegalStateException("observe(...) must be called after evolve(...)");

    current.observationCount = b.getDimension();
    
    RealMatrix WG  = C.weigh(G);
    RealVector Wbo = C.weigh(b);
    //RealMatrix WG  = G.scalarMultiply(0.00000001);
    //RealVector Wbo = b.mapMultiply(0.00000001);
    
    if (steps.size()==0 || current.Rdiag==null) {               // normally the first step
      //if (current.observationCount <= current.stateDimension) { // not enough observations to warrant a QR facotrization
      //  current.Rdiag = WG;
      //  current.QTb   = Wbo;
      //} else {
      //Matrix.print(WG, "WG");

      QRDecomposition qr = new QRDecomposition(WG);
      //RealMatrix Q  = qr.getQ().getSubMatrix(0, current.observationCount-1, 0, current.stateDimension-1); // thin Q factor
      current.Rdiag = qr.getR().getSubMatrix(0, current.stateDimension-1,   0, current.stateDimension-1);
      current.QTb   = qr.getQ().transpose().operate(Wbo); 
      //}
      //Matrix.print(current.Rdiag, "newRdiag");
      //Matrix.print(Wbo.toArray(), "Wbo");
      //Matrix.print(current.QTb.toArray(), "QTb");
    } else {
      // stack Rdiag (from previous evolve) and WG, QTb (from previous evolve) and Wbo
      RealMatrix K = MatrixUtils.createRealMatrix(current.Rdiag.getRowDimension()+WG.getRowDimension(), current.stateDimension);
      K.setSubMatrix(current.Rdiag.getData(), 0,                               0);
      K.setSubMatrix(WG.getData()           , current.Rdiag.getRowDimension(), 0);
      
      //Matrix.print(current.Rdiag, "Rdiag");
      //Matrix.print(WG, "WG");
      //Matrix.print(K, "K");
      
      QRDecomposition qr = new QRDecomposition(K);
      RealMatrix Q  = qr.getQ().getSubMatrix(0,current.Rdiag.getRowDimension()+WG.getRowDimension()-1,0,current.stateDimension-1); // thin Q factor
      current.Rdiag = qr.getR().getSubMatrix(0,current.stateDimension-1,0,current.stateDimension-1);

      //Matrix.print(current.Rdiag, "newRdiag");
      
      RealVector k = current.QTb.append(Wbo); // current RHS  
      current.QTb   = qr.getQ().getSubMatrix(0, K.getRowDimension()-1, 0, K.getColumnDimension()-1).transpose().operate(k);       // transform; equivalent to Q.transpose().operate(k);
      //Matrix.print(Wbo.toArray(), "Wbo");
      //Matrix.print(current.QTb.toArray(), "QTb");
    }
    
    //System.out.printf("PSK o current qtb %d\n",current.QTb.getDimension());

    ///RealMatrix RdiagInverse = MatrixUtils.inverse(current.Rdiag);
    //current.estimatedCovariance = RdiagInverse.multiply(RdiagInverse.transpose());
    current.estimatedCovariance = new RealCovarianceMatrix(current.Rdiag,RealCovarianceMatrix.Representation.INVERSE_FACTOR);
    
    //Matrix.print(current.Rdiag, "R2");
    
    current.stepState = StepState.OBSERVED;

    steps.addLast(current);
  }

  
  @Override
  public void smooth() {
    Iterator<Step> i = steps.descendingIterator();
    
    // compute the estimates first
    
    RealVector v = null;
    while (i.hasNext()) {
      Step s = i.next();
      //Matrix.print(s.Rdiag, "Rdiag");
      
      if (v == null || s.Rsupdiag==null) v = s.QTb.copy();
      else                               v = s.QTb.subtract( s.Rsupdiag.operate(v) );
      
      MatrixUtils.solveUpperTriangularSystem(s.Rdiag, v);
      
      s.estimatedState = v;
 
      s.stepState = StepState.SMOOTHED;
      
      //Matrix.print(v.toArray(), "state");      
    }
    //return null;
    
    // now compute the diagonal blocks of the covariance matrices
    
    RealMatrix R = null;
    
    i = steps.descendingIterator();
    
    /*
     *  the last block row already has its covariance matrix computed, the
     *  covariance of the filtered estimate. We keep it and move on to the previous
     *  row. 
     */
    
    if (i.hasNext()) {
      Step s = i.next();
      
      R = s.Rdiag;
    }
    
    /*
     * Now process the other block rows.
     */
    
    while (i.hasNext()) {
      Step s = i.next();
      
      int n = R.getRowDimension();
      int m = s.Rdiag.getRowDimension();
      
      // build the nonzero part of the previous block column and factor
      RealMatrix P = MatrixUtils.createRealMatrix(n+m, n);
      P.setSubMatrix(s.Rsupdiag.getData(), 0, 0);
      P.setSubMatrix(R.getData()         , m, 0);
      QRDecomposition qr = new QRDecomposition(P);
      
      // build the corresponding block rows of the current block column and transform
      // the bottom n rows remain zero
      RealMatrix C = MatrixUtils.createRealMatrix(n+m, m);
      C.setSubMatrix(s.Rdiag.getData(),    0, 0); 
      RealMatrix QTC = qr.getQ().transpose().multiply(C);
     
      R = QTC.getSubMatrix(n, n+m-1, 0, m-1); // start row, end row, start column, end column
            
      s.estimatedCovariance = new RealCovarianceMatrix(R,RealCovarianceMatrix.Representation.INVERSE_FACTOR);; 
    }

  }
  
  @Override
  public CovarianceMatrix covariance() {
    return steps.getLast().estimatedCovariance;
  }
  
  @Override
  public void drop() {
    while (steps.size()>1) steps.removeFirst();
  }
  
  @Override
  public RealVector filtered() {
    
    Step step = steps.getLast();
    
    if (step.stepState != StepState.OBSERVED) throw new IllegalStateException("filtered() must be called after evolve(...) and observe(...)");
    
    //System.out.printf("PSK QTb %d Rdiag %d %d\n", steps.getLast().QTb.getDimension(), steps.getLast().Rdiag.getRowDimension(), steps.getLast().Rdiag.getColumnDimension());
    if (step.Rdiag.getRowDimension() < steps.getLast().Rdiag.getColumnDimension()) throw new DimensionMismatchException(steps.getLast().Rdiag.getRowDimension(),steps.getLast().Rdiag.getColumnDimension());
    
    //Matrix.print(steps.getLast().Rdiag, "Rdiag");
    
    RealVector solution = step.QTb.copy(); 
    try {
      MatrixUtils.solveUpperTriangularSystem(step.Rdiag, solution);
    } catch (MathArithmeticException mae) {
      System.out.printf("MathArithmeticException: %s\n", mae.getMessage());
      Matrix.print(step.Rdiag,"Rdiag (exception)");
      throw mae;
    }

    step.estimatedState = solution; // keep the filter's estimate, not sure if useful
    
    return solution;
  }
  
  @Override
  public RealVector predicted() {
    
    if (current.stepState != StepState.EVOLVED) throw new IllegalStateException("predicted() must be called after evolve(...) and before observe(...)");
    
    //System.out.printf("PSK QTb %d Rdiag %d %d\n", steps.getLast().QTb.getDimension(), steps.getLast().Rdiag.getRowDimension(), steps.getLast().Rdiag.getColumnDimension());
    if (current.Rdiag.getRowDimension() < current.Rdiag.getColumnDimension()) throw new DimensionMismatchException(current.Rdiag.getRowDimension(),current.Rdiag.getColumnDimension());
    
    //Matrix.print(steps.getLast().Rdiag, "Rdiag");
    
    RealVector solution = current.QTb.copy(); 
    try {
      MatrixUtils.solveUpperTriangularSystem(current.Rdiag, solution);
    } catch (MathArithmeticException mae) {
      System.out.printf("MathArithmeticException: %s\n", mae.getMessage());
      Matrix.print(current.Rdiag,"Rdiag (exception)");
      throw mae;
    }

    return solution;
  }
  
  @Override
  public RealVector smoothed(long i) {
    Step step_i = step(i);
    
    if (step_i.stepState != StepState.SMOOTHED) throw new IllegalStateException("step "+i+" has not been smoothed yet, call smooth()");
    
    return step_i.estimatedState;
  }

  @Override
  public CovarianceMatrix covariance(long i) {
    Step step_i = step(i);
    
    if (step_i.stepState != StepState.SMOOTHED) throw new IllegalStateException("step "+i+" has not been smoothed yet, call smooth()");
    
    return step_i.estimatedCovariance;
  }

  
  //@Override
  //public RealVector state() {
  //  return steps.getLast().estimatedState;
  //}

}
