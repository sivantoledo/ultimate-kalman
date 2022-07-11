package sivantoledo.kalman;

import org.apache.commons.math3.exception.DimensionMismatchException;
//import org.apache.commons.math3.exception.MathArithmeticException;
//import org.apache.commons.math3.linear.MatrixUtils;
//import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


/**
 * 
 * An Ultimate Kalman filter/smoother.
 * 
 * Each state is created by calling three methods in a particular order:
 *   advance
 *   evolve
 *   observe
 *   
 * In the first step, omit the call to advance.
 * 
 * @author Sivan Toledo
 *
 */
public interface UltimateKalman {

  /**
   * This method creates a new step in the filter. This must be called before
   * both evolve() and observe().
   * 
   * The dimension can vary between states.
   * 
   * The time stamp for this step will be NaN.
   * 
   * @param dimension the dimension of the state vector of the next state
   */
  default public void advance(int dimension) {
    advance(dimension, Double.NaN);
  }
  
  /**
   * This method creates a new step in the filter. This must be called before
   * both evolve() and observe().
   * 
   * The dimension can vary between states.
   * 
   * @param dimension the dimension of the state vector of the next state
   * @param timestamp a label for this step (not used internally)
   */
  void advance(int dimension, double timestamp);


  /**
   * Adds evolution equations from the previous step to the current step.
   * 
   * The matrix equation is H*u_i = F*u_{i-1} + b + e
   * 
   * where e represents Gaussian random errors with zero mean and covariance C. 
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
  void evolve(RealMatrix H, RealMatrix F, RealVector b, CovarianceMatrix C);

  /**
   * A simplified version with H=I.
   * 
   * @param F
   * @param b
   * @param C
   */
  void evolve(RealMatrix F, RealVector b, CovarianceMatrix C);


  /**
   * Adds observation equations to the current step,
   * 
   *   b = G*u_i + e
   *   
   * where e is a Gaussian random vector with zero mean and covariance C.
   * 
   * @param G
   * @param b
   * @param C
   */
  void observe(RealMatrix G, RealVector b, CovarianceMatrix C);

  /**
   * A simplified version for steps with no observations at all.
   */
  void observe();
  
  /**
   * The index of the current step. 
   * Indexing starts from zero. The method returns -1 if advance()
   * has never been called (so even step 0 is not in progress yet).
   * 
   * @return index of the current step
   */
  
  long currentIndex();
  
  /**
   * Index of the first (oldest) step that is still retained (has not been dropped).
   * 
   * @return index of the first (oldest) step that has not been dropped
   */
  long firstIndex();
  
  /**
   * The time stamp associated with a particular step.
   * If step i has been dropped or has not been reached,  
   * the method returns Double.NaN.
   * 
   * @param i index of the step
   * @return the time stamp associated with step i
   */
  double timestamp(long i);
  
  /**
   * Computes the predicted state vector of the current step.
   * Must be called after evolve() but before observe();
   * 
   * @return the predicted estimate of the state of the current step
   */

  RealVector predicted();

  /**
   * Computes the filtered state vector of the current state.
   * Must be called after evolve() and after observe();
   * 
   * @return the filtered estimate of the state of the current step
   */

  RealVector filtered();
    
  /**
   * The covariance matrix of either the current predicted or
   * the current filtered state estimate, depending on whether
   * observed has been called or not.  
   *
   * @return covariance matrix of the current estimate of the state vector
   */
  
  CovarianceMatrix covariance();
  
  /**
   * Computes the state vectors of all the steps that have not been dropped.
   */

  void smooth();

  /**
   * Drops from memory all but the current step.
   */
  
  void drop();
  
  
  /**
   * Computes the filtered state vector of the current state.
   * Must be called after evolve() and after observe();
   * 
   * @param i index of the step
   * @return the smoothed estimate of the state of step i
   */

  RealVector smoothed(long i);
    
  /**
   * The covariance matrix of the smoothed state estimate of step i.  
   *
   * @param i index of the step
   * @return the covariance matrix of the smoothed estimate of state i
   */
  
  CovarianceMatrix covariance(long i);


  /**
   * 
   * @return a copy of the object.
   */
  
  UltimateKalman copy();

}