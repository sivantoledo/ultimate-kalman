package sivantoledo.kalman.tests;

import java.util.Arrays;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.CovarianceMatrix;
import sivantoledo.kalman.DiagonalCovarianceMatrix;
//import sivantoledo.kalman.Matrix;
import sivantoledo.kalman.PaigeSaundersKalman;
//import sivantoledo.kalman.PaigeSaundersKalman.Step;

public class ConstantAcceleration {
  
  

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
      
    //RealMatrix positionObservation = MatrixUtils.createRealMatrix(new double[][] {
    //    { 1, 0, 0, 0 },
    //    { 0, 1, 0, 0 } });
    //RealVector positionObsStdDevs  = MatrixUtils.createRealVector(new double[] { 0.1, 0.1 });
    //CovarianceMatrix positionObsVariances = new DiagonalCovarianceMatrix( positionObsStdDevs.map( (s) -> (s*s) ) );
    
    PaigeSaundersKalman kalman = new PaigeSaundersKalman();

    kalman.advance(d);
    kalman.observe(initialObservationMatrix,initialStateExpectation,initialObsVariances);
    
    //kalman.drop();
    kalman.filtered();
    kalman.drop();
    //System.out.println( Arrays.toString( kalman.filter().toArray() ));
    
    for (int i=0; i<50; i++) {
      //System.out.println(i);

      kalman.advance(d);      
      kalman.evolve(constantAccelertion, transformedControl, transitionVariances);
      kalman.observe();

      kalman.drop();
      RealVector filtered = kalman.filtered();
      //kalman.covariance();

      System.out.println( Arrays.toString( filtered.toArray() ));
      //kalman.drop();
      
      //Matrix.print(kalman.state().toArray(), "rt_state");
      //Matrix.print(kalman.covariance(), "rt_cov");
     
    }
    
    kalman.smooth();
    
    //for (Step s: kalman.steps) Matrix.print(s.estimatedState.toArray(), "state");

    //for (Step s: kalman.steps) Matrix.print(s.estimatedCovariance, "cov");

    
  }


}
