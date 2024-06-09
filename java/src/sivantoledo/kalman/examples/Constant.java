package sivantoledo.kalman.examples;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.exception.InsufficientDataException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.CovarianceMatrix;
import sivantoledo.kalman.DiagonalCovarianceMatrix;
//import sivantoledo.kalman.Matrix;
import sivantoledo.kalman.PaigeSaundersKalman;
//import sivantoledo.kalman.PaigeSaundersKalman.Step;
import sivantoledo.kalman.UltimateKalman;

/**
 * Kalman filter example for a constant value.
 * 
 * @author Sivan Toledo
 *
 */

public class Constant {
  
  private static Random random = new Random(69978);
  
  public static void main(String[] args) {
    
    /*
     * Part one, evolve the system for n time steps.
     */
    
    int n = 250;
    double evolutionStdDev = 1;
    RealMatrix evolutionMatrix = MatrixUtils.createRealMatrix(new double[][] {{ 1 } });
    CovarianceMatrix evolutionVariances = new DiagonalCovarianceMatrix( new double[] { evolutionStdDev }, DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS );

    RealVector[] trajectory = new RealVector[n+1];
    trajectory[0] = MatrixUtils.createRealVector(new double[] { 0 });
    
    for (int i=1; i<=n; i++) {
      RealVector noise = MatrixUtils.createRealVector(new double[] { random.nextGaussian() }).mapMultiply(evolutionStdDev);
      trajectory[i] = evolutionMatrix.operate(trajectory[i-1]).add(noise);
    }
    
    run(trajectory, evolutionMatrix, evolutionVariances, "Constant");
        
    trajectory[0] = MatrixUtils.createRealVector(new double[] { 0 });
    
    RealVector shift = MatrixUtils.createRealVector(new double[] { 0.2 });
    
    for (int i=1; i<=n; i++) {
      trajectory[i] = trajectory[i-1].add( shift );
    }
    
    run(trajectory, evolutionMatrix, evolutionVariances, "Slope");

  }
  
  public static void run(RealVector[] trajectory, RealMatrix evolutionMatrix, CovarianceMatrix evolutionVariances, String fname) {
    int n = trajectory.length-1;
    
     double noiseStdDev     = 10;

    RealVector transformedControl =  MatrixUtils.createRealVector(new double[] { 0 });

    /*
     * generate observations.
     */
    
    RealVector[] observations = new RealVector[n+1];
    
    for (int i=0; i<=n; i++) {
      RealVector noise = MatrixUtils.createRealVector(new double[] { random.nextGaussian() }).mapMultiply(noiseStdDev); 
      observations[i] = trajectory[i].add(noise);
    }
    
    /*
     * Part three, estimate trajectory from observations.
     */
        
    RealMatrix observationMatrix = MatrixUtils.createRealMatrix(new double[][] { { 1 } });
    CovarianceMatrix observationVariances = new DiagonalCovarianceMatrix( new double[] { noiseStdDev }, DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS );
    
    UltimateKalman kalman = new PaigeSaundersKalman();
        
    RealVector[] filtered = new RealVector[n+1];
    RealVector[] filtstd   = new RealVector[n+1];
    
    for (int i=0; i<=n; i++) {
      kalman.advance(1);
      if (i>0) kalman.evolve(evolutionMatrix, transformedControl, evolutionVariances);
      if (i >= 0)
        kalman.observe(observationMatrix, observations[i], observationVariances);      
      else 
        kalman.observe();
      try { 
        filtered[i] = kalman.filtered();
        double cov = kalman.covariance().get().getEntry(0, 0);
        filtstd[i] = MatrixUtils.createRealVector(new double[] { Math.sqrt(cov) });
      } catch (InsufficientDataException isde) {}
      //System.out.printf("... %.17e %.17e\n", observations[i].getEntry(0),estimates[i].getEntry(0));
    }
    
    //System.out.printf("firstIndex=%d\n", three.firstIndex());
    //System.out.printf("currentIndex=%d\n", three.currentIndex());
    
    kalman.smooth();

    RealVector[] smoothed = new RealVector[n+1];
    RealVector[] smthstd  = new RealVector[n+1];

    for (int i=0; i<=n; i++) {
      smoothed[i] = kalman.smoothed(i);
      double cov = kalman.covariance().get().getEntry(0, 0);
      smthstd[i] = MatrixUtils.createRealVector(new double[] { Math.sqrt(cov) });
    }

    
    try {
      Matlab script = new Matlab("outputs/"+fname+".m");
      
      script.printIndexedMatrix(trajectory,"trajectory");
      script.printIndexedMatrix(observations,"observations");
      script.printIndexedMatrix(filtered,"filtered");
      script.printIndexedMatrix(filtstd,"filtstd");
      script.printIndexedMatrix(smoothed,"smoothed");
      script.printIndexedMatrix(smthstd,"smthstd");

      script.figure();
      script.printf("plot(trajectory(:,1),trajectory(:,2),'k-','LineWidth',1);\n");
      script.printf("plot(observations(:,1),observations(:,2),'r.');\n");
      script.printf("plot(filtered(:,1),filtered(:,2),'b-','LineWidth',1);\n");

      //script.printf("plot(filtered(:,1),filtered(:,2)+3*filtstd(:,2),'m-');\n");
      //script.printf("plot(filtered(:,1),filtered(:,2)-3*filtstd(:,2),'m-');\n");
      
      script.printf("xs = [filtered(:,1)', fliplr(filtered(:,1)')];\n");
      script.printf("ys = [filtered(:,2)'+3*filtstd(:,2)', fliplr(filtered(:,2)'-3*filtstd(:,2)')];\n");
      script.printf("fill(xs, ys,'k','FaceAlpha',0.2,'LineStyle','none');\n");

      //script.printf("plot(smoothed(:,1),smoothed(:,2),'g-');\n");
      script.printf("xlim([ min(observations(:,1)) max(observations(:,1)) ]);\n");    
      script.printf("ylim([ min(observations(:,2)) max(observations(:,2)) ]);\n");    
    
      script.save(fname+"Filtered");

      script.figure();
      script.printf("plot(trajectory(:,1),trajectory(:,2),'k-','LineWidth',1);\n");
      script.printf("plot(observations(:,1),observations(:,2),'r.');\n");
      script.printf("plot(smoothed(:,1),smoothed(:,2),'m-','LineWidth',1);\n");

      //script.printf("plot(filtered(:,1),filtered(:,2)+3*std(:,2),'m-');\n");
      //script.printf("plot(filtered(:,1),filtered(:,2)-3*std(:,2),'m-');\n");
      
      script.printf("xs = [smoothed(:,1)', fliplr(smoothed(:,1)')];\n");
      script.printf("ys = [smoothed(:,2)'+3*smthstd(:,2)', fliplr(smoothed(:,2)'-3*smthstd(:,2)')];\n");
      script.printf("fill(xs, ys,'k','FaceAlpha',0.2,'LineStyle','none');\n");

      //script.printf("plot(smoothed(:,1),smoothed(:,2),'g-');\n");
      script.printf("xlim([ min(observations(:,1)) max(observations(:,1)) ]);\n");    
      script.printf("ylim([ min(observations(:,2)) max(observations(:,2)) ]);\n");    
    
      script.save(fname+"Smoothed");
      
      
      //script.close();
      
    } catch (FileNotFoundException e) {
      System.err.println("File not found, exiting\n");
      System.exit(1);
    }    
  }
}
