package sivantoledo.kalman.tests;

import java.util.Arrays;
import java.util.Random;

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
 * Kalman filter example from the paper "A fresh look at the Kalman filter" by
 * Humpherys, Redd, and West, SIAM Review 54(4):801--823.
 * 
 * @author Sivan Toledo
 *
 */

public class FreshLook {
  
  public static void matlabPrintMatrix(RealVector[] rows, String label) {
    System.out.printf("%s = [",label);
    for (RealVector row: rows)
      if (row!=null) System.out.printf("\n%.3f %.3f",row.getEntry(0), row.getEntry(1));
    System.out.printf("]\n");
   
  }
  
  public static void main(String[] args) {
    
    double dt = 0.1;
    double g  = 9.8;
    double b  = 1e-4;
    
    double verticalAccel = -g;

    int d = 4;

    RealMatrix evolutionMatrix = MatrixUtils.createRealMatrix(new double[][] {
      { 1, 0,   dt,    0 }, 
      { 0, 1,    0,   dt },
      { 0, 0,  1-b,    0 },
      { 0, 0,    0,  1-b } });

    RealVector transformedControl =  MatrixUtils.createRealVector(new double[] { 0, 0, 0, -g*dt });

    /*
     * Part one, evolve the system for 1200 time steps.
     */
    
    Random r = new Random(69978);
    
    RealVector[] trajectory = new RealVector[1201];
    trajectory[0] = MatrixUtils.createRealVector(new double[] {0, 0, 300, 600 });
    
    for (int i=1; i<=1200; i++) {
      RealVector noise = MatrixUtils.createRealVector(new double[] { r.nextGaussian(), r.nextGaussian(), r.nextGaussian(), r.nextGaussian() })
                         .mapMultiply(Math.sqrt(0.1)); // variances of 0.1
      trajectory[i] = evolutionMatrix.operate(trajectory[i-1]).add(transformedControl).add(noise);
    }
            
    /*
     * Part two, generate observations for time steps 400 to 600.
     */
    
    RealVector[] observations = new RealVector[1201];
    
    for (int i=400; i<=600; i++) {
      RealVector noise = MatrixUtils.createRealVector(new double[] { r.nextGaussian(), r.nextGaussian() })
                         .mapMultiply(Math.sqrt(5000)); // variances of 0.1
      observations[i] = trajectory[i].getSubVector(0, 2).add(noise);
    }
    
    /*
     * Part three, estimate trajectory from observations.
     */

    
    RealMatrix initialObservationMatrix  = MatrixUtils.createRealIdentityMatrix(4);
    RealVector initialStateExpectation   = MatrixUtils.createRealVector(new double[] {
        observations[400].getEntry(0), 
        observations[400].getEntry(1), 
        (observations[410].getEntry(0) - observations[400].getEntry(0)) / (10*dt), 
        (observations[410].getEntry(1) - observations[400].getEntry(1)) / (10*dt) });
    CovarianceMatrix initialObsVariances = new DiagonalCovarianceMatrix( new double[] { 1e5, 1e5, 1e5, 1e5 }, 
        DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );

    CovarianceMatrix evolutionVariances = new DiagonalCovarianceMatrix( new double[] { 0.1, 0.1, 0.1, 0.1 }, 
        DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );
    
    RealMatrix observationMatrix = MatrixUtils.createRealMatrix(new double[][] {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 } });
    CovarianceMatrix observationVariances = new DiagonalCovarianceMatrix( new double[] { 5000, 5000 }, 
        DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );
    
    UltimateKalman three = new PaigeSaundersKalman();
    
    three.advance(4);
    three.observe(initialObservationMatrix, initialStateExpectation, initialObsVariances);
    
    RealVector[] filtered = new RealVector[1201];
    
    filtered[400] = three.filtered();
    
    for (int i=401; i<=600; i++) {
      three.advance(4);
      three.evolve(evolutionMatrix, transformedControl, evolutionVariances);
      three.observe(observationMatrix, observations[i], observationVariances);      
      filtered[i] = three.filtered();
    }
    
    three.smooth();
    
    RealVector[] smoothed = new RealVector[1201];
    
    for (int i=400-400; i<=600-400; i++) {
      smoothed[i+400] = three.smoothed(i);
    }
    
    System.out.printf("clear all; close all\n");

    matlabPrintMatrix(trajectory,"trajectory");
    matlabPrintMatrix(observations,"observations");
    matlabPrintMatrix(filtered,"filtered");
    matlabPrintMatrix(smoothed,"smoothed");
    
    System.out.printf("subplot(2,2,1)\n");
    System.out.printf("plot(trajectory(:,1),trajectory(:,2),'k-');\n");
    System.out.printf("hold on;\n");
    System.out.printf("plot(observations(:,1),observations(:,2),'r.');\n");
    System.out.printf("hold off;\n");
    
    System.out.printf("subplot(2,2,2)\n");
    System.out.printf("plot(trajectory(:,1),trajectory(:,2),'k-');\n");
    System.out.printf("hold on;\n");
    System.out.printf("plot(observations(:,1),observations(:,2),'r.');\n");
    System.out.printf("plot(filtered(:,1),filtered(:,2),'b-');\n");
    System.out.printf("plot(smoothed(:,1),smoothed(:,2),'g-');\n");
    System.out.printf("xlim([ min(observations(:,1)) max(observations(:,1)) ]);\n");    
    System.out.printf("ylim([ min(observations(:,2)) max(observations(:,2)) ]);\n");    
    System.out.printf("hold off;\n");
    
   }
}
