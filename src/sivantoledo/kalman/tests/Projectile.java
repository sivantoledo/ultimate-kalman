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

public class Projectile {
  
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
     * The initial state expectation is defined as Humpherys et al defined: position 
     * from the first observation, and velocity from the first 10 seconds of measurement.
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
    
    RealVector[] filtered = new RealVector[1300]; // stay on the safe side, we predict to hitting ground....
    
    filtered[400] = three.filtered();
    
    for (int i=401; i<=600; i++) {
      three.advance(4);
      three.evolve(evolutionMatrix, transformedControl, evolutionVariances);
      three.observe(observationMatrix, observations[i], observationVariances);      
      filtered[i] = three.filtered();
    }
    
    /*
     * Part four: continue to predict until the altitude drops below zero (the 
     * projectile hits ground).
     */
 
    for (int i=601; true ; i++) {
      three.advance(4);
      three.evolve(evolutionMatrix, transformedControl, evolutionVariances);
      three.observe();      
      filtered[i] = three.filtered();
      if (filtered[i].getEntry(1) <= 0) break;
    }
    
    /*
     * Part five: predict the point of departure. We simply smooth the entire
     * filter (we dropped nothing yet). We do not reverse the dynamical system.
     */
    
    three.smooth();
    
    RealVector[] smoothed = new RealVector[1201];
    
    for (int i=400-400; i<=600-400; i++) {
      smoothed[i+400] = three.smoothed(i);
    }
    
    /*
     * A more direct solution for parts 4 and 5.
     */
    
    UltimateKalman four = new PaigeSaundersKalman();
        
    for (int i=0; i<=1200; i++) {
      four.advance(4);
      if (i>0) four.evolve(evolutionMatrix, transformedControl, evolutionVariances);
      if (i>=400 && i<=600) four.observe(observationMatrix, observations[i], observationVariances); 
      else                  four.observe();
    }
    
    four.smooth();

    RealVector[] full = new RealVector[1201];

    for (int i=0; i<=1200; i++) {
      full[i] = four.smoothed(i);
    }
    
    System.out.printf("clear all; close all\n");

    matlabPrintMatrix(trajectory,"trajectory");
    matlabPrintMatrix(observations,"observations");
    matlabPrintMatrix(filtered,"filtered");
    matlabPrintMatrix(smoothed,"smoothed");
    matlabPrintMatrix(full,"full");
    
    //System.out.printf("subplot(2,2,1)\n");
    System.out.printf("figure\n");
    System.out.printf("plot(trajectory(:,1),trajectory(:,2),'k-');\n");
    System.out.printf("hold on;\n");
    System.out.printf("plot(filtered(:,1),filtered(:,2),'c+');\n");
    System.out.printf("plot(observations(:,1),observations(:,2),'r.');\n");
    System.out.printf("plot(filtered(:,1),filtered(:,2),'c-');\n");
    System.out.printf("plot(smoothed(:,1),smoothed(:,2),'g-');\n");
    System.out.printf("plot(full(:,1),full(:,2),'m-');\n");
    System.out.printf("hold off;\n");
    
    //System.out.printf("subplot(2,2,2)\n");
    System.out.printf("figure\n");
    System.out.printf("plot(trajectory(:,1),trajectory(:,2),'k-');\n");
    System.out.printf("hold on;\n");
    System.out.printf("plot(observations(:,1),observations(:,2),'r.');\n");
    System.out.printf("plot(filtered(:,1),filtered(:,2),'c+');\n");
    System.out.printf("plot(smoothed(:,1),smoothed(:,2),'g-');\n");
    System.out.printf("plot(full(:,1),full(:,2),'m-');\n");
    System.out.printf("xlim([ min(observations(:,1)) max(observations(:,1)) ]);\n");    
    System.out.printf("ylim([ min(observations(:,2)) max(observations(:,2)) ]);\n");    
    System.out.printf("hold off;\n");
    
   }
}
