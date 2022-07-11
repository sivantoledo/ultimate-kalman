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
 * Kalman filter example for a constant value.
 * 
 * @author Sivan Toledo
 *
 */

public class Constant {
  
  public static void matlabPrintMatrix(RealVector[] rows, String label) {
    System.out.printf("%s = [",label);
    int x = 1;
    for (RealVector row: rows)
      if (row!=null) System.out.printf("\n%d %.3f",x++, row.getEntry(0));
    System.out.printf("]\n");
   
  }
  
  public static void main(String[] args) {
    
    int n = 100;
    
    int d = 1;

    RealMatrix evolutionMatrix = MatrixUtils.createRealMatrix(new double[][] {{ 1 } });

    RealVector transformedControl =  MatrixUtils.createRealVector(new double[] { 0 });

    /*
     * Part one, evolve the system for n time steps.
     */
    
    Random r = new Random(69978);
    
    RealVector[] trajectory = new RealVector[n+1];
    trajectory[0] = MatrixUtils.createRealVector(new double[] { 0 });
    
    for (int i=1; i<=n; i++) {
      RealVector noise = MatrixUtils.createRealVector(new double[] { r.nextGaussian() }).mapMultiply(1);
      trajectory[i] = evolutionMatrix.operate(trajectory[i-1]).add(transformedControl).add(noise);
    }
            
    /*
     * Part two, generate observations.
     */
    
    RealVector[] observations = new RealVector[n+1];
    
    for (int i=0; i<=n; i++) {
      RealVector noise = MatrixUtils.createRealVector(new double[] { r.nextGaussian() }).mapMultiply(10); 
      observations[i] = trajectory[i].add(noise);
    }
    
    /*
     * Part three, estimate trajectory from observations.
     */
    
    CovarianceMatrix evolutionVariances = new DiagonalCovarianceMatrix( new double[] { 1 }, DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS );
    
    RealMatrix observationMatrix = MatrixUtils.createRealMatrix(new double[][] { { 1 } });
    CovarianceMatrix observationVariances = new DiagonalCovarianceMatrix( new double[] { 10 }, DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS );
    
    UltimateKalman three = new PaigeSaundersKalman();
        
    RealVector[] filtered = new RealVector[n+1];
    
    for (int i=0; i<=n; i++) {
      three.advance(1);
      if (i>0) three.evolve(evolutionMatrix, transformedControl, evolutionVariances);
      three.observe(observationMatrix, observations[i], observationVariances);      
      filtered[i] = three.filtered();
      //System.out.printf("... %.17e %.17e\n", observations[i].getEntry(0),estimates[i].getEntry(0));
    }
    
    //System.out.printf("firstIndex=%d\n", three.firstIndex());
    //System.out.printf("currentIndex=%d\n", three.currentIndex());
    
    three.smooth();

    RealVector[] smoothed = new RealVector[n+1];

    for (int i=0; i<=n; i++) {
      smoothed[i] = three.smoothed(i);
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
