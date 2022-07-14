package sivantoledo.kalman.tests;

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
 * Kalman filter example for the offsets of a groups of clocks in a distributed system.
 * 
 * @author Sivan Toledo
 *
 */

public class ClockOffsets {
    
  public static double[] constant(int n, double value) {
    double[] a = new double[ n ];
    for (int i=0; i<n; i++) a[i] = value;
    return a;
  }

  public static double[][] constant(int m, int n, double value) {
    double[][] a = new double[ m ][ n ];
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) a[i][j] = value;
    return a;
  }

  private static Random random = new Random(69978);
  //     Random random = new Random(System.currentTimeMillis());
  
  public static void main(String[] args) {
    
    //Random random = new Random(69978);
    
    int clockCount = 3;
    int stepCount  = 100;
    
    double evolutionStdDev     =  10e-9;
    double observationStdDev   = 100e-9;
    double initialOffsetStdDev =   1e-6;
    double delayStdDev         =  10e-6;
    
    /*
     * Simulation
     */
    
    double[][] t   = new double[ stepCount ][ clockCount ];
    double[]   tau = new double[ stepCount ];
    double  [] d   = new double             [ clockCount ];
    double[][] o   = new double[ stepCount ][ clockCount ];
    double  [] r   = new double             [ clockCount ]; // rate errors
    
    for (int i=0; i<stepCount;  i++) tau[i]  = i;
    for (int j=0; j<clockCount; j++) d[j]    = delayStdDev         * random.nextGaussian();
    for (int j=0; j<clockCount; j++) r[j]    = evolutionStdDev     * random.nextGaussian();
    for (int j=0; j<clockCount; j++) o[0][j] = initialOffsetStdDev * random.nextGaussian();
    //o[0][0] = 0; // constraint to be able to plot absolute offsets rather than relative ones

    for (int i=1; i<stepCount;  i++) {
      for (int j=0; j<clockCount; j++) {
        o[i][j] = o[i-1][j] + evolutionStdDev * random.nextGaussian();
        o[i][j] = o[i-1][j] + r[j];
      }
    }

    for (int i=0; i<stepCount;  i++) {
      for (int j=0; j<clockCount; j++) {
        if (i==0 || random.nextDouble()<=0.90) // the first step must have at least 2 receivers...
          t[i][j] = tau[i] + d[j] + o[i][j] + observationStdDev * random.nextGaussian();
        else 
          //t[i][j] = tau[i] + d[j] + o[i][j] + observationStdDev * random.nextGaussian();
          t[i][j] = Double.NaN; // the receiver did not hear the packet
      }
    }
    
    /*
     * Estimate using UltimateKalman
     */
    
    RealMatrix evolutionMatrix = MatrixUtils.createRealMatrix(clockCount,clockCount+1); // also tau[i]
    for (int j=0; j<clockCount; j++) evolutionMatrix.setEntry(j, j, 1);

    CovarianceMatrix evolutionVariances = new DiagonalCovarianceMatrix( constant(clockCount,evolutionStdDev), DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS );

    RealVector transformedControl =  MatrixUtils.createRealVector(constant(clockCount,0));
    
    
    UltimateKalman kalman = new PaigeSaundersKalman();
        
    double[][] filtered = new double[stepCount][];
    
    for (int i=0; i<stepCount; i++) {
      int dim = clockCount+1;
      kalman.advance(dim);
      if (i>0) kalman.evolve(evolutionMatrix, transformedControl, evolutionVariances);

      RealVector observations = null;
      RealMatrix observationMatrix = null;
      CovarianceMatrix observationVariances = null;
      
      int obsCount = 0;
      for (int j=0; j<clockCount; j++) if (!Double.isNaN(t[i][j])) obsCount++;
    

      if (obsCount >= 1) {
        if (i==0) {
          // the last entry is the extra constaint, equal to zero
          observationMatrix = MatrixUtils.createRealMatrix     (obsCount+1,clockCount+1);
          observationMatrix.setEntry(obsCount,0, 1); // coefficient of o[0][0]
          observations =  MatrixUtils.createRealVector(constant(obsCount+1,0));
        } else {
          observationMatrix = MatrixUtils.createRealMatrix     (obsCount,clockCount+1);
          observations =  MatrixUtils.createRealVector(constant(obsCount,0));
        }
        int obsIndex = 0;
        for (int j=0; j<clockCount; j++) {
          if (Double.isNaN(t[i][j])) continue;
          observations.setEntry(obsIndex, t[i][j]-d[j]);
          observationMatrix.setEntry(obsIndex,j         , 1); // coefficient of o[i][j]
          observationMatrix.setEntry(obsIndex,clockCount, 1); // coefficient of tau[i]   
          obsIndex++;
        }

        observationVariances = new DiagonalCovarianceMatrix( constant(observationMatrix.getRowDimension(),observationStdDev), DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS );        
      
        System.out.printf("step %d, %d receivers, obs %d by %d\n", i,clockCount, observationMatrix.getRowDimension(), observationMatrix.getColumnDimension());
      } else {        
        System.out.printf("step %d, %d receivers, obsCount %d\n", i,clockCount, obsCount);
      }
      
      if (obsCount<=0) 
        kalman.observe(); // can't estimate tau_i from one constraint...
      else
        kalman.observe(observationMatrix, observations, observationVariances); 
      try {
        filtered[i] = kalman.filtered().toArray();
      } catch (InsufficientDataException e) {}
    }
    
    //kalman.smooth();

    //for (int i=0; i<stepCount; i++) {
      //filtered[i] = kalman.smoothed(i).toArray();
    //}

    try {
      Matlab script = new Matlab("matlab/ClockOffsets.m");
      script.printMatrix(o,"o");
      script.printVector(tau,"tau");
      script.printMatrix(t,"t");
      script.printMatrix(filtered,"filtered");
      
      script.figure();
    
      script.printf("plot(tau,o-o(:,1));\n");
      script.printf("set(gca,'ColorOrderIndex',1);\n");
      script.printf("plot(tau,filtered(:,1:%d)-filtered(:,1),'--');\n",clockCount); // don't plot th eestimates of tau....
      
      script.save("ClockOffsets");
      script.close();
      
    } catch (FileNotFoundException e) {
      System.err.println("File not found, exiting\n");
      System.exit(1);
    }

  }
}
