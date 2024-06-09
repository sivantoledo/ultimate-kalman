package sivantoledo.kalman.examples;

import java.util.Random;

import org.apache.commons.math3.exception.InsufficientDataException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.CovarianceMatrix;
import sivantoledo.kalman.DiagonalCovarianceMatrix;
import sivantoledo.kalman.UltimateKalman;

public class Rotation {
    
  /**
   * A simple rotating point.
   * 
   * The argument gives the dimension of the observations, and should
   * be between 1 and 6.
   */
  
  public static void rotation(long seed, int GRowDim) {
    
    int k = 16;
    GRowDim = 2;
    
    Random random = new Random(seed);

    RealMatrix evolutionErrors = MatrixUtils.createRealMatrix(new double[][] {
      { -0.343003152130103,-0.766711794483284,-0.016814112314737, 0.684339759945504,-1.401783282955619,-1.521660304521858,-0.127785244107286, 0.602860572524585,-0.139677982915557, 0.407768714902350, 0.397539533883833,-0.317539749169638,-0.779285825610984,-1.935513755513929, 0.678730596165904 },
      { 1.666349045016822, 2.635481573310387, 0.304155468427342, 0.055808274805755,-1.360112379179931, 1.054743814037827,-1.410338023439304,-0.456929290517258,-0.983310072206319, 0.242994841538368,-0.175692485792199,-1.101615186229668,-1.762205119649466, 1.526915548584107,-2.277161011565906  }
    });

    RealMatrix observationErrors = MatrixUtils.createRealMatrix(new double[][] {
      { -1.428567988496096, 0.913205695955837,-1.576872295738796,-1.888336147279610, 1.116853507009928, 1.615888145666843,-0.102585012191329,-0.192732954692481, 0.160906008337421,-0.024849020282298,-1.001561909251739,-0.314462113181954,0.276865687293751, 0.175430340572582, 0.746792737753047, 1.648965874319728 },
      { -1.114618464565160, 0.976371425014641, 0.204080086636545, 0.736193913185726, 0.743379272133998,-1.666530392059792, 0.622727541956653, 0.794595441386172, 0.539084689771962,-2.548385761079745,-1.161623730001803, 1.066876935479899,1.748562141782206, 0.362976707912966, 0.842263598054067, 1.725578381396231 }
    });
    

    double alpha = 2 * Math.PI / 16; // angle
    
    RealMatrix F = MatrixUtils.createRealMatrix(new double[][] {
      {   Math.cos(alpha), -Math.sin(alpha) },
      {   Math.sin(alpha),  Math.cos(alpha) }
    });

    RealMatrix G = MatrixUtils.createRealMatrix(new double[][] {
      { 1, 0 },
      { 0, 1 },
      { 1, 1 },
      { 2, 1 },
      { 1, 2 },
      { 3, 1 }
    }).getSubMatrix(0, GRowDim-1, 0, 1);
    
    double evolutionStdDev   = 1e-3;
    double observationStdDev = 1e-1;
    
    CovarianceMatrix K = new DiagonalCovarianceMatrix(F.getRowDimension(), 
        evolutionStdDev, 
        DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS);

    CovarianceMatrix C = new DiagonalCovarianceMatrix(G.getRowDimension(),
        observationStdDev,
        DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS);

    RealVector[] states       = new RealVector[k];

    states[0] = MatrixUtils.createRealVector(new double[] { 1, 0 });
    
    for (int i=1; i<k; i++) {
      states[i] = F.operate(states[i-1]).add(evolutionErrors.getColumnVector(i-1).mapMultiply(evolutionStdDev));
    }
    
    printMatrix(states, "states", "%.4f");
    
    RealVector[] observations = new RealVector[k];

    for (int i=0; i<k; i++) {
      observations[i] = G.operate(states[i]).add(observationErrors.getColumnVector(i).mapMultiply(observationStdDev));
    }

    printMatrix(observations, "observations", "%.4f");
    
    RealVector zero = MatrixUtils.createRealVector(new double[] { 0, 0 });

    RealVector[] predicted = new RealVector[k];

    UltimateKalman kalman = new UltimateKalman();

    kalman.evolve(2);
    kalman.observe(G, observations[0], C);
    predicted[0] = kalman.estimate();
    
    for (int i=1; i<k; i++) {
      kalman.evolve(2, F, zero, K); // simplified interface without H
      kalman.observe(); // no observations in prediction mode of this example
      predicted[i] = kalman.estimate();
    }

    printMatrix(predicted, "predicted", "%.4f");
    
    kalman.rollback(1);
    
    System.out.printf("earliest->latest %d->%d (after rollback)\n",
                      kalman.earliest(), kalman.latest());
    
    kalman.observe(G, observations[1], C);
    
    RealVector[] filtered = new RealVector[k];

    filtered[0] = kalman.estimate(0);
    filtered[1] = kalman.estimate(1);
    
    for (int i=2; i<k; i++) {
      kalman.evolve(2, F, zero, K); // simplified interface without H
      kalman.observe(G, observations[i], C);
      filtered[i] = kalman.estimate();
    }
   
    printMatrix(filtered, "filtered", "%.4f");

    RealVector[] smoothed = new RealVector[k];

    kalman.smooth();
    
    for (int i=0; i<k; i++) {
      smoothed[i] = kalman.estimate(i);
    }
   
    printMatrix(smoothed, "smoothed", "%.4f");
    
    CovarianceMatrix E = kalman.covariance(0);
    printMatrix(E.get(), "Covariance of smoothed estimate of state 0","%.2e");
    printMatrix(E.weigh(MatrixUtils.createRealIdentityMatrix(2)),"W such that W^T*W=cov^-1","%.2e");
    System.out.printf("std deviation of first state coordinate %.2e (of first observation %.2e)\n",
                     Math.sqrt(E.get().getEntry(0, 0)), observationStdDev);
    

    //Simulation sim = new Simulation(1000, null, F, G, null, 1e-3, 0.1);
    //sim.simulate(k, MatrixUtils.createRealVector(new double[] { 1, 0 }));

    //Matlab script = new Matlab(System.out);
    
    //script.printMatrix(sim.states, "states");
    //script.printMatrix(sim.observations, "obs");
    /*
    RealVector[] filtered = new RealVector[k];
    RealVector[] smoothed = new RealVector[k];
    
    UltimateKalman kalman = new UltimateKalman();
    for (int i=0; i<k; i++) {
      kalman.advance(F.getRowDimension());
      if (i>0) kalman.evolve(sim.F, sim.be, sim.Ce);
      kalman.observe(sim.G, sim.observations[i], sim.Co);
      try {
        filtered[i] = kalman.filtered();
      } catch (InsufficientDataException isde) {
        System.err.printf("%% insufficient data for filtering in step %d\n", i);
      }
    }
    
    kalman.smooth();

    for (int i=0; i<k; i++) {
      smoothed[i] = kalman.smoothed(i);
    }

    script.printMatrix(filtered, "filtered");
    script.printMatrix(smoothed, "smoothed");
    script.figure();
    script.printf("plot(states(:,1),states(:,2),'k-','LineWidth',1);\n");
    script.printf("if size(obs,2)==2; plot(obs(:,1),obs(:,2),'r-'); end;\n");
    script.printf("plot(filtered(:,1),filtered(:,2),'b-','LineWidth',1);\n");
    script.printf("plot(smoothed(:,1),smoothed(:,2),'m-','LineWidth',1);\n");
    script.printf("xlim([-1.25 1.25]);\n");
    script.printf("ylim([-1.25 1.25]);\n");
    //script.printf("daspect([1 1 1]);\n");
    //if (GRowDim < 2) script.save("Underdetermined");
    //if (GRowDim > 2) script.save("Overdetermined");
  */
  }
  
  /*
  public void printMatrix(RealVector[] a, String label) {
    ps.printf("%s = [",label);
    for (RealVector row: a) {
      if (row==null) continue;
      ps.printf("\n");
      for (int i=0; i<row.getDimension(); i++) {
        ps.printf(" %.17e",row.getEntry(i));
      }
    }
    ps.printf("];\n");
  }

  public void printVector(RealVector row, String label) {
    ps.printf("%s = [",label);
    for (int i=0; i<row.getDimension(); i++) {
        ps.printf(" %.17e",row.getEntry(i));
    }
    ps.printf("];\n");
  }
*/
  public static void printMatrix(RealMatrix a, String label, String format) {
    printMatrix(a.getData(), label, format);
  }
  
  public static void printMatrix(double[][] a, String label, String format) {
    int n = 0;
    for (double[] row: a) if (row!=null) n = row.length;
    
    System.out.printf("%s = [",label);
    for (double[] row: a) {
      System.out.printf("\n");
      if (row==null) {
        for (int j=0; j<n; j++) System.out.printf(" NaN");
      } else {
        for (double d: row) {
          System.out.printf(" "+format,d);
        }
      }
    }
    System.out.printf("];\n");
  }

  // this assumes that the array holds column vectors of the same size
  public static void printMatrix(RealVector[] a, String label, String format) {
    int n = a.length;
    
    System.out.printf("%s = [",label);
    for (int i=0; i<a[0].getDimension(); i++) {
      System.out.printf("\n");
      for (int j=0; j<n; j++) {
        System.out.printf(" "+format,a[j].getEntry(i));
      }
    }
    System.out.printf("];\n");
  }

  
  public static void main(String[] args) {
    //dynamic();
    //simple(1);
    rotation(5,1);
  }

}
