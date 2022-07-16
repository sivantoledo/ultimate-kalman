package sivantoledo.kalman.tests;

import java.util.Random;

import org.apache.commons.math3.exception.InsufficientDataException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.CovarianceMatrix;
import sivantoledo.kalman.DiagonalCovarianceMatrix;
import sivantoledo.kalman.PaigeSaundersKalman;
import sivantoledo.kalman.UltimateKalman;

public class Rotation {
    
  /**
   * A simple rotating point.
   * 
   * The argument gives the dimension of the observations, and should
   * be between 1 and 6.
   */
  
  public static void rotation(long seed, int GRowDim) {
    
    Random random = new Random(seed);
    
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
    
    Simulation sim = new Simulation(1000, null, F, G, null, 1e-3, 0.1);
    
    int k = 16;
    
    sim.simulate(k, MatrixUtils.createRealVector(new double[] { 1, 0 }));

    Matlab script = new Matlab(System.out);
    
    script.printMatrix(sim.states, "states");
    script.printMatrix(sim.observations, "obs");
    
    RealVector[] filtered = new RealVector[k];
    RealVector[] smoothed = new RealVector[k];
    
    UltimateKalman kalman = new PaigeSaundersKalman();
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

  }
  
  public static void main(String[] args) {
    //dynamic();
    //simple(1);
    rotation(5,1);
  }

}
