package sivantoledo.kalman.examples;

import java.util.Random;

import org.apache.commons.math3.exception.InsufficientDataException;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.CovarianceMatrix;
import sivantoledo.kalman.DiagonalCovarianceMatrix;
import sivantoledo.kalman.PaigeSaundersKalman;
import sivantoledo.kalman.UltimateKalman;

public class VaryingVariances {
  
  private static Random random = new Random(69978);
  
  
  /**
   * Tracking a constant, but with a middle observation that is more accurate
   */
  
  public static void varyingVariance() {
    
    RealMatrix F = MatrixUtils.createRealMatrix(new double[][] {
      { 1 }
    });

    RealMatrix G = MatrixUtils.createRealMatrix(new double[][] {
      { 1 }
    });
    
    int k=20;
    int mid = k/2;
    
    Simulation sim = new Simulation(k, null, F, G, null, 1e-3, 1);
    
    // simulate for 10 steps with these parameters
    sim.simulate(mid, MatrixUtils.createRealVector(new double[] { 0 }));
   
    double normal = sim.observationStdDev;
    double tight = 0.1;
    sim.observationStdDev = tight;
    sim.simulate(1);

    sim.observationStdDev = normal;
    sim.simulate(k-mid-1);
    
    Matlab script = new Matlab(System.out);
    
    script.printMatrix(sim.states, "states");
    script.printMatrix(sim.observations, "obs");
    
    RealVector[] filtered = new RealVector[k];
    RealVector[] smoothed = new RealVector[k];
    RealVector[] filtvar  = new RealVector[k];
    RealVector[] smthvar  = new RealVector[k];
    
    UltimateKalman kalman = new PaigeSaundersKalman();
    for (int i=0; i<k; i++) {
      kalman.advance(F.getRowDimension());
      if (i>0) kalman.evolve(sim.F, sim.be, sim.Ce);
      if (i==mid) kalman.observe(sim.G, sim.observations[i], new DiagonalCovarianceMatrix(1,tight, DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS));
      else        kalman.observe(sim.G, sim.observations[i], new DiagonalCovarianceMatrix(1,normal,DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS));
      try {
        filtered[i] = kalman.filtered();
        filtvar[i]  = kalman.covariance().get().getRowVector(0);
      } catch (InsufficientDataException isde) {
        System.err.printf("%% insufficient data for filtering in step %d\n", i);
      }
    }
    
    kalman.smooth();

    for (int i=0; i<k; i++) {
      smoothed[i] = kalman.smoothed(i);
      smthvar[i]  = kalman.covariance(i).get().getRowVector(0);
    }

    script.printMatrix(filtered, "filtered");
    script.printMatrix(smoothed, "smoothed");
    script.printMatrix(filtvar , "filtvar");
    script.printMatrix(smthvar , "smthvar");
    script.figure();
    
    script.printf("xs = [0:%d, fliplr(0:%d)];\n",k-1,k-1);
    script.printf("ys = [filtered'+3*(filtvar.^0.5)', fliplr(filtered'-3*(filtvar.^0.5)')];\n");
    script.printf("fill(xs, ys,'k','FaceAlpha',0.2,'LineStyle','none');\n");

    //script.printf("plot(smoothed(:,1),smoothed(:,2),'g-');\n");
    script.printf("xlim([ 0 %d ]);\n",k-1);    
    //script.printf("ylim([ min(obs) max(obs) ]);\n");    

    script.printf("plot(0:%d,states(:,1),'k-','LineWidth',1);\n",k-1);
    script.printf("plot(0:%d,obs(:,1),'r.');\n",k-1);
    script.printf("plot(0:%d,filtered(:,1),'b-','LineWidth',1);\n",k-1);
    script.printf("plot(0:%d,smoothed(:,1),'m-','LineWidth',1);\n",k-1);
    script.save("VaryingVariances");
  }

  public static void main(String[] args) {
    //dynamic();
    //simple(1);
    varyingVariance();
  }

}
