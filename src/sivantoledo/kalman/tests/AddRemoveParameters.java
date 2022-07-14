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

public class AddRemoveParameters {
  
  private static Random random = new Random(69978);
  
  
  /**
   * An example that shows how to add and drop parameters from the model.
   */
  
  public static void addRemove() {
    
    Matlab script = new Matlab(System.out);
    
    UltimateKalman kalman = new PaigeSaundersKalman();
    
    CovarianceMatrix C_F = new DiagonalCovarianceMatrix( 1, 1.0, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );
    RealMatrix F = MatrixUtils.createRealMatrix(new double[][] {
      { 1 } });
    RealVector b_e = MatrixUtils.createRealVector(new double[] { 0 });

   
    RealMatrix G = MatrixUtils.createRealMatrix(new double[][] {
        { 1 } });
    CovarianceMatrix C_G = new DiagonalCovarianceMatrix( 1, 1.0, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );

    RealVector u_hat; 
    RealVector b_o;
    
    // 1
    
    kalman.advance(1);
    b_o = MatrixUtils.createRealVector(new double[] { 10+random.nextGaussian() });
    kalman.observe(G, b_o, C_G);
    u_hat = kalman.filtered();
    
    script.printVector(b_o,"o1");
    script.printVector(u_hat,"u1");

    // 2

    kalman.advance(1);
    kalman.evolve(F, b_e, C_F);
    b_o = MatrixUtils.createRealVector(new double[] { 10+random.nextGaussian() });
    kalman.observe(G, b_o, C_G);
    u_hat = kalman.filtered();
    
    script.printVector(b_o,"o2");
    script.printVector(u_hat,"u2");
    
    // 3

    kalman.advance(2);
    kalman.evolve(F, b_e, C_F); // evolution only on the existing parameter
    b_o = MatrixUtils.createRealVector(new double[] { 10+random.nextGaussian(), 20+random.nextGaussian() });
    
    G = MatrixUtils.createRealIdentityMatrix(2);
    C_G = new DiagonalCovarianceMatrix( 2, 1.0, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );

    kalman.observe(G, b_o, C_G);
    u_hat = kalman.filtered();
    
    script.printVector(b_o,"o3");
    script.printVector(u_hat,"u3");

    // 4

    kalman.advance(2);
    
    b_e = MatrixUtils.createRealVector(new double[] { 0, 0 });
    F = MatrixUtils.createRealIdentityMatrix(2);
    C_F = new DiagonalCovarianceMatrix( 2, 1.0, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );
    kalman.evolve(F, b_e, C_F); // evolution only on the existing parameter

    b_o = MatrixUtils.createRealVector(new double[] { 10+random.nextGaussian(), 20+random.nextGaussian() });
    kalman.observe(G, b_o, C_G);

    u_hat = kalman.filtered();
    
    script.printVector(b_o,"o4");
    script.printVector(u_hat,"u4");
    
    // 5

    kalman.advance(1);
    
    //b_e = MatrixUtils.createRealVector(new double[] { 0, 0 });
    //F = MatrixUtils.createRealIdentityMatrix(2);
    //C_F = new DiagonalCovarianceMatrix( new double[] { 1, 1 }, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );
    RealMatrix H = MatrixUtils.createRealMatrix(new double[][] {
      { 0 }, 
      { 1 }, 
      });
    
    kalman.evolve(H, F, b_e, C_F); // evolution only on the existing parameter

    G = MatrixUtils.createRealIdentityMatrix(1);
    C_G = new DiagonalCovarianceMatrix( 1, 1.0, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );

    b_o = MatrixUtils.createRealVector(new double[] { 20+random.nextGaussian() });
    kalman.observe(G, b_o, C_G);

    u_hat = kalman.filtered();
    
    script.printVector(b_o,"o5");
    script.printVector(u_hat,"u5");

    // 6

    kalman.advance(1);

    C_F = new DiagonalCovarianceMatrix( 1, 1.0, DiagonalCovarianceMatrix.Representation.DIAGONAL_VARIANCES );
    F = MatrixUtils.createRealMatrix(new double[][] {
      { 1 } });
    b_e = MatrixUtils.createRealVector(new double[] { 0 });
    kalman.evolve(F, b_e, C_F);
    
    b_o = MatrixUtils.createRealVector(new double[] { 20+random.nextGaussian() });
    kalman.observe(G, b_o, C_G);
    
    u_hat = kalman.filtered();
    
    script.printVector(b_o,"o6");
    script.printVector(u_hat,"u6");

  }
  
  /**
   * A simple rotating point.
   * 
   * The argument gives the dimension of the observations, and should
   * be between 1 and 6.
   */
  
  public static void simple(int GRowDim) {
    
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
    
    int k = 17;
    
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
    script.printf("plot(states(:,1),states(:,2),'k-');\n");
    script.printf("if size(obs,2)==2; plot(obs(:,1),obs(:,2),'r-'); end;\n");
    script.printf("plot(filtered(:,1),filtered(:,2),'b-');\n");
    script.printf("plot(smoothed(:,1),smoothed(:,2),'m-');\n");

  }
  
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
