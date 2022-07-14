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
  
  public static void main(String[] args) {
    addRemove();
  }

}
