package sivantoledo.kalman.examples;

import java.util.Random;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import sivantoledo.kalman.CovarianceMatrix;
import sivantoledo.kalman.DiagonalCovarianceMatrix;

/**
 * This class simulates a dynamical system
 * 
 * @author Sivan Toledo
 *
 */
public class Simulation {
  
  public Random random = new Random();

  public RealMatrix H;
  public RealMatrix F;
  public RealMatrix G;
  public RealVector be;
  public double evolutionStdDev;
  public double observationStdDev;
  
  public RealVector[] states;
  public RealVector[] observations;
  
  public CovarianceMatrix Ce;
  public CovarianceMatrix Co;
  
  public int lastStep = -1; // no steps yet
  
  public Simulation(int maxSteps, RealMatrix H, RealMatrix F, RealMatrix G, RealVector be, double evolutionStdDev, double observationStdDev) {
    this.H = H;
    this.F = F;
    this.G = G;
    this.be = be;
    this.evolutionStdDev   = evolutionStdDev;
    this.observationStdDev = observationStdDev;  
    
    if (this.be == null) {
      double[] z = new double[ F.getRowDimension() ];
      this.be = MatrixUtils.createRealVector(z);
    }
    
    states = new RealVector[ maxSteps ];
    observations = new RealVector[ maxSteps ];
    
    Ce = new DiagonalCovarianceMatrix(F.getRowDimension(), 
                                      evolutionStdDev, 
                                      DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS);

    Co = new DiagonalCovarianceMatrix(G.getRowDimension(),
                                      observationStdDev,
                                      DiagonalCovarianceMatrix.Representation.DIAGONAL_STANDARD_DEVIATIONS);
  }
  
  public RealVector noise(int dim, double std) {
    double[] a = new double[ dim ];
    for (int i=0; i<dim; i++) a[i] = std*random.nextGaussian();
    return MatrixUtils.createRealVector(a);
  }
  
  public void simulate(int k, RealVector initialState) {
    assert(lastStep == -1);
    states[0] = initialState.copy();
    observations[0] = G.operate(states[0]).add( noise(G.getRowDimension(), observationStdDev) );
    //System.out.printf("sim 0\n");
    lastStep++;
    simulate(k-1);
  }
  
  public void simulate(int k) {
    for (int i=lastStep+1; i<=lastStep+k; i++) {
      //System.out.printf("sim %d\n",i);
      states[i] = F.operate(states[i-1]).add(be).add( noise(F.getRowDimension(), evolutionStdDev) );
      observations[i] = G.operate(states[i]).add( noise(G.getRowDimension(), observationStdDev) );
    }
    lastStep += k;
  }

  
  public CovarianceMatrix evolutionCov()   { return Ce; }
  public CovarianceMatrix observationCov() { return Co; }
}
