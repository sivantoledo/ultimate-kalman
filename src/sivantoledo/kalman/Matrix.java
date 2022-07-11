package sivantoledo.kalman;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
//import org.apache.logging.log4j.Level;
//import org.apache.logging.log4j.LogManager;
//import org.apache.logging.log4j.Logger;

public class Matrix {
  
  //private final static Logger log = LogManager.getLogger();

  public static double[] confidenceEllipse(double[][] cov, double fraction) {
    double[][] V = new double[2][2];
    double  [] D = new double[2];
    eig2x2(cov,V,D);
    
    double angle;
    double[] half_axes_angle = new double[3];
    if (D[0] > D[1]) {
      angle = Math.atan(V[0][0]/V[1][0]);
      half_axes_angle[0] = Math.sqrt(5.991*D[0]);
      half_axes_angle[1] = Math.sqrt(5.991*D[1]);
    } else {
      angle = Math.atan(V[0][1]/V[1][1]);
      half_axes_angle[0] = Math.sqrt(5.991*D[1]);
      half_axes_angle[1] = Math.sqrt(5.991*D[0]);
    }
    
    double angle_degree_north = 360*(angle/(2*Math.PI)); 
    if (angle_degree_north < 0)   angle_degree_north += 360;
    if (angle_degree_north > 180) angle_degree_north -= 180; // ellipse is symmetric
    half_axes_angle[2] = angle_degree_north;
    
    //System.out.printf("angle %.2f half-axes=%.1f %.1f\n",half_axes_angle[2],half_axes_angle[0],half_axes_angle[1]);
    
    return half_axes_angle;
  }
  
  public static void eig2x2(double[][] A, double[][]V, double[] D) {
    /*
    double y = (A[0][0] - A[1][1])/2;
    double p = A[0][1];
    double p_squared = p*p;
    double d = Math.abs(y) + Math.sqrt(y*y + p_squared);
    double r = Math.sqrt(d*d + p_squared);
    double c = d/r;
    double s = p/r;
    double t = p_squared/r;
    if (y<0) { s=-s; t=-t; }
    
    V[0][0] =  c;
    V[0][1] =  s;
    V[1][0] = -s;
    V[1][1] =  c;
    
    //System.out.printf("c=%.3f s=%.3f\n", c,s);
    //print(V,"V");
    
    //print(MatrixUtils.createRealMatrix(A).multiply(MatrixUtils.createRealMatrix(V)),"A*V");
    D[0] = c*c*A[0][0] - 2*s*c*A[0][1] + s*s*A[1][1];
    D[1] = s*s*A[0][0] + 2*s*c*A[0][1] + c*c*A[1][1];
    
    //print(D,"D");
    log.printf(Level.TRACE, "eig2x2 => [%.1f %.1f] [%.1f %.1f]", D[0]   ,D[1]   ,V[0][0],V[1][0]);
     */    
    
    try {
      EigenDecomposition eig = new EigenDecomposition(MatrixUtils.createRealMatrix(A));
      D[0] = eig.getD().getEntry(0, 0);
      D[1] = eig.getD().getEntry(1, 1);
      V[0][0] = eig.getV().getEntry(0, 0);
      V[0][1] = eig.getV().getEntry(0, 1);
      V[1][0] = eig.getV().getEntry(1, 0);
      V[1][1] = eig.getV().getEntry(1, 1);
      //log.printf(Level.TRACE, "eig2x2 => [%.1f %.1f] [%.1f %.1f]", D[0]   ,D[1]   ,V[0][0],V[1][0]);
    } catch (Exception e) { // the excpetion I have seen is MaxCountExceededException
      D[0] = D[1] = V[0][0] = V[0][1] = V[1][0] = V[1][1] = Double.NaN;      
    }
  }

  public static void print(double[] v, String label) {
    System.out.printf("%s = [ ",label);
    for (int i=0; i<v.length; i++) System.out.printf("%.2e ",v[i]);
    System.out.printf("]\n");
  }

  public static boolean isIdentical(double[] stored_v, double[] new_v) {
    if (stored_v == null || new_v == null) return false;
    if (stored_v.length != new_v.length)   return false;
    for (int i=0; i<new_v.length; i++) 
      if (new_v[i] != stored_v[i]) return false;
    return true;
  }
  
  public static void copy(double[] dest, double[] src) {
    for (int i=0; i<dest.length; i++) dest[i] = src[i];
  }
  
  public static double norm(double[] input) {
    double sum_of_squares = 0;
    for (int i=0; i<input.length; i++) sum_of_squares += input[i]*input[i];
    return Math.sqrt(sum_of_squares);
  }

  public static double normMax(double[] input) {
    double max = Double.NEGATIVE_INFINITY;
    for (int i=0; i<input.length; i++) {
      double a = Math.abs(input[i]);
      if (a > max) max = a;
    }
    return max;
  }
  
  public static double dotProduct(double[] x, double[] y) {
    double d = 0;
    for (int i=0; i<x.length; i++) d += x[i]*y[i];
    return d;
  }

  public static void print(RealMatrix A, String label) {
    System.out.printf("%s = [\n",label);
    for (int d=0; d<A.getRowDimension(); d++) {
      System.out.printf("    ");
      for (int i=0; i<A.getColumnDimension(); i++) {
        System.out.printf("%.2e ", A.getEntry(d,i));
      }
      System.out.printf("\n");
    }
    System.out.printf("  ]\n");
  }

  public static void print(double[][] A, String label) {
    System.out.printf("%s = [\n",label);
    for (int d=0; d<A.length; d++) {
      System.out.printf("    ");
      for (int i=0; i<A[d].length; i++) {
        System.out.printf("%.2e ", A[d][i]);
      }
      System.out.printf("\n");
    }
    System.out.printf("  ]\n");
  }
  
  public static String toString(double[][] A, String format) {
    StringBuilder s = new StringBuilder();
    s.append('[');
    for (int d=0; d<A.length; d++) {
      s.append('[');
      for (int i=0; i<A[d].length; i++) {
        s.append(String.format(format, A[d][i]));
      }
      s.append(']');
    }
    s.append(']');
    return s.toString();
  }

  
  public static void main(String[] args) {
    double[][] A = { { 4, -1 },
                     { -1, 100 } };
    double[][] V = new double[2][2];
    double[]   d = new double   [2];
    
    //eig2x2(A,V,d);
    
    confidenceEllipse(A,0.95);
  }
}
