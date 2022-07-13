package sivantoledo.kalman.tests;

import java.io.FileNotFoundException;
import java.io.PrintStream;

import org.apache.commons.math3.linear.RealVector;

public class Matlab {
  
  private final PrintStream ps;
  
  public Matlab(String fname) throws FileNotFoundException {
    this(new PrintStream(fname));
  }

  public Matlab(PrintStream ps) {
    this.ps = ps;
    //ps.printf("addpath 'g:\\My Drive\\Software\\matlab\\export_fig' -end;");
    ps.printf("clear all; close all;\n");
  }
  
  public void printf(String f, Object... p) {
    ps.printf(f, p);
  }
  
  public void figure() {
    ps.printf("figure\n");    
    ps.printf("axis square\n");  
    //ps.printf("set(gca,'FontName','Helvetica');\n");  
    //ps.printf("set(gca,'FontSize',16);\n");  
    //ps.printf("daspect([1 1 1]);\n");  
    ps.printf("set(gca,'Box','on');;\n");
    ps.printf("hold on;\n");  
  }
  
  public void close() {
    //ps.printf("set(gca,'xtick',[],'ytick',[]);\n");
    //ps.printf("set(gca,'xticklabel',[],'yticklabel',[]);\n");
    ps.printf("hold off;\n");
  }
  
  public void save(String fname) {
    //ps.printf("set(gca,'xtick',[],'ytick',[]);\n");
    //ps.printf("set(gca,'xticklabel',[],'yticklabel',[]);\n");
    ps.printf("hold off;\n");
    ps.printf("exportgraphics(gca,'%s.pdf');\n",fname);
  }

  public void printIndexedMatrix(RealVector[] rows, String label) {
    ps.printf("%s = [",label);
    int x = 1;
    for (RealVector row: rows)
      if (row!=null) ps.printf("\n%d %.3f",x++, row.getEntry(0));
      else           ps.printf("\n%d %.3f",x++, Double.NaN);
    ps.printf("];\n");
  
  }

  public void printVector(double[] a, String label) {
    ps.printf("%s = [",label);
    for (double d: a)
      ps.printf("\n%.17e",d);
    ps.printf("];\n");
  }

  public void printMatrix(double[][] a, String label) {
    ps.printf("%s = [",label);
    for (double[] row: a) {
      ps.printf("\n");
      for (double d: row) {
        ps.printf(" %.17e",d);
      }
    }
    ps.printf("];\n");
  }

}
