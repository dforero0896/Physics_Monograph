//g++ -o gsl_interpolation.o gsl_interpolation.cpp `gsl-config --cflags --libs`           
//./gsl_interpolation.o > interp.dat                                                      
//graph -T ps < interp.dat > interp.ps     
#include <iostream>
using namespace std;
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
int main(){
  int i;
  // Any function.
  double x[10], y[10];
  printf ("#m=0,S=2\n");
  for(i=0;i<10;i++){
    x[i]=0.2*i;
    y[i]=x[i]*x[i]+cos(x[i]);
    printf ("%g %g\n", x[i], y[i]);
  }
  printf ("#m=1,S=0\n");
  {
  //Define accelerator object
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  //Define spline object and allocate necessary memory for CUBIC spline
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 10);
  //Initialize spline with x and y arrays and number of points
  gsl_spline_init(spline, x, y, 10);
  //Evaluate spline
  double x_int, y_int;
  for(x_int=x[0]; x_int<x[9];x_int+=0.01){
    y_int=gsl_spline_eval(spline, x_int, acc);
    printf ("%g %g\n", x_int, y_int);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
  return 0;
}
