/*
g++ -fopenmp -o prob_weight.o earth_simul.cpp prob_weight.cpp `gsl-config --cflags --libs`
*/

#include "earth_simul.h"
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
gsl_spline *spectrum_cspline;
gsl_interp_accel *accel;
//This program expects Emin, Emax, Steps_in_energy, dist, bse_model, i, k

int main(int argc, char const *argv[]) {
  if(argc!=8){
    cout << "ERROR\nInvalid number of arguments passed." << endl;
    return 1;
  }
  string dist, bse_model;
  dist=argv[4];
  bse_model=argv[5];
  float E_min, E_max;
  int Steps;
  E_min = atof(argv[1]);
  E_max = atof(argv[2]);
  Steps = atoi(argv[3]);
  if(E_min<0.0005||E_max>4.5005){
    cout << "ERROR\nChoose energy values in [0.0005, 4.5005] MeV" << endl;
    return 1;
  }
  int i,k;
  i=atoi(argv[6]);
  k=atoi(argv[7]);
  vector<float> Energy_linspace=linspace(E_min, E_max, Steps);
  Planet *earth = new Planet();
  earth->Planet::initializeCoords(0);
  earth->Planet::initializeDensity();
  earth->Planet::initializeAbundanceCrust();
  earth->Planet::initializeAbundanceMantle(dist, bse_model);
  float factorU, factorTh;
  factorU=(6.)*(earth->asArray[i][k].abundanceU*1e-9)*(0.9927)*(4.916*1e-18*1e-6)*(earth->asArray[i][k].massDensity*1e-3)*earth->asArray[i][k].volume/(238.051*1.661e-27);
  factorTh=(4.)*(earth->asArray[i][k].abundanceTh*1e-9)*(1.)*(1.563*1e-18*1e-6)*(earth->asArray[i][k].massDensity*1e-3)*earth->asArray[i][k].volume/(232.038*1.661e-27);
  float factor_norm = factorU+factorTh;
  if(factor_norm==0){
    factor_norm++;
  }
  string filenameU = "../Models/AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt";
  string filenameTh = "../Models/AntineutrinoSpectrum_all/AntineutrinoSpectrum_232Th.knt";
  double spectrumU[4500][2];
  read_file_into_2D_array(filenameU, spectrumU);
  double spectrumTh[4500][2];
  read_file_into_2D_array(filenameTh, spectrumTh);
  accel = gsl_interp_accel_alloc(); //Allocate memory for interpolation acceleration.
  spectrum_cspline = gsl_spline_alloc(gsl_interp_cspline, 4500); //Allocate memory for spline object for cubic spline in array of length 4500.
  double e_array[4500], y_arrTh[4500];
  split_array(spectrumTh, e_array, 0);
  double total_spectrum[4500];
  for(int n=0 ;n<4500;n++){
    total_spectrum[n]=(factorU*spectrumU[n][1]+factorTh*spectrumTh[n][1])/(factor_norm);
  }
  gsl_spline_init(spectrum_cspline, e_array, total_spectrum, 4500); //Initialize spline object.
  vector <float> weights;
  weights.reserve(Steps);
  double energy;
  float this_weight;
  float total_weight=0;
  for(int n=0;n<Steps;n++){
    energy=Energy_linspace[n];
    this_weight = gsl_spline_eval(spectrum_cspline, energy, accel);
    weights.push_back(this_weight);
    total_weight+=this_weight;
  }
  float comp=0;
  ofstream outfile;
  outfile.open("prob_weight.dat");
  if(total_weight==0){
    total_weight++;
  }
  for(int n=0;n<Steps;n++){
    weights[n]/=total_weight;
    outfile<<weights[n]<< endl;
    comp+=weights[n];
  }
  outfile.close();

  cout << comp << endl;
/*
0 1 1
0 2 1
0 3 1
0 4 1
0 5 1
0 6 1
0 7 1
0 8 1
0 9 1
0 10 1
0 11 1
0 12 1
0 13 1
0 14 1
0 15 1
0 16 1

*/
  return 0;
}
