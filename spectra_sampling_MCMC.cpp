/*
g++ -o spectra_sampling_MCMC.o spectra_sampling_MCMC.cpp `gsl-config --cflags --libs`
*/

#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include<cmath>

const gsl_rng_type *Gen_Type;
gsl_rng *Gen;
int i, N;
gsl_interp_accel *acc;
gsl_spline *spline;




//"../../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt"
void spectrum_array2D(string filename, float to_fill[4500][2]){
  ifstream infile(filename.c_str());

  string x, y;
  int i=0;
  //vector<float> x_vec;
  //vector<float> y_vec;
  while(infile >> y >> x){
    i++; //i corresponds to iteration number
    if(i>=12){
      istringstream x_str(x);
      istringstream y_str(y);
      double x_num, y_num;
      x_str>>x_num ;
      y_str>>y_num;
      //y_vec.push_back(y_num);
      to_fill[i-12][0]=x_num;
      to_fill[i-12][1]=y_num;
      //x_vec.push_back(x_num);
      //cout << y  << endl;
      x_str.clear();
      y_str.clear();
      //cout << y_num << endl;

    }

  }
}
vector<double> MH_spectrum_sampling(float spectrum_array[4500][2], int num_to_gen){
  vector <double> markov_chain;
  double initial_value= 5*gsl_rng_uniform_pos(Gen);
  markov_chain.push_back(initial_value);
  for(int i=0;i<num_to_gen;i++){
    double possible_jump=gsl_ran_gaussian(Gen, 0.1)+markov_chain[i];
    double criteria = gsl_spline_eval(spline, possible_jump, acc)/gsl_spline_eval(spline, markov_chain[i], acc);
    if(criteria>=1.){
      markov_chain.push_back(abs(possible_jump));
    }
    else{
      double other_random = gsl_rng_uniform_pos(Gen);
      if(other_random<=criteria){
        markov_chain.push_back(possible_jump);
      }
      else{
        markov_chain.push_back(markov_chain[i]);
      }
    }
  }
  return markov_chain;
}
double *split_array(float to_split[4500][2], int comp){
  double to_return[4500];
  for(int k=0;k<4500;k++){
    to_return[k]=to_split[k][i];
  }
  return to_return;
}
int main(int argc, char const *argv[]) {
  float U_238[4500][2];
  spectrum_array2D("../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt", U_238);
  gsl_rng_env_setup();
  Gen_Type = gsl_rng_taus;
  Gen = gsl_rng_alloc(Gen_Type);
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, 4500);
  double *x_arr = split_array(U_238, 0);
  double *y_arr = split_array(U_238, 1);
  cout << x_arr[30] << endl;
  gsl_spline_init(spline, x_arr, y_arr, 4500);
  vector <double> rand_sampl = MH_spectrum_sampling(U_238, 100);

  return 0;
}
