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
#include <sys/time.h>

const gsl_rng_type *Gen_Type;
gsl_rng *Gen;
int i;
const int N=5000;
gsl_interp_accel *acc;
gsl_spline *spline;
struct timeval t;



//"../../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt"
void spectrum_array2D(string filename, double to_fill[4500][2]){
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
      x_str>>x_num;
      y_str>>y_num;
      //y_vec.push_back(y_num);
      to_fill[i-12][0]=x_num/1000;
      to_fill[i-12][1]=y_num;
      //x_vec.push_back(x_num);
      //cout << y  << endl;
      x_str.clear();
      y_str.clear();

    }

  }
}
void MH_spectrum_sampling(double spectrum_array[4500][2],double markov_chain[N]){
  double initial_value= 3.5*gsl_rng_uniform_pos(Gen);
  markov_chain[0]=initial_value;
  cout << initial_value << endl;
  for(int i=0;i<N-1;i++){
    double possible_jump=abs(gsl_ran_gaussian(Gen, 0.1)+markov_chain[i]);
    if(possible_jump<spectrum_array[0][0]){possible_jump=spectrum_array[0][0];}
    double criteria = gsl_spline_eval(spline, possible_jump, acc)/gsl_spline_eval(spline, markov_chain[i], acc);
    if(criteria>=1.){
      cout << abs(possible_jump) << endl;
      markov_chain[i+1]=abs(possible_jump);
    }
    else{
      double other_random = gsl_rng_uniform_pos(Gen);
      if(other_random<=criteria){
        cout << possible_jump << endl;
        markov_chain[i+1]=possible_jump;
      }
      else{
        cout << markov_chain[i] << endl;
        markov_chain[i+1]=markov_chain[i];
      }
    }
  }
}
void split_array(double to_split[4500][2], double to_return[4500], int comp){
  for(int k=0;k<4500;k++){
    to_return[k]=to_split[k][comp];
  }
}
int main(int argc, char const *argv[]) {
  double U_238[4500][2];
  spectrum_array2D("../AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt", U_238);
  Gen_Type = gsl_rng_taus;
  Gen = gsl_rng_alloc(Gen_Type);
  gsl_rng_env_setup();
  gettimeofday(&t,NULL);
  int seed = (t.tv_sec * 1000) + (t.tv_usec / 1000);
  gsl_rng_set(Gen, seed);
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, 4500);
  double x_arr[4500], y_arr[4500];
  split_array(U_238, x_arr, 0);
  split_array(U_238, y_arr, 1);
  gsl_spline_init(spline, x_arr, y_arr, 4500);
  //cout << gsl_spline_eval(spline, 3.5, acc) << endl;
  double rand_sampl[N];
  MH_spectrum_sampling(U_238, rand_sampl);

  return 0;
}
