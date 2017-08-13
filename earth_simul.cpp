//g++ -o earth_simul.o earth_simul.cpp `gsl-config --cflags --libs`


#include<iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string>
using namespace std;
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>



const int N = 1000;
const int PREM_len = 187;
const int totalNeutrinos=10000000;
const int path_resolution = 100; //km
int counter =0;
const float R_earth = 6371.;
float dz = 2*R_earth/N;
float dx = dz;
const float PI = 3.1415962589793238;
float abundanceU_BSE[3]={12., 20., 35.};
float abundanceTh_BSE[3]={43., 80., 140.};
//format = {cosm, geoch, geodyn}
const gsl_rng_type *Gen_Type; //Type of random number generator to use.
gsl_rng *Gen; //The actual generator object.
gsl_spline *spectrum_spline;
gsl_interp_accel *acc;


void read_file_into_2D_array(string filename, double to_fill[4500][2]){
  //Fills the array to_fill with 4500 elements for energy and neutrino/energy.

  //Open the file filename
  ifstream spectrum_file(filename.c_str());
  //Strings to save the input as there's header.
  string x, y;
  int it=0; //Count the lines with header.
  //Read file line by line.
  while(spectrum_file >> y >> x){
    it++; //corresponds to iteration number.
    //When header has been read.
    if(it>=12 && it<4512){ //Avoid stack smashing!
      //Convert read strings in stringstreams.
      istringstream x_str(x);
      istringstream y_str(y);
      double x_num, y_num;
      //Cast stringstream into doubles.
      x_str >> x_num;
      y_str >> y_num;
      //Fill array.
      to_fill[it-12][0]=x_num/1e3; //MeV
      to_fill[it-12][1]=y_num*1e3; //1/MeV
      //Clear the stringstreams to avoid problems in next iteration.
      x_str.clear();
      y_str.clear();
    }
  }
}
void split_array(double to_split[4500][2], double to_return[4500], int comp){
  for(int k=0;k<4500;k++){
    to_return[k]=to_split[k][comp];
  }
}
vector<double> mh_sampling(double spectrum_array[4500][2], int len){
  vector<double> markov_chain;
  markov_chain.reserve(len);
  double initial_value = ((4.3-0.0005)*gsl_rng_uniform(Gen)) + 0.0005; //Uniformly random sample in [0.0005, 4.5).
  markov_chain.push_back(initial_value);
  //cout << initial_value << endl;

  for(int m=0;m<len;m++){
    double possible_jump;
    do {
      possible_jump = gsl_ran_gaussian(Gen, 0.1) + markov_chain[m];
    } while(possible_jump <= 0.0005 || possible_jump>=4.5);
    double criteria = gsl_spline_eval(spectrum_spline, possible_jump, acc)/gsl_spline_eval(spectrum_spline, markov_chain[m], acc);
    if(criteria>=1.){
      markov_chain.push_back(possible_jump);
    //  cout << possible_jump << endl;
    }
    else{
      double other_random = gsl_rng_uniform(Gen);
      if(other_random<=criteria){
        markov_chain.push_back(possible_jump);
      //  cout << possible_jump << endl;
      }
      else{
        markov_chain.push_back(markov_chain[m]);
        //cout << markov_chain[i] << endl;
      }
    }
  }
  return markov_chain;
}

vector<float> linspace(float min, float max, int arraylen){
  vector<float> array;
  array.reserve(arraylen);
	float step=(max-min)/arraylen;
	for(int i=0; i<=arraylen; i++){
		array.push_back(min+i*step);
	}
  return array;
}
float the_r(float x, float z, float t, char component){
    float the_r_z=R_earth-z;
    float the_r_x = x;
    float the_r_mag = sqrt(the_r_x*the_r_x + the_r_z*the_r_z);
    if (component=='x'){
      return x-t*300000*the_r_x/the_r_mag;
    }
    else if(component=='z'){
      return z+t*300000*the_r_z/the_r_mag;
    }
}

float get_r(float x, float y){
  return sqrt(x*x + y*y);
}


float density_polynomials(float radius){
    float x = radius/6371.;
    //inner core
    if( radius<= 1221.5){
        return 13.0885-8.8381*x*x;}
    //outer core
    else if (radius<=3480){
        return 12.5815-1.2638*x-3.6426*x*x-5.5281*x*x*x;}
    //lower mantle
    else if (radius <= 5701){
        return 7.9565-6.4761*x+5.5283*x*x-3.0807*x*x*x;}
    //transition zone
    else if (radius <= 5771){
        return 5.3197-1.4836*x;}
    else if (radius <= 5971){
        return 11.2494-8.0298*x;}
    else if (radius <= 6151){
        return 7.1089-3.8045*x;}
    //lvz + lid
    else if (radius <= 6346.6){
        return 2.6910+0.6924*x;}
    //crust
    else if (radius <= 6356){
        return 2.9;}
    else if (radius <= 6368){
        return 2.6;}
    //ocean
    else if (radius <= 6371){
        return 1.020;}
}
vector<float> calculateMantleAbundances(float c_mass, float m_mass, float t_mass, float abundance_isot[3], float crust_abundance_isot){
  vector <float> abundances;
  abundances.reserve(3);
  for(int i =0;i<3;i++){
    abundances.push_back((t_mass*abundance_isot[i] - c_mass*crust_abundance_isot)/m_mass);
  }
  return abundances;
}
vector<float> copy_vector(vector<float> to_copy){
  vector<float> copy;
  copy.reserve(10);
  for(int n=0;n<10;n++){
    copy.push_back(to_copy[n]);
  }
  return copy;
}
/*
vector< vector<float> > import_model(string filename){
  float rad, depth, density, Vpv, Vph, Vsv, Vsh, eta, Q_mu, Q_kappa;
  string line;
  ifstream infile(filename.c_str());
  int i = 0;
  vector< vector<float> > model_matrix;
  model_matrix.reserve(199);
  vector<float> last_row;
  last_row.reserve(10);
  while(getline(infile, line) && i<=199){
    istringstream splittable_line(line);
    string field;
    vector<float> row;
    row.reserve(10);
    while(getline(splittable_line, field, ',')){
      istringstream field_ss(field);
      float field_num;
      field_ss >> field_num;
      row.push_back(field_num);
    }
    if(i==0){
      model_matrix.push_back(row);
      last_row = copy_vector(row);
      counter++;
    }
    else if(row[0]<last_row[0]){
      model_matrix.push_back(row);
      last_row = copy_vector(row);
      counter++;
    }
    splittable_line.clear();
    i++;
  }
  return model_matrix;
}

void split_array(vector< vector<float> > to_split, float container[PREM_len], int comp, bool invert){
  if(invert){
    for(int i=0;i<PREM_len;i++){
      container[PREM_len-1-i]=to_split[i][comp];
    }
  }
  else{
    for(int i=0;i<PREM_len;i++){
      container[i]=to_split[i][comp];
    }  }
}
*/

class RingNode{
  public:
    float x;
    float z;
    float r;
    float getRadius(){
      r = sqrt(x*x + z*z);
      return r;
    }
    bool isEarth;
    bool isSE;
    bool isCrust;
    bool isMantle;
    float mass;
    float massDensity;
    float solidAngle;
    float getSolidAngle(){
      solidAngle = 2*PI*(x*dx/((R_earth-z)*(R_earth-z) + x*x));
      return solidAngle;
    }
    float volume;
    float getVolume(){
      volume = 2*PI*x*dx*dz;
      return volume;
    }
    float abundanceU;
    float abundanceTh;
    float neutrinoFlux;
    float neutrinoThFlux;
    float neutrinoUFlux;
    float relativeNeutrinoTh;
    float relativeNeutrinoU;
    float relativeNeutrino;
    float neutrinosProduced;
    float neutrinosProducedU;
    float neutrinosProducedTh;
    float slope;
    vector<float> path;
    float pathLen;
    float distanceToDetector;
    vector<double> allowedEnergies;


};

class Planet{
  public:
    RingNode asArray[N/2][N];
    float totalMass;
    float crustMass;
    float mantleMass;
    float totalFlux;

    void initializeCoords(){
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          asArray[i][k].x=i*dx;
          asArray[i][k].z=-6371. + k*dz;
          float r = asArray[i][k].getRadius();
          if(r<6371){
            asArray[i][k].isEarth=1;
            float dk=1e3-k;
            float di=-i;
            if(i!=0){asArray[i][k].slope = dk/di;}
            if(r>3480){
              asArray[i][k].isSE=1;
              if(r>R_earth-32.5){
                asArray[i][k].isCrust=1;
                asArray[i][k].isMantle=0;
              }
              else{
                asArray[i][k].isCrust=0;
                asArray[i][k].isMantle=1;
              }
            }
            else{
              asArray[i][k].isSE=0;
            }
          }
          else{asArray[i][k].isEarth=0;}
          float dummy_sa = asArray[i][k].getSolidAngle();
          float dummy_vol = asArray[i][k].getVolume();
        }
      }
    }
    void initializeDensity(){
      /*
      vector< vector<float> > PREM_complete;
      PREM_complete = ::import_model("../Models/PREM_1s.csv");
      float radiusArray[PREM_len], densityArray[PREM_len];
      ::split_array(PREM_complete, radiusArray, 0, 1);
      ::split_array(PREM_complete, densityArray, 2, 1);
      gsl_interp_accel *acc = gsl_interp_accel_alloc();
      gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, PREM_len);
      gsl_spline_init(spline, radiusArray, densityArray, PREM_len);
      */
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          if(asArray[i][k].isEarth){
            //asArray[i][k].massDensity=gsl_spline_eval(spline, asArray[i][k].r, acc);
            asArray[i][k].massDensity=density_polynomials(asArray[i][k].r);
            asArray[i][k].mass=asArray[i][k].massDensity*1e3*asArray[i][k].volume*1e9;
            totalMass+=asArray[i][k].mass;
            if(asArray[i][k].isCrust){
              crustMass+=asArray[i][k].mass;
            }
            else if(asArray[i][k].isMantle){
              mantleMass+=asArray[i][k].mass;
            }
          }
          else{asArray[i][k].massDensity=-1;}
        }
      }
      //gsl_spline_free(spline);
      //gsl_interp_accel_free(acc);
    }
    void initializeAbundanceCrust(){
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          if(asArray[i][k].isCrust){
            asArray[i][k].abundanceU=1.31; //ppb
            asArray[i][k].abundanceTh=5.61; //ppb
          }
          else{
            asArray[i][k].abundanceU=0; //%
            asArray[i][k].abundanceTh=0; //%
          }
        }
      }
    }
    void initializeAbundanceMantle(string key, string bse_model){
      int model;
      if(bse_model=="cosmo"){model=0;}
      else if(bse_model=="geoch"){model=1;}
      else if(bse_model=="geodyn"){model=2;}
      if(key=="unif"){
        for(int i =0 ; i<N/2;i++){
          for(int k = 0;k<N;k++){
            if(asArray[i][k].isMantle){
              asArray[i][k].abundanceU=calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceU_BSE, 1.31)[model]; //ppb
              asArray[i][k].abundanceTh=calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceTh_BSE, 5.61)[model]; //ppb
            }
          }
        }
      }
      else if(key=="two_layer"){
        float bulk_mantle_U = calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceU_BSE, 1.31)[model];
        float bulk_mantle_Th = calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceTh_BSE, 5.61)[model];
        float mantleMass_fraction = 0.1*mantleMass;
        float mass_count=0;
        float limit_rad;
        int n=0;
        do {
          mass_count+=4*PI*(asArray[n][500].r)*(asArray[n][500].r)*dx*1e3*1e9*asArray[n][500].massDensity;
          limit_rad=asArray[n][500].r;
          n++;
        } while(mass_count<=mantleMass_fraction);
        for(int i =0 ; i<N/2;i++){
          for(int k = 0;k<N;k++){
            if(asArray[i][k].isMantle && asArray[i][k].r>limit_rad){
              asArray[i][k].abundanceU=4.7; //ppb
              asArray[i][k].abundanceTh=13.7; //ppb
            }
            else if(asArray[i][k].isMantle && asArray[i][k].r<=limit_rad){
              asArray[i][k].abundanceU=(-4.7*(mantleMass-mass_count)+bulk_mantle_U*mantleMass)/mass_count; //ppb
              asArray[i][k].abundanceTh=(-13.7*(mantleMass-mass_count)+bulk_mantle_Th*mantleMass)/mass_count; //ppb
            }
          }
        }

      }
    }
    void initializeFluxes(){
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          if(asArray[i][k].isSE){
            asArray[i][k].neutrinoUFlux=(asArray[i][k].abundanceU*1e-6)*(asArray[i][k].massDensity*1e3)*(6)*(4.916*1e-18)*(0.9927)/(238.051*1.661e-27)*(asArray[i][k].volume*1e9)*asArray[i][k].solidAngle;
            asArray[i][k].neutrinoThFlux=(asArray[i][k].abundanceTh*1e-6)*(asArray[i][k].massDensity*1e3)*(4)*(1.563*1e-18)*(1)/(232.038*1.661e-27)*(asArray[i][k].volume*1e9)*asArray[i][k].solidAngle;
            asArray[i][k].neutrinoFlux=  asArray[i][k].neutrinoUFlux+asArray[i][k].neutrinoThFlux;
            asArray[i][k].relativeNeutrinoU=asArray[i][k].neutrinoUFlux/asArray[i][k].neutrinoFlux;
            asArray[i][k].relativeNeutrinoTh=asArray[i][k].neutrinoThFlux/asArray[i][k].neutrinoFlux;
            totalFlux+=asArray[i][k].neutrinoFlux;
          }
          else{
            asArray[i][k].neutrinoUFlux=0;
            asArray[i][k].neutrinoThFlux=0;
            asArray[i][k].neutrinoFlux= 0;
            asArray[i][k].relativeNeutrinoU=0;
            asArray[i][k].relativeNeutrinoTh=0;
          }
        }
      }
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          asArray[i][k].relativeNeutrino=asArray[i][k].neutrinoFlux/totalFlux;
          asArray[i][k].neutrinosProduced=  asArray[i][k].relativeNeutrino*totalNeutrinos;
          asArray[i][k].neutrinosProducedU=asArray[i][k].neutrinosProduced*asArray[i][k].relativeNeutrinoU;
          asArray[i][k].neutrinosProducedTh=asArray[i][k].neutrinosProduced*asArray[i][k].relativeNeutrinoTh;
        }
      }
    }
    void initializePaths(){
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          if(asArray[i][k].isSE){
            float z, x;
            z=asArray[i][k].z;
            x=asArray[i][k].x;
            vector<float> path;
            float path_len = get_r(R_earth-z,x );
            asArray[i][k].distanceToDetector=path_len;
            float element_num = roundf(path_len/path_resolution);
            asArray[i][k].pathLen=element_num;
            path.reserve(element_num);
            vector<float> times;
            times=linspace(0, path_len/300000, int(element_num));
            for(int n=0;n<element_num;n++){
              float t = times[n];
              path.push_back(density_polynomials(get_r(the_r(x, z, t, 'x'), the_r(x, z, t, 'y'))));
            }
            asArray[i][k].path=path;
          }
        }
      }
    }
    void initializeEnergySamples(string isotope){
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          string file;
          int num_neut;
          if(isotope == "238U"){
            file = "../Models/AntineutrinoSpectrum_all/AntineutrinoSpectrum_238U.knt";
            num_neut = int(roundf(asArray[i][k].neutrinosProducedU));
          }
          else if(isotope == "232Th"){
            file = "../Models/AntineutrinoSpectrum_all/AntineutrinoSpectrum_232Th.knt";
            num_neut = int(roundf(asArray[i][k].neutrinosProducedTh));
          }

          struct timeval time;
          gettimeofday(&time,NULL);
          double isot_spectrum[4500][2];
          read_file_into_2D_array(file, isot_spectrum);
          unsigned long int seed = (time.tv_sec * 1000) + (time.tv_usec / 1000)+i+k;
          gsl_rng_env_setup(); //Setup environment variables.
          Gen_Type = gsl_rng_taus; //The fastest random number generator.
          Gen = gsl_rng_alloc(Gen_Type); //Allocate necessary memory, initialize generator object.
          gsl_rng_set(Gen, seed); //Seed the generator
          acc = gsl_interp_accel_alloc(); //Allocate memory for interṕolation acceleration.
          spectrum_spline = gsl_spline_alloc(gsl_interp_cspline, 4500); //Allocate memory for spline object for cubic spline in array of length 4500.
          double x_arr[4500], y_arr[4500];
          split_array(isot_spectrum, x_arr, 0);
          split_array(isot_spectrum, y_arr, 1);
          gsl_spline_init(spectrum_spline, x_arr, y_arr, 4500); //Initialize spline object.
          vector<double> rand_sampl;
          rand_sampl = mh_sampling(isot_spectrum, int(asArray[i][k].pathLen));
          asArray[i][k].allowedEnergies=rand_sampl;
          gsl_rng_free(Gen);
        }
        }
      }

    void initialize(string key, string bse_model){
      initializeCoords();
      initializeDensity();
      initializeFluxes();
      initializePaths();
      initializeEnergySamples("232Th");
    }
};


int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  earth->initialize("two_layer", "geodyn");
//  cout << earth->totalMass << endl;

  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      cout << earth->asArray[i][k].distanceToDetector  << ',' ;
      }
      cout << 0 << endl;
    }

  delete earth;
//cout << earth->asArray[200][500].invPath[0] << endl;
  return 0;
}
