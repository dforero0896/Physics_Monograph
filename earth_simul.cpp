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



const int N = 1000;
const int PREM_len = 187;
int counter =0;
double R_earth = 6371.;
double dz = 2*R_earth/N;
double dx = dz;
const double PI = 3.1415962589793238;

vector<double> copy_vector(vector<double> to_copy){
  vector<double> copy;
  copy.reserve(10);
  for(int n=0;n<10;n++){
    copy.push_back(to_copy[n]);
  }
  return copy;
}
vector< vector<double> > import_model(string filename){
  float rad, depth, density, Vpv, Vph, Vsv, Vsh, eta, Q_mu, Q_kappa;
  string line;
  ifstream infile(filename.c_str());
  int i = 0;
  vector< vector<double> > model_matrix;
  model_matrix.reserve(199);
  vector<double> last_row;
  last_row.reserve(10);
  while(getline(infile, line) && i<=199){
    istringstream splittable_line(line);
    string field;
    vector<double> row;
    row.reserve(10);
    while(getline(splittable_line, field, ',')){
      istringstream field_ss(field);
      double field_num;
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
void split_array(vector< vector<double> > to_split, double container[PREM_len], int comp, bool invert){
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


class RingNode{
  public:
    double x;
    double z;
    double r;
    double getRadius(){
      r = sqrt(x*x + z*z);
      return r;
    }
    bool isEarth;
    double massDensity;
    double solidAngle;
    double getSolidAngle(){
      solidAngle = 2*PI*(x*dx/((R_earth-z)*(R_earth-z) + x*x));
      return solidAngle;
    }
    double volume;
    double getVolume(){
      volume = 2*PI*x*dx*dz;
      return volume;
    }
    double abundanceU;
    double abundanceTh;
    double neutrinoFlux;
    double neutrinoThFlux;
    double neutrinoUFlux;
    double relativeNeutrinoTh;
    double relativeNeutrinoU;
    double relativeNeutrino;
};
class Planet{
  public:
    RingNode asArray[N/2][N];
    void initializeCoords(){
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          asArray[i][k].x=i*dx;
          asArray[i][k].z=-6371. + k*dz;
          double r = asArray[i][k].getRadius();
          if(r<6371){asArray[i][k].isEarth=1;}
          else{asArray[i][k].isEarth=0;}
        }
      }
    }
    void initializeDensity(){
      vector< vector<double> > PREM_complete;
      PREM_complete = ::import_model("../Models/PREM_1s.csv");
      double radiusArray[PREM_len], densityArray[PREM_len];
      ::split_array(PREM_complete, radiusArray, 0, 1);
      ::split_array(PREM_complete, densityArray, 2, 1);
      gsl_interp_accel *acc = gsl_interp_accel_alloc();
      gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, PREM_len);
      gsl_spline_init(spline, radiusArray, densityArray, PREM_len);
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          if(asArray[i][k].isEarth){
            asArray[i][k].massDensity=gsl_spline_eval(spline, asArray[i][k].r, acc);
          }
          else{asArray[i][k].massDensity=-1;}
        }
      }
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
    }
    void initialize(){
      initializeCoords();
      initializeDensity();
    }
};
int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  RingNode node;
  earth->initialize();
  for(int i =0 ; i<N/2;i++){
    for(int k = 0;k<N;k++){
      cout << earth->asArray[i][k].massDensity  << ',' ;
      }
      cout << 0 << endl;
    }
  return 0;
}
