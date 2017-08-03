//g++ -o earth_simul.o earth_simul.cpp


#include<iostream>
using namespace std;
#include <cmath>
const int N = 1000;
double R_earth = 6371.;
double dz = 2*R_earth/N;
double dx = dz;
const double PI = 3.1415962589793238;
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
    RingNode asArray[N][N];
    void initializeCoords(){
      for(int i =0 ; i<N;i++){
        for(int k = 0;k<N;k++){
          asArray[i][k].x=-6371. + i*dx;
          asArray[i][k].z=-6371. + k*dz;
          double r = asArray[i][k].getRadius();
          if(r<6371){asArray[i][k].isEarth=1;}
          else{asArray[i][k].isEarth=0;}
        }
      }


    }
};
int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  RingNode node;
  earth->initializeCoords();
  cout << earth->asArray[500][500].isEarth << ','<< earth->asArray[0][999].z<< endl;
  return 0;
}
