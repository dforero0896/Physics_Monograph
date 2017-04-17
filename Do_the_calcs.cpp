#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
using namespace std;

/*
Parameters to oscillations program should be Energy Electron density (A) L
*/
const int npoints = 100;
//The density profile

double LengthsArray[npoints];
double DensityArray[npoints];
double Energies[npoints];

/*
The mean
matter density of the mantle and the core were chosen to
be ρ mantle = 4.5 g/cm 3 (A mantle  1.7 · 10 −13 eV, L mantle =
2885 km) and ρ core = 11.5 g/cm 3 (A core  4.4 · 10 −13 eV,
L core = 6972 km), respectively. Parameter values: h = 0, θ 1 =
45 ◦ , θ 2 = 5 ◦ , θ 3 = 45 ◦ , ∆m 2 = 0, and ∆M 2 = 3.2 · 10 −3 eV 2
*/
int main(int argc, char* argv[]){
clock_t t1,t2;
t1=clock();

int i;
double totalLen = 2885. + 6972.; //km
double interval = totalLen/(10.*npoints);
double energia = 1e10 - 1e9;
double energySpacing = energia/npoints;
int k;
for(k=1;k<npoints;k++){
	double actualEnergy=k*energySpacing;
	for(i=0;i<npoints;i++){
		if(LengthsArray[i]<2885.){
		DensityArray[i]=1.7e-13; //eV
		}
		else{
		DensityArray[i]=4.4e-13; //eV
		}	
		std::stringstream stream;
		stream <<"./oscillations " << " " <<  actualEnergy  << " " << DensityArray[i] << " " << LengthsArray[i] ;
		system(stream.str().c_str());

	}


}


t2=clock();
float diff ((float)t2-(float)t1);
cout<<"Time elapsed " << diff/CLOCKS_PER_SEC * 1E3 <<  "ms" << endl;


return 0;
}
