#include <iostream>
using namespace std;
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <complex>
#include <cmath>
#include <gsl/gsl_sf_trig.h>
/* This program should be compiled as:
g++ -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
*/
//Fermi coupling constant
long double G_f=1.16637*pow(10, -5); //GeV^-2
long double A;
double m_N=939.57 ; //MeV

//Define matrix elements
gsl_complex Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3;
gsl_complex E1, E2, E3;
//Define squared mass differences.
float dM32 = 3.2 * pow(10, -3); //eV^2
float dm21 = 0; //eV^2

//Define vacuum mixing angles.
int theta1 = 45; //Degrees
int theta2 = 5; //Degrees
int theta3 = 45; //Degrees
double delta = 0;


//Distance
float L=1;
int main(){

double density( double lon);
//Calculate matrix elements.
gsl_complex phaseCKM = gsl_complex_polar(1, delta); 
Ue1=gsl_complex_rect(gsl_sf_cos(theta2)*gsl_sf_cos(theta3), 0);
Ue2=gsl_complex_rect(gsl_sf_sin(theta3)*gsl_sf_cos(theta2), 0);
Ue3=gsl_complex_mul(gsl_complex_rect(gsl_sf_sin(theta2), 0), phaseCKM);
Umu1=gsl_complex_add(gsl_complex_rect(gsl_sf_sin(theta3)*gsl_sf_cos(theta1), 0), gsl_complex_mul(gsl_complex_rect(-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_cos(theta3), 0), phaseCKM));
Umu2=gsl_complex_add(gsl_complex_rect(gsl_sf_cos(theta1)*gsl_sf_cos(theta3),0), gsl_complex_mul(gsl_complex_rect(-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_cos(theta3), 0), phaseCKM)   );
Umu3=gsl_complex_rect(gsl_sf_sin(theta1)*gsl_sf_cos(theta2), 0);
Ut1=gsl_complex_add(gsl_complex_rect(gsl_sf_sin(theta1)*gsl_sf_sin(theta3), 0), gsl_complex_mul(gsl_complex_rect(-gsl_sf_sin(theta2)*gsl_sf_cos(theta1)*gsl_sf_cos(theta3), 0), phaseCKM));
Ut2=gsl_complex_add(gsl_complex_rect(-gsl_sf_sin(theta1)*gsl_sf_cos(theta3), 0), gsl_complex_mul(gsl_complex_rect(-gsl_sf_sin(theta2)*gsl_sf_sin(theta3)*gsl_sf_cos(theta1), 0), phaseCKM));
Ut3=gsl_complex_rect(gsl_sf_cos(theta1)*gsl_sf_cos(theta2), 0);

cout << "Matrix elements calculated, phase = " << delta << endl;

//Build CKM matrix.
gsl_matrix_complex *CKM = gsl_matrix_complex_alloc(3, 3);
gsl_matrix_complex_set (CKM, 0, 0, Ue1);
gsl_matrix_complex_set (CKM, 0, 1, Ue2);
gsl_matrix_complex_set (CKM, 0, 2, Ue3);
gsl_matrix_complex_set (CKM, 1, 0, Umu1);
gsl_matrix_complex_set (CKM, 1, 1, Umu2);
gsl_matrix_complex_set (CKM, 1, 2, Umu3);
gsl_matrix_complex_set (CKM, 2, 0, Ut1);
gsl_matrix_complex_set (CKM, 2, 1, Ut2);
gsl_matrix_complex_set (CKM, 2, 2, Ut3); 
cout << "CKM matrix built..." << endl;


//Invert CKM matrix.
gsl_matrix_complex *ICKM = gsl_matrix_complex_alloc(3, 3);
gsl_matrix_complex_transpose_memcpy(ICKM, CKM);
cout << "CKM matrix inverted (transposed)..." << endl;

//Build Hamiltonian in mass basis
gsl_matrix_complex *H_m = gsl_matrix_complex_alloc(3, 3);
gsl_matrix_complex_set_zero(H_m);
//Matter density A
A=(1/sqrt(2))*G_f*(1/(m_N*pow(10, -3)))*density(L);

cout << "The matter density is " << A << endl;
//Build matrix T
gsl_matrix_complex *T = gsl_matrix_complex_alloc(3, 3);

return 0;
}

double density ( double lon){

return 1;

}
