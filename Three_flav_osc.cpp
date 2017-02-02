#include <iostream>
#include <time.h>
using namespace std;
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <complex>
#include <cmath>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
/* This program should be compiled as:
g++ -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
*/
//Fermi coupling constant
long double G_f=1.16637*pow(10, -5); //GeV^-2
long double A;
double m_N=939.57 ; //MeV

//Define matrix elements
double Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3;
double E12, E13, E21, E23, E31, E32;
//Define squared mass differences.
double dM32 = 3.2E-3; //eV^2
double dm21 = 0.0; //eV^2

//Define vacuum mixing angles.
int theta1 = 45; //Degrees
int theta2 = 5; //Degrees
int theta3 = 45; //Degrees
double delta = 0;


//Distance
float L=1;

//Calculate energy differences.
double energies ( double massq, double neut_energy);


//Density profile.
double density( double lon);




int main(){
clock_t t1,t2;
t1=clock();
E21=energies( dm21, 1E9);
E32=energies( dM32, 1e9);

E12=-E21;
E23=-E32;
E31=-E12-E23;
E13=-E31;


//Calculate matrix elements.
Ue1 = gsl_sf_cos(theta2)*gsl_sf_cos(theta3);
Ue2 = gsl_sf_sin(theta3)*gsl_sf_cos(theta2);
Ue3 = gsl_sf_sin(theta2);
Umu1=-gsl_sf_sin(theta3)*gsl_sf_cos(theta1)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_cos(theta3);
Umu2=gsl_sf_cos(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_sin(theta3);
Umu3=gsl_sf_sin(theta1)*gsl_sf_cos(theta2);
Ut1=gsl_sf_sin(theta1)*gsl_sf_sin(theta3)-gsl_sf_sin(theta2)*gsl_sf_cos(theta2)*gsl_sf_cos(theta3);
Ut2=-gsl_sf_sin(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta2)*gsl_sf_sin(theta3)*gsl_sf_cos(theta1);
Ut3=gsl_sf_cos(theta1)*gsl_sf_cos(theta2);
/*
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
*/
cout << "Matrix elements calculated" << endl;

//Build CKM matrix.

gsl_matrix *CKM = gsl_matrix_alloc(3, 3);

gsl_matrix_set (CKM, 0, 0, Ue1);
gsl_matrix_set (CKM, 0, 1, Ue2);
gsl_matrix_set (CKM, 0, 2, Ue3);
gsl_matrix_set (CKM, 1, 0, Umu1);
gsl_matrix_set (CKM, 1, 1, Umu2);
gsl_matrix_set (CKM, 1, 2, Umu3);
gsl_matrix_set (CKM, 2, 0, Ut1);
gsl_matrix_set (CKM, 2, 1, Ut2);
gsl_matrix_set (CKM, 2, 2, Ut3); 

cout << "CKM matrix built..." << endl;


//Invert CKM matrix.
gsl_matrix *ICKM = gsl_matrix_alloc(3, 3);
gsl_matrix_transpose_memcpy(ICKM, CKM);
cout << "CKM matrix inverted (transposed)..." << endl;


//Build perturbation in flavor basis
gsl_matrix *V_f = gsl_matrix_alloc(3, 3);
gsl_matrix_set (V_f, 0, 0, 1.0); 
gsl_matrix_set (V_f, 0, 1, 0); 
gsl_matrix_set (V_f, 0, 2, 0); 
gsl_matrix_set (V_f, 1, 0, 0); 
gsl_matrix_set (V_f, 1, 1, 0); 
gsl_matrix_set (V_f, 1, 2, 0); 
gsl_matrix_set (V_f, 2, 0, 0); 
gsl_matrix_set (V_f, 2, 1, 0); 
gsl_matrix_set (V_f, 2, 2, 0); 
//Build perturbation in mass basis
//V_m=U^-1 V_f U

gsl_matrix *V_m = gsl_matrix_alloc(3, 3);
gsl_matrix *first = gsl_matrix_alloc(3, 3);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ICKM, V_f, 0.0, first);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, first, CKM, 0.0, V_m); 

gsl_matrix_free(first);

//Matter density A
A=1.0;//(1/sqrt(2))*G_f*(1/(m_N*1E-3))*density(L);
gsl_matrix_scale(V_m, A);
cout << "The matter density is " << A << endl;

//Build matrix T
gsl_matrix *T = gsl_matrix_alloc(3, 3);
gsl_matrix *traceHm = gsl_matrix_alloc(3, 3);
//Matrix traceHm is defined as H_m - 1/3trH_m.
gsl_matrix_set(traceHm, 0, 0, (1/3)*(E12+E13-A));
gsl_matrix_set(traceHm, 0, 1, 0);
gsl_matrix_set(traceHm, 0, 2, 0);
gsl_matrix_set(traceHm, 1, 0, 0);
gsl_matrix_set(traceHm, 1, 1, (1/3)*(E21+E23-A));
gsl_matrix_set(traceHm, 1, 2, 0);
gsl_matrix_set(traceHm, 2, 0, 0);
gsl_matrix_set(traceHm, 2, 1, 0);
gsl_matrix_set(traceHm, 2, 2, (1/3)*(E31+E32-A));

cout.precision(30);
cout << gsl_matrix_fprintf(stdout, first, "%f");


gsl_matrix_add(V_m, traceHm);
gsl_matrix *otherT = gsl_matrix_alloc(3, 3);
gsl_matrix_memcpy(V_m, T);
gsl_matrix_memcpy(V_m, otherT);

//gsl_matrix_free(V_m);
//gsl_matrix_free(traceHm);


//Calculate the eigenvalues of T
//Calculate c_i
//c_0 = -detT

gsl_permutation *p = gsl_permutation_alloc(3);
int signum;
gsl_linalg_LU_decomp(otherT, p, &signum);

//cout << gsl_matrix_fprintf(stdout, otherT, "%g");

double c_0 = gsl_linalg_LU_det(otherT, signum);
cout << "c_0= " << c_0 << endl;



t2=clock();
float diff ((float)t2-(float)t1);
cout<<"Time elapsed " << diff << endl;


return 0;
}

double density ( double lon){

return 1.0;

}

double energies ( double massq, double neut_energy){

double Eij = massq/(2*neut_energy);
return Eij;
}
