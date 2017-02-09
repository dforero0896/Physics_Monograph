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
long double G_f=1.16637E-5; //GeV^-2
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

//Eigenvalues
gsl_complex lam1, lam2, lam3;

//Matrices.
gsl_matrix CKM;

void toFlavor(const gsl_matrix* toTransform, gsl_matrix *destiny, const gsl_matrix *CKM);
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



//Matter density A
A=1E-13;//(1/sqrt(2))*G_f*(1/(m_N*1E-3))*density(L);

cout << "The matter density is " << A << endl;

//Build matrix T
gsl_matrix *T = gsl_matrix_alloc(3, 3);

gsl_matrix_set(T, 0, 0, A*pow(Ue1,2)+(1./3)*(E12+E13-A));
gsl_matrix_set(T, 0, 1, A*Ue1*Ue2);
gsl_matrix_set(T, 0, 2, A*Ue1*Ue3);
gsl_matrix_set(T, 1, 0, A*Ue1*Ue2);
gsl_matrix_set(T, 1, 1, A*pow(Ue2,2)+(1./3)*(E21+E23-A));
gsl_matrix_set(T, 1, 2, A*Ue2*Ue3);
gsl_matrix_set(T, 2, 0, A*Ue1*Ue3);
gsl_matrix_set(T, 2, 1, A*Ue2*Ue3);
gsl_matrix_set(T, 2, 2, A*pow(Ue3,2)+(1./3)*(E31+E32-A));
cout << "Matrix T calculated..." << endl;

//Calculate c0.
long double c0=gsl_matrix_get(T,0, 0)*gsl_matrix_get(T,1, 1)*gsl_matrix_get(T,2, 2)-gsl_matrix_get(T,0,0)*gsl_matrix_get(T,1,2)*gsl_matrix_get(T,2,1)-gsl_matrix_get(T,0,1)*gsl_matrix_get(T,1,0)*gsl_matrix_get(T,2,2)+gsl_matrix_get(T,0,1)*gsl_matrix_get(T,2,0)*gsl_matrix_get(T,1,2)+gsl_matrix_get(T,0,2)*gsl_matrix_get(T,1,0)*gsl_matrix_get(T,2,1)-gsl_matrix_get(T,0,2)*gsl_matrix_get(T,2,0)*gsl_matrix_get(T,1,1);
c0 *= -1.;



//Calculate c1.
long double c1=gsl_matrix_get(T,0,0)*gsl_matrix_get(T,1,1)+gsl_matrix_get(T,0,0)*gsl_matrix_get(T,2,2)+gsl_matrix_get(T,1,1)*gsl_matrix_get(T,2,2)-gsl_matrix_get(T,1,2)*gsl_matrix_get(T,2,1)-gsl_matrix_get(T,0,1)*gsl_matrix_get(T,1,0)-gsl_matrix_get(T,0,2)*gsl_matrix_get(T,2,0);

cout << "Constants c_i calculated..." << endl;

//Calculate eigenvalues.
gsl_complex atanArg = gsl_complex_rect((1./c0)*sqrt(-pow(c0, 2)-(4./27.)*pow(c1, 3)), 0);
gsl_complex atanVal=gsl_complex_mul_real(gsl_complex_arctan(atanArg), 1./3.);
gsl_complex half = gsl_complex_rect(2*sqrt((-1./3.)*c1), 0);
gsl_complex s1Ps2 = gsl_complex_mul(half, gsl_complex_cos(atanVal));
gsl_complex dummy_s1Ms2 = gsl_complex_mul(half, gsl_complex_sin(atanVal));
gsl_complex s1Ms2 = gsl_complex_mul(dummy_s1Ms2, gsl_complex_rect(0., 1.));

gsl_complex lam1 = gsl_complex_add(gsl_complex_mul_real(s1Ps2, -1./2.), gsl_complex_mul_real(gsl_complex_mul(gsl_complex_rect(0., 1.), s1Ms2), sqrt(3.)/2.));

gsl_complex lam2 = gsl_complex_sub(gsl_complex_mul_real(s1Ps2, -1./2.), gsl_complex_mul_real(gsl_complex_mul(gsl_complex_rect(0., 1.), s1Ms2), sqrt(3.)/2.));

gsl_complex lam3 = s1Ps2;

cout << GSL_REAL(lam2) << endl;
cout << "Eigenvalues calculated..." << endl;
// /*
//Build matrix T**2
gsl_matrix *Tsq = gsl_matrix_alloc(3, 3);
gsl_matrix_set(Tsq, 0, 0,(1./3.)*(pow(A, 2)*(pow(Ue1, 2) + (1./3.)) + 2.*A*(pow(Ue1, 2)-(1./3.))*(E12+E13)+(1./3.)*pow((E12+E13), 2)));
gsl_matrix_set(Tsq, 0, 1, (1./3.)*Ue1*Ue2*A*(A+E13+E23));
gsl_matrix_set(Tsq, 0, 2, (1./3.)*Ue1*Ue3*A*(A+E12+E32));
gsl_matrix_set(Tsq, 1, 0, (1./3.)*Ue1*Ue2*A*(A+E13+E23));
gsl_matrix_set(Tsq, 1, 1, (1./3.)*(pow(A, 2)*(pow(Ue2, 2) + (1./3.)) + 2.*A*(pow(Ue2, 2)-(1./3.))*(E21+E23)+(1./3.)*pow((E21+E23), 2)));
gsl_matrix_set(Tsq, 1, 2, (1./3.)*Ue2*Ue3*A*(A+E21+E31));
gsl_matrix_set(Tsq, 2, 0, (1./3.)*Ue1*Ue3*A*(A+E12+E32));
gsl_matrix_set(Tsq, 2, 1, (1./3.)*Ue2*Ue3*A*(A+E21+E31));
gsl_matrix_set(Tsq, 2, 2, (1./3.)*(pow(A, 2)*(pow(Ue3, 2) + (1./3.)) + 2.*A*(pow(Ue3, 2)-(1./3.))*(E31+E32)+(1./3.)*pow((E31+E32), 2)));



cout << gsl_matrix_get(Tsq, 0, 2) << endl;
cout << "Matrix T**2 calculated..." << endl;

//Transform matrices to flavor basis
gsl_matrix *TFsq = gsl_matrix_alloc(3, 3);
toFlavor(Tsq, TFsq, CKM);
gsl_matrix *TF = gsl_matrix_alloc(3, 3);
toFlavor(T, TF, CKM);
int i, k;
for ( i =0;i<3; i++){
for ( k=0; k<3; k++){
cout.precision(30);
cout << gsl_matrix_get(TFsq, i, k)<< endl;

}
}
cout << "T and T**2 have been converted to flavor basis..." <<endl;
//*/
//Build operator.


/*
gsl_matrix *That = gsl_matrix_alloc(3, 3);
gsl_matrix *first = gsl_matrix_alloc(3, 3);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ICKM, T, 0.0, first);
gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, first, CKM, 0.0, That); 


cout.precision(30);
cout << gsl_matrix_fprintf(stdout, T_hat.matrix, "%f");
*/

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

void toFlavor (const gsl_matrix *toTransform, gsl_matrix *destiny, const gsl_matrix *CKM){
int alpha, beta, a, b;
long double sum;
for (alpha=0;alpha<3;alpha++){
	for (beta=0;beta<3; beta++){
sum=0.;
		for (a=0;a<3;a++){
			for (b=0;b<3;b++){
sum += gsl_matrix_get(CKM, alpha, a)*gsl_matrix_get(CKM, beta, b)*gsl_matrix_get(toTransform, a, b);

			}
		}
gsl_matrix_set(destiny, alpha, beta, sum);
	}
}
}
