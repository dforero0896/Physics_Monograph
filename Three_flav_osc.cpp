#include <iostream>
#include <time.h>
using namespace std;
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <omp.h>
/* This program should be compiled as:
g++ -o oscillations Three_flav_osc.cpp `gsl-config --cflags --libs`
Or using the makefile attached.
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


//Distance
float L;
double neutrinoEnergy;
//Trace
long double traceHm;
//Calculate energy differences.
double energies ( double massq, double neut_energy);


//densityStep profile.
void densityStep ( double *fill , double *dist, int arraysize);

//Eigenvalues
gsl_complex lam1, lam2, lam3;
long double eigenVals[3];
long double c0, c1;
//Matrices.
gsl_matrix CKM, T, TF, Tsq, TFsq;
gsl_matrix Ident( gsl_matrix *matrix);
gsl_matrix *Itty = gsl_matrix_alloc(3, 3);


//Functions
void toFlavor(const gsl_matrix* toTransform, gsl_matrix *destiny, const gsl_matrix *CKM);
void scaleToOther(gsl_matrix *toScale, double scaleFactor, gsl_matrix *result);
gsl_complex real2comp(double real);
void sumTerms(int index, gsl_matrix *result, gsl_matrix *term3, gsl_matrix *term2, gsl_matrix *term1, gsl_matrix *sum, double eigVal,  gsl_matrix *Tmat, gsl_matrix *TmatSQ);
void addMatrices(gsl_matrix *mat1, gsl_matrix *mat2, gsl_matrix *mat3,gsl_matrix *result);
void addComplexMatrices(gsl_matrix_complex *mat1, gsl_matrix_complex *mat2, gsl_matrix_complex *mat3,gsl_matrix_complex *result);
void copyToComplexMat(gsl_matrix *initial, gsl_matrix_complex *result);
void initializeMatrixZero(gsl_matrix *initial);
void calculateOperator(double nEnergy, double dist, double dens, gsl_matrix_complex *container);
void copyToComplexMatFromComplex(gsl_matrix_complex *initial, gsl_matrix_complex *result);
void linspace(double array[], double max, double min, int arraylen);
gsl_matrix_complex IdentComp(gsl_matrix_complex *matrix);





int main(int argc, char* argv[]){


omp_set_num_threads(4);
int N=100;
double EnergyLins[N];
linspace(EnergyLins, 1e10, 1e9, N);
double spatialArr[10000];
linspace(spatialArr, 2885.+6972., 0, 10000);
double DensityStep[10000];
densityStep(DensityStep, spatialArr, 10000);

double Probabilities[N];

int i, k;

//#pragma omp parallel for private(i, k)
//Iterate over the energies
for(i=0;i<N;i++){
	double energy = EnergyLins[i];
	gsl_matrix_complex *OpProduct = gsl_matrix_complex_alloc(3, 3);
	IdentComp(OpProduct);
//	#pragma omp parallel for private(k)
	for(k=0;k<10000;k++){
		calculateOperator(energy, spatialArr[k], DensityStep[k], OpProduct);
	}
	gsl_matrix_complex *OpProductH = gsl_matrix_complex_alloc(3, 3);
	gsl_matrix_complex *Probs = gsl_matrix_complex_alloc(3, 3);
	copyToComplexMatFromComplex(OpProduct, OpProductH);
	gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, real2comp(1.), OpProduct, OpProductH, real2comp(0.), Probs);
	Probabilities[i]=GSL_REAL(gsl_matrix_complex_get(Probs, 0, 1));
	cout<<energy<<"," << GSL_REAL(gsl_matrix_complex_get(Probs, 0, 1)) << endl;
	gsl_matrix_complex_free(OpProduct);
	gsl_matrix_complex_free(OpProductH);
	gsl_matrix_complex_free(Probs);
	
	
	
}

	










return 0;
}
/*----------------------------------------------------------------------------------------------------------------------------------- 	*/

void densityStep ( double *fill , double *dist, int arraysize){

int i;
for(i=0;i<arraysize;i++){
	if(dist[i] < 2885.){fill[i]=1.7e-13;}
	else {fill[i] = 4.4e-13;}
}

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

gsl_matrix Ident(gsl_matrix *matrix){
int i, k;
	for(i=0;i<3;i++){
		for(k=0;k<3;k++){
			if(i==k){
				gsl_matrix_set(matrix, i, k, 1);
			}
			else{
				gsl_matrix_set(matrix, i, k,0);
			}
		}
	}
	return *matrix;
}

gsl_matrix_complex IdentComp(gsl_matrix_complex *matrix){
int i, k;
	for(i=0;i<3;i++){
		for(k=0;k<3;k++){
			if(i==k){
				gsl_matrix_complex_set(matrix, i, k, real2comp(1.));
			}
			else{
				gsl_matrix_complex_set(matrix, i, k,real2comp(0.));
			}
		}
	}
	return *matrix;
}

void scaleToOther(gsl_matrix *toScale, double scaleFactor, gsl_matrix *result){
int i, k;
for (i=0;i<3;i++){
	for (k=0;k<3;k++){
		long double matrixElem;
		matrixElem=scaleFactor*gsl_matrix_get(toScale, i, k);
		gsl_matrix_set(result, i, k, matrixElem);
		}
	}
}



gsl_complex real2comp(double real){
return gsl_complex_rect(real, 0);
}

void addMatrices(gsl_matrix *mat1, gsl_matrix *mat2, gsl_matrix *mat3,gsl_matrix *result){

int i, k;
for (i=0;i<3;i++){
	for (k=0;k<3;k++){
		double elem1 = gsl_matrix_get(mat1, i, k);
		double elem2 = gsl_matrix_get(mat2, i, k);
		double elem3 = gsl_matrix_get(mat3, i, k);
		gsl_matrix_set(result, i, k, elem1 + elem2 + elem3);
		}
	}
}

void addComplexMatrices(gsl_matrix_complex *mat1, gsl_matrix_complex *mat2, gsl_matrix_complex *mat3,gsl_matrix_complex *result){

int i, k;
for (i=0;i<3;i++){
	for (k=0;k<3;k++){
		gsl_complex elem1 = gsl_matrix_complex_get(mat1, i, k);
		gsl_complex elem2 = gsl_matrix_complex_get(mat2, i, k);
		gsl_complex elem3 = gsl_matrix_complex_get(mat3, i, k);
		gsl_complex sum = gsl_complex_add(gsl_complex_add(elem1, elem2), elem3);
		gsl_matrix_complex_set(result, i, k, sum);
		}
	}
}


void sumTerms(int index, gsl_matrix *result, gsl_matrix *term3, gsl_matrix *term2, gsl_matrix *term1, gsl_matrix *sum, double eigVal,  gsl_matrix *Tmat, gsl_matrix *TmatSQ){

	double extFactor = 1./(3.*eigVal*eigVal + c1);
	scaleToOther(Tmat, eigVal, term2);
	scaleToOther(Itty, eigVal*eigVal + c1, term1);
	addMatrices(term1, term2, TmatSQ, sum);
	scaleToOther(sum, extFactor, result);
	gsl_matrix_free(term3);
	gsl_matrix_free(term2);
	gsl_matrix_free(term1);
	gsl_matrix_free(sum);


}

void copyToComplexMat(gsl_matrix *initial, gsl_matrix_complex *result){
int i, k;
gsl_complex elem;
for (i=0;i<3;i++){
	for (k=0;k<3;k++){
		elem =real2comp(gsl_matrix_get(initial, i, k));
		gsl_matrix_complex_set(result, i, k, elem);
		}
	}

}

void copyToComplexMatFromComplex(gsl_matrix_complex *initial, gsl_matrix_complex *result){
int i, k;
gsl_complex elem;
for (i=0;i<3;i++){
	for (k=0;k<3;k++){
		elem =gsl_matrix_complex_get(initial, i, k);
		gsl_matrix_complex_set(result, i, k, elem);
		}
	}

}

void initializeMatrixZero(gsl_matrix *initial){
int i, k;
	for(i=0;i<3;i++){
		for(k=0;k<3;k++){
			gsl_matrix_set(initial, i, k, 0.);
		}
	}
}


void calculateOperator(double nEnergy, double dist, double dens, gsl_matrix_complex *container){


neutrinoEnergy=nEnergy;
L=dist;


//Energy differences.
E21=energies(dm21, neutrinoEnergy);
E32=energies(dM32, neutrinoEnergy);

E12=-E21;
E23=-E32;
E31=-E12-E23;
E13=-E31;

traceHm=dM32+0.5*(dm21/E21 + E21)+dens;
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

//cout << "Matrix elements calculated" << endl;

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

//cout << "CKM matrix built..." << endl;


//Invert CKM matrix.
gsl_matrix *ICKM = gsl_matrix_alloc(3, 3);
gsl_matrix_transpose_memcpy(ICKM, CKM);
//cout << "CKM matrix inverted (transposed)..." << endl;



//Matter densityStep A
A=dens;//(1/sqrt(2))*G_f*(1/(m_N*1E-3))*densityStep(L);

//cout << "The matter densityStep is " << A << endl;

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
//cout << "Matrix T calculated..." << endl;

//Calculate c0.
c0=gsl_matrix_get(T,0, 0)*gsl_matrix_get(T,1, 1)*gsl_matrix_get(T,2, 2)-gsl_matrix_get(T,0,0)*gsl_matrix_get(T,1,2)*gsl_matrix_get(T,2,1)-gsl_matrix_get(T,0,1)*gsl_matrix_get(T,1,0)*gsl_matrix_get(T,2,2)+gsl_matrix_get(T,0,1)*gsl_matrix_get(T,2,0)*gsl_matrix_get(T,1,2)+gsl_matrix_get(T,0,2)*gsl_matrix_get(T,1,0)*gsl_matrix_get(T,2,1)-gsl_matrix_get(T,0,2)*gsl_matrix_get(T,2,0)*gsl_matrix_get(T,1,1);
c0 *= -1.;



//Calculate c1.
c1=gsl_matrix_get(T,0,0)*gsl_matrix_get(T,1,1)+gsl_matrix_get(T,0,0)*gsl_matrix_get(T,2,2)+gsl_matrix_get(T,1,1)*gsl_matrix_get(T,2,2)-gsl_matrix_get(T,1,2)*gsl_matrix_get(T,2,1)-gsl_matrix_get(T,0,1)*gsl_matrix_get(T,1,0)-gsl_matrix_get(T,0,2)*gsl_matrix_get(T,2,0);

//cout << "Constants c_i calculated..." << endl;

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

//cout << GSL_REAL(lam2) << endl;
//cout << "Eigenvalues calculated..." << endl;
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



//cout << gsl_matrix_get(Tsq, 0, 2) << endl;
//cout << "Matrix T**2 calculated..." << endl;

//Transform matrices to flavor basis
gsl_matrix *TFsq = gsl_matrix_alloc(3, 3);
toFlavor(Tsq, TFsq, CKM);
gsl_matrix *TF = gsl_matrix_alloc(3, 3);
toFlavor(T, TF, CKM);

//cout << "T and T**2 have been converted to flavor basis..." <<endl;



//Build operator.

gsl_complex phi = gsl_complex_polar(1., -L*traceHm/3);
Ident(Itty);
eigenVals[0]=GSL_REAL(lam1);
eigenVals[1]=GSL_REAL(lam2);
eigenVals[2]=GSL_REAL(lam3);


gsl_matrix *term1_= gsl_matrix_alloc(3,3);
gsl_matrix *term2_= gsl_matrix_alloc(3,3);
gsl_matrix *term3_= gsl_matrix_alloc(3,3);
gsl_matrix *sum_= gsl_matrix_alloc(3,3);
gsl_matrix *projOp1= gsl_matrix_alloc(3,3);
gsl_matrix_complex *Eig1Term = gsl_matrix_complex_alloc(3, 3);
sumTerms(0, projOp1, term3_, term2_, term1_, sum_, eigenVals[0], TF, TFsq); //This is the P1
copyToComplexMat(projOp1, Eig1Term); 
gsl_matrix_complex_scale(Eig1Term, gsl_complex_polar(1, -L*eigenVals[0]));

gsl_matrix *term1_2= gsl_matrix_alloc(3,3);
gsl_matrix *term2_2= gsl_matrix_alloc(3,3);
gsl_matrix *term3_2= gsl_matrix_alloc(3,3);
gsl_matrix *sum_2= gsl_matrix_alloc(3,3);
gsl_matrix *projOp2= gsl_matrix_alloc(3,3);
gsl_matrix_complex *Eig2Term = gsl_matrix_complex_alloc(3, 3);
sumTerms(0, projOp2, term3_2, term2_2, term1_2, sum_2, eigenVals[1], TF, TFsq);  //This is the P2
copyToComplexMat(projOp2, Eig2Term);
gsl_matrix_complex_scale(Eig2Term, gsl_complex_polar(1, -L*eigenVals[1]));

gsl_matrix *term1_3= gsl_matrix_alloc(3,3);
gsl_matrix *term2_3= gsl_matrix_alloc(3,3);
gsl_matrix *term3_3= gsl_matrix_alloc(3,3);
gsl_matrix *sum_3= gsl_matrix_alloc(3,3);
gsl_matrix *projOp3= gsl_matrix_alloc(3,3);
gsl_matrix_complex *Eig3Term = gsl_matrix_complex_alloc(3, 3);
sumTerms(0, projOp3, term3_3, term2_3, term1_3, sum_3, eigenVals[2], TF, TFsq);  //This is the P3
copyToComplexMat(projOp3, Eig3Term);
gsl_matrix_complex_scale(Eig3Term, gsl_complex_polar(1, -L*eigenVals[2]));


gsl_matrix_complex *FlavorOp = gsl_matrix_complex_alloc(3, 3);

addComplexMatrices(Eig1Term, Eig2Term, Eig3Term, FlavorOp); //FlavorOp is the operator in flavor basis.

gsl_matrix *ProjectionOperators[3];

ProjectionOperators[0]=projOp1;
ProjectionOperators[1]=projOp2;
ProjectionOperators[2]=projOp3;



//cout << "Operator in flavor basis ready!" << endl;

//Build probability amplitude.

int a, b;
gsl_matrix *ProbAmps = gsl_matrix_alloc(3, 3);
initializeMatrixZero(ProbAmps);
gsl_matrix *thisIter = gsl_matrix_alloc(3, 3);
for(a=0;a<3;a++){
	for(b=a;b<3;b++){
		initializeMatrixZero(thisIter);
		double sinX = gsl_sf_sin((eigenVals[a]-eigenVals[b])*L/2);
		double sinSqX=sinX*sinX;
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., ProjectionOperators[a], ProjectionOperators[b], 0., thisIter);
		gsl_matrix_scale(thisIter, -4*sinSqX);
		gsl_matrix_add(ProbAmps, thisIter);
}
}

gsl_matrix_add(ProbAmps, Itty);

gsl_matrix_complex *copyPrevOp = gsl_matrix_complex_alloc(3, 3);
copyToComplexMatFromComplex(container, copyPrevOp);
gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, real2comp(1.), FlavorOp, copyPrevOp, real2comp(0.), container);

gsl_matrix_complex_free(copyPrevOp);
}


void linspace(double array[], double max, double min, int arraylen){
int i;
double step=(max-min)/arraylen;
	for(i=0; i<=arraylen; i++){
		array[i]=min+i*step;
		}
} 
