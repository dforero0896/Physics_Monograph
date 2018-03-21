/*
g++ -fopenmp -o uandino.o earth_simul.cpp uandino.cpp `gsl-config --cflags --libs`
*/
#include "uandino.h"
double dM32 = 2.53685e-3; //ev^2
double dm21 = 7.37e-5; //ev^2
//CKM Elements
double Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3;
gsl_matrix *CKM;
float path_resol=path_resolution;//km
const double PI= 3.1415926589793238;

double longitude_units_conversion(double lon_in_km){
	return lon_in_km*1e3/(1.972e-7);
}
double deg2rad(double deg){
	return deg*PI/180.;
}

void fill_real_matrix(gsl_matrix *empty, double elem_11, double elem_12, double elem_13, double elem_21, double elem_22, double elem_23, double elem_31, double elem_32, double elem_33){
  gsl_matrix_set(empty, 0, 0, elem_11);
  gsl_matrix_set(empty, 0, 1, elem_12);
  gsl_matrix_set(empty, 0, 2, elem_13);
  gsl_matrix_set(empty, 1, 0, elem_21);
  gsl_matrix_set(empty, 1, 1, elem_22);
  gsl_matrix_set(empty, 1, 2, elem_23);
  gsl_matrix_set(empty, 2, 0, elem_31);
  gsl_matrix_set(empty, 2, 1, elem_32);
  gsl_matrix_set(empty, 2, 2, elem_33);
}
void print_real_matrix(gsl_matrix *to_print){
  for(int i=0;i<3;i++){
    cout << gsl_matrix_get(to_print, i,0) << "," << gsl_matrix_get(to_print, i, 1) << "," << gsl_matrix_get(to_print, i, 2) << endl;

  }
}
string print_complex_number(gsl_complex to_print){
  stringstream real;
  stringstream imag;
  real << GSL_REAL(to_print);
  imag << GSL_IMAG(to_print);
  string result= real.str() + "+ j" + imag.str();
  //cout << result << endl;
  return result;

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
gsl_matrix generate_real_identity(gsl_matrix *matrix){
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
gsl_matrix scale_real_matrix(gsl_matrix *to_scale, double factor){
  gsl_matrix *result =gsl_matrix_alloc(3,3);
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      double element=gsl_matrix_get(to_scale, i, k);
      element*=factor;
      gsl_matrix_set(result, i,k,element);
    }
  }
  return *result;
}
gsl_matrix add_real_matrices(gsl_matrix *term_1, gsl_matrix *term_2){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      double element_1=gsl_matrix_get(term_1, i, k);
      double element_2=gsl_matrix_get(term_2, i, k);
      element_1+=element_2;
      gsl_matrix_set(term_1, i,k,element_1);
    }
  }
  return *term_1;
}
void copy_to_complex_from_real(gsl_matrix *real, gsl_matrix_complex *container){
  /*Copies matrix real matrix "real" into complex matric "container".*/
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_matrix_complex_set(container, i,k, gsl_complex_rect(gsl_matrix_get(real, i, k), 0));
    }
  }
}
void copy_to_real_from_real(gsl_matrix *real, gsl_matrix *container){
  /*Copies matrix real matrix "real" into complex matric "container".*/
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_matrix_set(container, i,k, gsl_matrix_get(real, i, k));
    }
  }
}
gsl_matrix_complex scale_complex_matrix(gsl_matrix_complex *to_scale, gsl_complex complex_factor, double real_factor){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_complex element=gsl_matrix_complex_get(to_scale, i, k);
      element = gsl_complex_mul(element, complex_factor);
      element = gsl_complex_mul_real(element, real_factor);
      gsl_matrix_complex_set(to_scale, i,k,element);
    }
  }
  return *to_scale;}
void print_complex_matrix(gsl_matrix_complex *to_print){
  for(int i=0;i<3;i++){
    cout << print_complex_number(gsl_matrix_complex_get(to_print, i,0)) << "," << print_complex_number(gsl_matrix_complex_get(to_print, i,1)) << "," << print_complex_number(gsl_matrix_complex_get(to_print, i,2)) << endl;

  }
}
void fill_complex_matrix(gsl_matrix_complex *empty, gsl_complex elem_11, gsl_complex elem_12, gsl_complex elem_13, gsl_complex elem_21, gsl_complex elem_22, gsl_complex elem_23, gsl_complex elem_31, gsl_complex elem_32, gsl_complex elem_33){
  gsl_matrix_complex_set(empty, 0, 0, elem_11);
  gsl_matrix_complex_set(empty, 0, 1, elem_12);
  gsl_matrix_complex_set(empty, 0, 2, elem_13);
  gsl_matrix_complex_set(empty, 1, 0, elem_21);
  gsl_matrix_complex_set(empty, 1, 1, elem_22);
  gsl_matrix_complex_set(empty, 1, 2, elem_23);
  gsl_matrix_complex_set(empty, 2, 0, elem_31);
  gsl_matrix_complex_set(empty, 2, 1, elem_32);
  gsl_matrix_complex_set(empty, 2, 2, elem_33);
}
gsl_matrix_complex copy_to_complex_from_complex(gsl_matrix_complex *complex, gsl_matrix_complex *container){
  for(int i=0;i<3;i++){
    for(int k=0; k<3; k++){
      gsl_matrix_complex_set(container, i,k, gsl_matrix_complex_get(complex, i, k));
    }
  }
}
void calculateOperator(double neutrinoEnergy, double A, double L, gsl_matrix_complex* matrix){
	/*Calculates the time-evolution operator for a(n) (anti)neutrino with energy "neutrinoEnergy" passing through a potential "A" and a distance "L".*/
  double E21=dm21/(2*neutrinoEnergy);
  long double E32=dM32/(2*neutrinoEnergy);
  //long double E32 = -(neutrinoEnergy + E21) + GSL_REAL(gsl_complex_sqrt_real((neutrinoEnergy+E21)*(neutrinoEnergy+E21) + dM32));
  //cout << E32 << "_______" << E32p << endl;
  double E12=-E21;
  double E23=-E32;
  double E31=E12-E23;
  double E13=-E31;
  //Elements of the Tmatrix in mass basis
  double T_11=A*Ue1*Ue1-(1./3)*A+(1./3)*(E12+E13);
  double T_12=A*Ue1*Ue2;
  double T_13=A*Ue1*Ue3;
  double T_21=T_12;
  double T_22=A*Ue2*Ue2-(1./3)*A+(1./3)*(E21+E23);
  double T_23=A*Ue2*Ue3;
  double T_31=T_13;
  double T_32=T_23;
  double T_33=A*Ue3*Ue3-(1./3)*A+(1./3)*(E31+E32);
  gsl_matrix *T_mass_mat = gsl_matrix_alloc(3,3); //Allocate memory for 3x3 matrix.
  fill_real_matrix(T_mass_mat, T_11, T_12, T_13, T_21, T_22, T_23, T_31, T_32, T_33);
  //cout << "T matrix is:"<< endl;
  //print_real_matrix(T_mass_mat);


  //Elements of the T**2 matrix in mass basis
  double T_sq_11=(1./3)*(A*A*(Ue1*Ue1+(1./3))+2*A*(Ue1*Ue1-(1./3))*(E12+E13)+(1./3)*(E12+E13)*(E12+E13));
  double T_sq_12=(1./3)*Ue1*Ue2*A*(A+E13+E23);
  double T_sq_13=(1./3)*Ue1*Ue3*A*(A+E12+E32);
  double T_sq_21=T_sq_12;
  double T_sq_22=(1./3)*(A*A*(Ue2*Ue2+(1./3))+2*A*(Ue2*Ue2-(1./3))*(E21+E23)+(1./3)*(E21+E23)*(E21+E23));
  double T_sq_23=(1./3)*Ue2*Ue3*A*(A+E21+E31);
  double T_sq_31=T_sq_13;
  double T_sq_32=T_sq_23;
  double T_sq_33=(1./3)*(A*A*(Ue3*Ue3+(1./3))+2*A*(Ue3*Ue3-(1./3))*(E31+E32)+(1./3)*(E31+E32)*(E31+E32));
  gsl_matrix *T_sq_mass_mat = gsl_matrix_alloc(3, 3);
  fill_real_matrix(T_sq_mass_mat, T_sq_11, T_sq_12, T_sq_13, T_sq_21, T_sq_22, T_sq_23, T_sq_31, T_sq_32, T_sq_33);
  //cout << "T**2 matrix is:"<< endl;
  //print_real_matrix(T_sq_mass_mat);


  //T matrix in flavor basis
  gsl_matrix *T_flav_mat = gsl_matrix_alloc(3, 3);
  toFlavor(T_mass_mat, T_flav_mat, CKM);
  //cout << "T matrix in flavor basis is:"<< endl;
  //print_real_matrix(T_flav_mat);


  //T**2 matrix in flavor basis
  gsl_matrix *T_sq_flav_mat = gsl_matrix_alloc(3, 3);
  toFlavor(T_sq_mass_mat, T_sq_flav_mat, CKM);
  //cout << "T**2 matrix in flavor basis is:"<< endl;
  //print_real_matrix(T_sq_flav_mat);


  //Get rid of T and T**2 in mass basis as they are no longer useful
  gsl_matrix_free(T_mass_mat);//Free allocated memory.
  gsl_matrix_free(T_sq_mass_mat);


  //Calculate c's
  //double c1=T_11*T_22-T_12*T_21+T_11*T_33-T_13*T_31+T_22*T_33-T_23*T_32;
  double c1 = -A*A/3 + (A/(6*neutrinoEnergy))*(Ue1*Ue1*(dM32+2*dm21)+Ue2*Ue2*(dM32-dm21)-Ue3*Ue3*(2*dM32+dm21))-(1./(12*neutrinoEnergy*neutrinoEnergy))*(dM32*dM32+dm21*dm21+dM32*dm21);

  //double c0=-(T_11*T_22*T_33-T_11*T_23*T_32-T_12*T_21*T_33+T_12*T_31*T_23+T_13*T_21*T_32-T_13*T_31*T_22);
  double c0 = (-2./27)*A*A*A+(A*A/(18*neutrinoEnergy))*(Ue1*Ue1*(dM32+2*dm21) + Ue2*Ue2*(dM32-dm21) - Ue3*Ue3*(2*dM32+dm21))+(A/(36*neutrinoEnergy*neutrinoEnergy))*(Ue1*Ue1*(2*dM32+dm21)*(dM32-dm21)+Ue2*Ue2*(2*dM32+dm21)*(dM32+2*dm21)-Ue3*Ue3*(dM32+2*dm21)*(dM32-dm21)-(dM32*dM32+dm21*dm21+dM32*dm21))-(1./(216*neutrinoEnergy*neutrinoEnergy*neutrinoEnergy))*(2*dM32+dm21)*(dM32+2*dm21)*(dM32-dm21);

  //Calculate eigenvalues
  long double q=c1/3;
  long double r=-0.5*c0;
  gsl_complex atanArg = gsl_complex_rect((1./c0)*sqrt(-c0*c0-(4./27.)*c1*c1*c1), 0);
  gsl_complex atanVal=gsl_complex_mul_real(gsl_complex_arctan(atanArg), 1./3.);
  gsl_complex half = gsl_complex_rect(2*sqrt((-1./3.)*c1), 0);
  //s1 + s2
  gsl_complex s1Ps2 = gsl_complex_mul(half, gsl_complex_cos(atanVal));
  gsl_complex dummy_s1Ms2 = gsl_complex_mul(half, gsl_complex_sin(atanVal));
  //s1 - s2
  gsl_complex s1Ms2 = gsl_complex_mul(dummy_s1Ms2, gsl_complex_rect(0., 1.));

  gsl_complex lam1 = gsl_complex_add(gsl_complex_mul_real(s1Ps2, -1./2.), gsl_complex_mul_real(gsl_complex_mul(gsl_complex_rect(0., 1.), s1Ms2), sqrt(3.)/2.));
  //double lam1 = -sqrt(-(1./3)*c1)*gsl_sf_cos(GSL_REAL(atanVal))+sqrt(-c1)*gsl_sf_sin(GSL_REAL(atanVal));

  gsl_complex lam2 = gsl_complex_sub(gsl_complex_mul_real(s1Ps2, -1./2.), gsl_complex_mul_real(gsl_complex_mul(gsl_complex_rect(0., 1.), s1Ms2), sqrt(3.)/2.));
  //double lam2 = -sqrt(-(1./3)*c1)*gsl_sf_cos(GSL_REAL(atanVal))-sqrt(-c1)*gsl_sf_sin(GSL_REAL(atanVal));
  gsl_complex lam3 = s1Ps2;
  //double lam3 = 2*sqrt((-1./3)*c1)*gsl_sf_cos(GSL_REAL(atanVal));

  //cout << "Eigenvals" << endl;
  //print_complex_number(lam1);
  //print_complex_number(lam2);
  //print_complex_number(lam3);
  //print_complex_number(gsl_complex_mul_real(gsl_complex_sub(lam3,lam2), 2*neutrinoEnergy));


  //Calculate Operator Eq 46 Pre print
  double trace_hamiltonian=0.5*E21+E32+3*neutrinoEnergy+A;
  gsl_complex phi_phase = gsl_complex_polar(1., -L*trace_hamiltonian/3);//gsl_complex_polar(1., -L*trace_hamiltonian/3);

  //print_complex_number(phi_phase);


  gsl_complex eigenvalues[3]={lam1, lam2, lam3};
  gsl_matrix_complex *evol_operator = gsl_matrix_complex_alloc(3, 3);
  gsl_complex complex_zero=gsl_complex_rect(0,0);
  fill_complex_matrix(evol_operator, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero, complex_zero);


  for(int n=0;n<3;n++){
    gsl_matrix *I_term = gsl_matrix_alloc(3, 3);
    generate_real_identity(I_term);
    //cout << "Created Identity"<< endl;
    //print_real_matrix(I_term);
    gsl_matrix_scale(I_term, GSL_REAL(eigenvalues[n])*GSL_REAL(eigenvalues[n])+c1);
    //cout << "Scaled Identity"<< endl;
    //print_real_matrix(I_term);
    gsl_matrix *T_term = gsl_matrix_alloc(3, 3);
    copy_to_real_from_real(T_flav_mat,T_term);
    gsl_matrix_scale(T_term, GSL_REAL(eigenvalues[n]));
    gsl_matrix_add(I_term, T_term);
    //Add function stores result in first matrix given.
    gsl_matrix_add(I_term, T_sq_flav_mat);
    //I_term stores the addition of all 3 terms.

    gsl_matrix_scale(I_term, 1./(3*GSL_REAL(eigenvalues[n])*GSL_REAL(eigenvalues[n])+c1));
    gsl_matrix_complex *sum_term = gsl_matrix_complex_alloc(3, 3);
    copy_to_complex_from_real(I_term, sum_term);
    gsl_matrix_complex_scale(sum_term, gsl_complex_polar(1,-L*GSL_REAL(eigenvalues[n]) ));
    gsl_matrix_complex_scale(sum_term, phi_phase);

    gsl_matrix_complex_add(evol_operator, sum_term);
    //print_complex_matrix(evol_operator);

    gsl_matrix_free(I_term);
    gsl_matrix_free(T_term);
    gsl_matrix_complex_free(sum_term);
  }
  copy_to_complex_from_complex(evol_operator, matrix);
  gsl_matrix_complex_free(evol_operator);
  gsl_matrix_free(T_flav_mat);
  gsl_matrix_free(T_sq_flav_mat);
}
void calculateProbabilities(vector<float> path, int N, int Steps, float E_min, float E_max){
  /*Writes a file with energies and all three probabilities*/

  double theta1 = 0.7222;
  double theta2 = 0.1468;
  double theta3 = 0.5764;
  //*/
  Ue1 = gsl_sf_cos(theta2)*gsl_sf_cos(theta3);
  Ue2 = gsl_sf_sin(theta3)*gsl_sf_cos(theta2);
  Ue3 = gsl_sf_sin(theta2);
  Umu1=-gsl_sf_sin(theta3)*gsl_sf_cos(theta1)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_cos(theta3);
  Umu2=gsl_sf_cos(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta1)*gsl_sf_sin(theta2)*gsl_sf_sin(theta3);
  Umu3=gsl_sf_sin(theta1)*gsl_sf_cos(theta2);
  Ut1=gsl_sf_sin(theta1)*gsl_sf_sin(theta3)-gsl_sf_sin(theta2)*gsl_sf_cos(theta1)*gsl_sf_cos(theta3);
  Ut2=-gsl_sf_sin(theta1)*gsl_sf_cos(theta3)-gsl_sf_sin(theta2)*gsl_sf_sin(theta3)*gsl_sf_cos(theta1);
  Ut3=gsl_sf_cos(theta1)*gsl_sf_cos(theta2);


  CKM=gsl_matrix_alloc(3, 3);
  fill_real_matrix(CKM, Ue1, Ue2, Ue3, Umu1, Umu2, Umu3, Ut1, Ut2, Ut3);
  //Define spatial limits for the Earth in km.
  vector<double> EnergyLins = linspace(E_min*1e6, E_max*1e6, N);
	//omp_set_num_threads(4);//Number of threads to use.
	int i,k;

	long double Probabilities[N][3];//Array to save probabilities.
	//double Probabilities[N];
	#pragma omp parallel for private(i)
	for(i=0;i<N;i++){//For each energy...
	cout << "using " << omp_get_num_threads() << " cores" << endl;
	  long double energy=EnergyLins[i];
	  gsl_matrix *Id =gsl_matrix_alloc(3, 3);
    //A matrix to save the product of operators.
	  gsl_matrix_complex *operator_product = gsl_matrix_complex_alloc(3, 3);
	  generate_real_identity(Id);
	  copy_to_complex_from_real(Id, operator_product);
    gsl_matrix_free(Id);
		//#pragma omp parallel for private(k)
    for(k=0;k<Steps;k++){
	    double density=path[k]; //eV
      //double density = density_to_potential(sun_rho(coord),0);
      //Increase coordinate value.
      //coord += step_len;
      //A matrix to store the operator for this step.
	    gsl_matrix_complex *iter_operator = gsl_matrix_complex_alloc(3,3);
	    calculateOperator(energy, density, longitude_units_conversion(path_resolution),iter_operator);
	    gsl_matrix_complex *operator_product_copy = gsl_matrix_complex_alloc(3,3);
      //Copy operator product so far.
	    copy_to_complex_from_complex(operator_product, operator_product_copy);
      //Multiply the operator for this step and the copy. Store them in the matrix containing the whole product.
	    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1., 0), iter_operator, operator_product_copy, gsl_complex_rect(0., 0.),operator_product);

      //Free memory.
      gsl_matrix_complex_free(operator_product_copy);
	    gsl_matrix_complex_free(iter_operator);



	  }

		//Probabilities[i] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,1));
		Probabilities[i][0] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 0,0));
		Probabilities[i][1] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 1,0));
		Probabilities[i][2] = gsl_complex_abs2(gsl_matrix_complex_get(operator_product, 2,0));
		//cout << EnergyLins[i] << "," << Probabilities[i][0] << "," << Probabilities[i][1] << "," << Probabilities[i][2] << endl;
    gsl_matrix_complex_free(operator_product);

	}
  //Write file with probabilities
  ofstream myfile;
  myfile.open ("raw_probs.csv");

	for(i=0;i<N;i++){
		myfile << EnergyLins[i] << "," << Probabilities[i][0] << "," << Probabilities[i][1] << "," << Probabilities[i][2] << endl;
		//cout << EnergyLins[i] << "," << Probabilities[i] << endl;
	}
 //Write file with potential values.
  myfile.close();
  ofstream potentialfile;
  potentialfile.open("potentialTest.csv");
  for(k=0;k<Steps;k+=1000){
    //coord = coord_init + k*step_len;
    potentialfile << path[k] << endl;
  }
  potentialfile.close();
}
