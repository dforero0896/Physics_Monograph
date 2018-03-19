  /*
g++ -fopenmp -o earth_simul.o earth_simul.cpp `gsl-config --cflags --libs`

 */
#include "earth_simul.h"

const int N_x=500;
const int N_z=1000;
 //km
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

//RingNode asArray[N_x][N_z];
void import_probability(string filename, float prob_matrix[500][1000]){
  string line;
  ifstream infile(filename.c_str());
  int i = 0; //iteration/line
  string i_str;
  string k_str;
  string prob_str;
  while(infile >> i_str >> k_str >> prob_str){
    istringstream i_ss(i_str);
    istringstream k_ss(k_str);
    istringstream prob_ss(prob_str);
    int i_int, k_int;
    float prob_num;
    i_ss >> i_int;
    k_ss >> k_int;
    prob_ss >> prob_num;
    prob_matrix[i_int][k_int]=prob_num;
    i_ss.clear();
    k_ss.clear();
    prob_ss.clear();
  }
}
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

float *retrieve_energies(string filename){
  ifstream energy_repo(filename.c_str());
  float *energy_repo_arr = new float[10000000];
  string value_str;
  int ind=0;
  while (energy_repo >> value_str && ind <10000000){
    istringstream value_ss(value_str);
    float value_num;
    value_ss >> value_num;
    energy_repo_arr[ind]=value_num;
    ind++;
  }
  return energy_repo_arr;
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
    float the_r_z=6371-z;
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
float density_to_potential(float dty, bool antineutrino){
  float to_return = (1/sqrt(2))*dty*1e-3*8.96189e-47*1e9   /1.672e-27;
  if(antineutrino){
    return -1*to_return;
  }
  else{
    return to_return;
  }
}
    float RingNode::getRadius(){
      r = sqrt(x*x + z*z);
      return r;
    }
    float RingNode::getSolidAngle(){
      solidAngle = (1./(4.*PI))*(1./((R_earth-z)*(R_earth-z) + x*x));
      return solidAngle;
    }
    float volume;
    float RingNode::getVolume(){
      volume = 2*PI*x*dx*dz;
      return volume;
    }
    void Planet::initializeCoords(bool expo){
      cout << "Initializing Coordinates" << endl;
      ofstream export_file;
      if(expo){
        //In this case, writes in "planet_coords.csv", which pairs (i,k) are Silicate Earth.
        export_file.open("planet_coords.csv");
      }
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          asArray[i][k].x=i*dx;
          asArray[i][k].z=-6371. + k*dz;
          float r = asArray[i][k].getRadius();
          if(r<6371){
            asArray[i][k].isEarth=1;
            float dk=1e3-k;
            float di=-i;
            if(r>3480){
              asArray[i][k].isSE=1;
              if(expo){
                export_file << i << ' ' << k <<' ' << 1 << endl;
              }
              if(r>(R_earth-42.5)){
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
      if(expo){
        export_file.close();
      }
    }
    void Planet::initializeDensity(){
      cout << "Initializing Density" << endl;
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
    void Planet::initializeAbundanceCrust(){
      cout << "Initializing Abundances in the Crust" << endl;
      for(int i =0 ; i<N/2;i++){
        for(int k = 0;k<N;k++){
          if(asArray[i][k].isCrust){
            asArray[i][k].abundanceU=453.193965399; //ppb
            asArray[i][k].abundanceTh=1940.64130183; //ppb
          }
          else{
            asArray[i][k].abundanceU=0.0001; //%
            asArray[i][k].abundanceTh=0.0001; //%
          }
        }
      }
    }
    void Planet::initializeAbundanceMantle(string key, string bse_model){
      cout << "Initializing Abundances in the Mantle" << endl;
      int model;
      if(bse_model=="cosmo"){model=0;}
      else if(bse_model=="geoch"){model=1;}
      else if(bse_model=="geodyn"){model=2;}
      if(key=="unif"){
        for(int i =0 ; i<N/2;i++){
          for(int k = 0;k<N;k++){
            if(asArray[i][k].isMantle){
              asArray[i][k].abundanceU=calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceU_BSE,453.193965399 )[model]; //ppb
              asArray[i][k].abundanceTh=calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceTh_BSE, 1940.64130183)[model]; //ppb
            }
          }
        }
      }
      else if(key=="two_layer"){
        float bulk_mantle_U = calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceU_BSE, 453.193965399)[model];
        float bulk_mantle_Th = calculateMantleAbundances(crustMass, mantleMass, crustMass+mantleMass, abundanceTh_BSE, 1940.64130183)[model];
        float mantleMass_fraction = 0.1*mantleMass;
        float mass_count=0;
        float limit_rad;
        int n=0;

        do {
          if(asArray[n][500].isMantle){
            mass_count+=4*PI*(asArray[n][500].r)*(asArray[n][500].r)*dx*1e3*1e9*asArray[n][500].massDensity;
            limit_rad=asArray[n][500].r;
          }
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

    void Planet::initializeFluxes(bool oscillated, string hpe_dist, string bse_model){
      totalFlux=0;
      totalUFlux=0;
      totalThFlux=0;
      cout << "Initializing Fluxes" << endl;
      float prob_matrix[500][1000];
      for(int n=0;n<N/2;n++){
        for(int m = 0;m<N;m++){
          prob_matrix[n][m]=0;
        }
      }
      import_probability("probability_planet"+hpe_dist+"_"+bse_model+".csv", prob_matrix);
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          float prob =1.;
          if(oscillated){
            prob = prob_matrix[i][k];
            asArray[i][k].meanSurvProb=prob;
            //cout << prob << " , " << i << " , "<< k << endl;
          }
          if(asArray[i][k].isMantle){//flux from mantle and crust
            asArray[i][k].neutrinoUFlux=prob*(6.)*(asArray[i][k].abundanceU*1e-9)*(0.9927)*(4.916*1e-18*1e-6)*(asArray[i][k].massDensity*1e-3)*asArray[i][k].volume*asArray[i][k].solidAngle*1e5/(238.051*1.661e-27);
            asArray[i][k].neutrinoThFlux=prob*(4.)*(asArray[i][k].abundanceTh*1e-9)*(1.)*(1.563*1e-18*1e-6)*(asArray[i][k].massDensity*1e-3)*asArray[i][k].volume*asArray[i][k].solidAngle*1e5/(232.038*1.661e-27);
            asArray[i][k].neutrinoFlux=asArray[i][k].neutrinoUFlux+asArray[i][k].neutrinoThFlux;
            asArray[i][k].relativeNeutrinoU=asArray[i][k].neutrinoUFlux/asArray[i][k].neutrinoFlux;
            asArray[i][k].relativeNeutrinoTh=asArray[i][k].neutrinoThFlux/asArray[i][k].neutrinoFlux;
          }
          else{
            asArray[i][k].neutrinoUFlux=0;
            asArray[i][k].neutrinoThFlux=0;
            asArray[i][k].neutrinoFlux= 0;
            asArray[i][k].relativeNeutrinoU=0;
            asArray[i][k].relativeNeutrinoTh=0;
          }
          totalFlux+=asArray[i][k].neutrinoFlux;
          totalUFlux+=asArray[i][k].neutrinoUFlux;
          totalThFlux+=asArray[i][k].neutrinoThFlux;
        }
      }
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          asArray[i][k].relativeNeutrino=asArray[i][k].neutrinoFlux/totalFlux;
          asArray[i][k].neutrinosProduced=  roundf(asArray[i][k].relativeNeutrino*totalNeutrinos);
          asArray[i][k].neutrinosProducedU=roundf(asArray[i][k].neutrinosProduced*asArray[i][k].relativeNeutrinoU);
          asArray[i][k].neutrinosProducedTh=roundf(asArray[i][k].neutrinosProduced*asArray[i][k].relativeNeutrinoTh);
          totalNeut+=asArray[i][k].neutrinosProduced;
        }
      }
    }
    void Planet::initializePaths(bool all, int i_o, int k_o){
      cout << "Initializing Potential (density) paths" << endl;
      omp_set_num_threads(4);
      int i, k;
      if(all){
        for(i=0;i<N/2;i++){
          for(k=0;k<N;k++){
            if(asArray[i][k].isSE){
              float z, x;
              z=asArray[i][k].z;
              x=asArray[i][k].x;
              vector<float> path;
              float path_len = get_r(R_earth-z,x );
              asArray[i][k].distanceToDetector=path_len;
              float element_num = roundf(path_len/path_resolution);
              asArray[i][k].pathLen=element_num;
              path.reserve(int(element_num));
              vector<float> times;
              times=linspace(0, path_len/300000, int(element_num));
              for(int n=0;n<int(element_num);n++){
                float t = times[n];
                float dty = density_polynomials(get_r(the_r(x, z, t, 'x'), the_r(x, z, t, 'z')));
                path.push_back(density_to_potential(dty, 1)); //eV
                }
              asArray[i][k].path=path;
            }
          }
        }
      }
        else{
          i=i_o;
          k=k_o;
          if(asArray[i][k].isEarth){
            float z, x;
            z=asArray[i][k].z;
            x=asArray[i][k].x;
            vector<float> path;
            float path_len = get_r(R_earth-z,x );
            asArray[i][k].distanceToDetector=path_len;
            float element_num = roundf(path_len/path_resolution);
            asArray[i][k].pathLen=element_num;
            path.reserve(int(element_num));
            vector<float> times;
            times=linspace(0, path_len/300000, int(element_num));
            for(int n=0;n<int(element_num);n++){
              float t = times[n];
              float dty = density_polynomials(get_r(the_r(x, z, t, 'x'), the_r(x, z, t, 'z')));
              path.push_back(density_to_potential(dty, 1)); //eV
              }
            asArray[i][k].path=path;
          }
          else{
            cout << "ERROR\nPaths cannot be initialized outside the planet."<<endl;
          }
        }
      }

    void Planet::initializeEnergySamples(){
      cout << "Initializing Energies for Neutrinos" << endl;
      float *uran_energy_repo;
      float *thor_energy_repo;
      uran_energy_repo = retrieve_energies("energy_repo_238U.knt");
      thor_energy_repo = retrieve_energies("energy_repo_232Th.knt");
      //cout << thor_energy_repo[9999239] << ',' << uran_energy_repo[9935793] << endl;
      int uran_i=0;
      int thor_i=0;
      int i, k, n, m;
      for(i=0;i<N/2;i++){
        for(k=0;k<N;k++){
          if(asArray[i][k].isSE){
            vector <float> path_U;
            vector <float> path_Th;
            int U_len = int(asArray[i][k].neutrinosProducedU);
            if (U_len>0){
              path_U.reserve(U_len);
              asArray[i][k].allowedEnergiesU.reserve(U_len);
              for(n=0;n<U_len;n++){
                path_U.push_back(1e6*uran_energy_repo[uran_i]);
                uran_i++;
              }
            }
            int Th_len = int(asArray[i][k].neutrinosProducedTh);
            if(Th_len>0){
              asArray[i][k].allowedEnergiesTh.reserve(Th_len);
              path_Th.reserve(Th_len);
              for(m=0;m<Th_len;m++){
                path_Th.push_back(1e6*thor_energy_repo[thor_i]);
                thor_i++;
              }
            }
            asArray[i][k].allowedEnergiesU=path_U;
            asArray[i][k].allowedEnergiesTh=path_Th;
          }
        }
      }
    }

    /*
    void simulateProbabilities(){
      cout << "Simulating Oscillations" << endl;
      for(int i=0;i<N/2;i++){
        for(int k=0;k<N;k++){
          int U_len = int(asArray[i][k].neutrinosProducedU);
          int Th_len = int(asArray[i][k].neutrinosProducedTh);
    //      vector <double> probabilitiesU;
  //        vector <double> probabilitiesTh;
          int n, m;
          if (U_len>0){
//            probabilitiesU.reserve(U_len);
            asArray[i][k].probabilitiesU.reserve(U_len);
            for(n=0;n<U_len;n++){
              asArray[i][k].probabilitiesU.push_back(calculateProbability(int(asArray[i][k].pathLen), asArray[i][k].path, asArray[i][k].allowedEnergiesU[n]));
            }

          }
          if(Th_len>0){
            asArray[i][k].probabilitiesTh.reserve(Th_len);
      //      probabilitiesTh.reserve(Th_len);
            for(m=0;m<Th_len;m++){
              asArray[i][k].probabilitiesTh.push_back(calculateProbability(int(asArray[i][k].pathLen), asArray[i][k].path, asArray[i][k].allowedEnergiesTh[n]));
            }

          }
        }
      }
    }
*/

    void Planet::initialize(string key, string bse_model){
      cout << "Building Planet" << endl;
      Planet::initializeCoords(0);
      Planet::initializeDensity();
      Planet::initializeAbundanceCrust();
      Planet::initializeAbundanceMantle(key, bse_model);
      Planet::initializeFluxes(0, key, bse_model);
      //Planet::initializePaths(1, 0, 0);
      //Planet::initializeEnergySamples();
      //simulateProbabilities();
      cout << "Done" << endl;
    }



/*

int main(int argc, char const *argv[]) {
  Planet *earth = new Planet();

  earth->initialize("two_layer", "geodyn");
  cout << "total flux " << earth->totalFlux << endl;
  ofstream outfile;
  outfile.open("earth_simul_plots.csv");
  for(int k=0;k<N;k++){
    for(int i =0 ; i<N/2;i++){
      outfile << earth->asArray[i][k].neutrinoFlux << ',' ;
      }
      outfile << 0 << endl;
    }
  outfile.close();

int xtest, ytest;
xtest=3;
ytest=998;
  ofstream test_path_file;
  test_path_file.open("test_path.csv");
  int test_N = int(earth->asArray[xtest][ytest].pathLen);
  for(int step=0;step<test_N;step++){
    test_path_file << earth->asArray[xtest][ytest].path[step] << endl;
  }
  test_path_file.close();
  cout << "an energy " << (earth->asArray[xtest][ytest].allowedEnergiesU[30]) << endl;
  //cout << "the prob " << calculateProbability(test_N, earth->asArray[xtest][ytest].path, (earth->asArray[xtest][ytest].allowedEnergiesU[30])) << endl;
  //calculateProbabilitiesFunctionEnergy(test_N, earth->asArray[xtest][ytest].path);
  cout << "total mass " << earth->totalMass << endl;
  delete earth;
  return 0;
}
*/
