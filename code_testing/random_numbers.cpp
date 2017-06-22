//g++ -o random_numbers.o random_numbers.cpp `gsl-config --cflags --libs`
#include <iostream>
using namespace std;
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
 int main(int argc, char const *argv[]) {
   //Define a variable Gen_Type for the type of random number generator
   const gsl_rng_type *Gen_Type;
   gsl_rng *Gen;
   int i, N;
   //cout << "How many numbers do you want?" << endl;
   //cin >> N;
   N=10000;
   //Setup environment variables for generator
   gsl_rng_env_setup();
   //Set type as taus because it is the fastest
   Gen_Type = gsl_rng_taus;
   //Allocate necessary memory space
   Gen = gsl_rng_alloc(Gen_Type);
   for(i=0;i<N;i++){
     //Returns a positive number uniformly distributed in the range [0,1)
     //double rand_num = gsl_rng_uniform(Gen);
     //Return values that follow a normal distribution.
     double rand_num = gsl_ran_gaussian(Gen, 0.1);
     cout << rand_num << endl;

   }
   //Free allocated memory
   gsl_rng_free(Gen);
   return 0;
 }
