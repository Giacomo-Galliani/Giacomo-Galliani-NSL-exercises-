/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;

double f(double x); //function to evaluate with a uniform distribution

double f(double x){

   return (M_PI/2)*cos(M_PI*x/2); 
}

double f2(double x){ //function to evaluate with importance sampling

   return (((M_PI/4)*cos(M_PI*x/2.))/(1-x)); 
}

 
int main (int argc, char *argv[]){

   //initialize rnd
   
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();


   //Evaluation of the integral using a uniform distribution

   int M = 100000; //number of steps
   int N = 100; //number of blocks
   int L = int (M/N); //lenght of the blocks 

   double mean[N]; //prepare the vector of the acutal mean of each block 
   

   //fill the vector 
   for (int i=0; i<N; i++){
      mean[i] =0; 
      for (int j=0; j<L; j++){
         mean[i]+= f(rnd.Rannyu()); 
      }
      mean[i] /= double(L); 
   }

   
   //prepare and fill the progressive mean vector 
   double dynmean[N]; 
   for (int i=0; i<N; i++){

      dynmean[i] =0; 
      for (int j=0; j<i+1; j++){
         dynmean[i] += mean[j]; 
      }

      dynmean[i] /= double(i+1); 
   }

   
   //prepare and fill the progressive variance vector 
   double dynvar[N]; 
   dynvar[0] =0; 
   for (int i=1; i<N; i++){
      dynvar[i] =0; 
      for (int j=0; j<i+1; j++){
         dynvar[i] += pow(mean[j], 2); 
      }

      dynvar[i] /= double(i+1); 
      dynvar[i] -= pow(dynmean[i], 2); 

   }


   //print to file the results
   ofstream fileout1;
   ofstream fileout2; 
   fileout1.open("mean_01.txt"); 
   fileout2.open("err_01.txt"); 

   fileout1 << dynmean[0] << " "; 
   fileout2 << dynvar[0]<< " "; 

   for (int i=1; i<N; i++){
      fileout1 << dynmean[i] << " "; 
      fileout2 << sqrt(dynvar[i]/double(i)) << " "; 
   }


   fileout1.close(); 
   fileout2.close(); 

   fileout1.open("integers.txt"); 
   fileout1 << N << " " << L; 
   fileout1.close(); 

   //evaluation of the integral with importance sampling technique 

   //fill the mean vector with the actual mean of each block
   for (int i=0; i<N; i++){
      mean[i]=0; 
      for (int j=0; j<L; j++){
         mean[i]+= f2(rnd.Taylor())/double(L); //rnd.Taylor() gives the correctly sampled distribution 
      }
   }

   //fill the vector of progressive means 
   for (int i=0; i<N; i++){
      dynmean[i] =0; 
      for (int j=0; j<i+1; j++){
         dynmean[i] += mean[j]; 
      }

      dynmean[i] /= double(i+1); 
   }

   //fill the vector of progressive variance 
   dynvar[0] =0; 
   for (int i=1; i<N; i++){
      dynvar[i] =0; 
      for (int j=0; j<i+1; j++){
         dynvar[i] += pow(mean[j], 2); 
      }

      dynvar[i] /= double(i+1); 
      dynvar[i] -= pow(dynmean[i], 2); 

   }

   //print to file the results
   fileout1.open("mean_02.txt"); 
   fileout2.open("err_02.txt"); 

   fileout1 << dynmean[0] << " "; 
   fileout2 << dynvar[0]<< " "; 

   for (int i=1; i<N; i++){
      fileout1 << dynmean[i] << " "; 
      fileout2 << sqrt(dynvar[i]/double(i)) << " "; 
   }

   fileout1.close(); 
   fileout2.close(); 


   
   //random walk in a lattice 
   double x[3]; 
   //starting in the (0,0,0) position 
   for (int i=0; i<3; i++){
      x[i]=0; 
   }
   double a=1; //setting the step 
   double rand; 
   double mean2[N*N]; 
   
   //setting to zero the accumulators 
   for (int i=0; i<N*N; i++){ 
      mean2[i]=0; 
   }
   for (int i=0; i<N; i++){
      mean[i]=0; 
      dynvar[i]=0; 
   }
   

   for (int i=0; i<N; i++){ //loop over blocks 
      for(int j=0; j<L; j++){ //loop in each block 
         for (int l=0; l<N; l++){ //loop over the steps 

         //select the direction of the step and perform the step 
         rand = rnd.Rannyu(-3, 3); 
         if (rand >0)
            x[static_cast<int>(rand)] +=a; 
         else 
            x[static_cast<int>(-rand)] -=a;

         //saves the distance from the origin 
         mean2[l+ i*N] += sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])/double(L); 
         }

         //reset starting point 
         for (int i=0; i<3; i++){
            x[i]=0; 
         }

      }
   }


   //fill the vector mean with the means in each block 
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         mean[i] += mean2[i +j*N]/double(N); 
      }
   }

   //fill the vector dynvar with the progressive variance
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         dynvar[i] += mean2[i+j*N]*mean2[i+j*N]/double(N); 
      }

      dynvar[i] -= mean[i]*mean[i]; 
   }

   //print to file the results 
   fileout1.open("mean_03.txt"); 
   fileout2.open("err_03.txt");

   fileout1 << mean[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1 << mean[i] << " "; 
      fileout2 << sqrt(dynvar[i]/double(i)) << " ";  
   }

   fileout1.close(); 
   fileout2.close(); 

   //also print the mean squared with the appropriate error (for statistical purposes)

   fileout1.open("mean_04.txt"); 
   fileout2.open("err_04.txt"); 

   fileout1 << mean[0]*mean[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1 << mean[i]*mean[i] << " "; 
      fileout2 << mean[i]*sqrt(2*dynvar[i]/double(i)) << " ";  
   }

   fileout1.close(); 
   fileout2.close(); 

   //random walk in the continuum 

   //set to zero the accumulators 
   for (int i=0; i<N*N; i++){ 
      mean2[i]=0; 
   }
   for (int i=0; i<N; i++){
      mean[i]=0; 
      dynvar[i]=0; 
   }

   double phi; //new variable needed for the angle with the z-axis  


   for (int i=0; i<N; i++){ //loop over blocks 
      for(int j=0; j<L; j++){ //loop in the block 
         for (int l=0; l<N; l++){ //loop over the steps 

         //select the direction (two angles)
         rand = 2*M_PI*rnd.Rannyu(); 
         phi = rnd.Sin(); 
         
         //perform the step 
         x[0] += a*sin(phi)*cos(rand); 
         x[1] += a*sin(phi)*sin(rand); 
         x[2] += a*cos(phi); 

         //save the distance from the origin 
         mean2[l+ i*N] += sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])/double(L); 
         }

         //reset the starting point 
         for (int i=0; i<3; i++){
            x[i]=0; 
         }

      }
   }

   //fill the mean vector 
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         mean[i] += mean2[i +j*N]/double(N); 
      }
   }

   //fill the dynvar vector 
   for (int i=0; i<N; i++){
      for (int j=0; j<N; j++){
         dynvar[i] += mean2[i+j*N]*mean2[i+j*N]/double(N); 
      }

      dynvar[i] -= mean[i]*mean[i]; 
   }

   //print to file the results 
   fileout1.open("mean_05.txt"); 
   fileout2.open("err_05.txt");

   fileout1 << mean[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1 << mean[i] << " "; 
      fileout2 << sqrt(dynvar[i]/double(i)) << " ";  
   }

   fileout1.close(); 
   fileout2.close(); 

   //also print the mean squared with the appropriate error (for statistical purposes) 

   fileout1.open("mean_06.txt"); 
   fileout2.open("err_06.txt"); 

   fileout1 << mean[0]*mean[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1 << mean[i]*mean[i] << " "; 
      fileout2 << mean[i]*sqrt(2*dynvar[i]/double(i)) << " ";  
   }

   fileout1.close(); 
   fileout2.close(); 


 
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/