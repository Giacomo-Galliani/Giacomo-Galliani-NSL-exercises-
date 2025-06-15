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
  

   //parameters and starting values 
   double S=100; 
   double T= 1; 
   double K=100; 
   double r =0.1; 
   double sigma= 0.25; 
   int L=100; //number of subintervals we divide [0, T]

   //direct computation

   //number of steps and blocks
   int M = 10e04; 
   int N=100; 
   int R = int (M/N); 


   //filling call[] and put[] as the vectors of the mean of each block 
   double call[N]; 
   double put[N]; 
   double appo; 
   for (int i=0; i<N; i++){ //loop over blocks 
      call[i]=0; 
      put[i]=0; 
      for (int j=0; j<R; j++){ //loop in the block 
         appo = S*exp((r- 0.5*pow(sigma, 2))*T + sigma*rnd.Gauss(0,1)*sqrt(T)) - K;
         if (appo < 0)
            put[i] -= exp(-r*T)*appo; 
         else 
            call[i] += exp(-r*T)*appo; 

      }
      put[i] /= double(R); 
      call[i] /= double(R); 
   }


   //fill the dyncall and dyncallvar vectors as the progressive mean and variance of the call in function of the number of blocks 
   double dyncall[N]; 
   
   for (int i=0; i<N; i++){
      dyncall[i] =0; 
      for (int j=0; j<i+1; j++){
         dyncall[i] += call[j]; 
      }
      dyncall[i] /= double(i+1); 
   }

   double dyncallvar[N]; 
   dyncallvar[0]=0; 

   for(int i=1; i<N; i++){
      dyncallvar[i]=0; 
      for(int j=0; j<i+1; j++){
         dyncallvar[i] += call[j]*call[j]; 
      }

      dyncallvar[i] /= double(i+1); 

      dyncallvar[i] -= dyncall[i]*dyncall[i]; 
   }

   //fill the dynput and dynputvar vectors as the progressive mean and variance of the put in function of the number of blocks 
   double dynput[N]; 
   double dynputvar[N]; 
   
   for (int i=0; i<N; i++){
      dynput[i] =0; 
      for (int j=0; j<i+1; j++){
         dynput[i] += put[j]; 
      }
      dynput[i] /= double(i+1); 
   }

   dynputvar[0]=0; 

   for(int i=1; i<N; i++){
      dynputvar[i]=0; 
      for(int j=0; j<i+1; j++){
         dynputvar[i] += put[j]*put[j]; 
      }
      dynputvar[i] /= double(i+1); 

      dynputvar[i] -= dynput[i]*dynput[i]; 
   }   



   //print the results 
   ofstream fileout1; 
   ofstream fileout2; 

   fileout1.open("mean_01.txt");
   fileout2.open("err_01.txt"); 

   fileout1 << dyncall[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1  << dyncall[i] << " ";
      fileout2 << sqrt(dyncallvar[i]/double(i)) << " "; 
   }

   fileout1.close(); 
   fileout2.close(); 

   fileout1.open("integers.txt"); 
   fileout1 << N << " " << R; 
   fileout1.close(); 

   fileout1.open("mean_02.txt");
   fileout2.open("err_02.txt"); 

   fileout1 << dynput[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1  << dynput[i] << " ";
      fileout2 << sqrt(dynputvar[i]/double(i)) << " "; 
   }

   fileout1.close(); 
   fileout2.close();


   //computation using subintervals  


   //fill the new call[] and put[] vectors 
   for (int i=0; i<N; i++){ //loop over blocks 
      call[i] =0; 
      put[i] =0; 

      for (int j=0; j<R; j++){ //loop in the block 
         appo = S*exp((r- 0.5*pow(sigma, 2))*double(T/L) + sigma*rnd.Gauss(0,1)*sqrt(double(T/L))); 
         for (int k=1; k<L; k++){ //loop in the subintervals 
            appo *=exp((r- 0.5*pow(sigma, 2))*double(T/L) + sigma*rnd.Gauss(0,1)*sqrt(double(T/L)));  
         }

         if (appo - K< 0)
            put[i] += exp(-r*T)*(K-appo); 
         else 
            call[i] += exp(-r*T)*(appo - K); 

      }
      put[i] /= double(R); 
      call[i] /= double(R); 

   }

   
   //fill the dyncall and dyncallvar vectors as the progressive mean and variance of the call in function of the number of blocks 
   for (int i=0; i<N; i++){
      dyncall[i] =0; 
      for (int j=0; j<i+1; j++){
         dyncall[i] += call[j]; 
      }
      dyncall[i] /= double(i+1); 
   }


   dyncallvar[0]=0; 

   for(int i=1; i<N; i++){
      dyncallvar[i]=0; 
      for(int j=0; j<i+1; j++){
         dyncallvar[i] += call[j]*call[j]; 
      }
      dyncallvar[i] /= double(i+1); 

      dyncallvar[i] -= dyncall[i]*dyncall[i]; 
   }

   //fill the dynput and dynputvar vectors as the progressive mean and variance of the put in function of the number of blocks 
   for (int i=0; i<N; i++){
      dynput[i] =0; 
      for (int j=0; j<i+1; j++){
         dynput[i] += put[j]; 
      }
      dynput[i] /= double(i+1); 
   }

   dynputvar[0]=0; 

   for(int i=1; i<N; i++){
      dynputvar[i]=0; 
      for(int j=0; j<i+1; j++){
         dynputvar[i] += put[j]*put[j]; 
      }
      dynputvar[i] /= double(i+1); 

      dynputvar[i] -= dynput[i]*dynput[i]; 
   }   


   //print the results 
   fileout1.open("mean_03.txt");
   fileout2.open("err_03.txt"); 

   fileout1 << dyncall[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1  << dyncall[i] << " ";
      fileout2 << sqrt(dyncallvar[i]/double(i)) << " "; 
   }

   fileout1.close(); 
   fileout2.close(); 


   fileout1.open("mean_04.txt");
   fileout2.open("err_04.txt"); 

   fileout1 << dynput[0] << " "; 
   fileout2 << 0 << " "; 

   for (int i=1; i<N; i++){
      fileout1  << dynput[i] << " ";
      fileout2 << sqrt(dynputvar[i]/double(i)) << " "; 
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
