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
#include <cmath>
#include "random.h"

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

   //Evaluation of the integral 


   int M = 100000; //number of steps 
   int N = 100; //number of blocks
   int L = int (M/N); //length of the blocks

   //initializing and filling the vector of the means of each block 
   double mean[N];

   for (int i=0; i<N; i++){
      mean[i] =0; 
      for (int j=0; j<L; j++){
         mean[i] += rnd.Rannyu();
      }

      mean[i] /= double(L); 
   }
   
   //initializing and filling the vector of progressive means 
   double dynmean[N]; 

   for (int i=0; i<N; i++){

      dynmean[i] = 0; 
      for (int j=0; j<i+1; j++){
         dynmean[i] += mean[j]; 
      }

      dynmean[i] /= double(i+1); 

   }


   //initializing and filling the dynvar vector of the progressive variance 
   double dynvar[N]; 
   dynvar[0] = 0; 

   for (int i=1; i<N; i++){
      dynvar[i] = 0; 
      for (int j=0; j<i+1; j++){
         dynvar[i]+= pow(mean[j], 2);
      }
      dynvar[i] /= double(i+1); 
      dynvar[i]-= pow(dynmean[i], 2);
   }

   
   //print to file the results
   ofstream fileout1; 
   ofstream fileout2; 


   fileout1.open("mean_01.txt"); 
   fileout2.open("err_01.txt"); 

   fileout1 << dynmean[0]; 
   fileout2 << 0; 

   for(int i=1; i<N; i++){
      fileout1 << " "<< dynmean[i]; 
      fileout2 << " " << sqrt(dynvar[i]/double(i)); 
   }


   fileout1.close();
   fileout2.close(); 

 
   fileout1.open("integers.txt"); 
   fileout1 << N << " " << L; 

   fileout1.close(); 

   //Evaluation of the variance of the integral 

   //filling the vectors the means of each block 
   for (int i=0; i<N; i++){
      mean[i] =0; 
      for (int j=0; j<L; j++){
         mean[i] += pow(rnd.Rannyu()-0.5, 2);
      }

      mean[i] /= double (L);  
   }
   
   //filling the vector of the progressive means 
   for (int i=0; i<N; i++){

      dynmean[i] = 0; 
      for (int j=0; j<i+1; j++){
         dynmean[i] += mean[j]; 
      }

      dynmean[i] /= double(i+1); 

   }

   //filling the vector of the progressive variance
   for (int i=1; i<N; i++){

      dynvar[i] = 0; 
      for (int j=0; j<i+1; j++){
         dynvar[i]+= pow(mean[j], 2);
      }

      dynvar[i] /= double(i+1); 
      dynvar[i]-= pow(dynmean[i], 2);

   }
   
   //printing on file the results 

   fileout1.open("mean_02.txt"); 
   fileout2.open("err_02.txt"); 

   fileout1 << dynmean[0]; 
   fileout2 << 0; 

   for(int i=1; i<N; i++){
      fileout1 << " "<< dynmean[i]; 
      fileout2 << " " << sqrt(dynvar[i]/double(i)); 
   }


   fileout1.close();
   fileout2.close(); 

   //computation of the pearson test on each block 

   M = 100; //number of blocks 
   L = 10E04; //length of the blocks 
   int hits[M]; //histogram 
   double chi[M]; //chi values for each block 
   //setting to zero the counters
   for (int i=0; i<M; i++){
      hits[i] =0; 
      chi[i]=0; 
   }


   for (int i=0; i<M; i++){ //loop over the blocks 
      
      //randomly filling the istogram 
      for (int j=0; j<L; j++){
         hits[static_cast<int>(M*rnd.Rannyu())]++; 
      }

      //compute the chi for the i-block 
      for (int j=0; j<M; j++){
         chi[i] += M*pow(hits[j]-  static_cast<double>(L)/M, 2)/static_cast<double>(L); 
      }

      //reset the instogram for next block 
      for (int i=0; i<M; i++){
         hits[i] =0; 
      }
   }

   //print to file the results 
   fileout1.open("chi.txt"); 
   for (int i=0; i<M; i++){
      fileout1 << chi[i] << " "; 
   }

   //simulating one trhow of a dice 

   fileout1.open("histo11.txt");
 
   for (int i=0; i<M; i++){ //normal dice 
      fileout1 << static_cast<int>(6*rnd.Rannyu()) +1<< " "; //print to file 
   }

   fileout1.close(); 

   double rand1;
   fileout1.open("histo12.txt");  
   for (int i=0; i<M; i++){ //exponential dice 
      rand1 =rnd.Exp(1);  
      if (rand1 < 6)
      fileout1 << static_cast<int>(rand1) +1<< " "; //print to file 
   }

   fileout1.close();

   fileout1.open("histo13.txt"); 
   for (int i=0; i<M; i++){
      rand1 =abs(rnd.Lorentz(0,1)); //Lorentzian dice 
      if (rand1 < 6)
      fileout1 << static_cast<int>(rand1) +1<< " "; //print to file 
   }

   fileout1.close(); 

   //now simulating N=2 throws of a dice

   fileout1.open("histo21.txt");
 
   for (int i=0; i<M; i++){//normal dice 
      fileout1 << 0.5*(static_cast<int>(6*rnd.Rannyu()+1)+ static_cast<int>(6*rnd.Rannyu()+1)) << " "; //print to file 
   }

   fileout1.close(); 

   double rand2;
   fileout1.open("histo22.txt");  

   for (int i=0; i<M; i++){ //exponential dice 
      do { rand1 =static_cast<int>(rnd.Exp(1))+1; 
      } while(rand1 > 6);
      do { rand2 =static_cast<int>(rnd.Exp(1))+1; 
      } while(rand2 > 6); 
      fileout1 << 0.5*(rand1 + rand2) << " ";  //print to file 
   }

   fileout1.close();

   fileout1.open("histo23.txt"); 

   for (int i=0; i<M; i++){ //Lorentzian dice 
      do { rand1 =abs(static_cast<int>(rnd.Lorentz(0,1)))+1; 
      } while(rand1 > 6);
      do { rand2 =abs(static_cast<int>(rnd.Lorentz(0,1)))+1; 
      } while(rand2 > 6); 
      fileout1 << 0.5*(rand1 + rand2) << " ";  //print to file 
   }

   fileout1.close(); 

   //simulating N=10 throws of a dice 

   fileout1.open("histo31.txt");
   double sum=0; 
 
   for (int i=0; i<M; i++){ //normal dice 
      for (int j=0; j<10; j++){
         sum += (static_cast<int>(6*rnd.Rannyu())+1)/double(10); 
      }

      fileout1 << sum << " "; //print to file 
      sum=0; 
   }

   fileout1.close(); 
   fileout1.open("histo32.txt");  

   for (int i=0; i<M; i++){ //exponential dice 
      
      for (int j=0; j<10; j++){
         do {
            rand1 = static_cast<int>(rnd.Exp(1)+1); 
         } while (rand1 >6); 
         sum += rand1/double(10); 
      }
      
      fileout1 << sum << " "; //print to file
      sum=0; 
   }

   fileout1.close();

   fileout1.open("histo33.txt"); 

   for (int i=0; i<M; i++){ //lortentzian dice 
      
      for (int j=0; j<10; j++){
         do {
            rand1 = static_cast<int>(abs(rnd.Lorentz(0,1))+1); 
         } while (rand1 >6); 
         sum += rand1/double(10); 
      }
      
      fileout1 << sum << " "; //print to file 
      sum=0; 
   }

   fileout1.close(); 

   //simulating N=1000 throws of a dice 

   fileout1.open("histo41.txt");
 
   for (int i=0; i<M; i++){ //normal dice 
      for (int j=0; j<100; j++){
         sum += static_cast<int>(6*rnd.Rannyu()+1)/double(100); 
      }

      fileout1 << sum << " "; //print to file 
      sum=0; 
   }

   fileout1.close(); 
   fileout1.open("histo42.txt");  

   for (int i=0; i<M; i++){ //exponential dice 
      
      for (int j=0; j<100; j++){
         do {
            rand1 = static_cast<int>(rnd.Exp(1)+1); 
         } while (rand1 >6); 
         sum += rand1/double(100); 
      }
      
      fileout1 << sum << " "; //print to file 
      sum=0; 
   }

   fileout1.close();

   fileout1.open("histo43.txt"); 

   for (int i=0; i<M; i++){ //Lortentzian dice 
      
      for (int j=0; j<100; j++){
         do {
            rand1 = static_cast<int>(abs(rnd.Lorentz(0,1))+1); 
         } while (rand1 >6); 
         sum += rand1/double(100); 
      }
      
      fileout1 << sum << " "; //print to file 
      sum=0; 
   }

   fileout1.close(); 


   //simulation of the Buffon esperiment 

   double d = 1; //distance of the lines 
   double l = 0.5; //lenght of the stick 
   M = 100000; //number of throws 
   N = 100; //number of blocks 
   L = int (M/N); //lenght of the blocks


   
   /*for (int i=0; i<N; i++){
      mean[i] =0; 
      for (int j=0; j<L; j++){

         if (d*rnd.Rannyu() + l*sqrt((1-pow(rnd.Rannyu(),2))) > d)
            mean[i]++; 
      }

      mean[i] /= double(L); 

      mean[i] = 2*l/(d*mean[i]); 
   }*/

   //simulate experiment 
   for (int i=0; i<N; i++){
      mean[i] =0; //reset accumulators 
      for (int j=0; j<L; j++){
         if(d*rnd.Rannyu() + l*sin(M_PI*rnd.Rannyu()) > d) //check if the needle intersects the lines 
            mean[i] ++; 
      }
      mean[i] /= double(L); 
      mean[i] =2*l/(d*mean[i]); 
   }
 
   //compute progressive mean and error 
   for (int i=0; i<N; i++){

      dynmean[i] = 0; 
      for (int j=0; j<i+1; j++){
         dynmean[i] += mean[j]; 
      }

      dynmean[i] /= double(i+1); 

   }

   dynvar[0] = 0; 
   for (int i=1; i<N; i++){

      dynvar[i] = 0; 
      for (int j=0; j<i+1; j++){
         dynvar[i]+= pow(mean[j], 2);
      }

      dynvar[i] /= double(i+1); 
      dynvar[i]-= pow(dynmean[i], 2);
   }


   //print to file results 
   fileout1.open("mean_03.txt"); 
   fileout2.open("err_03.txt"); 

   for (int i=0; i<N; i++){
      fileout1 << dynmean[i] << " "; 
      fileout2 << dynvar[i] << " "; 
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
