/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Gallià
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "GA.h"

using namespace std;
 
int main (int argc, char* argv[]){

   //checks simulation type 
   int sim_type; 
   if (argc !=2){
      do {
         cout << "Please enter type of simulation: " << endl; 
         cout << "type 0 for cities on a circumference, type 1 for cities inside a square" << endl; 
         cin >> sim_type; 
      } while (sim_type <0 || sim_type >1); 
   }
    
   if (argc ==2){
      sim_type = std::atoi(argv[1]);
      if (sim_type ==0)
         cout << "perfoming TSP with cities on a circumference"<< endl; 
      else if (sim_type ==1)
         cout << "performing TSP with cities inside a square" << endl; 
      else {
         do{
            cout << "Please enter type of simulation: " << endl; 
            cout << "type 0 for cities on a circumference, type 1 for cities inside a square" << endl; 
            cin >> sim_type; 
         } while (sim_type <0 || sim_type >1);
      }  

   }

   ofstream config; 
   Random Rnd; 
   
   //initialize Rnd
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
            Rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   Rnd.SaveSeed();

   
   //generate cities position 
   config.open("config.xy"); 
   config << 34 << " " << 400 << endl; 
   double theta; 
   if (sim_type ==0){
      for (int i=0; i<34; i++){

         theta = 2*M_PI*Rnd.Rannyu(); 
         config << cos(theta) << " " << sin(theta) << endl; 
      }
   }

   if (sim_type ==1){
      for (int i=0; i<34; i++){
         config << Rnd.Rannyu() << " " << Rnd.Rannyu() << endl; 
      }
   }
 

   //initialize GenAlg
   GA GenAlg; 
   GenAlg.initialize(sim_type); 
   GenAlg.check(); 
   

   ofstream results; 
   if (sim_type==0)
      results.open("lenght_0.dat"); 
   if (sim_type==1)
      results.open("lenght_1.dat"); 

   //performs generations and print lenght results 
   for (int i=0; i<400;  i++){
      if (i%10 ==0 || i==399){ 
         results << GenAlg.get_nsteps() << " " << GenAlg.get_min_lenght() << " " << GenAlg.get_half_mean_lenght() << endl; }
      GenAlg.generation(); 
   }

   //prints file with cities in the best order 
   GenAlg.sort_cities(); 


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
