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
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include "random.h"
#include <cmath>
#include "SA.h"
#include <vector>

using namespace std;

int main (int argc, char *argv[]){

   //initialize SA and prepare output files 
   SA SimAnn; 
   SimAnn.initialize(); 
   ofstream energy; 
   energy.open("energy.dat"); 
   ofstream parameters; 
   parameters.open("parameters.dat"); 
   ofstream pro_energy; 
   pro_energy.open("pro_energy.dat"); 


   //set beta (1/temp) of the simulation
   vector <double> beta;  
   beta.push_back(1); 
   beta.push_back(10); 
   beta.push_back(20); 
   beta.push_back(50); 
   beta.push_back(100);
   beta.push_back(200); 
   beta.push_back(500); 


   //simulation
   for (unsigned int i=0; i<beta.size()+1; i++){ //loop over temperatures 

      SimAnn.equilibrate(); 

      //print results 
      energy << SimAnn.get_beta() << "  " << SimAnn.get_energy() << " " << SimAnn.get_error() << endl;  
      parameters << SimAnn.get_mu() << " " << SimAnn.get_sigma() << endl; 

      //prepare parameters for next temperature 
      SimAnn.set_beta(beta[i]); 
      SimAnn.set_step(i);
      SimAnn.reset_ratio2(); 
   }

   energy.close(); 
   //print psi2 distribution 
   ofstream histo; 
   histo.open("histo.dat"); 

   double x = SimAnn.get_mu(); 

   for (int i=0; i<1000000; i++){
      x = SimAnn.move(x, SimAnn.get_mu(), SimAnn.get_sigma());
      histo << x << " ";  
   }


   //print progressive value of GS energy 
   vector <double> progressive_H = SimAnn.get_progressive_energy(); 
   pro_energy << "#BLOCK:      H_AVE:      ERROR:" << endl; 

   for (unsigned int i=0; i<progressive_H.size()/2; i++){
      pro_energy << i+1 << setw(12) << progressive_H[i] << setw(12) << progressive_H[progressive_H.size()/2 + i] << endl; 
   }

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
