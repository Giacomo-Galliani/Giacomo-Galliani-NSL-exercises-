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
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  //initialize SYS
  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  
  //parameters for T dependence
  double tstart = 0.5; 
  double deltat = 0.15; 
  
  for (int l=0; l<11; l++){ //loop over temperatures 

    SYS.set_temp(tstart + l*deltat); 
    
  
    //let the system equilibrate before starting to take measures 
    for (int i=0; i<200; i++){
      SYS.step(); 
    }

    //actual measure 
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }

    //print results and reset accumulators 
    SYS.save_T_results(l); 
    SYS.new_measurement(); 
  }
  
  SYS.finalize();
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
