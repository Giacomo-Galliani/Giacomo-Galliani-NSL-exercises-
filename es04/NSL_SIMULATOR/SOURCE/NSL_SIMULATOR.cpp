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

  int nconf = 1;
  System SYS;

  //initialize sitem, reads half box fcc configuration
  SYS.initialize();
  SYS.initialize_properties();

  SYS.read_configuration_half_fcc(); 
  

  //prepare a file with velocities configuration to read 
  ofstream config; 
  config.open("../INPUT/CONFIG/configvelocities.xyz"); 
  config << SYS.get_npart() << endl; 


  
  config << "velocità_deltiformi" << endl; 

  for (int i=0; i<SYS.get_npart()/(2.*double(SYS.get_ndim())); i++){
    for (int l=0; l<2*SYS.get_ndim(); l++){
      if (l==0)
        config << "LJV "  <<sqrt(SYS.get_temp()*3) << " " << 0 << " " << 0 << endl; 
      if (l==1)
        config << "LJV "  << 0  << " " << sqrt(SYS.get_temp()*3) << " " << 0 << endl;
      if (l==2)
        config << "LJV "  << 0  << " " << 0 << " " << sqrt(SYS.get_temp()*3) << endl;
      if (l==3)
        config << "LJV "  <<-sqrt(SYS.get_temp()*3) << " " << 0 << " " << 0 << endl; 
      if (l==4)
        config << "LJV "  << 0  << " " << -sqrt(SYS.get_temp()*3) << " " << 0 << endl;
      if (l==5)
        config << "LJV "  << 0  << " " << 0 << " " << -sqrt(SYS.get_temp()*3) << endl;

    }
  }

  //read velocities 
  SYS.read_configuration_velocities();

  //start of the simulation

  SYS.block_reset(0);
  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if (i==0){
        if(j%50 == 0 && nconf <5){
          SYS.write_pofv(nconf); 
          nconf++;
        }
      }
    }
     
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  
  //invert velocities and tries to recover initial configuration
  SYS.invert_velocities(); 
  SYS.new_measure(); 


   for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(j%50 == 0){
        //SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
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

