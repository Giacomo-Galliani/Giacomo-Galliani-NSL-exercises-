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
#include "GA.h"
#include "mpi.h"


using namespace std;
 
int main (int argc, char *argv[]){

    MPI_Init(&argc, &argv); 
    //set sim type 
    int sim_type; 

    if (argc > 1)
        sim_type =  std::atoi(argv[1]);
    else 
        sim_type =0; 
    
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 

   
    //initialize GenAlg
    GA GenAlg; 
    GenAlg.initialize(); 
    GenAlg.check(); 

    //prepare output files 
    if (rank==0){
        if (sim_type==0){
            ofstream results; 
            results.open("lenght.dat"); 
            results.close(); 
        }
        else{
            ofstream results; 
            results.open("lenght_comp.dat"); 
            results.close();
        }
    }


    //perform generations and save results 
    for (int i=0; i<1000; i++){
        GenAlg.generation(); 
        if (i%5 ==0){ 
            GenAlg.paste_min_lenght(sim_type);    
        }
        if (sim_type == 0)
            if(i%3 ==0) GenAlg.migration(); //migrate every 3 generations 

    }


    //get sorted cities 
    if (sim_type ==0)
    GenAlg.sort_cities(); 



   MPI_Finalize(); 

   return 0;
}