#ifndef __SA__
#define __SA__

#include <iostream>
#include <cstdlib>
#include <fstream>
#include "random.h"
#include <vector>

using namespace std; 


double psi (double, double, double); //wave function
double psi_2(double, double, double); //modulus square of the wave function 
double KE(double, double, double); //local kinetic energy 
double V(double); //potential energy 


class SA {

private: 

Random _rnd; //random generator 
int _ntrial=0; //number of metropolis trial moves (in psi2 evaulation)
int _naccepted=0; //number of metropolis accepted moves (in psi2 evaluation)
int _ntrial2 =0; //number of metropolis trial moves (in simulated anniling)
int _naccepted2 =0; //number of metropolis accepted moves (in simulated anniling)
int _nblocks; //number of blocks 
int _nsteps; //numer of total steps 
int _block_length; //lenght of the blocks 
double _sigma; //wave function current parameter
double _mu; //wave function current parameter
double _beta; //inverse temperature  
double _mu_new; //wave function proposed parameter
double _sigma_new; //wave function proposed parameter 
vector<double> _H_current; //current hamiltionian value and error 
vector<double> _H_new; //poposed hamiltionian value and error 
double _step; //maximum step to move parameters 

public: 

void initialize(); //initialize algorithm 
double move(double, double, double); //given a position x, gives xnew distribuited according to psi2
bool metro(double, double, double, double); //metropolis for psi2 distribution
double get_ratio(); //get acceptance ratio for psi2 
double get_ratio2(); //get acceptance ratio for simulated anniling
vector<double> get_hamiltonian(double, double); //gets hamiltionian and uncertanty given mu and sigma  
void equilibrate(); //equilibrate system after change of temperature 
void move_parameters(); //try to move parameters 
double get_sigma(); //get sigma of the current wave function
double get_mu(); //get mu of the current wave function 
void set_beta(double);
double get_beta(); 
void set_mu(double); 
void set_sigma(double); 
bool metro2(); //metropolis for simulated anniling 
double get_energy(); //get hamiltionian for the current wace function
void update_energy(double, double); //update hamiltionian and uncertanty for new parameters 
void change_energy(double, double); //compute hamiltionian and uncertanty for new parameters 
double get_energy_new(); //get energy of the proposed new parameters 
void reset_ratio2(); //sets to zero the metropolis simulated anniling ratio 
double get_error(); //get uncertanty of the hamiltonian of the current wave function
void set_step(int); //tune step to fix the acceptance rate
double get_step(); //gets current step 
vector <double> get_progressive_energy(); //get progressive mean values and uncertanties of the hamiltonian 

}; 



















#endif 