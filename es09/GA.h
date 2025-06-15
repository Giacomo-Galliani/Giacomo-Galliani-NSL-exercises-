#ifndef __GA__
#define __GA__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "trajectory.h"


using namespace std;
using namespace arma;

class GA {

private: 

int _sim_type; 
int _n_cities; 
int _n_trajectories; 
int _n_elite; 
int _n_steps=0; 
Random _rnd; 
arma :: field <trajectory> _trajectory; 
vec _cities_position_x;
vec _cities_position_y; 

public: 


void initialize(int); //creats starting population of n trajectories of m cities   
void crossover(int, int); //crossover function 
void order();             //orders the trajectories in function of the lenght 
int selection(int, int);          //selects an individual in a certain range
void generation();        //generates new individuals 
void mutation();         //mutates the individuals with random probabilities  
void check();             //Checks that all the indiduals are valid 
double get_min_lenght();  //Gets the min lenght of the trajectories   
int get_nsteps();         //Gets number of the generations
double get_lenght(int);   //Gets the best trajectory length 
double get_half_mean_lenght(); //Gets the mean of the half shourtest lenghts
trajectory get_trajectory(int); //Gets the i-th trajectory of the field 
void sort_cities();        //Sortes citities depending on the best trajectory found 


}; 
#endif 