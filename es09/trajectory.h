#ifndef __TRAJECTORY__
#define __TRAJECTORY__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "random.h"


using namespace std;
using namespace arma;

class trajectory {

private: 


Random _rnd; 
int _ncities; //number of cities 
vec _trajectory; //vector of the number of cities 
vec _cities_position_x; //x component of the cities 
vec _cities_position_y; //y component of the cities 

public: 

void check(); //Checks that the trajectory is an allowed one
void initialize_rnd();  //initialize _rnd

//set-functions  
void set_ncities(int);      //set _ncities  
void set_component(int,int);    //set the n-th component of _trajectory 
void set_trajectory(vec);      //sets _trajectory 
void set_x_position(vec);     //sets the cities x position
void set_y_position(vec);     //sets the cities y position

//Get-functions
int get_ncities();          //Get the number of cities 
int get_city(int); //Get the i-th component of the vec trajectory
vec get_trajectory();       //Get the vec _trajectory
int get_position(int);     //Gets the position of the given component
double get_lenght();        //Gets the lenght of the trajectory (the cost function)

//Mutation functions 
void pair_permutation(int, int);    //swaps a random pair 
void shift(int, int ,int);               //shifts for a random n a random m contiguos cities
void block_permutation(int, int, int);   //swaps two random blocks 
void order_inversion(int, int);      //inverts the order of a block 
}; 




#endif 