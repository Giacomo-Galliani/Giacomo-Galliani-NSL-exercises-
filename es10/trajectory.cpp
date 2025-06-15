#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include "trajectory.h"

using namespace std;
using namespace arma;

void trajectory :: check (){
    
    if (_trajectory.size() != _ncities){
        cerr << "traiettorie inizializzate male" << endl; 
        exit(EXIT_FAILURE); 
    }
    
    
    int sum=0; 
    for (int i=0; i<_ncities; i++){
        sum+= _trajectory(i); 
    }

    if (2*sum != _ncities*(_ncities +1)){
        cerr << "Generated a not allowed trajectory" << endl; 
        exit(EXIT_FAILURE); 
    }
    else 
        return;
}


void trajectory :: initialize_rnd(){

    //initialize _rnd    
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
                _rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    _rnd.SaveSeed();

}

void trajectory :: set_ncities(int n){
    _ncities =n; 
    _trajectory.resize(_ncities); 
    return; 
}

void trajectory :: set_component(int index, int value){
    _trajectory(index) = value; 
}

void trajectory :: set_trajectory(vec trajectory){
    _trajectory = trajectory; 
    return; 
}

void trajectory :: set_x_position(vec x){

    _cities_position_x = x; 
    return; 
}

void trajectory :: set_y_position(vec y){

    _cities_position_y = y; 
    return; 
}

int trajectory :: get_ncities(){
    return _ncities; 
}

int trajectory :: get_city(int dim){
    return _trajectory(dim); 
}

vec trajectory :: get_trajectory(){
    return _trajectory; 
}


int trajectory::get_position(int n) {
    for (int i = 0; i < _ncities;  i++) {
        if (_trajectory(i) == n)
            return i;
    }

    std::cerr << "Errore: elemento " << n << " non trovato nella traiettoria." << std::endl;
    exit(EXIT_FAILURE);
}



double trajectory :: get_lenght(){

    double lenght=0; 
    for (int i=0; i<_ncities - 1; i++){  
        lenght += pow((_cities_position_x(_trajectory(i)-1)-_cities_position_x(_trajectory(i+1)-1)),2) + pow((_cities_position_y(_trajectory(i)-1)-_cities_position_y(_trajectory(i+1)-1)),2);  
    }

    lenght += pow((_cities_position_x(_trajectory(_ncities-1)-1)-_cities_position_x(_trajectory(0)-1)),2) + pow((_cities_position_y(_trajectory(_ncities-1)-1)-_cities_position_y(_trajectory(0)-1)),2);

    return lenght; 
}

void trajectory :: pair_permutation(int l, int m){

    std :: swap (_trajectory(l), _trajectory(m)); 

    return; 
}


void trajectory :: shift (int n, int m, int l){
    
    //prepare index 
    uvec indx(_ncities); 
    int k=0; 
    //first part stays the same 
    for (int i=0; i<l; i++){
        indx(k++) = i; 
    }

    //push back 
    for (int i=l+m; i<l+m+n; i++){
        indx(k++) = i; 
    }

    //push forward 
    for (int i=l; i<l+m; i++){
        indx(k++) = i; 
    }

    //last part stays the same 
    if (l+n+m < _ncities){
        for (int i=l+n+m; i<_ncities; i++){
            indx(k++) = i; 
        }
    }

    //sorts according to index 
    _trajectory = _trajectory.elem(indx); 
    return; 
}

void trajectory :: block_permutation(int m, int l, int k){

    if (_ncities < 4){
        return; }


    for (int i=0; i<m; i++){
        swap(_trajectory(k+i), _trajectory(l+i)); 
    }

    return;     
}

void trajectory :: order_inversion(int m, int l){

        
    for (int i=l; i<l+(m-1)/2; i++){
            swap(_trajectory(i), _trajectory(2*l + m -1 -i)); 
        }
    
    return;
}
