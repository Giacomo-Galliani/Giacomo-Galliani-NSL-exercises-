#include <cmath>
#include <cstdlib>
#include <string>
#include "GA.h"

using namespace std;
using namespace arma;


int get_position(vec, int); 

int get_position(vec v, int n){
 for (int i = 0; i < v.size();  i++) {
        if (v(i) == n)
            return i;
    }

    std::cerr << "Errore: elemento " << n << " non trovato nella traiettoria." << std::endl;
    exit(EXIT_FAILURE);

}

void GA :: initialize(int sim_type){
  
    _sim_type = sim_type; 

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
    
    
    ifstream config; 
    config.open("config.xy"); 

    if (config.is_open() == true){
        config >> _n_cities; 
        config >> _n_trajectories;  

        _n_elite = static_cast<int>(_n_trajectories/20); //elite members 

        
        _cities_position_x.resize(_n_cities); 
        _cities_position_y.resize(_n_cities); 

        _trajectory.set_size(_n_trajectories); 
    
        //create starting population
        for (int i=0; i<_n_trajectories; i++){
            _trajectory(i).set_ncities(_n_cities); 
        }

        //set cities position 
        for (int i=0; i<_n_cities; i++){
            config >> _cities_position_x(i); 
            config >> _cities_position_y(i); 
        }

        for (int i=0; i<_n_trajectories; i++){
            _trajectory(i).set_x_position(_cities_position_x); 
            _trajectory(i).set_y_position(_cities_position_y);  
        }

    }
    else 
        exit(EXIT_FAILURE); 

    //randomize starting population
    vec one; 
    one.resize(_n_cities); 
    for (int i=0; i<_n_cities; i++){
        one(i) = i+1; 
    }

    for (int i=0; i<_n_trajectories; i++){
        _trajectory(i).set_trajectory(one); 
    }

    for (int i=0; i<_n_trajectories; i++){
        for (int j=0; j<_rnd.Rannyu(1, _n_cities); j++){
            _trajectory(i).pair_permutation(int(_rnd.Rannyu(1, _n_cities)), int(_rnd.Rannyu(1, _n_cities))); 
        }

    }

}

void GA :: crossover (int j, int l){

    int k = int(_rnd.Rannyu(1, _n_cities)); 

    trajectory son, daugther, temp; 

    son.set_ncities(_n_cities); 
    daugther.set_ncities(_n_cities);
    temp.set_ncities(_n_cities);  

    son = _trajectory(j); 
    daugther = _trajectory(l); 

    int appo, count=0;
    

    //ciclo su son
    temp = son; 

    for (int i=0; i<_n_cities; i++){
        appo = get_position(son.get_trajectory(), daugther.get_city(i)); 
        if (appo >k -1){
           
             if (k+count > _n_cities || appo > _n_cities){
                cerr << "indici fuori limiti!"; 
                exit(EXIT_FAILURE); 
             }
            else
            temp.set_component(k+count++, daugther.get_city(i)); 
        }
    }

    temp.check(); 
    son = temp; 
    count =0; 

    //ciclo su daughter
    temp = daugther; 


    for (int i=0; i<_n_cities; i++){
        appo = get_position(daugther.get_trajectory(), son.get_city(i)); 

        if (appo >k-1){
            if (k+count > _n_cities){
                cerr << "indici fuori limiti!"; 
                exit(EXIT_FAILURE); 
            }
            else 
            temp.set_component(k+count++, son.get_city(i)); 
        }
    }

    temp.check(); 
    daugther = temp; 
    count =0; 


    _trajectory(_n_trajectories -1 - int(selection(0, _n_trajectories - _n_elite))) = son;
    _trajectory(_n_trajectories -1 - int(selection(0, _n_trajectories - _n_elite))) = daugther; 

    return; 
}

void GA :: order(){

    for (int i=0; i<_trajectory.n_elem; i++){
        for (int j=0; j<i; j++){
            if (_trajectory(j).get_lenght() > _trajectory(i).get_lenght())
                swap(_trajectory(j), _trajectory(i)); 
        }  
    }
}


int GA :: selection(int min, int max){
    return int((max-min)*pow(_rnd.Rannyu(),8))+ min; 
}


void GA :: generation(){ 

    for (int i=0; i<_n_trajectories/2; i++){
        if (_rnd.Rannyu()> 0.3)
            crossover(selection(0, _n_trajectories), selection(0, _n_trajectories)); 
        }   
    
    _n_steps++; 

    mutation(); 

    order(); 

    check(); 

    return; 
}

void GA :: mutation(){
    for (int i=_n_elite; i<_n_trajectories; i++){
    if (_rnd.Rannyu() > 0.90)
        _trajectory(i).pair_permutation(int(_rnd.Rannyu(1, _n_cities)), int(_rnd.Rannyu(1, _n_cities)));

    if (_rnd.Rannyu() > 0.90){
        int l,m,n; 

        do{
        n = int(_rnd.Rannyu(1, _n_cities-1)); 
        m = int(_rnd.Rannyu(1, _n_cities-2)); 
        l = int(_rnd.Rannyu(2, _n_cities-2)); 
        } while (n+m+l > _n_cities);

        _trajectory(i).shift(m,n,l); 

        }

    if (_rnd.Rannyu() > 0.90){
        int m,l,k; 
        m = int(_rnd.Rannyu(1, _n_cities/2));
        l = int(_rnd.Rannyu(2+ m, _n_cities -m)); 
        k = int(_rnd.Rannyu(2, l));


        _trajectory(i).block_permutation(m,l,k);


        }

    if (_rnd.Rannyu() > 0.90){

        int m,l; 
        
        do {
        m = int(_rnd.Rannyu(1, _n_cities)); //dimension of the block 
        l = int(_rnd.Rannyu(1, _n_cities));  //index of the block 
        } while(m+l>_n_cities || l<0 || m<0);
        
        
        _trajectory(i).order_inversion(m,l);


        }
    }
}


void GA :: check(){

    for (int i=0; i<_n_trajectories; i++){
        _trajectory(i).check(); 
    }
}


double GA :: get_min_lenght(){
    return _trajectory(0).get_lenght(); 
}



int GA :: get_nsteps(){
    return _n_steps; 

}


double GA :: get_lenght(int l){
    return _trajectory(l).get_lenght(); 
}


double GA :: get_half_mean_lenght(){
    double mean =0; 

    for (int i=0; i<_n_trajectories/2; i++){
        mean += _trajectory(i).get_lenght(); 
    }

    return 2*mean/double(_n_trajectories); 

}

trajectory GA :: get_trajectory(int i){
    return _trajectory(i); 
}

void GA :: sort_cities(){
    
    uvec indx(_n_cities); 
    
    for (int i=0; i<_n_cities; i++){
      indx(i) = get_trajectory(0).get_city(i)-1; 
    }
    
    _cities_position_x = _cities_position_x.elem(indx);
    _cities_position_y = _cities_position_y.elem(indx); 


    if (_sim_type ==0){
        ofstream cities("sorted_cities_0.dat"); 
        for (int i=0; i<_n_cities; i++){
            cities << _cities_position_x(i) << " " << _cities_position_y(i) << endl; 
        }
    }
     if (_sim_type ==1){
        ofstream cities("sorted_cities_1.dat"); 
        for (int i=0; i<_n_cities; i++){
            cities << _cities_position_x(i) << " " << _cities_position_y(i) << endl; 
        }
    }

}
