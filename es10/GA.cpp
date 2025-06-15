#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include "GA.h"
#include "mpi.h"

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

void GA :: initialize(){
  
   
    int rank, size; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
   
    //initialize _rnd 
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for (int i=0; i<10*rank; i++){ //sets differents seeds for different populations 
            Primes >> p1; 
        }
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
        
    //reads data from file   
    ifstream config; 
    config.open("cap_prov_ita.dat"); 

    if (config.is_open() == true){
        _n_cities =110; 
        _n_trajectories = 4*_n_cities;  

        _n_elite = 0;  

        _cities_position_x.resize(_n_cities); 
        _cities_position_y.resize(_n_cities); 

        _trajectory.set_size(_n_trajectories); 
    
        for (int i=0; i<_n_trajectories; i++){
            _trajectory(i).set_ncities(_n_cities); 
            _trajectory(i).initialize_rnd(); 
        }

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

    //prepare starting population
    vec one; 
    one.resize(_n_cities); 
    for (int i=0; i<_n_cities; i++){
        one(i) = i+1; 
    }

    for (int i=0; i<_n_trajectories; i++){
        _trajectory(i).set_trajectory(one); 
    }

    for (int i=0; i<_n_trajectories; i++){
        for (int j=0; j< _n_cities; j++){
            _trajectory(i).pair_permutation(int(_rnd.Rannyu(1, _n_cities)), int(_rnd.Rannyu(1, _n_cities))); 
        }

    }
}

void GA :: crossover (int j, int l){


    int k = int(_rnd.Rannyu(1, _n_cities)); 

    //prepare new individuals 
    trajectory son, daugther, temp; 

    son.set_ncities(_n_cities); 
    daugther.set_ncities(_n_cities);
    temp.set_ncities(_n_cities);  

    son = _trajectory(j); 
    daugther = _trajectory(l); 

    int appo, count=0;
    

    //son loop
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

    //daughter loop
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


    _trajectory(_n_trajectories -1 - int(selection(0, _n_trajectories - _n_elite))) = son; //replaces random individual from the worsts 
    _trajectory(_n_trajectories -1 - int(selection(0, _n_trajectories - _n_elite))) = daugther; //replaces random individual from the worsts 

    return; 
}

void GA :: order(){

    
    for (int i=0; i<_n_trajectories; i++){
        for (int j=0; j<i; j++){
            if (_trajectory(j).get_lenght() > _trajectory(i).get_lenght())
                swap(_trajectory(j), _trajectory(i)); 
        }  
    }
}


int GA :: selection(int min, int max){
    return int((max-min)*pow(_rnd.Rannyu(),4))+ min; 
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
    if (_rnd.Rannyu() > 0.90){
        _trajectory(i).pair_permutation(int(_rnd.Rannyu(1, _n_cities)), int(_rnd.Rannyu(1, _n_cities)));
    }

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


void GA :: paste_min_lenght(int sim_type){


    int size, rank; 

    double lengths = _trajectory(0).get_lenght(); //element to reduce (best length for each population)
    double length; 

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 

    MPI_Reduce(&lengths, &length, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); //give lenght the best length found 

    //print to file 
    if (rank==0){
        if (sim_type==0){
        ofstream results; 
        results.open("lenght.dat", ios:: app); 
        results << _n_steps << " " << length << endl; 
        }
        else{
            ofstream results; 
            results.open("lenght_comp.dat", ios:: app); 
            results << _n_steps << " " << length << endl; 
        }

        
    }
    
    return; 
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

trajectory& GA :: get_trajectory(int i){
    return _trajectory(i); 
}


void GA :: sort_cities(){
    

    uvec indx(_n_cities); //index of the order of cities 


    double lengths[2];  
    double length[2]; 

    lengths[0] = get_trajectory(0).get_lenght(); 
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    lengths[1] = rank; 

    MPI_Reduce(&lengths, &length, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, 0 , MPI_COMM_WORLD); //find witch generation has the best lenght 


    if(rank==0){
        
        swap_best(0, int(length[1])); //gives to 0 the best trajectory 

        //fill index 
        for (int i=0; i<_n_cities; i++){
            indx(i) = get_trajectory(0).get_city(i)-1; 
        }
        
        //sort cities 
        _cities_position_x = _cities_position_x.elem(indx);
        _cities_position_y = _cities_position_y.elem(indx); 

        //print results 
        ofstream cities("sorted_cities.dat"); 
            for (int i=0; i<_n_cities; i++){
                cities << _cities_position_x(i) << " " << _cities_position_y(i) << endl; 
            }
    }
}


void GA :: migration(){

    int size, rank; 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 


    //create a "trajectory" a of the indexes of the swapping 
    trajectory a; 
    a.set_ncities(size); 
    for (int i=0; i<size; i++){
        a.set_component(i, i); 
    }

    if (rank==0){
        for (int i=0; i<size; i++){
        a.pair_permutation(int(_rnd.Rannyu(0,size)),int(_rnd.Rannyu(0,size))); //use the pair_permutation function to create a random swapping that involves all populations 
        }
    }

    int b[size]; 
    for (int i=0; i<size; i++){
        b[i] = a.get_city(i); 
    }

    MPI_Bcast(b, size, MPI_INT, 0, MPI_COMM_WORLD); //broadcast so all processes know who should be swapped 

    //swap best of the selected pairs 
    for (int i=0; i<size/2; i++){ 
        if (rank == b[i] || rank==b[size/2 +i]) swap_best(b[i], b[size/2 +i]); 
    }
}

void GA :: swap_best(int l, int j){

     
    MPI_Status stat1, stat2;  

    if (l==j) //no swap needed 
        return; 
    
    int rank; 
    int tag1=1000 + l*100 + j, tag2=2000 + l*100 +j; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    int x[_n_cities], y[_n_cities]; 


    if (rank == l){
        for (int i=0; i<_n_cities; i++){
            x[i] = get_trajectory(0).get_city(i); //converts l best trajectory to double 
        }
        }

    if (rank == j){
        for (int i=0; i<_n_cities; i++){
            y[i] = get_trajectory(0).get_city(i); //converts j best trajecctory to double 
        }
        }


    //sends the message as double array 
    if (rank ==l){
        MPI_Send(&x[0], _n_cities, MPI_INT, j, tag1, MPI_COMM_WORLD); 
        MPI_Recv(&y[0], _n_cities, MPI_INT, j, tag2, MPI_COMM_WORLD, &stat2); 
    }

    if (rank ==j){ 
        MPI_Recv(&x[0], _n_cities, MPI_INT, l, tag1, MPI_COMM_WORLD, &stat1); 
        MPI_Send(&y[0], _n_cities, MPI_INT, l, tag2, MPI_COMM_WORLD);
    }


    //converts back to trajectory after the message is received 
    if (rank ==l){
        for (int i=0; i<_n_cities; i++){
            get_trajectory(0).set_component(i, y[i]); 
        }
    }

    if (rank ==j){
        for (int i=0; i<_n_cities; i++){
            get_trajectory(0).set_component(i,x[i]); 
        }
    }
}