#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include "SA.h"

using namespace std;


double psi(double mu, double sigma, double x){
   return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma)); 
}

double psi_2(double mu, double sigma, double x){
   return (exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma)))*(exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma)));
}


double KE(double mu, double sigma, double x){
    return -1/(2*sigma*sigma)*(((x-mu)/(sigma)*(x-mu)/(sigma)*exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + (x+mu)/(sigma)*(x+mu)/(sigma)*exp(-(x+mu)*(x+mu)/(2*sigma*sigma)))/(exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2*sigma*sigma))) -1); 
}



double V(double x){
    return x*x*x*x - (5./2.)*x*x; 
}


void SA :: initialize(){

    //initiali<e _rnd
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


    //set starting parameters, nsteps 
    _nblocks = 100; 
    _nsteps = 100000; 
    _block_length = int(_nsteps/_nblocks); 
    _mu = 2*_rnd.Rannyu(); 
    _sigma = _rnd.Rannyu(); 
    _beta = 0.1; 
    _step = 2.0; 

    //computes starting hamiltonian
    _H_current = get_hamiltonian(_mu, _sigma); 
}



double SA :: get_ratio(){
    return _naccepted/double(_ntrial); 
}



double SA :: move(double x, double mu, double sigma){

    _ntrial++; 
    double xnew = x + _rnd.Rannyu(-3*_sigma, +3*_sigma); 
    if (metro(x, xnew, mu, sigma)){
        _naccepted++; 
        return xnew;
    }
    else 
        return x; 
    
}

bool SA :: metro(double x, double xnew,  double mu, double sigma){ // Metropolis algorithm
    bool decision = false;
    double acceptance;


    acceptance = psi_2(mu, sigma, xnew)/psi_2(mu, sigma, x); 

    if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
    return decision;
}

vector<double>  SA :: get_hamiltonian(double mu, double sigma){

    double H=0; 
    double err =0; 
    vector<double> energy; 

    //prepare accumulators and starting position
    double sum[_nblocks];
    double x = _rnd.Rannyu(mu-sigma, mu+sigma); 

    for (int i=0; i<_nblocks; i++){
        sum[i] =0; 
    }

    //fill accumulators 
    for (int i=0; i<_nblocks; i++){
        for (int j=0; j<_block_length; j++){
            x = move(x, mu, sigma);  
            sum[i] += V(x) + KE(mu, sigma, x); 
        }
        sum[i] /= double(_block_length); 
    }

    //compute mean
    for (int i=0; i<_nblocks; i++){
        H += sum[i]; 
    }
        
    H /= double(_nblocks); 

    //compute statistical error 
    for (int i=0; i<_nblocks; i++){
        err += sum[i]*sum[i]; 
    }

    err /= double(_nblocks); 
    err -= H*H; 

    err = sqrt(err/double(_nblocks-1)); 

    energy.push_back(H); 
    energy.push_back(err); 

    return energy; 
}

void SA :: equilibrate(){

    for (int i=0; i<6000; i++){
        move_parameters(); 
    }
    return; 
}


double SA :: get_sigma(){
    return _sigma; 
}

double SA :: get_mu(){
    return _mu; 
}

void SA :: set_beta(double beta){
    _beta = beta; 
    return; 
}

double SA :: get_beta(){
    return _beta; 
}

void SA :: set_mu(double mu){
    _mu = mu; 
    return; 
}

void SA :: set_sigma(double sigma){
    _sigma = sigma; 
    return; 
}


double SA :: get_ratio2(){
    return _naccepted2/double(_ntrial2); 
}


double SA :: get_energy(){
    return _H_current[0]; 
}


void SA :: update_energy(double mu_new, double sigma_new){

    _mu = mu_new; 
    _sigma = sigma_new; 
    _H_current = get_hamiltonian(_mu, _sigma); 
    return; 
}

void SA :: change_energy(double mu_new, double sigma_new){
    _mu_new = mu_new; 
    _sigma_new = sigma_new; 

    _H_new = get_hamiltonian(_mu_new, _sigma_new);
    return; 
}

double SA :: get_energy_new(){
    return _H_new[0]; 
}

void SA :: move_parameters(){


    double sigma_new, mu_new; 
    _ntrial2++; 

    //generate new parameters 
    sigma_new = _sigma + _rnd.Rannyu(-_step/2.,_step/2.); 
    sigma_new = abs(sigma_new); 

    mu_new = _mu + _rnd.Rannyu(-_step,_step);

    //compute new parameters energy 
    change_energy(mu_new, sigma_new); 

    if(metro2()){
        update_energy(_mu_new, _sigma_new); 
        _naccepted2++; 
    }

    return; 
}

bool SA :: metro2(){

    bool decision = false; 

    //compute acceptance 
    double delta_E = get_energy_new() - get_energy(); 
    double acceptance = exp(-_beta*delta_E); 


    double a = _rnd.Rannyu(); 

    
    if (a < acceptance) { //metropolis acceptance step 
        decision = true; 
     
    }
 
    return decision; 
}

void SA :: reset_ratio2(){
    _ntrial2 =0; 
    _naccepted2 =0; 
}

double SA :: get_error(){
    return _H_current[1]; 
}

void SA :: set_step(int l){

    _step *= 1./3.; 
    return; 
}

double SA :: get_step(){
    return _step; 
}


vector <double> SA :: get_progressive_energy(){

    //prepare accumulators and starting position 
    vector<double> blocks, results, err; 
    for (int i=0; i<_nblocks; i++){
        results.push_back(0);  
        blocks.push_back(0); 
        err.push_back(0); 
    }

    double x = _rnd.Rannyu(_mu - _sigma, _mu + _sigma); 

    //fill accumulators 
    for (int i=0; i<_nblocks; i++){
        for (int j=0; j<_block_length; j++){
            x = move(x, _mu, _sigma); 
            blocks[i] += V(x) + KE(_mu, _sigma, x); 
        }
        blocks[i] /= double(_block_length); 
    }
 
    //compute progressive mean and error 
    for (int i=0; i<_nblocks; i++){
        for (int j=0; j<i+1; j++){
            results[i] += blocks[j]; 
            err[i] += blocks[j]*blocks[j]; 
        }
        results[i] /= double(i+1); 
        err[i] /= double(i+1); 
    }

    err[0] =0; 
    for (int i=1; i<_nblocks; i++){
        err[i] -= results[i]*results[i]; 
        err[i] = sqrt(err[i]/double(i)); 
    }

    for (int i=0; i<_nblocks; i++){
        results.push_back(err[i]); 
    }

    return results; 
}