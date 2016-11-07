#include "ising.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <random>
#include <armadillo>

//using std::cout; using std::endl; using std::setw; using std::setprecision;
using namespace std;
using namespace arma;

void Ising(int L, int MC_cycles)
{
    // might just be removed, im not sure if we ever will implement something here
    return;
}

void InitializeLattice(int L, mat &SpinMatrix, double &E, double &M)
{
    // Initialize spin matrix and magnetic moment
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            SpinMatrix(i,j) = 1;
            M += (double) SpinMatrix(i,j);
        }
    }
    // Initialize E
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            E -= (double) SpinMatrix(i,j) *
                    (SpinMatrix(periodic(i,L,-1),j) +
                     SpinMatrix(i,periodic(j,L,-1)));
        }
    }
    return;
}

void Metropolis(int L, int MCcycles, double T, vec &ExpectationValues)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // initializing stuff
    mat SpinMatrix = zeros<mat>(L,L);
    double E = 0.0; double M = 0.0;
    InitializeLattice(L,SpinMatrix,E,M);

    vec EnergyDifference = zeros<mat>(17);

    for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/T);
    // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
        // The sweep over the lattice, looping over all spin sites
        for(int x =0; x < L; x++) {
            for (int y= 0; y < L; y++){
                int ix = (int) (RandomNumberGenerator(gen)*(double)L);
                int iy = (int) (RandomNumberGenerator(gen)*(double)L);
                int deltaE =  2*SpinMatrix(ix,iy)*
                                (SpinMatrix( ix, periodic(iy,L,-1) )+
                                 SpinMatrix( periodic(ix,L,-1), iy )+
                                 SpinMatrix( ix, periodic(iy,L, 1) )+
                                 SpinMatrix( periodic(ix,L,1), iy) );
                if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
                    SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
                    M += (double) 2*SpinMatrix(ix,iy);
                    E += (double) deltaE;
                }
            }
        }
        // Update expectation values
        ExpectationValues(0) += E;
        ExpectationValues(1) += E*E;
        ExpectationValues(2) += M;
        ExpectationValues(3) += M*M;
        ExpectationValues(4) += fabs(M);
    }

    return;
}

void WriteToFile(int L, int MCcycles, double T, vec &ExpectationValues)
{
    // Currently only prints out
    double norm = 1.0/((double)MCcycles);
    double E_exp    = ExpectationValues(0)*norm;
    double E2_exp   = ExpectationValues(1)*norm;
    double M_exp    = ExpectationValues(2)*norm;
    double M2_exp   = ExpectationValues(3)*norm;
    double Mabs_exp = ExpectationValues(4)*norm;

    double E_variance = (E2_exp - E_exp*E_exp)/L/L;
    double M_variance = (M2_exp - Mabs_exp*Mabs_exp)/L/L;

    cout<<"T = "<<setprecision(8)<<T<<endl;
    cout<<"< E >   = "<<setprecision(8)<<E_exp<<endl;
    cout<<"<|M|>   = "<<setprecision(8)<<Mabs_exp<<endl;
    cout<<"sigma_E = "<<setprecision(8)<<E_variance<<endl;
    cout<<"sigma_M = "<<setprecision(8)<<M_variance<<endl;

    return;
}
