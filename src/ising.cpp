#include "ising.h"
#include "functions.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <random>
#include <armadillo>

using namespace std;
using namespace arma;

void InitializeLattice(int L, mat &SpinMatrix, double &E, double &M, string direction)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // Initialize spin matrix and magnetic moment
    if (direction == "up")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                SpinMatrix(i,j) = 1;
                M += (double) SpinMatrix(i,j);
            }
        }
    }
    else if (direction == "rand")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                double RD = RandomNumberGenerator(gen);
                if (RD >= 0.5)
                {
                    SpinMatrix(i,j) = 1;
                }
                else
                {
                    SpinMatrix(i,j) = -1;
                }
                M += (double) SpinMatrix(i,j);
            }
        }
    }
    else if (direction == "down")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                SpinMatrix(i,j) = -1;
                M += (double) SpinMatrix(i,j);
            }
        }
    }
    else
    {
       cout << "Direction must be set to either up, down or rand." << endl;
       return;
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

void Metropolis(int L, int MCcycles, double T, vec &ExpectationValues, char const *dir, string filename)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    //Because Valgrind and Eiriks CPU doesn't play along, I use this when running Valgrind
    //srand(100); //Set seed to 100
    //cout << (double) rand()/RAND_MAX << endl; // Random double [0,1]

    // initializing stuff
    mat SpinMatrix = zeros<mat>(L,L);
    double E = 0.0; double M = 0.0;

    InitializeLattice(L,SpinMatrix,E,M, dir);

    vec EnergyDifference = zeros<mat>(17);

    int r_counter = 0;
    int N_counter = 1000;

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
                if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) )
                {
                    SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
                    M += (double) 2*SpinMatrix(ix,iy);
                    E += (double) deltaE;
                    r_counter += 1;
                }
            }
        }
        // Update expectation values
        ExpectationValues(0) += E;
        ExpectationValues(1) += E*E;
        ExpectationValues(2) += M;
        ExpectationValues(3) += M*M;
        ExpectationValues(4) += fabs(M);

        if (cycles == N_counter)
        {
         WriteToFile(L, cycles, T, r_counter, ExpectationValues, filename);
         N_counter += MCcycles/1000;
        }
    }
    return;
}
void MetropolisD(int L, int MCcycles, double T, vec &ExpectationValues, char const *dir, string filename, int threshold)
{
        // Initialize random number generator
        std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
        std::mt19937_64 gen(rd());
        std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

        //Because Valgrind and Eiriks CPU doesn't play along, I use this when running Valgrind
        //srand(100); //Set seed to 100
        //cout << (double) rand()/RAND_MAX << endl; // Random double [0,1]

        // initializing stuff
        mat SpinMatrix = zeros<mat>(L,L);
        double E = 0.0; double M = 0.0; vec E_out = zeros<mat>(MCcycles/1000.0 - threshold/1000.0);

        InitializeLattice(L,SpinMatrix,E,M, dir);

        vec EnergyDifference = zeros<mat>(17);

        int r_counter = 0;
        int N_counter = 1000;

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
                    if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) )
                    {
                        SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
                        M += (double) 2*SpinMatrix(ix,iy);
                        E += (double) deltaE;
                        r_counter += 1;
                    }
                }
            }
            // Update expectation values
            ExpectationValues(0) += E;
            ExpectationValues(1) += E*E;
            ExpectationValues(2) += M;
            ExpectationValues(3) += M*M;
            ExpectationValues(4) += fabs(M);

            if (cycles == N_counter)
            {
                if (cycles >= threshold)
                {
                    E_out[cycles/1000.0 - threshold/1000.0] = E;
                }
             WriteToFile(L, cycles, T, r_counter, ExpectationValues, filename);
             N_counter += MCcycles/1000;
            }
        }

    E_handler(E_out, "bla", MCcycles/1000.0 - threshold/1000.0, T, L, dir);
    return;
}
