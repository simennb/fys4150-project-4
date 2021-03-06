#include "ising.h"
#include "functions.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <random>
//#include <armadillo>

using namespace std;
//using namespace arma;

void InitializeLattice(int L, double **SpinMatrix, double &E, double &M, string direction, int my_rank)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd()- my_rank);  // to get
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // Initialize spin matrix and magnetic moment
    if (direction == "up")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                SpinMatrix[i][j] = 1;
                M += (double) SpinMatrix[i][j];
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
                    SpinMatrix[i][j] = 1;
                }
                else
                {
                    SpinMatrix[i][j] = -1;
                }
                M += (double) SpinMatrix[i][j];
            }
        }
    }
    else if (direction == "down")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                SpinMatrix[i][j] = -1;
                M += (double) SpinMatrix[i][j];
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
            E -= (double) SpinMatrix[i][j] *
                    (SpinMatrix[periodic(i,L,-1)][j] +
                     SpinMatrix[i][periodic(j,L,-1)]);
        }
    }
    return;
}

void Metropolis(int L, int MCcycles, double T, double *ExpectationValues, char const *dir, ofstream &m_file, int my_rank)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd()- my_rank);  // unique seed
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    //Because Valgrind and Eiriks CPU doesn't play along, I use this when running Valgrind
    //srand(100); //Set seed to 100
    //cout << (double) rand()/RAND_MAX << endl; // Random double [0,1]

    // initializing stuff
//    mat SpinMatrix = zeros<mat>(L,L);

    double **SpinMatrix = new double *[L];
    for (int i=0; i<L; i++) SpinMatrix[i] = new double [L];

    double E = 0.0; double M = 0.0;

    InitializeLattice(L,SpinMatrix,E,M, dir, my_rank);

//    vec EnergyDifference = zeros<mat>(17);
    double *EnergyDifference = new double[17];
    for (int i=0; i<17; i++) EnergyDifference[i] = 0;

    int r_counter = 0;
    int N_counter = 100;

    for( int de =-8; de <= 8; de+=4) EnergyDifference[de+8] = exp(-de/T);

    // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
        // The sweep over the lattice, looping over all spin sites
        for(int x =0; x < L; x++) {
            for (int y= 0; y < L; y++){
                int ix = (int) (RandomNumberGenerator(gen)*(double)L);
                int iy = (int) (RandomNumberGenerator(gen)*(double)L);
                int deltaE =  2*SpinMatrix[ix][iy]*
                                (SpinMatrix[ ix][ periodic(iy,L,-1) ]+
                                 SpinMatrix[ periodic(ix,L,-1)][ iy ]+
                                 SpinMatrix[ ix][ periodic(iy,L, 1) ]+
                                 SpinMatrix[ periodic(ix,L,1)][ iy] );
                if ( RandomNumberGenerator(gen) <= EnergyDifference[deltaE+8] )
                {
                    SpinMatrix[ix][iy] *= -1.0;  // flip one spin and accept new spin config
                    M += (double) 2*SpinMatrix[ix][iy];
                    E += (double) deltaE;
                    r_counter += 1;
                }
            }
        }
        // Update expectation values
        ExpectationValues[0] += E;
        ExpectationValues[1] += E*E;
        ExpectationValues[2] += M;
        ExpectationValues[3] += M*M;
        ExpectationValues[4] += fabs(M);

        if (cycles == N_counter && my_rank == 0)
        {
            WriteToFile(L, cycles, T, r_counter, ExpectationValues, m_file);
            N_counter += MCcycles/10000;
        }
    }

    // Freeing memory
    for (int i = 0; i < L; i++)delete [] SpinMatrix[i];
    delete [] SpinMatrix;
    delete [] EnergyDifference;

    return;
}
void MetropolisD(int L, int MCcycles, double T, double *ExpectationValues, char const *dir, ofstream &m_file, int threshold, int my_rank)
{
        // Initialize random number generator
        std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
        std::mt19937_64 gen(rd()- my_rank);  // unique seed
        std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

        //Because Valgrind and Eiriks CPU doesn't play along, I use this when running Valgrind
        //srand(100); //Set seed to 100
        //cout << (double) rand()/RAND_MAX << endl; // Random double [0,1]

        // initializing stuff
//        mat SpinMatrix = zeros<mat>(L,L);
        double **SpinMatrix = new double *[L];
        for (int i=0; i<L; i++) SpinMatrix[i] = new double [L];

        double E = 0.0; double M = 0.0;// vec E_out = zeros<mat>(MCcycles/1000.0 - threshold/1000.0);
        int length = (int)(MCcycles/100.0 - threshold/100.0);
        double *E_out = new double [length];

        InitializeLattice(L,SpinMatrix,E,M, dir, my_rank);

//        vec EnergyDifference = zeros<mat>(17);
        double *EnergyDifference = new double[17];
        for (int i=0; i<17; i++) EnergyDifference[i] = 0;

        int r_counter = 0;
        int N_counter = 100;

        for( int de =-8; de <= 8; de+=4) EnergyDifference[de+8] = exp(-de/T);

        // Start Monte Carlo cycles
        for (int cycles = 1; cycles <= MCcycles; cycles++){
            // The sweep over the lattice, looping over all spin sites
            for(int x =0; x < L; x++) {
                for (int y= 0; y < L; y++){
                    int ix = (int) (RandomNumberGenerator(gen)*(double)L);
                    int iy = (int) (RandomNumberGenerator(gen)*(double)L);
                    int deltaE =  2*SpinMatrix[ix][iy]*
                                    (SpinMatrix[ ix][ periodic(iy,L,-1) ]+
                                     SpinMatrix[ periodic(ix,L,-1)][ iy ]+
                                     SpinMatrix[ ix][ periodic(iy,L, 1) ]+
                                     SpinMatrix[ periodic(ix,L,1)][ iy] );
                    if ( RandomNumberGenerator(gen) <= EnergyDifference[deltaE+8] )
                    {
                        SpinMatrix[ix][iy] *= -1.0;  // flip one spin and accept new spin config
                        M += (double) 2*SpinMatrix[ix][iy];
                        E += (double) deltaE;
                        r_counter += 1;
                    }
                }
            }
            // Update expectation values
            ExpectationValues[0] += E;
            ExpectationValues[1] += E*E;
            ExpectationValues[2] += M;
            ExpectationValues[3] += M*M;
            ExpectationValues[4] += fabs(M);

            if (cycles == N_counter && my_rank == 0)
            {
                if (cycles >= threshold)
                {
                    E_out[(int)(cycles/100.0 - threshold/100.0)] = E;

                }
             WriteToFile(L, cycles, T, r_counter, ExpectationValues, m_file);
             N_counter += MCcycles/10000;
            }
        }

    if (my_rank == 0)
    {
        E_handler(E_out, "bla", length, T, L, dir);
    }


    // Freeing memory
    for (int i = 0; i < L; i++)delete [] SpinMatrix[i];
    delete [] SpinMatrix;
    delete [] EnergyDifference;
    delete [] E_out;

    return;
}

void MetropolisE(int L, int MCcycles, double T, double **ExpectationValues, char const *dir, int my_rank, int temp_index)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd() - my_rank);  // unique seed
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    //Because Valgrind and Eiriks CPU doesn't play along, I use this when running Valgrind
    //srand(100); //Set seed to 100
    //cout << (double) rand()/RAND_MAX << endl; // Random double [0,1]

    // initializing stuff
//    mat SpinMatrix = zeros<mat>(L,L);
    double **SpinMatrix = new double *[L];
    for (int i=0; i<L; i++) SpinMatrix[i] = new double [L];

    double E = 0.0; double M = 0.0;

    InitializeLattice(L,SpinMatrix,E,M, dir, my_rank);

//    vec EnergyDifference = zeros<mat>(17);
    double *EnergyDifference = new double[17];
    for (int i=0; i<17; i++) EnergyDifference[i] = 0;

    int r_counter = 0;
    int N_counter = 1000;

    for( int de =-8; de <= 8; de+=4) EnergyDifference[de+8] = exp(-de/T);

    // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
        // The sweep over the lattice, looping over all spin sites
        for(int x =0; x < L; x++) {
            for (int y= 0; y < L; y++){
                int ix = (int) (RandomNumberGenerator(gen)*(double)L);
                int iy = (int) (RandomNumberGenerator(gen)*(double)L);
                int deltaE =  2*SpinMatrix[ix][iy]*
                                (SpinMatrix[ ix][ periodic(iy,L,-1) ]+
                                 SpinMatrix[ periodic(ix,L,-1)][ iy ]+
                                 SpinMatrix[ ix][ periodic(iy,L, 1) ]+
                                 SpinMatrix[ periodic(ix,L,1)][ iy] );
                if ( RandomNumberGenerator(gen) <= EnergyDifference[deltaE+8] )
                {
                    SpinMatrix[ix][iy] *= -1.0;  // flip one spin and accept new spin config
                    M += (double) 2*SpinMatrix[ix][iy];
                    E += (double) deltaE;
                    r_counter += 1;
                }
            }
        }
        // Update expectation values
        ExpectationValues[temp_index][0] += E;
        ExpectationValues[temp_index][1] += E*E;
        ExpectationValues[temp_index][2] += M;
        ExpectationValues[temp_index][3] += M*M;
        ExpectationValues[temp_index][4] += fabs(M);
    }

    // Freeing memory
    for (int i = 0; i < L; i++)delete [] SpinMatrix[i];
    delete [] SpinMatrix;
    delete [] EnergyDifference;

    return;
}
