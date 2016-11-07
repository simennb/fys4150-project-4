#ifndef ISING_H
#define ISING_H

#include <armadillo>

void Ising(int L, int MC_cycles);

void InitializeLattice(int L, arma::mat &SpinMatrix, double &E, double &M);

void Metropolis(int L, int MCcycles, double T, arma::vec &ExpectationValues);

void WriteToFile(int L, int MCcycles, double T, arma::vec &ExpectationValues);

// Function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

#endif // ISING_H
