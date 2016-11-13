#ifndef ISING_H
#define ISING_H

#include <armadillo>
#include <string>
#include <fstream>

void InitializeLattice(int L, arma::mat &SpinMatrix, double &E, double &M, std::string direction);

void Metropolis(int L, int MCcycles, double T, arma::vec &ExpectationValues, const char *dir, std::string filename);

// Function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);

}

#endif // ISING_H
