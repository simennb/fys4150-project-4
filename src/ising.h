#ifndef ISING_H
#define ISING_H

#include <armadillo>
#include <string>
#include <fstream>

void Ising(int L, int MC_cycles);

void InitializeLattice(int L, arma::mat &SpinMatrix, double &E, double &M, std::__cxx11::string direction);

void Metropolis(int L, int MCcycles, double T, arma::vec &ExpectationValues, const char *dir, std::__cxx11::string filename);

void WriteToFile(int L, int MCcycles, double T, int r_counter, arma::vec E_out, arma::vec &ExpectationValues, std::__cxx11::string filename);

void E_handler(arma::vec E_out, std::__cxx11::string filename, int length, double T, int L);

// Function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);

}

#endif // ISING_H
