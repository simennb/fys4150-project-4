#ifndef ISING_H
#define ISING_H

#include <string>
#include <fstream>

void InitializeLattice(int L, double **SpinMatrix, double &E, double &M, std::string direction, int my_rank);

void Metropolis(int L, int MCcycles, double T, double *ExpectationValues, char const *dir, std::ofstream &m_file, int my_rank);

void MetropolisD(int L, int MCcycles, double T, double *ExpectationValues, char const *dir, std::ofstream &m_file, int threshold, int my_rank);

void MetropolisE(int L, int MCcycles, double T, double *ExpectationValues, char const *dir, int my_rank);

// Function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
    return (i+limit+add) % (limit);

}

#endif // ISING_H
