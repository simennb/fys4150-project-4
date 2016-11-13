#ifndef FUNCTIONS
#define FUNCTIONS

#include<fstream>

std::string to_scieni(int number, int precision);
std::string to_scienf(double number, int precision);
std::string to_fixi(int number, int precision);
std::string to_fixf(double number, int precision);

void InitializeFile(std::string filename,std::ofstream &m_file);

void WriteToFile(int L, int MCcycles, double T, int r_counter, arma::vec &ExpectationValues, std::ofstream &m_file);

void E_handler(arma::vec E_out, std::string filename, int length, double T, int L, std::string stringdir);

#endif // FUNCTIONS

