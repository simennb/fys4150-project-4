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

string to_scieni(int number, int precision)
{
    //makes number a string on scientific form
    stringstream ss;
    ss << scientific << setprecision(precision) << (double) number;
    return ss.str();
}

string to_scienf(double number, int precision)
{
    //makes number a string on scientific form
    stringstream ss;
    ss << scientific << setprecision(precision) << number;
    return ss.str();
}
