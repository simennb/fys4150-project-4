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
    //makes int number a string on scientific form
    stringstream ss;
    ss << scientific << setprecision(precision) << (double) number;
    return ss.str();
}

string to_scienf(double number, int precision)
{
    //makes float number a string on scientific form
    stringstream ss;
    ss << scientific << setprecision(precision) << number;
    return ss.str();
}

string to_fixi(int number, int precision)
{
    //makes int number a string on a fixed form
    stringstream ss;
    ss << fixed << setprecision(precision) << (double) number;
    return ss.str();
}

string to_fixf(double number, int precision)
{
    //makes int number a string on a fixed form
    stringstream ss;
    ss << fixed << setprecision(precision) << number;
    return ss.str();
}
