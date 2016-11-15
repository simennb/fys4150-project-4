#include "ising.h"
#include "functions.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <random>
#include <sstream>
//#include <armadillo>

using namespace std;
//using namespace arma;

int max_value(double *array, int n)
{
    double max_curr = array[0];
    for (int i=0; i<n; i++){
        if (array[i] > max_curr) max_curr = array[i];
    }
    return (int)max_curr;
}

int min_value(double *array, int n)
{
    double min_curr = array[0];
    for (int i=0; i<n; i++){
        if (array[i] < min_curr) min_curr = array[i];
    }
    return (int)min_curr;
}

string time_to_print(double time)
{
    stringstream ss;
    int hours = ((int)time)/3600;
    time -= hours*60.0;
    int minutes = ((int)time)/60;
    double seconds = time - minutes*60.0;
    ss <<hours<<":"<<minutes<<":"<<setprecision(3)<<seconds;
    return ss.str();
}

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

void InitializeFile(string filename, ofstream &m_file)
{
    m_file.open(filename.c_str(), ofstream::out);
    if (!m_file.good()) {
    cout << "Error opening file " << filename << ". Aborting!" << endl;
    exit(1);
    }
    m_file<<" MCcycles   "<<"T        "<<"<E>         "<<"<|M|>       "<<"sigma_E      "<<"sigma_M     "<<"accepted"<<endl;
}

void WriteToFile(int L, int MCcycles, double T, int r_counter, double *ExpectationValues, ofstream &m_file)
{
    double norm = 1.0/((double)MCcycles);
    double E_exp    = ExpectationValues[0]*norm;
    double E2_exp   = ExpectationValues[1]*norm;
    double M_exp    = ExpectationValues[2]*norm;
    double M2_exp   = ExpectationValues[3]*norm;
    double Mabs_exp = ExpectationValues[4]*norm;

    double E_variance = (E2_exp - E_exp*E_exp)/L/L;
    double M_variance = (M2_exp - Mabs_exp*Mabs_exp)/L/L;

    // Writing to file
    m_file<<setw(8)<<MCcycles;
    m_file<<setw(13)<<setprecision(8)<<T;
    m_file<<setw(13)<<setprecision(8)<<E_exp;
    m_file<<setw(13)<<setprecision(8)<<Mabs_exp;
    m_file<<setw(13)<<setprecision(8)<<E_variance;
    m_file<<setw(13)<<setprecision(8)<<M_variance;
    m_file<<setw(13)<<setprecision(8)<<r_counter<<endl;
    return;
}

void E_handler(double *E_out, string filename, int length, double T, int L, string stringdir)
{
    int small_length = max_value(E_out,length) - min_value(E_out,length);
    int Emin = min_value(E_out,length); int Emax = max_value(E_out,length);
//    vec E_values = zeros<mat>(small_length + 1);
//    vec E_counter = zeros<mat>(small_length + 1);
    double *E_values = new double[small_length + 1];
    double *E_counter = new double[small_length + 1];

    for (int i = Emin; i <= Emax; i++)
    {
        E_values[i - Emin] = i;
    }
    cout << Emin << "   " << Emax << endl;

    for (int i = 0; i <= length; i++)
    {
        for (int j = i+1; j <= length; j++)  // start at i+1 to not count things twice
        {
            for (int k = 0; k <= small_length; k++)
            {
                if (E_out[i] == E_out[j] && E_out[i] == E_values[k])
                {
                    E_counter[k] += 1;
                }
            }
        }
    }

    ofstream E_file;
    E_file.open("../benchmarks/task_d/EnergyValues_T" + to_fixf(T,1) + "_L" + to_fixi(L,1) + "_dir_" + stringdir +".txt");
    E_file<< "E_counter = " << endl;
    E_file<< E_counter;
    E_file<<"E_values = " <<endl;
    E_file<<E_values;
    E_file.close();
    return;
}
