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


void WriteToFile(int L, int MCcycles, double T, int r_counter, vec &ExpectationValues, string filename)
{
    // Currently print outs and write to file
    double norm = 1.0/((double)MCcycles);
    double E_exp    = ExpectationValues(0)*norm;
    double E2_exp   = ExpectationValues(1)*norm;
    double M_exp    = ExpectationValues(2)*norm;
    double M2_exp   = ExpectationValues(3)*norm;
    double Mabs_exp = ExpectationValues(4)*norm;

    double E_variance = (E2_exp - E_exp*E_exp)/L/L;
    double M_variance = (M2_exp - Mabs_exp*Mabs_exp)/L/L;

    /*cout<<"MCcycles = "<<setprecision(8)<<MCcycles<<endl;
    cout<<"T = "<<setprecision(8)<<T<<endl;
    cout<<"< E >   = "<<setprecision(8)<<E_exp<<endl;
    cout<<"<|M|>   = "<<setprecision(8)<<Mabs_exp<<endl;
    cout<<"sigma_E = "<<setprecision(8)<<E_variance<<endl;
    cout<<"sigma_M = "<<setprecision(8)<<M_variance<<endl;
    cout<<"r_counter = "<<setprecision(8)<<r_counter<<endl;

    ofstream m_file;
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }
    */

    /*
    ofstream m_file;
    m_file.open(filename, ios::app);
    //ofile<<setiosflags(ios::showpoint | ios::uppercase);
    m_file<<"MCcycles = "<<setprecision(8)<<MCcycles<<endl;
    m_file<<"T = "<<setprecision(8)<<T<<endl;
    m_file<<"<E> = "<<setprecision(8)<<E_exp<<endl;
    m_file<<"<|M|> = "<<setprecision(8)<<Mabs_exp<<endl;
    m_file<<"sigma_E = "<<setprecision(8)<<E_variance<<endl;
    m_file<<"sigma_M = "<<setprecision(8)<<M_variance<<endl;
    m_file<<"r_counter = "<<setprecision(8)<<r_counter<<endl;
    //m_file<<"E_out = "<<setprecision(8)<<E_out<<endl;
    m_file << "  " << endl;
    */
    return;
}

void E_handler(vec E_out, string filename, int length, double T, int L, string stringdir)
{
    int small_length = E_out.max() - E_out.min();
    int Emin = E_out.min(); int Emax = E_out.max();
    vec E_values = zeros<mat>(small_length + 1);
    vec E_counter = zeros<mat>(small_length + 1);

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
