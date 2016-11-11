#include "ising.h"
#include "functions.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <random>
#include <armadillo>

//using std::cout; using std::endl; using std::setw; using std::setprecision;
using namespace std;
using namespace arma;

void Ising(int L, int MC_cycles)
{
    // might just be removed, im not sure if we ever will implement something here
    return;
}

void InitializeLattice(int L, mat &SpinMatrix, double &E, double &M, string direction)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // Initialize spin matrix and magnetic moment
    if (direction == "up")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                SpinMatrix(i,j) = 1;
                M += (double) SpinMatrix(i,j);
            }
        }
    }
    else if (direction == "rand")
    {
        for (int i=0; i<L; i++){
            double RD = RandomNumberGenerator(gen);
            for (int j=0; j<L; j++){
                if (RD >= 0.5)
                {
                    SpinMatrix(i,j) = 1;
                }
                else
                {
                    SpinMatrix(i,j) = -1;
                }
                M += (double) SpinMatrix(i,j);
            }
        }
    }
    else if (direction == "down")
    {
        for (int i=0; i<L; i++){
            for (int j=0; j<L; j++){
                SpinMatrix(i,j) = -1;
                M += (double) SpinMatrix(i,j);
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
            E -= (double) SpinMatrix(i,j) *
                    (SpinMatrix(periodic(i,L,-1),j) +
                     SpinMatrix(i,periodic(j,L,-1)));
        }
    }
    return;
}

void Metropolis(int L, int MCcycles, double T, vec &ExpectationValues, char const *dir, string filename)
{
    // Initialize random number generator
    std::random_device rd;  // rd() returns a number, so for parallelization add or subtract rank
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    //Because Valgrind and Eiriks CPU doesn't play along, I use this when running Valgrind
    //srand(100); //Set seed to 100
    //cout << (double) rand()/RAND_MAX << endl; // Random double [0,1]

    // initializing stuff
    mat SpinMatrix = zeros<mat>(L,L);
    int threshold = 300000;
    double E = 0.0; double M = 0.0; vec E_out = zeros<mat>(MCcycles/1000.0 - threshold/1000.0);

    InitializeLattice(L,SpinMatrix,E,M, dir);

    vec EnergyDifference = zeros<mat>(17);

    int r_counter = 0;
    int N_counter = 1000;

    for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/T);

    // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
        // The sweep over the lattice, looping over all spin sites
        for(int x =0; x < L; x++) {
            for (int y= 0; y < L; y++){
                int ix = (int) (RandomNumberGenerator(gen)*(double)L);
                int iy = (int) (RandomNumberGenerator(gen)*(double)L);
                int deltaE =  2*SpinMatrix(ix,iy)*
                                (SpinMatrix( ix, periodic(iy,L,-1) )+
                                 SpinMatrix( periodic(ix,L,-1), iy )+
                                 SpinMatrix( ix, periodic(iy,L, 1) )+
                                 SpinMatrix( periodic(ix,L,1), iy) );
                if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) )
                {
                    SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
                    M += (double) 2*SpinMatrix(ix,iy);
                    E += (double) deltaE;
                    r_counter += 1;
                }
            }
        }
        // Update expectation values
        ExpectationValues(0) += E;
        ExpectationValues(1) += E*E;
        ExpectationValues(2) += M;
        ExpectationValues(3) += M*M;
        ExpectationValues(4) += fabs(M);

        if (cycles == N_counter)
        {
            if (cycles >= threshold)
            {
                E_out[cycles/1000.0 - threshold/1000.0] = E;
            }
         //WriteToFile(L, cycles, T, r_counter, E_out, ExpectationValues, filename);
         N_counter += MCcycles/1000;
        }
    }

    E_handler(E_out, "bla", MCcycles/1000.0 - threshold/1000.0, T, L, dir);
    return;
}

void WriteToFile(int L, int MCcycles, double T, int r_counter, vec E_out, vec &ExpectationValues, string filename)
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

    //cout << E_values << endl;

    for (int i = Emin; i <= Emax; i++)
    {
        E_values[i - Emin] = i;
    }
    cout << Emin << "   " << Emax << endl;
    //cout << E_values << endl;
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
