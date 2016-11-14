#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
//#include <armadillo>
#include <time.h>
#include <mpi.h>
#include "ising.h"
#include "functions.h"

using namespace std;
//using namespace arma;

/* TODO:
 *  - fix command line arguments DONE
 *    - L, MC_cycles and task at the very least need to be read in DONE
 *    - might change a bit when parallelization is added, but thats a problem when we fix that ^^
 *  - structure for the different tasks, and actually what the tasks ask us to do     DONE
 *    - 2x InitializeLattice functions, one for only up spin start state, and one with random start state   DONE
 *  - parallelization with MPI!!!
 *  - fix so that WriteToFile actually writes to file   DONE
 *  - see if we change variable names a bit to make it less similar to Morten's example
 *  - compiler flags?
 *  - should time to see if we get improvements with parallelization (for a few runs at least), might not need to write to file
 *    just note it down somewhere before we add parallelization, and run with same parameters afterwards.
 *  - number of accepted configurations as func. of MC cycles, add counter to that in Metropolis(...), and write to file
 *
 *
 *                             !         !
 *                            ! !       ! !
 *                           ! . !     ! . !
 *                             ^^^^^^^^^^^
 *                           ^             ^
 *                         ^  (0)       (0)  ^
 *                        ^        ""         ^
 *                       ^   ***************    ^
 *                     ^   *                 *   ^
 *                    ^   *   /\   /\   /\    *    ^
 *                   ^   *                     *    ^           _____
 *                  ^   *   /\   /\   /\   /\   *    ^         /     \
 *                 ^   *                         *    ^        vvvvvvv  /|__/|
 *                 ^  *                           *   ^           I   /O,O   |
 *                 ^  *                           *   ^           I /_____   |      /|/|
 *                  ^ *                           *  ^           J|/^ ^ ^ \  |    /00  |    _//|
 *                   ^*                           * ^             |^ ^ ^ ^ |W|   |/^^\ |   /oo |
 *                    ^ *                        * ^               \m___m__|_|    \m_m_|   \mm_|
 *                    ^  *                      *  ^
 *                      ^  *       ) (         * ^
 *                          ^^^^^^^^ ^^^^^^^^^
 */


int main(int argc, char *argv[])
{
    if (argc < 5 || argc > 7)
    {
        cout << "Usage: " << argv[0] << ", task, L, Monte Carlo Cycles, T and spin direction as command line argument" << endl;
        cout << "Task can be set to b, c, d, e" << endl;
        cout << "L is the dimension of the spin lattice" << endl;
        cout << "Monte Carlo Cycles is the number of Monte Carlo Cycles" << endl;
        cout << "T is temperature" << endl;
        cout << "Spin configuration is optional, set to up by default." << endl;
        cout << "Spin configuration can either be set to up, down or rand to set all spins to be either up, down or random" << endl;
        exit(1);
    }
    // Initialize MPI
    int numprocs, my_rank;
    double time_start, time_end, total_time;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();

    // Initialize some values
    int L = atoi(argv[2]);
    int MC_cycles = atoi(argv[3]);  //1000000;
    double T = atof(argv[4]);

    char const *stringdir = "up";   //Standard direction

    if (argc == 6)
    {
        if (!(strcmp(argv[5],"up")==0 || strcmp(argv[5],"down")==0 || strcmp(argv[5], "rand") == 0))
        {
            cout<<"Spin direction not specified correctly, use 'up', 'down' or 'rand'. Set to 'up' by default."<<endl;
        }
        else
        {
            stringdir = argv[5];
        }
    }


    if (strcmp(argv[1], "b") == 0)
    {
        // Creating filename and initializing file
        ofstream m_file;
        if (my_rank == 0){
            string filename = "../benchmarks/task_b/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T, 1) +".xyz";
            InitializeFile(filename, m_file);
        }
        //vec totExpectationValues = zeros<mat>(5);
        //vec ExpectationValues = zeros<mat>(5);
        double *totExpectationValues = new double[5];
        double *ExpectationValues = new double[5];

        Metropolis(L, MC_cycles, T, ExpectationValues, stringdir, m_file, my_rank);

        for (int i=0; i<5; i++){
            MPI_Reduce(&ExpectationValues[i],&totExpectationValues[i],1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        time_end = MPI_Wtime();
        total_time = time_end - time_start;

        if (my_rank == 0){
            cout<<"Hey, my rank is "<<my_rank<<endl;
            cout<<"MC cycles = "<<MC_cycles<<endl;
            cout<<"local_E = "<<ExpectationValues[0]/MC_cycles<<endl;
            cout<<"total_E = "<<totExpectationValues[0]/MC_cycles/numprocs<<endl;
            cout<<"on number of processors = "<<numprocs<<endl;
            cout<<"time = "<<total_time<<endl;
        }
        //WriteToFile(L, MC_cycles, T, ExpectationValues, filename);

    }

    else if (strcmp(argv[1], "c") == 0)
    {
        //Creating filename and initializing file
        string filename = "../benchmarks/task_c/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T, 1) +".xyz";
        ofstream m_file;
        InitializeFile(filename, m_file);

        // Matrix to fill with expectation values
        double *totExpectationValues = new double[5];
        double *ExpectationValues = new double[5];
//        vec ExpectationValues = zeros<mat>(5);

        Metropolis(L, MC_cycles, T, ExpectationValues, stringdir, m_file, my_rank);

        //WriteToFile(L, MC_cycles, T, ExpectationValues, filename);

        m_file.close();
    }

    else if (strcmp(argv[1], "d") == 0)
    {
        //Creating filename and initializing file
        string filename = "../benchmarks/task_d/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T, 1) +".xyz";
        ofstream m_file;
        InitializeFile(filename, m_file);

        // Matrix to fill with expectation values
//        vec ExpectationValues = zeros<mat>(5);
        double *totExpectationValues = new double[5];
        double *ExpectationValues = new double[5];

        int thresholdT24 = 300000;
        int thresholdT1 = 100000;
        int threshold = 0;

        if (L != 20)
        {
            cout << "L must be 20. The program will now stop." << endl;
            exit(1);
        }

        if (T == 2.4 && MC_cycles <= thresholdT24 + 1000)
        {
            cout << "MCCycles must be bigger than " << (thresholdT24 + (int) 1000) << " when T = 2.4. The program will now stop." << endl;
            exit(1);
        }
        else if (T == 2.4 && MC_cycles > thresholdT24+1000)
        {
            threshold = thresholdT24;
        }
        else if (T == 1.0 && MC_cycles <= thresholdT1)
        {
            cout << "MCCycles must be bigger than " << (thresholdT1+1000) << " when T = 1.0. The program will now stop." << endl;
            exit(1);
        }
        else
        {
            threshold = thresholdT1;
        }

        MetropolisD(L, MC_cycles, T, ExpectationValues, stringdir, m_file, threshold, my_rank);

        m_file.close();
    }

    else if (strcmp(argv[1], "e") == 0)
    {
        int N = T;
        double Tstart = 2.0;
        double Tstop = 2.3;

        //vec Temp = linspace<vec>(Tstart, Tstop, N);
        // I'll make my own linspace, with mexican hookers. (wait what?)
        double *T = new double[N];
        double dt = (Tstop-Tstart)/(N-1);
        for (int i=0; i<N; i++){
            T[i] = Tstart + i*dt;
        }
        //double dt2 = Temp[1] - Temp[0];
        //cout<<dt<<" "<<dt2<<endl;

//      vec ExpectationValues = zeros<mat>(5);
        double **totExpectationValues = new double *[N];
        double **ExpectationValues = new double *[N];
        for (int i=0; i<N; i++){
            totExpectationValues[i] = new double [5];
            ExpectationValues[i] = new double [5];
            for (int j=0; j<5; j++){
                totExpectationValues[i][j] = 0.0;
                ExpectationValues[i][j] = 0.0;
            }
        }

        /* To put in after we are finished testing
        if (dt > 0.05)
        {
            cout << "Your dt is too big. It should be 0.05 or smaller, yours is " << dt << endl;
            exit(0);
        }
        */

        cout << N << endl;
        for (int i = 0; i < N; i++)
        {
            cout << T[i] << endl;
            //Creating filename and initializing file
            string filename = "../benchmarks/task_e/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T[i], 1) +".xyz";
            ofstream m_file;
            InitializeFile(filename, m_file);

            MetropolisE(L, MC_cycles, T[i], ExpectationValues[i], stringdir, my_rank);

            m_file.close();
        }

    }

    else
    {
        cout << "Wrong. Try again" << endl;
    }
    MPI_Finalize ();
    return 0;
}
