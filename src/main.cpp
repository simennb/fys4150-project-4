#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <armadillo>
#include <time.h>
#include "ising.h"
#include "functions.h"

using namespace std;
using namespace arma;

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
        cout << "Task can be set to b, c" << endl;
        cout << "L is the dimension of the spin lattice" << endl;
        cout << "T is temperature" << endl;
        cout << "Monte Carlo Cycles is the number of Monte Carlo Cycles" << endl;
        cout << "Spin configuration is optional, set to up by default." << endl;
        cout << "Spin configuration can either be set to up, down or rand to set all spins to be either up, down or random" << endl;
        exit(1);
    }
    // initialize some values
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
    // Task b)
        //string filename = "../benchmarks/task_b/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T, 1) +".xyz";
        vec ExpectationValues = zeros<mat>(5);
        string filename = "Tull";
        Metropolis(L, MC_cycles, T, ExpectationValues, stringdir, filename);

        //WriteToFile(L, MC_cycles, T, ExpectationValues, filename);

    }

    else if (strcmp(argv[1], "c") == 0)
    {
        //Creating filename
        //string filename = "../benchmarks/task_c/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T, 1) +".xyz";
        vec ExpectationValues = zeros<mat>(5);
        string filename = "Tull";

        ofstream m_file;
        m_file.open(filename);

        Metropolis(L, MC_cycles, T, ExpectationValues, stringdir, filename);

        //WriteToFile(L, MC_cycles, T, ExpectationValues, filename);

        m_file.close();
    }

    else if (strcmp(argv[1], "d") == 0)
    {
        //Creating filename
        string filename = "../benchmarks/task_d/eigenvalues_MC"+to_scieni(MC_cycles, 1) + "_dim"+to_string(L)+"_dir" + stringdir + "_T" + to_fixf(T, 1) +".xyz";
        vec ExpectationValues = zeros<mat>(5);

        ofstream m_file;
        m_file.open(filename);

        Metropolis(L, MC_cycles, T, ExpectationValues, stringdir, filename);

        //WriteToFile(L, MC_cycles, T, ExpectationValues, filename);

        m_file.close();
    }

    else
    {
        cout << "Wrong. Try again" << endl;
    }
    return 0;
}
