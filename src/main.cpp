#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <armadillo>
#include <time.h>
#include "ising.h"

using namespace std;
using namespace arma;

/* TODO:
 *  - fix command line arguments
 *    - L, MC_cycles and task at the very least need to be read in
 *    - might change a bit when parallelization is added, but thats a problem when we fix that ^^
 *  - structure for the different tasks, and actually what the tasks ask us to do
 *    - 2x InitializeLattice functions, one for only up spin start state, and one with random start state
 *  - parallelization with MPI!!!
 *  - fix so that WriteToFile actually writes to file
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
 *                   ^   *                     *    ^
 *                  ^   *   /\   /\   /\   /\   *    ^
 *                 ^   *                         *    ^
 *                 ^  *                           *   ^
 *                 ^  *                           *   ^
 *                  ^ *                           *  ^
 *                   ^*                           * ^
 *                    ^ *                        * ^
 *                    ^  *                      *  ^
 *                      ^  *       ) (         * ^
 *                          ^^^^^^^^ ^^^^^^^^^
 */


int main(int argc, char *argv[])
{
    // Task b)
    int L = 2;
    double T = 1.0;
    int MC_cycles = 1000000;

    vec ExpectationValues = zeros<mat>(5);

    Metropolis(L, MC_cycles, T, ExpectationValues);

    WriteToFile(L, MC_cycles, T, ExpectationValues);

    return 0;
}
