#include <iostream>
#include <tuple>
#include <mpi.h>
#include "parser.h"
#include "function.h"

using namespace std;


tuple<int, int, int> factor_number(int N){
    /* 
    Retruns a tuple with block numbers
    along axes in the following order: x, y, z.
    */
    // find maximum factor of N
    int i = 1, j = 1;
    for (int i = N - 1; i > 1; i--){
        if (N % i == 0) {
            break;
        }
    }
    if (i == 1){
        return make_tuple(N, 1 ,1);
    }
    j = N / i;
    return make_tuple(i, j, N / i / j);
}


int main(int argc, char** argv){
    int rank, world; // rank - process id, world - number of processes
    int Nx, Ny, Nz; // number of dots on the grid on each direction
    int nx, ny, nz; // number of dots in each block
    double Lx, Ly, Lz; // grid length
    int T; // number of timesteps for calculation
    tuple<int, int, int> tmp;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    // reading configuration
    ConfigParser* parser;
    try{
        parser = new ConfigParser(argc, argv);
    }
    catch (int e){
        return e;
    }
    cout << "Finished reading paramerers in process " << rank << endl;

    // filling variables and parameters
    tmp = factor_number(world);
    nx = get<0>(tmp);
    ny = get<1>(tmp);
    nz = get<2>(tmp);
    Function function = Function(parser);

    if (rank == 0){
        cout << nx << " " << ny << " " << nz << endl;
    }

    // cout << rank << " " << function(0.5, 0.75, 0.25, 1) << endl;

    MPI_Finalize();
    return 0;
}
