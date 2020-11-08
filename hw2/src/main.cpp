#include <iostream>
#include <tuple>
#include <mpi.h>
#include "parser.h"
#include "function.h"

using namespace std;

tuple<int, int, int> get_block_positions(int N, int nx, int ny, int nz){
    return make_tuple(N % nx, N / nx % ny, N / (nx * ny));
}

tuple<int, int, int> factor_number(int N){
    /* 
    Retruns a tuple with block numbers
    along axes in the following order: x, y, z.
    */
    // find maximum factor of N
    int i = 1, j = 1;
    for (i = N - 1; i >= 1; i--){
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
    int i, j, k; // counters for loops
    int rank, world; // rank - process id, world - number of processes
    int Nx, Ny, Nz; // number of dots on the grid on each direction
    int nx, ny, nz; // number of blocks in each direction
    int bx, by, bz; // number of dots in each block
    unsigned long B;
    double Lx, Ly, Lz, T, tau; // grid length & time
    double hx, hy, hz; // block length
    int block_pos_x, block_pos_y, block_pos_z;
    double x, y, z;
    tuple<int, int, int> tmp;


    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    // reading configuration

    Lx = atof(argv[1]); Ly = atof(argv[2]); Lz = atof(argv[3]);
    Nx = atoi(argv[4]); Ny = atoi(argv[5]); Nz = atoi(argv[6]);
    T = atof(argv[7]); tau = T / 20;

    ConfigParser* parser;
    try{
        parser = new ConfigParser(argc, argv);
    }
    catch (int e){
        return e;
    }
    

    // filling variables and parameters
    tmp = factor_number(world);
    nx = get<0>(tmp); bx = Nx / nx; hx = Lx / Nx;
    ny = get<1>(tmp); by = Ny / ny; hy = Ly / Ny;
    nz = get<2>(tmp); bz = Nz / nz; hz = Lz / Nz;
    if (rank == 0){
        cout << "AAAAAA" << endl;
        cout << bz << " " << Nz << " " << nz << endl;
    }
    B = bx * by * bz;
    //TODO: Add 2 from sides
    double grid_0[bx][by][bz], grid_1[bx][by][bz], grid_2[bx][by][bz];
    // Function function = Function(parser);
    Function function = Function(Lx, Ly, Lz);
    if (rank == 0)
        cout << "Function " << function(0.5, 0.75, 0.25, 1) << endl;
    
    if (rank == 0){
        cout << "Block numbers" << endl;
        cout << nx << " " << ny << " " << nz << endl;
    }

    tmp = get_block_positions(rank, nx, ny, nz);
    block_pos_x = get<0>(tmp);
    block_pos_y = get<1>(tmp);
    block_pos_z = get<2>(tmp);

    // Creating u0
    //TODO: can I optimize this?
    // optimize operations by excluding hx and pre-calculating block size
    if (rank == 0){
        cout << bx << " " << by << " " << bz << endl;
        cout << block_pos_z << endl;
    }

    for (i = 0; i < bx; i ++){
        for (j = 0; j < by; j ++){
            for (k = 0; k < bz; k++){
                x = i * hx + block_pos_x * bx * hx;
                y = j * hy + block_pos_y * by * hy;
                z = k * hz + block_pos_z * bz * hz;
                grid_0[i][j][k] = function.phi(x, y, z);
                // laplace phi
                grid_1[i][j][k] = tau * tau / 2

                grid_1[i][j][k] += grid_0[i][j][k];
            }
        }
    }



    MPI_Finalize();
    return 0;
}
