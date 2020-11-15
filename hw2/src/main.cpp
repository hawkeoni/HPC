#include <iostream>
#include <mpi.h>
#include "function.h"

using namespace std;

int* get_block_positions(int N, int nx, int ny, int nz){
    int *res = new int[3];
    res[0] = N % nx; res[1] = N / nx % ny; res[2] = N / (nx * ny);
    return res;
    // return make_tuple(N % nx, N / nx % ny, N / (nx * ny));
}

int* factor_number(int N){
    /* 
    Retruns a tuple with block numbers
    along axes in the following order: x, y, z.
    */
    // find maximum factor of N
    int *res = new int[3];
    int i = 1, j = 1;
    for (i = N - 1; i >= 1; i--){
        if (N % i == 0) {
            break;
        }
    }
    if (i == 1){
        res[0] = N; res[1] = 1; res[2] = 1;
    }
    else{
        res[0] = i;
        j = N / i;
        res[1] = j;
        res[2] = N / i / 5;
        // return make_tuple(i, j, N / i / j);
    }
    return res;
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
    int *tmp;


    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    // reading configuration

    Lx = atof(argv[1]); Ly = atof(argv[2]); Lz = atof(argv[3]);
    Nx = atoi(argv[4]); Ny = atoi(argv[5]); Nz = atoi(argv[6]);
    T = atof(argv[7]); tau = T / 20;
    

    // filling variables and parameters
    tmp = factor_number(world);

    nx = tmp[0]; bx = Nx / nx; hx = Lx / Nx;
    ny = tmp[1]; by = Ny / ny; hy = Ly / Ny;
    nz = tmp[2]; bz = Nz / nz; hz = Lz / Nz;
    if (rank == 0){
        cout << "AAAAAA" << endl;
        cout << bx << " " << by << " " << bz << endl;
    }
    B = bx * by * bz;
    //TODO: Add 2 from sides
    cout << "before array" << endl;
    double ***grid_0, ***grid_1, ***grid_2;
    grid_0 = new double**[bx];
    grid_1 = new double**[bx];
    for (i = 0; i < bx; i++){
        grid_0[i] =  new double*[by];
        grid_1[i] =  new double*[by];
        for (j = 0; j < by; j++){
            grid_0[i][j] = new double[bz];
            grid_1[i][j] = new double[bz];
        }
    }
    // double grid_0[bx][by][bz], grid_1[bx][by][bz], grid_2[bx][by][bz];
    cout << "after array" << endl;
    // Function function = Function(parser);
    Function function = Function(Lx, Ly, Lz);
    if (rank == 0)
        cout << "Function " << function(0.5, 0.75, 0.25, 1) << endl;
    
    if (rank == 0){
        cout << "Block numbers" << endl;
        cout << nx << " " << ny << " " << nz << endl;
    }

    tmp = get_block_positions(rank, nx, ny, nz);

    block_pos_x = tmp[0];
    block_pos_y = tmp[1];
    block_pos_z = tmp[2];

    // Creating u0
    //TODO: can I optimize this?
    // optimize operations by excluding hx and pre-calculating block size
    if (rank == 0){
        cout << "weird if" << endl;
        cout << bx << " " << by << " " << bz << endl;
        cout << block_pos_z << endl;
    }
    cout << rank << " before loop1" << endl; 
    for (int p = 0; p < 300; p++)
    for (i = 0; i < bx; i ++){
        for (j = 0; j < by; j ++){
            for (k = 0; k < bz; k++){
                x = i * hx + block_pos_x * bx * hx;
                y = j * hy + block_pos_y * by * hy;
                z = k * hz + block_pos_z * bz * hz;
                grid_0[i][j][k] = function.phi(x, y, z);
                // laplace phi
                grid_1[i][j][k] = tau * tau / 2;
                grid_1[i][j][k] += grid_0[i][j][k];
            }
        }
    }
    // double maxdiff = 100;
    // for (i = 0; i < bx; i ++){
    //     for (j = 0; j < by; j ++){
    //         for (k = 0; k < bz; k++){
    //             x = i * hx + block_pos_x * bx * hx;
    //             y = j * hy + block_pos_y * by * hy;
    //             z = k * hz + block_pos_z * bz * hz;
    //             // maxdiff = max(maxdiff, abs(grid_0[i][j][k] - function(x, y, z, 0)));
    //         }
    //     }
    // }
    // cout << "Process " << rank << " maxdiff = " << maxdiff << endl; 


    MPI_Finalize();
    return 0;
}
