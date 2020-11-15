#include <stdio.h>
#include <mpi.h>
#include "function.h"

using namespace std;



int i, j, k, p; // counters for loops
int Nx, Ny, Nz; // number of dots on the grid on each direction
int nx, ny, nz; // number of blocks in each direction
int bx, by, bz; // number of dots in each block
unsigned long B;
double Lx, Ly, Lz, T, tau; // grid length & time
double hx, hy, hz; // lengths between dots on all axes
int block_pos_x, block_pos_y, block_pos_z;
double block_x_len, block_y_len, block_z_len; // block lengths
double x, y, z;
int *tmp;
int K = 20;
bool debug = true;
double uijk, laplace;
double ***grid_0, ***grid_1, ***grid_2;
double *distances;
int *rankptr, *worldptr;
Function *u_analytical;

int MAIN_PROCESS = 0;


void calculate_error(double ***grid, double t, int step){
    if (debug)
        printf("Calculating error in process %d\n", *rankptr);
    double distance = -1;
    for (i = 1; i < bx + 2; i++){
        for (j = 1; j < by + 2; j++){
            for (k = 1; k < bz + 2; k++){
                x = (i - 1) * hx + block_pos_x * block_x_len;
                y = (j - 1) * hy + block_pos_y * block_y_len;
                z = (k - 1) * hz + block_pos_z * block_z_len;
                distance = max(distance, abs(grid[i][j][k] - (*u_analytical)(x, y, z, t)));
            }
        }
    }
    if (debug)
        printf("Process %d, error %f\n", *rankptr, distance);

    MPI_Gather(&distance, 1, MPI_DOUBLE, distances, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (*rankptr == 0){
        
        for (i = 0; i < *worldptr; i++){
            printf("%d\n", i);
            distance = max(distance, distances[i]);
        }
        printf("===========================================================\n");
        printf("At step %d and time %f maximum error is %f\n", step, t, distance);
        printf("===========================================================\n");
    }
}

int* get_block_positions(int N){
    /*
    There are N blocks, specifacally nx * ny * nz.
    This function calculates block indices by its number.
    */
    int *res = new int[3];
    res[0] = N % nx; res[1] = N / nx % ny; res[2] = N / (nx * ny);
    return res;
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
    int rank, world; // rank - process id, world - number of processes

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    // reading configuration

    Lx = atof(argv[1]); Ly = atof(argv[2]); Lz = atof(argv[3]);
    Nx = atoi(argv[4]); Ny = atoi(argv[5]); Nz = atoi(argv[6]);
    T = atof(argv[7]); tau = T / K;

    tmp = factor_number(world);

    nx = tmp[0]; bx = Nx / nx; hx = Lx / Nx;
    ny = tmp[1]; by = Ny / ny; hy = Ly / Ny;
    nz = tmp[2]; bz = Nz / nz; hz = Lz / Nz;

    rankptr = &rank; worldptr = &world;

    u_analytical = new Function(Lx, Ly, Lz);
    distances = new double[*worldptr];

    // input read
    if (rank == 0){
        printf("===========================================================\n");
        printf("Number of processes %d\n", world);
        printf("Lx = %f, Ly = %f, Lz = %f\n", Lx, Ly, Lz);
        printf("Number of dots across axes: %d, %d, %d\n", Nx, Ny, Nz);
        printf("T = %f, K = %d, tau = %f\n", T, K, tau);
        printf("Number of blocks across axes: %d, %d, %d\n", nx, ny, nz);
        printf("Number of dots in a block: %d, %d, %d\n", bx, by, bz);
        printf("Distance between dots in a block: %f, %f, %f\n", hx, hy, hz);
        printf("===========================================================\n");
    }

    B = bx * by * bz;

    // creating grids
    grid_0 = new double**[bx + 2];
    grid_1 = new double**[bx + 2];
    grid_2 = new double**[bx + 2];
    for (i = 0; i < bx + 2; i++){
        grid_0[i] = new double*[by + 2];
        grid_1[i] = new double*[by + 2];
        grid_2[i] = new double*[by + 2];
        for (j = 0; j < by + 2; j++){
            grid_0[i][j] = new double[bz + 2];
            grid_1[i][j] = new double[bz + 2];
            grid_2[i][j] = new double[bz + 2];
        }
    }


    tmp = get_block_positions(rank);
    block_pos_x = tmp[0];
    block_pos_y = tmp[1];
    block_pos_z = tmp[2];
    block_x_len = bx * hx; block_y_len = by * hy; block_z_len = bz * hz;


    if (rank == 0 && debug){
        printf("Preparing u_0 and u_1\n");    
    }
    
    for (i = 1; i < bx + 1; i ++){
        for (j = 1; j < by + 1; j ++){
            for (k = 1; k < bz + 1; k++){
                // dot x, y, z = block offset + dot offset inside the block
                // u0
                x = (i - 1) * hx + block_pos_x * block_x_len;
                y = (j - 1) * hy + block_pos_y * block_y_len;
                z = (k - 1) * hz + block_pos_z * block_z_len;
                grid_0[i][j][k] = u_analytical->phi(x, y, z);
                // u1
                uijk = grid_0[i][j][k];
                laplace = 0.;
                laplace += (u_analytical->phi(x - hx, y, z) - 2 * uijk + u_analytical->phi(x + hx, y, z)) / (hx * hx);
                laplace += (u_analytical->phi(x, y - hy, z) - 2 * uijk + u_analytical->phi(x, y + hy, z)) / (hy * hy);
                laplace += (u_analytical->phi(x, y, z - hz) - 2 * uijk + u_analytical->phi(x, y, z + hz)) / (hz * hz);
                grid_1[i][j][k] = uijk + tau * tau / 2 * laplace;
            }
        }
    }
    if (rank == 0 && debug){
        printf("Finished with u_0 and u_1\n");
    }
    // calculate errors
    calculate_error(grid_0, 0, 0);
    calculate_error(grid_1, tau, 1);


    //step
    // for (p = 2; p < K; p++){

    // }

    MPI_Finalize();
    return 0;
}
