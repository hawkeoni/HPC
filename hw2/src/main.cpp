#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <algorithm>
#include "function.h"

using namespace std;



int i, j, k, p; // counters for loops
int Nx, Ny, Nz; // number of dots on the grid on each direction
int nx, ny, nz; // number of blocks in each direction
int bx, by, bz; // number of dots in each block
double Lx, Ly, Lz, T, tau; // grid length & time
double hx, hy, hz; // lengths between dots on all axes
int block_pos_x, block_pos_y, block_pos_z; // block index on axes
double block_x_len, block_y_len, block_z_len; // block lengths
double x, y, z;
double time_start, time_end; // time variables for speed calculation
int *tmp;
int K = 20;
bool debug = false;
double uijk, laplace;
double ***grid_0, ***grid_1, ***grid_2, ***tmpptr;
double *distances;
double *xleft_from, *xright_from, *yleft_from, *yright_from, *zleft_from, *zright_from; // where it is coming from in the source block
double *xleft_to, *xright_to, *yleft_to, *yright_to, *zleft_to, *zright_to; // where it is coming to
int *rankptr, *worldptr;
Function *u_analytical;
MPI_Request xleft_request_from, xright_request_from, yleft_request_from, yright_request_from, zleft_request_from, zright_request_from;
MPI_Request xleft_request_to, xright_request_to, yleft_request_to, yright_request_to, zleft_request_to, zright_request_to;
MPI_Status status;
int MAIN_PROCESS = 0;


void calculate_error(double ***grid, double t, int step){
    double local_distance;
    if (debug)
        printf("Calculating error in process %d\n", *rankptr);
    double distance = -1;
    for (i = 1; i < bx + 1; i++){
        for (j = 1; j < by + 1; j++){
            for (k = 1; k < bz + 1; k++){
                x = (i - 1) * hx + block_pos_x * block_x_len;
                y = (j - 1) * hy + block_pos_y * block_y_len;
                z = (k - 1) * hz + block_pos_z * block_z_len;
                local_distance = abs(grid[i][j][k] - (*u_analytical)(x, y, z, t));
                distance = max(distance, local_distance);
            }
        }
    }
    if (debug)
        printf("Process %d, error %f\n", *rankptr, distance);

    MPI_Gather(&distance, 1, MPI_DOUBLE, distances, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (*rankptr == 0){
        for (i = 0; i < *worldptr; i++){
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

int get_rank_by_block(int x_offset, int y_offset, int z_offset){
    /*
    Return block given its offsets. A reverse function
    for get_block_positions.
    */
    return x_offset + y_offset * nx + z_offset * nx * ny;
}
 

int* factor_number(int N){
    /* 
    Retruns a tuple with block numbers
    along axes in the following order: x, y, z.
    */
    // find maximum factor of N
    int *res = new int[3];
    res[0] = 1; res[1] = 1; res[2] = 1;
    int pos = 0;
    int i = 2;
    while (N > 1){
        if (N % i == 0){
            res[pos] *= i;
            pos = (pos + 1) % 3;
            N /= i;
        }
        else i++;
    }
    // printf("%d %d %d\n", res[0], res[1], res[2]);
    return res;
}


void step(){
    
    // Sending xleft
    if (block_pos_x > 0){
        p = 0;
        for (j = 1; j < by + 1; j++)
            for (k = 1; k < bz + 1; k++)
                xleft_from[p++] = grid_1[1][j][k];
        MPI_Isend(xleft_from, by * bz, MPI_DOUBLE, *rankptr - 1, 0, MPI_COMM_WORLD, &xleft_request_from);
    }
    if (block_pos_x < nx - 1)
        MPI_Irecv(xleft_to, by * bz, MPI_DOUBLE, *rankptr + 1, 0, MPI_COMM_WORLD, &xleft_request_to);
    // Sending xright
    if (block_pos_x < nx - 1){
        p = 0;
        for (j = 1; j < by + 1; j++)
            for (k = 1; k < bz + 1; k++)
                xright_from[p++] = grid_1[bx][j][k];
        MPI_Isend(xright_from, by * bz, MPI_DOUBLE, *rankptr + 1, 0, MPI_COMM_WORLD, &xright_request_from);
    }
    if (block_pos_x > 0)
        MPI_Irecv(xright_to, by * bz, MPI_DOUBLE, *rankptr - 1, 0, MPI_COMM_WORLD, &xright_request_to);

    // Sending yleft
    if (block_pos_y > 0){
        p = 0;
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++)
                yleft_from[p++] = grid_1[i][1][k];
        MPI_Isend(yleft_from, bx * bz, MPI_DOUBLE, *rankptr - nx, 0, MPI_COMM_WORLD, &yleft_request_from);
    }
    if (block_pos_y < ny - 1)
        MPI_Irecv(yleft_to, bx * bz, MPI_DOUBLE, *rankptr + nx, 0, MPI_COMM_WORLD, &yleft_request_to);

    // Sending yright
    if (block_pos_y < ny - 1){
        p = 0;
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++)
                yright_from[p++] = grid_1[i][by][k];
        MPI_Isend(yright_from, bx * bz, MPI_DOUBLE, *rankptr + nx, 0, MPI_COMM_WORLD, &yright_request_from);
    }
    if (block_pos_y > 0)
        MPI_Irecv(yright_to, bx * bz, MPI_DOUBLE, *rankptr - nx, 0, MPI_COMM_WORLD, &yright_request_to);

    // Sending zleft
    p = 0;
    for (i = 1; i < bx + 1; i++)
        for (j = 1; j < by + 1; j++)
            zleft_from[p++] = grid_1[i][j][2];
    if (block_pos_z > 0) 
        MPI_Isend(zleft_from, bx * by, MPI_DOUBLE, *rankptr - nx * ny, 0, MPI_COMM_WORLD, &zleft_request_from);
    else if (block_pos_z == 0)
        MPI_Isend(zleft_from, bx * by, MPI_DOUBLE, get_rank_by_block(block_pos_x, block_pos_y, nz - 1), 0, MPI_COMM_WORLD, &zleft_request_from);
    if (block_pos_z < nz - 1)
        MPI_Irecv(zleft_to, bx * by, MPI_DOUBLE, *rankptr + nx * ny, 0, MPI_COMM_WORLD, &zleft_request_to);
    else if (block_pos_z == nz - 1)
        MPI_Irecv(zleft_to, bx * by, MPI_DOUBLE, get_rank_by_block(block_pos_x, block_pos_y, 0), 0, MPI_COMM_WORLD, &zleft_request_to);

    // Sending zright
    p = 0;
    for (i = 1; i < bx + 1; i++)
        for (j = 1; j < by + 1; j++)
            zright_from[p++] = grid_1[i][j][bz - 1];
    if (block_pos_z < nz - 1)
        MPI_Isend(zright_from, bx * by, MPI_DOUBLE, *rankptr + nx * ny, 0, MPI_COMM_WORLD, &zright_request_from);
    else if (block_pos_z == nz - 1)
        MPI_Isend(zright_from, bx * by, MPI_DOUBLE, get_rank_by_block(block_pos_x, block_pos_y, 0), 0, MPI_COMM_WORLD, &zright_request_from);
    if (block_pos_z > 0)
        MPI_Irecv(zright_to, bx * by, MPI_DOUBLE, *rankptr - nx * ny, 0, MPI_COMM_WORLD, &zright_request_to);
    else if (block_pos_z == 0)
        MPI_Irecv(zright_to, bx * by, MPI_DOUBLE, get_rank_by_block(block_pos_x, block_pos_y, nz - 1), 0, MPI_COMM_WORLD, &zright_request_to);

    if (debug) printf("Done sending in proc %d\n", *rankptr);


    /*
    from 2 to bx - 1, because boundary values have not been
    recieved yet.
    */

    for (i = 2; i < bx; i++){
        for (j = 2; j < by; j++){
            for (k = 2; k < bz; k++){
                grid_2[i][j][k] = 2 * grid_1[i][j][k] - grid_0[i][j][k];
                uijk = grid_1[i][j][k];
                laplace = 0;
                laplace += (grid_1[i - 1][j][k] - 2 * uijk + grid_1[i + 1][j][k]) / (hx * hx);
                laplace += (grid_1[i][j - 1][k] - 2 * uijk + grid_1[i][j + 1][k]) / (hy * hy);
                laplace += (grid_1[i][j][k - 1] - 2 * uijk + grid_1[i][j][k + 1]) / (hz * hz);
                grid_2[i][j][k] += tau * tau * laplace;
            }
        }
    }

    if (debug) printf("Finished main loop\n");

    // Wait xleft
    if (block_pos_x < nx - 1){ 
        MPI_Wait(&xleft_request_to, &status);
        for (j = 1; j < by + 1; j++){
            for (k = 1; k < bz + 1; k++){
                grid_1[bx + 1][j][k] = xleft_to[(k - 1) + bz * (j - 1)];
                // j = by, k = bz -> index = (bz - 1) + bz * (by - 1) =
                // bz * by - 1 *
            }
        }
    }

    // Wait xright
    if (block_pos_x > 0){
        MPI_Wait(&xright_request_to, &status);
        for (j = 1; j < by + 1; j++){
            for (k = 1; k < bz + 1; k++){
            //    printf("Proc %d, j = %d, k = %d\n", *rankptr, j, k);
                grid_1[0][j][k] = xright_to[(k - 1) + bz * (j - 1)];
            }
        }
    }

    // Wait yleft
    if (block_pos_y < ny - 1){
        MPI_Wait(&yleft_request_to, &status);
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++){
                grid_1[i][by + 1][k] = yleft_to[(k - 1) + bz * (i - 1)];
            }
    }

    // Wait yright
    if (block_pos_y > 0){
        MPI_Wait(&yright_request_to, &status);
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++){
                grid_1[i][0][k] = yright_to[(k - 1) + bz * (i - 1)];
            }
    }

    // Wait zleft
    MPI_Wait(&zleft_request_to, &status);
    for (i = 1; i < bx + 1; i++)
        for (j = 1; j < by + 1; j++)
            grid_1[i][j][bz + 1] = zleft_to[(j - 1) + (i - 1) * by];

    // Wait zright
    MPI_Wait(&zright_request_to, &status);
    for (i = 1; i < bx + 1; i++)
        for (j = 1; j < by + 1; j++)
            grid_1[i][j][0] = zright_to[(j - 1) + (i - 1) * by];


    if (debug) printf("Finished recv\n");


    // updating boundary values
    for (i = 1; i <= bx; i += bx - 1)
        for (j = 1; j < by + 1; j++)
            for (k = 1; k < bz + 1; k++){
                if ((i == 1 && block_pos_x == 0) || (j == 1 && block_pos_y == 0) || (i == bx && block_pos_x == nx - 1) || (j == by && block_pos_y == ny - 1)) {
                    grid_2[i][j][k] = 0;
                    break;
                }
                grid_2[i][j][k] = 2 * grid_1[i][j][k] - grid_0[i][j][k];
                uijk = grid_1[i][j][k];
                laplace = 0;
                laplace += (grid_1[i - 1][j][k] - 2 * uijk + grid_1[i + 1][j][k]) / (hx * hx);
                laplace += (grid_1[i][j - 1][k] - 2 * uijk + grid_1[i][j + 1][k]) / (hy * hy);
                laplace += (grid_1[i][j][k - 1] - 2 * uijk + grid_1[i][j][k + 1]) / (hz * hz);
                grid_2[i][j][k] += tau * tau * laplace;
            }
    if (debug) printf("Finished boundary on i\n");
    for (j = 1; j <= by; j+= by - 1)
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++){
                if ((i == 1 && block_pos_x == 0) || (j == 1 && block_pos_y == 0) || (i == bx && block_pos_x == nx - 1) || (j == by && block_pos_y == ny - 1)) {
                    grid_2[i][j][k] = 0;
                    break;
                }
                grid_2[i][j][k] = 2 * grid_1[i][j][k] - grid_0[i][j][k];
                uijk = grid_1[i][j][k];
                laplace = 0;
                laplace += (grid_1[i - 1][j][k] - 2 * uijk + grid_1[i + 1][j][k]) / (hx * hx);
                laplace += (grid_1[i][j - 1][k] - 2 * uijk + grid_1[i][j + 1][k]) / (hy * hy);
                laplace += (grid_1[i][j][k - 1] - 2 * uijk + grid_1[i][j][k + 1]) / (hz * hz);
                grid_2[i][j][k] += tau * tau * laplace;
            }
    if (debug) printf("Finished boundary on j\n");

    for (k = 1; k <= bz; k += bz - 1)
        for (i = 1; i < bx + 1; i++)
            for (j = 1; j < by + 1; j++){
                if ((i == 1 && block_pos_x == 0) || (j == 1 && block_pos_y == 0) || (i == bx && block_pos_x == nx - 1) || (j == by && block_pos_y == ny - 1)) {
                    grid_2[i][j][k] = 0;
                    break;
                }
                grid_2[i][j][k] = 2 * grid_1[i][j][k] - grid_0[i][j][k];
                uijk = grid_1[i][j][k];
                laplace = 0;
                laplace += (grid_1[i - 1][j][k] - 2 * uijk + grid_1[i + 1][j][k]) / (hx * hx);
                laplace += (grid_1[i][j - 1][k] - 2 * uijk + grid_1[i][j + 1][k]) / (hy * hy);
                laplace += (grid_1[i][j][k - 1] - 2 * uijk + grid_1[i][j][k + 1]) / (hz * hz);
                grid_2[i][j][k] += tau * tau * laplace;

            }
    if (debug) printf("Finished boundary on k in %d\n", *rankptr);


}


int main(int argc, char** argv){
    int rank, world; // rank - process id, world - number of processes
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    // reading configuration
    Lx = atof(argv[1]); Ly = atof(argv[2]); Lz = atof(argv[3]);
    Nx = atoi(argv[4]); Ny = atoi(argv[5]); Nz = atoi(argv[6]);
    tmp = factor_number(world);
  
    nx = tmp[0]; nx = atoi(argv[7]); bx = Nx / nx; hx = Lx / (Nx - 1);
    ny = tmp[1]; ny = atoi(argv[8]); by = Ny / ny; hy = Ly / (Ny - 1);
    nz = tmp[2]; nz = atoi(argv[9]); bz = Nz / nz; hz = Lz / (Nz - 1);
    rankptr = &rank; worldptr = &world;
    // Deprecated
    // T = atof(argv[7]); tau = T / K;
    // New
    tau = min(min(hx, hy), hz) / 2; T = K * tau;

    u_analytical = new Function(Lx, Ly, Lz);
    distances = new double[*worldptr];
    memset(distances, 0, *worldptr);

    xleft_to = new double[by * bz]; xright_to = new double[by * bz];
    yleft_to = new double[bx * bz]; yright_to = new double[bx * bz];
    zleft_to = new double[bx * by]; zright_to = new double[bx * by];
    xleft_from = new double[by * bz]; xright_from = new double[by * bz];
    yleft_from = new double[bx * bz]; yright_from = new double[bx * bz];
    zleft_from = new double[bx * by]; zright_from = new double[bx * by];

    memset(xleft_to, 0, by * bz); memset(xright_to, 0, by * bz);
    memset(yleft_to, 0, bx * bz); memset(yright_to, 0, bx * bz);
    memset(zleft_to, 0, by * bz); memset(zright_to, 0, by * bz);
    memset(xleft_from, 0, by * bz); memset(xright_from, 0, by * bz);
    memset(yleft_from, 0, bx * bz); memset(yright_from, 0, bx * bz);
    memset(zleft_from, 0, by * bz); memset(zright_from, 0, by * bz);

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
            memset(grid_0[i][j], 0, bz + 2);
            grid_1[i][j] = new double[bz + 2];
            memset(grid_1[i][j], 0, bz + 2);
            grid_2[i][j] = new double[bz + 2];
            memset(grid_2[i][j], 0, bz + 2);
        }
    }

    tmp = get_block_positions(rank);
    block_pos_x = tmp[0];
    block_pos_y = tmp[1];
    block_pos_z = tmp[2];
    //block_x_len = (bx - 1) * hx; block_y_len = (by - 1) * hy; block_z_len = (bz - 1) * hz;
    block_x_len = (bx) * hx; block_y_len = (by) * hy; block_z_len = (bz) * hz;

    if (rank == 0 && debug){
        printf("Preparing u_0 and u_1\n");    
    }
    time_start = MPI_Wtime();
    for (i = 1; i < bx + 1; i ++){
        for (j = 1; j < by + 1; j ++){
            for (k = 1; k < bz + 1; k++){
                if ((i == 1 && block_pos_x == 0) || (j == 1 && block_pos_y == 0) || (i == bx && block_pos_x == nx - 1) || (j == by && block_pos_y == ny - 1)) {
                    grid_0[i][j][k] = 0;
                    grid_1[i][j][k] = 0;
                    break;
                }
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
                //if (i == 3 || j == 3) printf("Rankd %d x = %f y = %f z = %f value = %f\n", rank, x, y, z, grid_1[i][j][k]);
            }
        }
    }
    if (rank == 0 && debug){
        printf("Finished with u_0 and u_1\n");
    }
    // calculate errors
    calculate_error(grid_0, 0, 0);
    calculate_error(grid_1, tau, 1);

    //steps
    for (int stepnum = 2; stepnum < K + 1; stepnum++){
        if (debug && rank == 0)
            printf("Step %d\n", stepnum);
        // TODO: remove this after debug is over
        MPI_Barrier(MPI_COMM_WORLD);
        step();
        calculate_error(grid_2, stepnum * tau, stepnum);
        tmpptr = grid_0;
        grid_0 = grid_1;
        grid_1 = grid_2;
        grid_2 = tmpptr;
    }
    time_end = MPI_Wtime();
    if (rank == 0)
        printf("The program ran for %f time\n", time_end - time_start);

    MPI_Finalize();
    return 0;
}
