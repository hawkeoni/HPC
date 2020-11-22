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
int *tmp;
int K = 20;
bool debug = true;
double uijk, laplace;
double ***grid_0, ***grid_1, ***grid_2, ***tmpptr;
double *distances;
double *xleft, *xright, *yleft, *yright, *zleft, *zright; // where is it coming from in the source block
int *rankptr, *worldptr;
Function *u_analytical;
MPI_Request xleft_request, xright_request, yleft_request, yright_request, zleft_request, zright_request;
MPI_Status status;
int MAIN_PROCESS = 0;


void calculate_error(double ***grid, double t, int step){
    if (debug)
        printf("Calculating error in process %d\n", *rankptr);
    double distance = -1;
    for (i = 1; i < bx + 1; i++){
        for (j = 1; j < by + 1; j++){
            for (k = 1; k < bz + 1; k++){
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
    int i = 1, j = 1;
    for (i = max(N - 1, 1); i >= 1; i--){
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
        res[2] = N / i / j;
    }
    return res;
}

// buffers for sending

void step(){
    
    //ISends here https://github.com/bhavikm/Open-MPI-examples/blob/master/mpi_isend.c
    // Send x, y, z to the neigbouring blocks if I am middle

    // Sending xleft
    if (block_pos_x > 0){
        p = 0;
        for (j = 1; j < by + 1; j++)
            for (k = 1; k < bz + 1; k++)
                xleft[p++] = grid_1[1][j][k];
        MPI_Isend(xleft, p + 1, MPI_DOUBLE, *rankptr - 1, 0, MPI_COMM_WORLD, &xleft_request);
    }
    if (block_pos_x < nx - 1)
        MPI_Irecv(xleft, by * bz, MPI_DOUBLE, *rankptr + 1, 0, MPI_COMM_WORLD, &xleft_request);

    // Sending xright
    if (block_pos_x < nx - 1){
        p = 0;
        for (j = 1; j < by + 1; j++)
            for (k = 1; k < bz + 1; k++)
                xright[p++] = grid_1[bx][j][k];
        MPI_Isend(xright, p + 1, MPI_DOUBLE, *rankptr + 1, 0, MPI_COMM_WORLD, &xright_request);
    }
    if (block_pos_x > 0)
        MPI_Irecv(xright, by * bz, MPI_DOUBLE, *rankptr - 1, 0, MPI_COMM_WORLD, &xright_request);

    // Sending yleft
    if (block_pos_y > 0){
        p = 0;
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++)
                yleft[p++] = grid_1[i][1][k];
        MPI_Isend(yleft, bx * bz, MPI_DOUBLE, *rankptr - nx, 0, MPI_COMM_WORLD, &yleft_request);
    }
    if (block_pos_y < ny - 1)
        MPI_Irecv(yleft, by * bz, MPI_DOUBLE, *rankptr + nx, 0, MPI_COMM_WORLD, &yleft_request);

    // Sending yright
    if (block_pos_y < ny - 1){
        p = 0;
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++)
                yright[p++] = grid_1[i][by][k];
        MPI_Isend(yright, bx * bz, MPI_DOUBLE, *rankptr + nx, 0, MPI_COMM_WORLD, &yright_request);
    }
    if (block_pos_y > 0)
        MPI_Irecv(yright, by * bz, MPI_DOUBLE, *rankptr - nx, 0, MPI_COMM_WORLD, &yright_request);

    // Sending zleft
    p = 0;
    for (i = 1; i < bx + 1; i++)
        for (j = 1; j < by + 1; j++)
            zleft[p++] = grid_1[i][j][1];
    if (block_pos_z > 0) 
        MPI_Isend(zleft, bx * by, MPI_DOUBLE, *rankptr - nx * ny, 0, MPI_COMM_WORLD, &zleft_request);
    else if (block_pos_z == 0)
        MPI_Isend(zleft, bx * by, MPI_DOUBLE, get_rank_by_block(block_pos_x, block_pos_y, nz - 1), 0, MPI_COMM_WORLD, &zleft_request);
    if (block_pos_z < nz - 1)
        MPI_Irecv(zleft, bx * by, MPI_DOUBLE, *rankptr + nx * ny, 0, MPI_COMM_WORLD, &zleft_request);
    else if (block_pos_z == nz - 1)
        MPI_Irecv(zleft, bx * by, MPI_DOUBLE, get_rank_by_block(block_pos_x, block_pos_y, 0), 0, MPI_COMM_WORLD, &zleft_request);

    // Sending zright
    if (block_pos_z < nz - 1){
        p = 0;
        for (i = 1; i < bx + 1; i++)
            for (j = 1; j < by + 1; j++)
                zright[p++] = grid_1[i][j][bz];
        MPI_Isend(zright, bx * by, MPI_DOUBLE, *rankptr + nx * ny, 0, MPI_COMM_WORLD, &zright_request);
    }
    if (block_pos_z > 0)
        MPI_Irecv(zright, bx * by, MPI_DOUBLE, *rankptr - nx * ny, 0, MPI_COMM_WORLD, &zright_request);

        


    /*
    from 2 to bx - 1, because boundary values have not been
    recieved yet.
    */

    for (i = 2; i < bx; i++){
        for (j = 2; j < by; j++){
            for (k = 2; k < bz; k++){
                grid_2[i][j][k] = 1 * grid_1[i][j][k] - grid_0[i][j][k];
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
    /* 
    Waiting for IRecvs
    Calculating boundary values
    */
   // Wait xleft
    if (block_pos_x < nx - 1){ 
        printf("Before wait\n");
        MPI_Wait(&xleft_request, &status);
        printf("After wait\n");
        for (j = 1; j < by + 1; j++){
            for (k = 1; k < bz + 1; k++){
                grid_1[bx + 1][j][k] = xleft[(k - 1) + bz * (j - 1)];
                if (*rankptr == 1) cout << *rankptr << j << " " << k << endl;
                /* j = by, k = bz -> index = (bz - 1) + bz * (by - 1) =
                 bz * by - 1 */
            }
        if (debug) printf("Wait xleft\n");
        }
    }

    // Wait xright
    if (block_pos_x > 0){
        MPI_Wait(&xright_request, &status);
        for (j = 1; j < by + 1; j++){
            for (k = 1; k < bz + 1; k++){
                grid_1[0][j][k] = xright[(k - 1) + bz * (j - 1)];
            }
        }
        if (debug) printf("Wait xright\n");
    }

    // Wait yleft
    if (block_pos_y < ny - 1){
        MPI_Wait(&yleft_request, &status);
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++)
                grid_1[i][by + 1][k] = yleft[(i - 1) * bx + (k - 1) * bz];
    }

    // Wait yright
    if (block_pos_y > 0){
        MPI_Wait(&yright_request, &status);
        for (i = 1; i < bx + 1; i++)
            for (k = 1; k < bz + 1; k++)
                grid_1[i][0][k] = yright[(i - 1) * bx + (k - 1) * bz];

    }

    // Wait zleft
    MPI_Wait(&zleft_request, &status);
    for (i = 1; i < bx + 1; i++)
        for (j = 1; j < by + 1; j++)
            grid_1[i][j][bz + 1] = zright[(i - 1) * bx + (j - 1) * by];

    // Wait zright
    if (block_pos_z > 0){
        MPI_Wait(&zright_request, &status);
        for (i = 1; i < bx + 1; i++)
            for (j = 1; j < by + 1; j++)
                grid_1[i][j][0] = zright[(i - 1) * bx + (j - 1) * by];
    }


    if (debug) printf("Finished recv\n");

    // Updating boundary values of the grid
    for (i = 1; i <= bx; i += bx - 1){
        for (j = 1; j <= by; j += by - 1){
            for (k = 1; k <= bz; k+= bz - 1){
                grid_2[i][j][k] = 1 * grid_1[i][j][k] - grid_0[i][j][k];
                uijk = grid_1[i][j][k];
                laplace = 0;
                laplace += (grid_1[i - 1][j][k] - 2 * uijk + grid_1[i + 1][j][k]) / (hx * hx);
                laplace += (grid_1[i][j - 1][k] - 2 * uijk + grid_1[i][j + 1][k]) / (hy * hy);
                laplace += (grid_1[i][j][k - 1] - 2 * uijk + grid_1[i][j][k + 1]) / (hz * hz);
                grid_2[i][j][k] += tau * tau * laplace;
            }
        }
    }

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

    xleft = new double[by * bz]; xright = new double[by * bz];
    yleft = new double[bx * bz]; yright = new double[bx * bz];
    zleft = new double[bx * by]; zright = new double[bx * by];

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


    //steps
    for (p = 2; p < K + 1; p++){
        if (rank == 0 && debug)
            printf("Step %d\n", p);
        step();
        calculate_error(grid_2, p * tau, p);
        tmpptr = grid_0;
        grid_0 = grid_1;
        grid_1 = grid_2;
        grid_2 = tmpptr;
    }

    MPI_Finalize();
    return 0;
}
