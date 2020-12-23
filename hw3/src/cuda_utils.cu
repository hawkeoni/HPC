__global__ void calculate_inner_grid(double* grid_0, double* grid_1, double* grid_2, int bx, int by, int bz){
    int N = (bx + 2) * (by + 2) * (bz * 2);
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    double uijk = grid_1[idx], laplace = 0.;
    i = N % (bx + 2);
    if (i < 2 && i >= bx) return;
    j = N / (bx + 2) % (by + 2);
    if (j < 2 && j >= by) return;
    k = N / ((bx + 2) * (by + 2));
    if (k < 2 && k >= bz) return;
    grid_2[idx] = 2 * grid_1[idx] - grid_0[idx];
}

__global__ void first_step(double* grid_0, double* grid_1, \
        int bx, int by, int bz, \
        double hx, double hy, double hz, \
        double block_x_len, double block_y_len, double block_z_len, \
        int Lx, int Ly, int Lz, \
        int Nx, int Ny, int Nz, \
        int nx, int ny, int nz, \
        int block_pos_x, int block_pos_y, int block_pos_z,
        double at, double t){
    int N = (bx + 2) * (by + 2) * (bz * 2);
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    double uijk = grid_1[idx], laplace = 0.;
    double x, y, z;
    i = N % (bx + 2);
    if (i < 1 && i > bx) return;
    j = N / (bx + 2) % (by + 2);
    if (j < 1 && j > by) return;
    k = N / ((bx + 2) * (by + 2));
    if (k < 1 && k > bz) return;
    x = (i - 1) * hx + block_pos_x * block_x_len + min(Nx % nx, block_pos_x) * hx;
    y = (j - 1) * hy + block_pos_y * block_y_len + min(Ny % ny, block_pos_y) * hy;
    z = (k - 1) * hz + block_pos_z * block_z_len + min(Nz % nz, block_pos_z) * hz;
    grid_0[idx] = sin(3.14 / Lx * x) * sin(3.14 / Ly * y) * sin(2 * 3.14 / Lz * z) * cos(at * t + 2 * 3.14);
}
