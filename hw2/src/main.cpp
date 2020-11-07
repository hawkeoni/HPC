#include <iostream>
#include <mpi.h>
#include "parser.h"
#include "function.h"

using namespace std;



int main(){
    // cout << (*parser)["Ly"] << endl;
    // cout << function.at << endl;
    // cout << function(1, 2, 3, 0) << endl;
    // cout << function.phi(1, 2, 3) << endl;
    ConfigParser* parser = new ConfigParser("configs/config_base.txt");
    Function function = Function(parser);
    MPI_Init(NULL, NULL);
    int rank;
    int world;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    cout << rank << " " << function(0.5, 0.75, 0.25, 1) << endl;
    MPI_Finalize();
    return 0;
}
