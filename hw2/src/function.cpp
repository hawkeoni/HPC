#include "function.h"

using namespace std;

Function::Function(ConfigParser *parser){
    Lx = (*parser)["Lx"];
    Ly = (*parser)["Ly"];
    Lz = (*parser)["Lz"];
    at = calculate_at();
}

Function::Function(unsigned int Lx, unsigned int Ly, unsigned int Lz){
    Lx = Lx;
    Ly = Ly;
    Lz = Lz;
    at = calculate_at();
}

double Function::calculate_at(){
    return PI * sqrt(1. / (Lx * Lx) + 1. / (Ly * Ly) + 4. / (Lz * Lz));

}

double Function::operator ()(double x, double y, double z, double t){
    return sin(PI / Lx * x) * sin(PI / Ly * y) * sin(2 * PI / Lz * z) * cos(at * t + 2 * PI);
}

double Function::phi(double x, double y, double z){
    return this->operator()(x, y, z, 0);
}
