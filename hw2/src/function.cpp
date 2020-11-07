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

float Function::calculate_at(){
    return PI * sqrt(1. / (Lx * Lx) + 1. / (Ly * Ly) + 4. / (Lz * Lz));

}

float Function::operator ()(float x, float y, float z, float t){
    return sin(PI / Lx * x) * sin(PI / Ly * y) * sin(2 * PI / Lz * z) * cos(at * t + 2 * PI);
}

float Function::phi(float x, float y, float z){
    return this->operator()(x, y, z, 0);
}
