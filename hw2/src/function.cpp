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
    return PI * sqrt(1. / (Lx * Lx) + 4. / (Ly * Ly) + 9. / (Lz * Lz));

}

float Function::operator ()(float x, float y, float z, float t){
    if (x == 0.0 || x == Lx) return 0;
    // TODO: missing derivative on y!!
    cout << "MISSING Y corner" << endl;
    if (z == 0 || z == Lz) return 0;
    return sin(PI / Lx * x) * sin(2 * PI / Ly * y) * sin(3 * PI / Lz * z) * cos(at * t);
}

float Function::phi(float x, float y, float z){
    return this->operator()(x, y, z, 0);
}