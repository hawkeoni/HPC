#include <iostream>
#include <cmath>
#include "parser.h"

#define PI 3.14159265

class Function{
private:
    double calculate_at();
public:
    unsigned int Lx, Ly, Lz;
    double at;
    Function(ConfigParser *parser);
    Function(unsigned int Lx, unsigned int Ly, unsigned int Lz);
    double operator ()(double x, double y, double z, double t);
    double phi(double x, double y, double z);
};