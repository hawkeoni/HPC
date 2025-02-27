#include <iostream>
#include <cmath>

#define PI 3.14159265

class Function{
private:
    double calculate_at();
public:
    double Lx, Ly, Lz;
    double at;
    Function(double Lx, double Ly, double Lz);
    double operator ()(double x, double y, double z, double t);
    double phi(double x, double y, double z);
};