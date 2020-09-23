#include <iostream>
#include <cmath>
#include "parser.h"

#define PI 3.14159265

class Function{
private:
    float calculate_at();
public:
    unsigned int Lx, Ly, Lz;
    float at;
    Function(ConfigParser *parser);
    Function(unsigned int Lx, unsigned int Ly, unsigned int Lz);
    float operator ()(float x, float y, float z, float t);
};