#include "parser.h"

using namespace std;

ConfigParser::ConfigParser(string filename){
        string line, varname;
        unsigned int value;
        size_t pos_0 = 0, pos_1 = 0;
        cout << "Reading configuration from " << filename << "." << endl;
        ifstream file(filename);
        while (getline(file, line)){
            pos_0 = line.find(" ");
            pos_1 = line.find(" ", pos_0 + 1);
            varname = line.substr(0, pos_0);
            value = atoi(line.substr(pos_1 + 1, line.size() - pos_1 - 1).c_str());
            this->parameters[varname] = value;
        }
    };

ConfigParser::ConfigParser(int argc, char **args){
    if (argc != 8){
        cout << "ConfigParser expectes 7 arguments, got " << argc - 1 << endl;
        cout << "Arguments should be Lx, Ly, Lz, Nx, Ny, Nz, T" << endl;
        throw 1;
    }
    const string varnames[] = {"Lx", "Ly", "Lz", "Nx", "Ny", "Nz", "T"};
    for (int i = 1; i < 4; i++){
        this->parameters[varnames[i - 1]] = atof(args[i]);
    }
    for (int i = 4; i < argc; i++){
        this->parameters[varnames[i - 1]] = atoi(args[i]);
    }
}

unsigned int ConfigParser::operator [](string key){
    return this->parameters[key];
}