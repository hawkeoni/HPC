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
            value = stoi(line.substr(pos_1 + 1, line.size() - pos_1 - 1));
            // cout << varname << "|" << value << endl;
            this->parameters[varname] = value;
        }
    };

unsigned int ConfigParser::operator [](string key){
    return this->parameters[key];
}