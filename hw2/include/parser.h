#include <string>
#include <iostream>
#include <fstream>
#include <map>


#ifndef PARSER_H
#define PARSER_H

class ConfigParser{
private:
    std::map<std::string, unsigned int> parameters;
public:
    ConfigParser(std::string filename);
    ConfigParser(int argc, char **args);
    unsigned int operator [](std::string key);
};



#endif

