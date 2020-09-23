#include <string>
#include <iostream>
#include <fstream>
#include <map>

class ConfigParser{
private:
    std::map<std::string, unsigned int> parameters;
public:
    ConfigParser(std::string filename);
    unsigned int operator [](std::string key);
};
