//
// Created by Iliya Dimov on 2020-09-23.
//

#include <iostream>
#include "parser.h"
#include "function.h"

using namespace std;



int main(){
    ConfigParser* parser = new ConfigParser("configs/config_base.txt");
    cout << (*parser)["Ly"] << endl;
    Function function = Function(parser);
    cout << function.at << endl;
    cout << function(1, 2, 3, 0) << endl;
    cout << function.phi(1, 2, 3) << endl;
}