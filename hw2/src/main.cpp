//
// Created by Iliya Dimov on 2020-09-23.
//

#include <iostream>
#include "parser.h"

using namespace std;



int main(){
    ConfigParser* parser = new ConfigParser("configs/config_base.txt");
    cout << (*parser)["Lx"] << endl;
}