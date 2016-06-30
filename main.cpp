/*
 * main.cpp
 *
 *  Created on: 2016/6/29
 *      Author: Tsukasa Fukunaga
 */

#include <string>
#include <stdlib.h>
#include <iostream>
#include "CapR.h"

using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 4){
    cout << "Error: The number of argument is invalid." << endl;
    return(0);
  }
  
  string input_file = argv[1];
  string output_file = argv[2];
  int constraint =  atoi(argv[3]);
  CapR capr(input_file, output_file, constraint);
  capr.Run();
  
  return(0);
}
