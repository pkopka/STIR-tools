#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include "copytree3.h"


using namespace std;

int main (int argc, char* argv[]) {

  if (argc==7) {
    setenv("ASCII_FILE_NAME", argv[1], 1);
    setenv("ROOT_EXP", argv[2], 1);
    setenv("ROOT_OUTPUT_FILE_NAME", argv[3], 1);
    setenv("MLN_OF_CON", argv[4], 1);
    setenv("SCANNER_LENGHT", argv[5], 1);
    setenv("NUMBER_OF_PARTS", argv[6], 1);
  }
  else {
    cout<<"ERROR input";
    return 0;
    
  }
  
  copyt3(); 

  return 0;
  
}

