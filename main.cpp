/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include <iostream>

int main(int argc, char *argv[]) {
  if(argc < 2) {
    std::cout << "Usage: <>" << std::endl;
  }

  DGFEMSpace1D(N, Nx);

  return 0;
}

