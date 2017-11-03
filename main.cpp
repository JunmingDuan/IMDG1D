/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include <iostream>
#include "DGFEMSpace1D.h"

int main(int argc, char *argv[]) {
  if(argc < 2) {
    std::cout << "Usage: <Np> <Nx> " << std::endl;
  }

  u_int K = atoi(argv[1]);
  u_int Nx = atoi(argv[2]);
  DGFEMSpace1D(K, Nx);

  return 0;
}

