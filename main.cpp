/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include <iostream>
#include "DGFEMSpace1D.h"

VEC<double> f0(double x, double t) {
  VEC<double> u(3);
  u[0] = pow(sin(x), 2);
  u[1] = 1;//pow(sin(x), 2);
  u[2] = x;//pow(sin(x), 2);
  return u;
}

int main(int argc, char *argv[]) {
  if(argc < 5) {
    std::cout << "Usage: <Nx> <xl> <xr> <t_end> " << std::endl;
    abort();
  }

  clock_t t1, t2;
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  double t_end = atof(argv[4]);
  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(5);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);
  std::cout << "Start to solve..." << std::endl;
  t1 = clock();
  Problem.run(t_end);
  t2 = clock();
  std::cout << "Time consuming: " << std::setw(8) << (t2-t1)/CLOCKS_PER_SEC << std::endl;

  return 0;
}

