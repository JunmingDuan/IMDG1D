/**
 * @file main.cpp
 * @brief Implicit DG for 1D conservation laws
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include <iostream>
#include "DGFEMSpace1D.h"

VEC<double> f0(const double x, const double t) {
  VEC<double> u(DIM);
  u[0] = 1;//pow(sin(x), 2);
  u[1] = 2*x;//pow(sin(x), 2);
  u[2] = x;//pow(sin(x), 2);
  return u;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  F = u;
  return F;
}

VEC<VEC<double> > f_prime(const VEC<double>& u) {
  VEC<VEC<double> > a(DIM);
  for(u_int i = 0; i < DIM; ++i)
    a[i].resize(DIM, 0);
  //f(u)=u;
  a[0][0] = 1;
  a[1][1] = 2;
  a[2][2] = 3;
  return a;
}

VEC<double> exact(const double x, const double t) {
  return (sin(4*x)-8*sin(2*x)+12*x)/32;
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
  Problem.BuildQuad(2);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);
  std::cout << "Start to solve..." << std::endl;
  t1 = clock();
  Problem.run(f_prime, t_end);
  t2 = clock();
  std::cout << "Time consumed: " << std::setw(8) << (t2-t1)/CLOCKS_PER_SEC << std::endl;

  return 0;
}

