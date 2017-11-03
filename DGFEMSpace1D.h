/**
 * @file DGFEMSpace1D.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H
#include "vvector.h"
#include "Quadrature.h"

/**
 * @brief dimension of the equation, 1 for scalar equation and 3 for Euler equations
 */
#define VEC vvector
const int DIM = 3;
typedef VEC<VEC<double> > bU;
typedef VEC<bU> SOL;
typedef int BM;
typedef std::vector<std::vector<double> > QUAD;

class DGFEMSpace1D {
  private:
    u_int K;
    u_int Nx;
    SOL sol;
    BM bml, bmr;

  public:
    DGFEMSpace1D(u_int K, u_int Nx) : K(K), Nx(Nx) {
      sol.resize(Nx);
      for(u_int i = 0; i < Nx; ++i) {
        sol[i].resize(DIM);
        for(u_int j = 0; j < DIM; ++j) {
          sol[i][j].resize(K);
        }
      }
    }
    void init();
    void run(double t_end);
};

#endif //DGFEMSPACE1D_H

