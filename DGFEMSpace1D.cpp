/**
 * @file DGFEMSpace1D.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include "DGFEMSpace1D.h"

DGFEMSpace1D::DGFEMSpace1D(u_int Nx, double xl, double xr)
  : Nx(Nx), xl(xl), xr(xr) {
  mesh.resize(Nx+1);
  h = (xr - xl)/Nx;
  for(u_int i = 0; i < Nx+1; ++i) {
    mesh[i] = i*h;
  }
  sol.resize(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    sol[i].resize(K);
    for(u_int k = 0; k < K; ++k) {
      sol[i][k].resize(DIM);
    }
  }
  //for(u_int i = 0; i < Nx+1; ++i) {
    //std::cout << mesh[i] << std::endl;
  //}
}

void DGFEMSpace1D::BuildQuad(u_int np) {
  std::vector<std::vector<double> > pw;
  pw = QUADINFO[0].LGL(np);
  TemQuad.set_np(np);
  TemQuad.set_points(pw[0]);
  TemQuad.set_weight(pw[1]);
  QUADINFO.resize(Nx);
  std::vector<double> gv(2), lv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    QUADINFO[i].set_np(np);
    gv[0] = mesh[i], gv[1] = mesh[i+1];
    QUADINFO[i].set_weight(pw[1]);
    std::vector<double> x(np);
    for(u_int j = 0; j < np; ++j) {
      local_to_global(pw[0][j], lv, gv, &x[j]);
    }
    QUADINFO[i].set_points(x);
    QUADINFO[i].set_jacobi( local_to_global_jacobian(lv, gv) );
  }
  for(u_int i = 0; i < Nx; ++i) {
    std::cout << i << "-th cell QUADINFO" << "\n";
    QUADINFO[i].print(std::cout);
  }
}

void DGFEMSpace1D::Projection(u_int cell, func f0, double t, bU& u) {
  //set to zero
  for(u_int i = 0; i < K; ++i) {
    for(u_int d = 0; d < DIM; ++d) {
      u[i][d] = 0;
    }
  }
  std::vector<double> x = TemQuad.points();
  std::vector<double> p = QUADINFO[cell].points();
  std::vector<double> w = QUADINFO[cell].weight();
  VEC<double> U;
  std::vector<double> V;
  double jab = QUADINFO[cell].l2g_jacobian();
  for(u_int k = 0; k < K; ++k) {
    double basis(0);
    for(u_int g = 0; g < x.size(); ++g) {
      U = f0(p[g], t);
      V = Poly(x[g]);
      u[k] += jab*U*V[k]*w[g];
      basis += jab*V[k]*V[k]*w[g];
    }
    u[k] /= basis;
  }
}

VEC<double> DGFEMSpace1D::Composition(u_int cell, double x, double t) {
  VEC<double> u(DIM);
  std::vector<double> lv(2), gv(2), V;
  lv[0] = -1, lv[1] = 1;
  gv[0] = mesh[cell], gv[1] = mesh[cell+1];
  double lp;
  global_to_local(x, lv, gv, &lp);
  V = Poly(lp);
  for(u_int k = 0; k < K; ++k) {
    u += sol[cell][k]*V[k];
  }

  return u;
}

void DGFEMSpace1D::init(func f0) {
  for(u_int i = 0; i < Nx; ++i) {
    Projection(i, f0, 0, sol[i]);
  }
  print_solution(std::cout);
  //std::cout << "Composition:" << std::endl;
  //std::cout << Composition(4, 0.8, 0) << std::endl;
  //std::cout << Composition(4, 0.9, 0) << std::endl;
  //std::cout << Composition(4, 1.0, 0) << std::endl;
}

double DGFEMSpace1D::cal_dt() {
  return 0.1;
}

int DGFEMSpace1D::forward_one_step(double dt, double* dtt) {
  return 0;
}

void DGFEMSpace1D::run(double t_end) {
  double t(0), dt(0), dtt(0);
  while (t < t_end) {
    dt = cal_dt();
    if(t+dt > t_end) dt = t_end-t;
    forward_one_step(dt, &dtt);
    t += dt;
  }
}

void DGFEMSpace1D::print_solution(std::ostream& os) {
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int k = 0; k < K; ++k) {
      os << sol[i][k] << "\n";
    }
    os << "\n";
  }
  os << std::endl;
}

