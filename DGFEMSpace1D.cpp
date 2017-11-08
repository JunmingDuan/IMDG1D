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
  A.resize(Nx*K*DIM, Nx*K*DIM);
  rhs.resize(Nx*K*DIM);
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
    QUADINFO[i].set_jacobi( local_to_global_jacobian(lv, gv),
        global_to_local_jacobian(lv, gv) );
  }
  //for(u_int i = 0; i < Nx; ++i) {
    //std::cout << i << "-th cell QUADINFO" << "\n";
    //QUADINFO[i].print(std::cout);
  //}
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

int DGFEMSpace1D::forward_one_step(const F FLUX, afunc g, double t, double dt, double* dtt) {
  double alpha = 1;
  Newton_iter(sol, FLUX, g, t, dt, alpha);
  return 0;
}

VEC<double> DGFEMSpace1D::LxF(const F FLUX, const VEC<double>& a, const VEC<double>& b, const double alpha) {
  VEC<double> lf(a.size());
  lf = 0.5*(FLUX(a)+FLUX(b) - alpha*(b-a));
  return lf;
}

void DGFEMSpace1D::Newton_iter(SOL& sol, const F FLUX, const afunc g, const double t, const double dt, const double alpha) {
  form_jacobian_rhs(sol, FLUX, g, t, dt, alpha);
  solve_leqn(A, rhs);
}

int kronecker(const int a, const int b) {
  if(a == b) return 1;
  else return 0;
}

EVEC DGFEMSpace1D::NLF(const F FLUX, const SOL& sol, const SOL& u, const double alpha, const double t, const double dt) {
  EVEC fk(Nx*K*DIM);
  int row;
  std::vector<double> x = TemQuad.points();
  std::vector<double> pnt, wei;
  std::vector<double> PolyVal, PGVal;
  std::vector<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {
    pnt = QUADINFO[i].points();
    wei = QUADINFO[i].weight();
    for(u_int k = 0; k < K; ++k) {
      VEC<double> fu(DIM), fu1(DIM);
      for(u_int d = 0; d < DIM; ++d) {
        row = i*(K*DIM) + k*DIM + d;//fixed row
        //time derivative
        fk[row] = (sol[i][k][d]-u[i][k][d])/(2*k+1);
        //element integral
        for(u_int g = 0; g < pnt.size(); ++g) {
          PGVal = PolyG(x[g]);
          fu = Composition(i,pnt[g],t);
          fk[row] += - dt/2 * fu[d] * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
        }
        //
        gv[0] = mesh[i], gv[1] = mesh[i+1];
        VEC<double> flux;
        //flux on the outer boundary
        PolyVal = Poly(lv[1]);
        fu = Composition(i,gv[1],t);
        if(i < Nx-1) {
          fu1 = Composition(i+1,gv[1],t);
        }
        else fu1 = fu;
        flux = LxF(FLUX, fu, fu1, alpha);
        fk[row] += flux[d] * dt/(gv[1]-gv[0]) * PolyVal[k];

        //flux on the inner boundary
        fu = Composition(i,gv[0],t);
        if(i > 1) {
          fu1 = Composition(i-1,gv[0],t);
        }
        else {
          for(u_int d = 0; d < DIM; ++d)
            fu1[d] = 0;
        }
        flux = LxF(FLUX, fu1, fu, alpha);
        fk[row] -= flux[d] * dt/(gv[1]-gv[0]) * PolyVal[k];

      }
    }
  }
  return fk;
}

/**
 * @brief form_jacobian_rhs
 *        jacobian matrix is a (Nx*K) dimensional square matrix.
 *        rhs is a (Nx*K) dimensional vector, we order the vector
 *        in space firstly, then in polynomial space and last in physical space.
 *        rhs[i*(K*DIM)+k*DIM+d], is the d-th physical variable of the k-th
 *        polynomial in the i-the cell.
 *
 * @param sol
 */
void DGFEMSpace1D::form_jacobian_rhs(SOL& sol, const F FLUX, afunc fp, const double t, const double dt, const double alpha) {
  int row, col;
  double val;
  std::vector<double> pnt, wei;
  std::vector< T > List;
  std::vector<double> x = TemQuad.points();
  std::vector<double> PolyVal, PGVal;
  std::vector<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  for(u_int i = 0; i < Nx; ++i) {//i-th cell
    //std::cout << "i: " << i << std::endl;
    pnt = QUADINFO[i].points();
    wei = QUADINFO[i].weight();
    for(u_int k = 0; k < K; ++k) {//k-th phi
      //std::cout << "k: " << k << std::endl;
      for(u_int d = 0; d < DIM; ++d) {//d-th equation
        //std::cout << "d: " << d << std::endl;
        row = i*(K*DIM) + k*DIM + d;//fixed row
        ////time derivative: u_{i,d}^(k)
        col = row;
        val = 1./(2*k+1);
        List.push_back( T(row, col, val) );

        //element integral: u^(q)_{i,d}, q=0..K-1
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m}
          //std::cout << "q: " << q << std::endl;
          for(u_int m = 0; m < DIM; ++m) {
            //std::cout << "m: " << m << std::endl;
            val = 0;
            col = i*(K*DIM) + q*DIM + m;
            for(u_int g = 0; g < pnt.size(); ++g) {//numerical integral
              //std::cout << "g: " << g << std::endl;
              VEC<VEC<double> > AF = fp(Composition(i, pnt[g], 0));
              PolyVal = Poly(x[g]);
              PGVal = PolyG(x[g]);
              //only the m-th component of U, i.e.,
              //U_{i,m}=\sum u^(q)_{i,m}*poly^(q) has component u^(q)_{i,m}
              double pd = AF[d][m] * PolyVal[q];
              //scale in the prime of the polynomial
              val += pd * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
            }
            val *= -dt/2;
            List.push_back( T(row, col, val) );
          }
        }
        //fLux on the outer boundary
        gv[0] = mesh[i], gv[1] = mesh[i+1];
        PolyVal = Poly(lv[1]);
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m} and u^(q)_{i+1,m}
          for(u_int m = 0; m < DIM; ++m) {
            //i-th cell
            col = i*(K*DIM) + q*DIM + m;
            VEC<VEC<double> > AF = fp(Composition(i, gv[1], 0));
            double pd = 0.5 * (AF[d][m]+alpha*kronecker(d,m)) * PolyVal[q];
            val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            List.push_back( T(row, col, val) );

            //(i+1)-th cell
            if(i < Nx-1) {
              col = (i+1)*(K*DIM) + q*DIM + m;
              AF = fp(Composition(i+1, gv[1], 0));
              pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * PolyVal[q];
              val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
              List.push_back( T(row, col, val) );
            }
            else {//outflow right boundary
              //(Nx)-th cell is the same as (Nx-1)-the cell
              //we plus the coefficients to the col of (Nx-1)-th cell
              col = (i)*(K*DIM) + q*DIM + m;//modified here
              AF = fp(Composition(i, gv[1], 0));
              pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * PolyVal[q];
              val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
              List.push_back( T(row, col, val) );
            }
          }
        }

        //fLux on the inner boundary
        PolyVal = Poly(lv[0]);
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m} and u^(q)_{i-1,m}
          for(u_int m = 0; m < DIM; ++m) {
            VEC<VEC<double> > AF;
            double pd;
            //(i-1)-th cell
            if(i > 0) {
              col = (i-1)*(K*DIM) + q*DIM + m;
              AF= fp(Composition(i-1, gv[1], 0));
              pd = 0.5 * (AF[d][m]+alpha*kronecker(d,m)) * PolyVal[q];
              val = - pd * dt/(gv[1]-gv[0]) * PolyVal[k];
              List.push_back( T(row, col, val) );
            }
            else {//zero left boundary
            }

            //i-th cell
            col = i*(K*DIM) + q*DIM + m;
            AF = fp(Composition(i, gv[1], 0));
            pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * PolyVal[q];
            val = - pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            List.push_back( T(row, col, val) );
          }
        }

        //RHS
        rhs = -NLF(FLUX, sol, sol, alpha, t, dt);

      }
    }
  }
  A.setFromTriplets(List.begin(), List.end());
  std::cout << setiosflags(std::ios::fixed) << std::setprecision(3)
    << std::setiosflags(std::ios::right) << std::setw(5)
    << std::setiosflags(std::ios::showpos) << std::setfill(' ') << A << std::endl;
  std::cout << A.rows() << std::endl;
  std::cout << A.cols() << std::endl;
  std::cout << A.nonZeros() << std::endl;
  std::cout << std::defaultfloat;
}

void DGFEMSpace1D::solve_leqn(MAT& A, EVEC& rhs) {
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
  solver.compute(A);
  EVEC u = solver.solve(rhs);
  std::cout << "iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;
}

void DGFEMSpace1D::run(F FLUX, afunc g, double t_end) {
  double t(0), dt(0), dtt(0);
  while (t < t_end) {
    dt = cal_dt();
    if(t+dt > t_end) dt = t_end-t;
    forward_one_step(FLUX, g, t, dt, &dtt);
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
  os << std::defaultfloat;
}

