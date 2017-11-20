/**
 * @file DGFEMSpace1D_GSL.cpp
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-19
 */

#include <stdio.h>
#include <stdlib.h>
#include "DGFEMSpace1D_GSL.h"

DGFEMSpace1D::DGFEMSpace1D(u_int Nx, double xl, double xr,
    double Nt_tol, double Nt_Ftol, double TOL)
  : Nx(Nx), xl(xl), xr(xr), Nt_tol(Nt_tol), Nt_Ftol(Nt_Ftol), TOL(TOL) {
  mesh.resize(Nx+1);
  h = (xr - xl)/Nx;
  for(u_int i = 0; i < Nx+1; ++i) {
    mesh[i] = i*h;
  }
  sol.resize(Nx);
  sol1.resize(Nx);
  for(u_int i = 0; i < Nx; ++i) {
    sol[i].resize(K);
    sol1[i].resize(K);
    for(u_int k = 0; k < K; ++k) {
      sol[i][k].resize(DIM);
      sol1[i][k].resize(DIM);
    }
  }
  A = gsl_spmatrix_alloc(Nx*K*DIM, Nx*K*DIM);
  rhs = gsl_vector_alloc(Nx*K*DIM);
  vec_u1 = gsl_vector_alloc(Nx*K*DIM);
  vec_u2 = gsl_vector_alloc(Nx*K*DIM);
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
}

void DGFEMSpace1D::Projection(u_int cell, func f0, double t, bU& u) {
  //set to zero
  for(u_int k = 0; k < K; ++k) {
    for(u_int d = 0; d < DIM; ++d) {
      u[k][d] = 0;
    }
  }
  std::vector<double> x = TemQuad.points();
  std::vector<double> p = QUADINFO[cell].points();
  std::vector<double> w = QUADINFO[cell].weight();
  VEC<double> U(DIM), u0(DIM);
  std::vector<double> V;
  double jab = QUADINFO[cell].l2g_jacobian();
  for(u_int k = 0; k < K; ++k) {
    double basis(0);
    for(u_int g = 0; g < x.size(); ++g) {
      U = f0(u0, p[g], t);
      V = Poly(x[g]);
      u[k] += jab*U*V[k]*w[g];
      basis += jab*V[k]*V[k]*w[g];
    }
    u[k] /= basis;
  }
}

VEC<double> DGFEMSpace1D::Composition(const SOL& sol, u_int cell, double x, double t) {
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
  //print_solution(std::cout);
}

double DGFEMSpace1D::cal_dt() {
  return h;
}

int DGFEMSpace1D::forward_one_step(const SOL& sol, const F FLUX, afunc g, func source,
    double t, double dt, double* dtt, SOL& sol_new) {
  double alpha = 1;
  Newton_iter(sol, FLUX, g, source, t, dt, alpha, sol_new);
  //std::cout << sol_new << std::endl;
  return 0;
}

VEC<double> DGFEMSpace1D::LxF(const F FLUX, const VEC<double>& a, const VEC<double>& b, const double alpha) {
  VEC<double> lf(a.size());
  lf = 0.5*(FLUX(a)+FLUX(b) - alpha*(b-a));
  return lf;
}

void DGFEMSpace1D::Newton_iter(const SOL& sol, const F FLUX, const afunc g, func source,
    const double t, const double dt, const double alpha, SOL& sol_new) {
  int Nt_ite(0);
  double Nt_err(1), Fval_norm(1);
  sol_new = sol;
  SOL2EVEC(sol_new, vec_u1);
  while (Nt_err > Nt_tol && Fval_norm > Nt_Ftol) {
    form_jacobian_rhs(sol_new, sol, FLUX, g, source, t, dt, alpha);
    Nt_err = gsl_blas_dnrm2(vec_u2);
    //std::cout << "=====sol^n,A,rhs,sol^{n+1}=====" << std::endl;
    //std::cout << "vec_u1:" << std::endl;
    //gsl_vector_fprintf(stdout, vec_u1, "%.6lf");
    //std::cout << "A:" << std::endl;
    //gsl_spmatrix_fprintf(stdout, A, "%.6lf");
    //std::cout << "rhs:" << std::endl;
    //gsl_vector_fprintf(stdout, rhs, "%.6lf");
    //std::cout << "rhs_norm:" << std::endl;
    //std::cout << gsl_blas_dnrm2(rhs) << std::endl;
    solve_leqn(A, rhs, vec_u2);
    gsl_vector_add(vec_u1, vec_u2);
    //std::cout << "vec_u1:" << std::endl;
    //gsl_vector_fprintf(stdout, vec_u1, "%.6lf");
    EVEC2SOL(sol_new, vec_u1);
    gsl_vector * tmp = gsl_vector_alloc(Nx*K*DIM);
    *tmp = NLF(FLUX, sol_new, sol, source, alpha, t, dt);
    Fval_norm = gsl_blas_dnrm2(tmp);
    Nt_ite++;
    std::cout << "Nt_ite: " << Nt_ite
      << ", Nt_err: " << Nt_err
      << ", Fval: " << Fval_norm
      << std::endl;
  }
}

int kronecker(const int a, const int b) {
  if(a == b) return 1;
  else return 0;
}

EVEC DGFEMSpace1D::NLF(const F FLUX, const SOL& sol, const SOL& soln, func source,
    const double alpha, const double t, const double dt) {
  EVEC * fk = gsl_vector_alloc(Nx*K*DIM);
  int row;
  std::vector<double> x = TemQuad.points();
  std::vector<double> pnt, wei;
  std::vector<double> PolyVal, PGVal;
  std::vector<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  gsl_vector_set_zero(fk);
  bU tmp_u(K);
  for(u_int k = 0; k < K; ++k) {
    tmp_u[k].resize(DIM);
  }
  for(u_int i = 0; i < Nx; ++i) {
    pnt = QUADINFO[i].points();
    wei = QUADINFO[i].weight();
    Projection(i, source, t, tmp_u);//projection of the source term
    for(u_int k = 0; k < K; ++k) {
      VEC<double> fu(DIM), fu1(DIM);
      VEC<double> U(DIM), U1(DIM);
      for(u_int d = 0; d < DIM; ++d) {
        row = i*(K*DIM) + k*DIM + d;//fixed row
        //time derivative
        gsl_vector_set(fk, row, (sol[i][k][d]-soln[i][k][d])/(2*k+1));
        //element integral
        for(u_int g = 0; g < pnt.size(); ++g) {
          PGVal = PolyG(x[g]);
          fu = FLUX(Composition(sol,i,pnt[g],t));
          //fk[row] += - dt/2 * fu[d] * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
          *(fk->data+row) += - dt/2 * fu[d] * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
        }
        //
        gv[0] = mesh[i], gv[1] = mesh[i+1];
        VEC<double> flux;
        //flux on the outer boundary
        PolyVal = Poly(lv[1]);
        U = Composition(sol,i,gv[1],t);
        if(i < Nx-1) {
          U1 = Composition(sol,i+1,gv[1],t);
        }
        else U1 = U;
        flux = LxF(FLUX, U, U1, alpha);
        *(fk->data+row) += flux[d] * dt/(gv[1]-gv[0]) * PolyVal[k];

        //flux on the inner boundary
        PolyVal = Poly(lv[0]);
        U = Composition(sol,i,gv[0],t);
        if(i > 0) {
          U1 = Composition(sol,i-1,gv[0],t);
        }
        else {
          for(u_int d = 0; d < DIM; ++d)
            U1[d] = 0;
        }
        flux = LxF(FLUX, U1, U, alpha);
        *(fk->data+row) -= flux[d] * dt/(gv[1]-gv[0]) * PolyVal[k];

        *(fk->data+row) -= tmp_u[k][d] * dt/(2*k+1);//source  term
      }
    }
  }
  return *fk;
}

/**
 * @brief form_jacobian_rhs
 *        jacobian matrix is a (Nx*K) dimensional square matrix.
 *        rhs is a (Nx*K) dimensional vector, we order the vector
 *        in space firstly, then in polynomial space and last in physical space.
 *        rhs[i*(K*DIM)+k*DIM+d], is the d-th physical variable of the k-th
 *        polynomial in the i-the cell.
 *
 * @param sol current sol
 *        soln sol at t^n
 */
void DGFEMSpace1D::form_jacobian_rhs(const SOL& sol, const SOL& soln, const F FLUX, afunc fp, func source,
    const double t, const double dt, const double alpha) {
  int row, col;
  double val;
  double * ptr;
  std::vector<double> pnt, wei;
  std::vector<double> x = TemQuad.points();
  std::vector<double> PolyVal, PGVal, LocPolyVal;
  std::vector<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  gsl_spmatrix_set_zero(A);
  for(u_int i = 0; i < Nx; ++i) {//i-th cell
    //std::cout << "i: " << i << std::endl;
    pnt = QUADINFO[i].points();
    wei = QUADINFO[i].weight();
    for(u_int k = 0; k < K; ++k) {//k-th phi
      //std::cout << "k: " << k << std::endl;
      for(u_int d = 0; d < DIM; ++d) {//d-th equation
        //std::cout << "d: " << d << std::endl;
        row = i*(K*DIM) + k*DIM + d;//fixed row
        //time derivative: u_{i,d}^(k)
        col = row;
        val = 1./(2*k+1);
        ptr = gsl_spmatrix_ptr(A, row, col);
        if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
        else *ptr += val;

        //element integral: u^(q)_{i,d}, q=0..K-1
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m}
          //std::cout << "q: " << q << std::endl;
          for(u_int m = 0; m < DIM; ++m) {
            //std::cout << "m: " << m << std::endl;
            val = 0;
            col = i*(K*DIM) + q*DIM + m;
            for(u_int g = 0; g < pnt.size(); ++g) {//numerical integral
              //std::cout << "g: " << g << std::endl;
              VEC<VEC<double> > AF = fp(Composition(sol,i, pnt[g], 0));
              PolyVal = Poly(x[g]);
              PGVal = PolyG(x[g]);
              //only the m-th component of U, i.e.,
              //U_{i,m}=\sum u^(q)_{i,m}*poly^(q) has component u^(q)_{i,m}
              double pd = AF[d][m] * PolyVal[q];
              //scale in the prime of the polynomial
              val += pd * (PGVal[k]*QUADINFO[i].g2l_jacobian()) * wei[g];
            }
            val *= -dt/2;
            //gsl_spmatrix_set(A, row, col, val);
            ptr = gsl_spmatrix_ptr(A, row, col);
            if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
            else *ptr += val;
          }
        }
        //fLux on the outer boundary
        gv[0] = mesh[i], gv[1] = mesh[i+1];
        PolyVal = Poly(lv[1]);
        for(u_int q = 0; q < K; ++q) {//derivative of u^(q)_{i,m} and u^(q)_{i+1,m}
          for(u_int m = 0; m < DIM; ++m) {
            //i-th cell
            LocPolyVal = Poly(lv[1]);
            col = i*(K*DIM) + q*DIM + m;
            VEC<VEC<double> > AF = fp(Composition(sol,i, gv[1], 0));
            double pd = 0.5 * (AF[d][m]+alpha*kronecker(d,m)) * LocPolyVal[q];
            val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            //List.push_back( T(row, col, val) );
            ptr = gsl_spmatrix_ptr(A, row, col);
            if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
            else *ptr += val;

            //(i+1)-th cell
            if(i < Nx-1) {
              LocPolyVal = Poly(lv[0]);
              col = (i+1)*(K*DIM) + q*DIM + m;
              AF = fp(Composition(sol,i+1, gv[1], 0));
              pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * LocPolyVal[q];
              val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
              //List.push_back( T(row, col, val) );
              ptr = gsl_spmatrix_ptr(A, row, col);
              if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
              else *ptr += val;
            }
            else {//outflow right boundary
              LocPolyVal = Poly(lv[1]);
              //(Nx)-th cell is the same as (Nx-1)-the cell
              //we plus the coefficients to the col of (Nx-1)-th cell
              col = (i)*(K*DIM) + q*DIM + m;//modified here
              AF = fp(Composition(sol,i, gv[1], 0));
              pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * LocPolyVal[q];
              val = pd * dt/(gv[1]-gv[0]) * PolyVal[k];
              //List.push_back( T(row, col, val) );
              ptr = gsl_spmatrix_ptr(A, row, col);
              if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
              else *ptr += val;
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
            LocPolyVal = Poly(lv[1]);
            if(i > 0) {
              col = (i-1)*(K*DIM) + q*DIM + m;
              AF= fp(Composition(sol,i-1, gv[0], 0));
              pd = 0.5 * (AF[d][m]+alpha*kronecker(d,m)) * LocPolyVal[q];
              val = - pd * dt/(gv[1]-gv[0]) * PolyVal[k];
              //List.push_back( T(row, col, val) );
              ptr = gsl_spmatrix_ptr(A, row, col);
              if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
              else *ptr += val;
            }
            else {//zero left boundary
            }

            //i-th cell
            LocPolyVal = Poly(lv[0]);
            col = i*(K*DIM) + q*DIM + m;
            AF = fp(Composition(sol,i, gv[0], 0));
            pd = 0.5 * (AF[d][m]-alpha*kronecker(d,m)) * LocPolyVal[q];
            val = - pd * dt/(gv[1]-gv[0]) * PolyVal[k];
            //List.push_back( T(row, col, val) );
            ptr = gsl_spmatrix_ptr(A, row, col);
            if(ptr == NULL) gsl_spmatrix_set(A, row, col, val);
            else *ptr += val;
          }
        }

      }
    }
  }
  //gsl_spmatrix_ccs(A);
  //RHS
  *rhs = NLF(FLUX, sol, soln, source, alpha, t, dt);
  gsl_vector_scale(rhs, -1);

}

void DGFEMSpace1D::solve_leqn(MAT *A, const EVEC *rhs, EVEC *u) {
  std::cout << "======solve_leqn by GSL GMRES======" << std::endl;
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  gsl_splinalg_itersolve *work =
    gsl_splinalg_itersolve_alloc(T, Nx*K*DIM, 0);
  size_t iter = 0;
  double residual, tol(1e-14);
  int status;
  /* initial guess u = 0 */
  gsl_vector_set_zero(u);
  /* solve the system A u = f */
  do
  {
    status = gsl_splinalg_itersolve_iterate(A, rhs, tol, u, work);
    /* print out residual norm ||A*u - f|| */
    residual = gsl_splinalg_itersolve_normr(work);
    //fprintf(stdout, "iter %zu residual = %.12e\n", iter++, residual);

    if (status == GSL_SUCCESS)
      fprintf(stdout, "Converged, residual%.6e\n", residual);
  }
  while (status == GSL_CONTINUE);
  std::cout << "===================================" << std::endl;
}

VEC<double> exact0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  U[0] = (sin(4*x)-8*sin(2*x)+12*x)/32;
  return U;
}

void DGFEMSpace1D::run(F FLUX, afunc g, func source, double t_end) {
  int ite(0), pp(1);
  double t(0), dt(0), dtt(0);
  VEC<double> err(DIM,1), tol(DIM,TOL);
  while ( err > tol ) {//|| pp == 0 ) {
    dt = cal_dt();
    forward_one_step(sol, FLUX, g, source, t, dt, &dtt, sol1);
    pp = judge_positivity(sol1);
    //while (pp == 0) {
      //std::cout << "Not positive!!" << std::endl;
      //dt = 0.5*dt;
      //forward_one_step(sol, FLUX, g, t, dt, &dtt, sol1);
      //pp = judge_positivity(sol1);
    //}
    err = cal_norm(sol, sol1, 2);
    sol = sol1;
    ite++;
    t += dt;
    std::cout << "ite: " << ite << ", dt: " << dt << ", t: " << t << ", err: ";
    std::cout << err << std::endl;
  }
  //use exact solution to form the lneq
  //SOL tmp1, tmp2;
  //for(u_int i = 0; i < Nx; ++i) {
    //tmp1.resize(Nx); tmp2.resize(Nx);
    //for(u_int k = 0; k < K; ++k) {
      //tmp1[i].resize(K); tmp2[i].resize(K);
      //for(u_int d = 0; d < DIM; ++d) {
        //tmp1[i][k].resize(DIM);
        //tmp2[i][k].resize(DIM);
      //}
    //}
  //}
  //for(u_int i = 0; i < Nx; ++i) {
    //Projection(i, exact0, 0, tmp1[i]);
    //tmp1[i] = tmp2[i];
  //}
  //form_jacobian_rhs(tmp1, tmp2, FLUX, g, source, t, dt, 1);
  //gsl_vector * tmp3 = gsl_vector_alloc(Nx*K*DIM);
  //*rhs = NLF(FLUX, tmp1, sol, source, 1, t, dt);
  //gsl_vector_fprintf(stdout, tmp3, "%.15e");
}

int DGFEMSpace1D::judge_positivity(const SOL& sol) {
  double eps(1e-13);
  std::vector<double> p;
  VEC<double> average(DIM);
  for(u_int i = 0; i < Nx; ++i) {
    p = QUADINFO[i].points();
    for(u_int g = 0; g < p.size(); ++g) {
      average += Composition(sol,i,p[g],0);
    }
    for(u_int d = 0; d < DIM; ++d) {
      if(average[d] < eps) return 0;
    }
  }
  return 1;
}

void DGFEMSpace1D::SOL2EVEC(const SOL& sol, EVEC *vec_u) {
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int k = 0; k < K; ++k) {
      for(u_int d = 0; d < DIM; ++d) {
        gsl_vector_set(vec_u, i*(K*DIM)+k*DIM+d, sol[i][k][d]);
      }
    }
  }
}

void DGFEMSpace1D::EVEC2SOL(SOL& sol, const EVEC *vec_u) {
  for(u_int i = 0; i < Nx; ++i) {
    for(u_int k = 0; k < K; ++k) {
      for(u_int d = 0; d < DIM; ++d) {
        sol[i][k][d] = gsl_vector_get(vec_u, i*(K*DIM)+k*DIM+d);
      }
    }
  }
}

VEC<double> DGFEMSpace1D::cal_norm(const SOL& s1, const SOL& s2, int n) {
  VEC<double> norm(DIM), tmp;
  double center;
  if(n == 2) {
    for(u_int i = 0; i < Nx; ++i) {
      center = 0.5*(mesh[i]+mesh[i+1]);
      tmp = Composition(s1,i,center,0)-Composition(s2,i,center,0);
      for(u_int d = 0; d < DIM; ++d) {
        norm[d] += pow(tmp[d], 2);
      }
    }
    for(u_int d = 0; d < DIM; ++d) {
      norm[d] = sqrt(norm[d]*h);
    }
    return norm;
  }
  else {
    std::cout << "Wrong norm choice!" << std::endl;
    return norm;
  }
}

void DGFEMSpace1D::print_solution(std::ostream& os) {
	os.precision(16);
	os << std::showpos;
  os.setf(std::ios::scientific);
  std::vector<double> x = TemQuad.points();
  for(u_int i = 0; i < Nx; ++i) {
    std::vector<double> w = QUADINFO[i].weight();
    std::vector<double> p = QUADINFO[i].points();
    for(u_int g = 0; g < x.size(); ++g) {
      os << p[g] << " "  << w[g] << " " << Composition(sol,i,p[g],0) << "\n";
    }
    os << "\n";
  }
  os << std::endl;
  os << std::defaultfloat;
  //os.setf(std::ios::floatfield);
}


