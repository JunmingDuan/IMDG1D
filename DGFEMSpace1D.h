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
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

#define VEC vvector
typedef VEC<VEC<double> > bU;
typedef VEC<bU> SOL;
typedef int BM;
typedef std::vector<std::vector<double> > QUAD;
typedef VEC<double> (*func)(const double x, const double t);
typedef VEC<VEC<double> > (*afunc)(const VEC<double>&);
typedef VEC<double> (*F)(const VEC<double>&);
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> MAT;
typedef Eigen::VectorXd EVEC;

/**
 * @brief dimension of the equation, 1 for scalar equation and 3 for Euler equations
 */
const u_int DIM = 3;

class DGFEMSpace1D {
  private:
    u_int Nx;
    double xl, xr;
    double h;
    VEC<double> mesh;
    TemplateQuadrature TemQuad;
    std::vector<Quadrature> QUADINFO;
    SOL sol;
    BM bml, bmr;
    MAT A;
    EVEC rhs;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr);
    void BuildQuad(u_int np);
    void Projection(u_int cell, func f0, double t, bU&);
    VEC<double> Composition(u_int cell, double x, double t);
    void init(func f0);
    double cal_dt();
    /**
     * @brief forward_one_step
     *
     * @param double dt
     *
     * @return 0, donnot change dt; 1, change dt to dtt
     */
    int forward_one_step(const F, afunc g, double t, double dt, double* dtt);
    /**
     * @brief Newton_iter
     *
     * @param sol solution at t^n
     * @param dt
     */
    VEC<double> LxF(const F, const VEC<double>&, const VEC<double>&, const double);
    void Newton_iter(SOL& sol, const F, const afunc g, const double, const double, const double);
    /**
     * @brief NLF nonlinear function
     *
     * @param sol
     *
     * @return RHS, i.e., F(sol)
     */
    EVEC NLF(const F, const SOL& sol, const SOL& para_u, const double alpha, const double t, const double dt);
    void form_jacobian_rhs(SOL& sol, const F, afunc, const double, const double, const double);
    void solve_leqn(MAT& A, EVEC& rhs);
    void run(F, afunc g, double t_end);
    void print_solution(std::ostream&);
};

#endif //DGFEMSPACE1D_H

