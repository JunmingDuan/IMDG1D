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

#define VEC vvector
typedef VEC<VEC<double> > bU;
typedef VEC<bU> SOL;
typedef int BM;
typedef std::vector<std::vector<double> > QUAD;
typedef VEC<double> (*func)(double x, double t);

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
    int forward_one_step(double dt, double* dtt);
    void run(double t_end);
    void print_solution(std::ostream&);
};

#endif //DGFEMSPACE1D_H

