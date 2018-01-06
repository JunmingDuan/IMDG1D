#ifndef PARA_H_
#define PARA_H_

#define VEC vvector
#define bU VEC<VEC<double>>
#define SOL VEC<bU>

//dimension of the equation, 1 for scalar equation and 3 for Euler equations
const u_int DIM = 3;
//number of basis function
const u_int K = 2;
//0 for ghost = 0, 1 for flux = 0, 2 for period BD
//ex9
const u_int BDL = 1; const u_int BDR = 1;
//positivity preserving limiter
const u_int PP_limiter = 1;
//parameters for Newton iteration
const int MaxNt_ite = 1e1;
const double Nt_tol = 1e-14;
const double Nt_Ftol = 1e-14;
//tol for linear equation solver
//const double tol = 1e-13;
const double tol = 1e-14;
//tol for steady solution solver
const double TOL = 1e-12;
//eps for scaling_limiter
const double EPS = 1e-13;
//GAMMA for Euler equations
const double GAMMA = 1.4;

#endif //PARA_H_

