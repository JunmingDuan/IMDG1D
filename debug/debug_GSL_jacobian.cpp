//#include <cstdlib>
#include "DGFEMSpace1D_GSL.h"

VEC<double> f0(const VEC<double>& u, double x, double t) {
  VEC<double> U(DIM);
  U[0] = x;
  return U;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  F[0] = 0.5*pow(u[0],2);
  //F[0] = u[0];
  return F;
}

VEC<double> source(const VEC<double>& u, double x, double t) {
  VEC<double> F(DIM);
  F[0] = sin(x/4.);
  return F;
}

VEC<VEC<double> > f_prime(const VEC<double>& u) {
  VEC<VEC<double> > a(DIM);
  for(u_int i = 0; i < DIM; ++i)
    a[i].resize(DIM, 0);
  a[0][0] = u[0];
  //a[0][0] = 1;
  return a;
}

//debug jacobian by differential F
//verifying the jacobian matrix
int main(int argc, char* argv[]) {
  if(argc < 4) {
    std::cout << "Usage: <Nx> <xl> <xr> " << std::endl;
    abort();
  }
  u_int Nx = atoi(argv[1]);
  double xl = atof(argv[2]);
  double xr = atof(argv[3]);
  double Nt_tol(1e-14), Nt_Ftol(1e-14), TOL(1e-15);

  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr, Nt_tol, Nt_Ftol, TOL);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(K+1);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);

  std::cout << "########## DEBUG jacobian #########" << std::endl;
  SOL x1, x2, xn;
  double eps(1e-8), a(0), dt(0.1);
  int i1, k1, d1, i2, k2, d2;
  for(u_int loop = 0; loop < 1e4; ++loop) {
    double haha;
    //std::cin >> haha;
    //i1 = 1, k1 = 0, d1 = 0;
    //i2 = 0, k2 = 1, d2 = 2;
    srand((unsigned)time(NULL));
    std::cin.get();
    a = rand()/(double)RAND_MAX;
    i1 = rand()%Nx, k1 = rand()%K, d1 = rand()%DIM;
    i2 = rand()%Nx, k2 = rand()%K, d2 = rand()%DIM;
    int row = i1*(K*DIM)+k1*DIM+d1;
    int col = i2*(K*DIM)+k2*DIM+d2;
    //initialize xn, x1, x2
    x1.resize(Nx); x2.resize(Nx); xn.resize(Nx);
    for(u_int i = 0; i < Nx; ++i) {
      x1[i].resize(K); x2[i].resize(K); xn[i].resize(K);
      for(u_int k = 0; k < K; ++k) {
        x1[i][k].resize(DIM); x2[i][k].resize(DIM); xn[i][k].resize(DIM);
        for(u_int d = 0; d < DIM; ++d) {
          xn[i][k][d] = rand()/(double)RAND_MAX;
          x2[i][k][d] = rand()/(double)RAND_MAX;
          x1[i][k][d] = x2[i][k][d];
        }
      }
    }
    x1[i2][k2][d2] = x2[i2][k2][d2]+eps;

    EVEC r1 = Problem.NLF(f, x1, xn, source, a, 0, dt);
    EVEC r2 = Problem.NLF(f, x2, xn, source, a, 0, dt);
    Problem.form_jacobian_rhs(x2, xn, f, f_prime, source, 0, dt, a);
    MAT B = Problem.get_A();

    double tmp = (r1.data[row]-r2.data[row])/eps - gsl_spmatrix_get(&B, row, col);
    std::cout << "DIFF-A_ij: " << tmp << std::endl;
    std::cout << "row,col: " << row << " " << col << std::endl;
    std::cout << "DIFF[row]: " << (r1.data[row]-r2.data[row])/eps << " ,A_ij: " << gsl_spmatrix_get(&B, row, col) << std::endl;
    if(fabs(tmp) > 1e-6) {
      std::cout << "row,col: " << row << " " << col << std::endl;
      std::cout << "i1,k1,d1: " << i1 << " " << k1 << " " << d1 << std::endl;
      std::cout << "i2,k2,d2: " << i2 << " " << k2 << " " << d2 << std::endl;
      //std::cout << "NLF(x1, xn): " << Problem.NLF(f, x1, xn, a, 0, dt)[row]
        //<< " NLF(x2, xn): " << Problem.NLF(f, x2, xn, a, 0, dt)[row] << std::endl;
      //std::cout << "DIFF[row]: " << DIFF[row] << " ,A_ij: " << Problem.get_A().coeffRef(row, col) << std::endl;
      std::cout << "xn: " << xn << "\n";
      std::cout << "x1: " << x1 << "\n";
      std::cout << "x2: " << x2 << std::endl;
      abort();
    }
  }

  std::cout << "########## DEBUG jacobian ##########" << std::endl;

  return 0;
}

