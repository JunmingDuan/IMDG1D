//#include <cstdlib>
#include "DGFEMSpace1D.h"

VEC<double> f0(const double x, const double t) {
  VEC<double> u(DIM);
  u[0] = 1;//pow(sin(x), 2);
  u[1] = 2*x;//pow(sin(x), 2);
  u[2] = 3*x;//pow(sin(x), 2);
  return u;
}

VEC<double> f(const VEC<double>& u) {
  VEC<double> F(DIM);
  F[0] = u[0];
  F[1] = 2*u[1];
  F[2] = 3*u[2];
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
  std::cout << "Set up problem..." << std::endl;
  DGFEMSpace1D Problem(Nx, xl, xr);
  std::cout << "Build quadrature info..." << std::endl;
  Problem.BuildQuad(2);
  std::cout << "Initialize..." << std::endl;
  Problem.init(f0);

  std::cout << "########## DEBUG jacobian #########" << std::endl;
  SOL x1, x2, xn;
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
  double eps(1e-6), a(2), dt(0.1);
  int i1, k1, d1, i2, k2, d2;
  for(u_int loop = 0; loop < 1e3; ++loop) {
    i1 = 1, k1 = 0, d1 = 1;
    i2 = 0, k2 = 1, d2 = 2;
    //srand((unsigned)time(NULL));
    //i1 = rand()%Nx, k1 = rand()%K, d1 = rand()%DIM;
    //i2 = rand()%Nx, k2 = rand()%K, d2 = rand()%DIM;
    int row = i1*(K*DIM)+k1*DIM+d1;
    int col = i2*(K*DIM)+k2*DIM+d2;
    x1[i2][k2][d2] = x2[i2][k2][d2]+eps;
    EVEC DIFF = (Problem.NLF(f, x1, xn, a, 0, dt) - Problem.NLF(f, x2, xn, a, 0, dt))/eps;
    Problem.form_jacobian_rhs(x1, f, f_prime, 0, dt, a);
    double tmp = DIFF[row] - Problem.get_A().coeffRef(row, col);
    std::cout << "DIFF-A_ij: " << tmp << std::endl;
    std::cout << "row,col: " << row << " " << col << std::endl;
    std::cout << "DIFF[row]: " << DIFF[row] << " ,A_ij: " << Problem.get_A().coeffRef(row, col) << std::endl;
    if(fabs(tmp) > 1e-6) {
      std::cout << "row,col: " << row << " " << col << std::endl;
      std::cout << "i1,k1,d1: " << i1 << " " << k1 << " " << d1 << std::endl;
      std::cout << "i2,k2,d2: " << i2 << " " << k2 << " " << d2 << std::endl;
      std::cout << "DIFF[row]: " << DIFF[row] << " ,A_ij: " << Problem.get_A().coeffRef(row, col) << std::endl;
      abort();
    }
  }

  std::cout << "########## DEBUG jacobian ##########" << std::endl;
  return 0;
}

