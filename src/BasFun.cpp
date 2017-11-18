#include "BasFun.h"

u_int K = 5;

std::vector<double> Poly(double x) {
  std::vector<double> poly(K);
  poly[0] = 1;
  poly[1] = x;
  poly[2] = (3*x*x-1)/2;
  poly[3] = (5*x*x*x-3*x)/2;
  poly[4] = (35*x*x*x*x-30*x*x+3)/8;
  return poly;
}

/**
 * @brief PG Gradient of polynomial P.
 *
 * @param x
 *
 * @return A vector of length DIM.
 */
std::vector<double> PolyG(double x) {
  std::vector<double> poly(K);
  poly[0] = 0;
  poly[1] = 1;
  poly[2] = 3*x;
  poly[3] = (15*x*x-3)/2;
  poly[4] = (140*x*x*x-60*x)/8;
  return poly;
}

