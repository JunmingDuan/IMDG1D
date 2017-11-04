/**
 * @file Quadrature.cpp
 * @brief Quadrature points and weights on [-1,1]
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include "Quadrature.h"

void Quadrature::set_jacobi(const double j) {
  local_to_global_jacobian = j;
}

double Quadrature::l2g_jacobian() {
  return local_to_global_jacobian;
}

void Quadrature::print(std::ostream& os) {
  os << "Number of quadrature points: " << np << std::endl;
  for(u_int i = 0; i < np; ++i) {
    os << pnt[i] << "\t";
  }
  os << "\n";
  for(u_int i = 0; i < np; ++i) {
    os << wei[i] << "\t";
  }
  os << "\n";
  os << local_to_global_jacobian;
  os << std::endl;
}

