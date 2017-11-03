/**
 * @file Quadrature.cc
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#include "Quadrature.h"

std::vector<std::vector<double> > Quad_LGL(u_int np) {
  std::vector<std::vector<double> > pw;
  pw.resize(2);
  pw[0].resize(np);
  pw[1].resize(np);
  switch (np) {
    case 2: {
              pw[0][0] = -1, pw[0][1] = 1;
              pw[1][0] = 1, pw[1][1] = 1;
              break;
            }
    case 3: {
              pw[0][0] = -1, pw[0][1] = 0, pw[0][2] = 1;
              pw[1][0] = 1./3, pw[1][1] = 4./3, pw[1][2] = 1./3;
              break;
            }
    case 4: {
              pw[0][0] = -1, pw[0][1] = 0, pw[0][2] = 1;
              pw[1][0] = 1./3, pw[1][1] = 4./3, pw[1][2] = 1./3;
              break;
            }
  }
  return pw;
}


