/**
 * @file Quadrature.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef QUADRTURE_H
#define QUADRTURE_H

class QuadratureInfo {
  private:
    std::vector<double> x;
    std::vector<double> pnt;
    std::vector<double> wei;
  public:
    std::vector<std::vector<double> > LGL(u_int np);

};

#endif //QUADRTURE_H

