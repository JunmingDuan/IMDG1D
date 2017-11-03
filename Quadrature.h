/**
 * @file Quadrature.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-03
 */

#ifndef QUADRTURE_H
#define QUADRTURE_H
#include <vector>

class Quadrature {
  private:
    std::vector<double> pnt;
    std::vector<double> wei;
  public:
    void local_to_global();

};

#endif //QUADRTURE_H

