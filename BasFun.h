/**
 * @file BasFun.h
 * @brief Basis function on [-1,1], Legendre polynomial not normalized.
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2017-11-04
 */

#ifndef BASFUN_H
#define BASFUN_H

double P0(double x) {
  return 1;
}
double P1(double x) {
  return x;
}
double P2(double x) {
  return (3*x*x-1)/2;
}
double P3(double x) {
  return (5*x*x*x-3*x)/2;
}
double P4(double x) {
  return (35*x*x*x*x-30*x*x+3)/8;
}

double gradient_0(double x) {
  return 0;
}
double gradient_1(double x) {
  return 1;
}
double gradient_2(double x) {
  return 3*x;
}
double gradient_3(double x) {
  return (15*x*x-3)/2;
}
double gradient_4(double x) {
  return (140*x*x*x-60*x)/8;
}

#endif //BASFUN_H

