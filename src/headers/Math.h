/*
 * This program is licensed granted by STATE UNIVERSITY OF CAMPINAS - UNICAMP ("University")
 * for use of MassCCS software ("the Software") through this website
 * https://github.com/cces-cepid/MassCCS (the "Website").
 *
 * By downloading the Software through the Website, you (the "License") are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about MassCCS please contact: 
 * skaf@unicamp.br (Munir S. Skaf)
 * guido@unicamp.br (Guido Araujo)
 * samuelcm@unicamp.br (Samuel Cajahuaringa)
 * danielzc@unicamp.br (Daniel L. Z. Caetano)
 * zanottol@unicamp.br (Leandro N. Zanotto)
 */

#ifndef MASSCCS_V1_MATH_H
#define MASSCCS_V1_MATH_H

#include <cmath>
#include <exception>
#include <omp.h>
#include <stdexcept>
#include <vector>

class Math {
public:
  template <typename T>
  static T dotProduct(std::vector<T> const &v1, std::vector<T> const &v2);

  template <typename T> static T vecModulus(std::vector<T> const &v1);

  template <typename T>
  __attribute__((vector, nothrow)) static T distance(std::vector<T> const &v1,
                                                     std::vector<T> const &v2);

  template <typename T>
  __attribute__((nothrow)) static T distance(std::vector<T> const *v1,
                                             std::vector<T> const *v2);

  template <typename T>
  static std::vector<T> crossProduct(std::vector<T> const &v1,
                                     std::vector<T> const &v2);

  template <typename T>
  static std::vector<T> addVec(std::vector<T> const &v1,
                               std::vector<T> const &v2);

  template <typename T>
  static std::vector<T> subVec(std::vector<T> const &v1,
                               std::vector<T> const &v2);

  template <typename T>
  static std::vector<T> normalize(std::vector<T> const &vector);
};

template <typename T>
T Math::dotProduct(std::vector<T> const &v1, std::vector<T> const &v2) {
  T ret{};

  for (unsigned int i = 0; i < v1.size(); i++)
    ret += v1.at(i) * v2.at(i);

  return ret;
}

template <typename T> T Math::vecModulus(std::vector<T> const &v1) {
  T ret{};
  for (double i : v1)
    ret += i * i;

  return sqrt(ret);
}

template <typename T>
std::vector<T> Math::crossProduct(std::vector<T> const &v1,
                                  std::vector<T> const &v2) {
  std::vector<T> ret{};

  // x
  ret.emplace_back(v1[1] * v2[2] - v1[2] * v2[1]);

  // y
  ret.emplace_back(v1[2] * v2[0] - v1[0] * v2[2]);

  // z
  ret.emplace_back(v1[0] * v2[1] - v1[1] * v2[0]);

  return ret;
}

#pragma omp declare simd
template <typename T>
T Math::distance(std::vector<T> const &v1, std::vector<T> const &v2) {
  double sum[3];
  double res = 0.0;

  for (unsigned int i = 0; i < 3; ++i)
    sum[i] = pow(v1[i] - v2[i], 2.0);

  for (unsigned int i = 0; i < 3; ++i)
    res += sum[i];

  return sqrt(res);
}

template <typename T>
std::vector<T> Math::subVec(std::vector<T> const &v1,
                            std::vector<T> const &v2) {
  std::vector<T> ret{};
  for (unsigned int i = 0; i < 3; ++i)
    ret.emplace_back(v1.at(i) - v2.at(i));

  return ret;
}

template <typename T>
std::vector<T> Math::addVec(std::vector<T> const &v1,
                            std::vector<T> const &v2) {
  std::vector<T> ret{};
  for (unsigned int i = 0; i < v1.size(); ++i)
    ret.emplace_back(v1.at(i) + v2.at(i));

  return ret;
}

template <typename T>
std::vector<T> Math::normalize(std::vector<T> const &vector) noexcept(false) {

  if (vector.size() == 0)
    throw std::range_error("MATH::Vector must have a valid size.");

  std::vector<T> ret{};

  T vecMod = vecModulus(vector);

  for (unsigned int i = 0; i < vector.size(); ++i)
    ret.emplace_back(vector.at(i) / vecMod);

  return ret;
}

#pragma omp declare simd
template <typename T>
T Math::distance(const std::vector<T> *v1, const std::vector<T> *v2) {
  double sum[3];
  double res = 0.0;

  for (unsigned int i = 0; i < 3; ++i)
    sum[i] = pow((*v1)[i] - (*v2)[i], 2.0);

  for (unsigned int i = 0; i < 3; ++i)
    res += sum[i];

  return sqrt(res);
}

#endif // MASSCCS_V1_MATH_H
