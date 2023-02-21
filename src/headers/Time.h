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

#ifndef MASSCCS_V1_TIME_H
#define MASSCCS_V1_TIME_H

#include <chrono>
#include <iostream>

class Time {

public:
  template <typename TFunc>
  static void RunAndMeasure(const char *title, TFunc func)
      __attribute__((nothrow));
};

template <typename TFunc>
void Time::RunAndMeasure(const char *title, TFunc func) {
  const auto start = std::chrono::steady_clock::now();
  func();
  const auto end = std::chrono::steady_clock::now();
  std::cout << title << ": "
            << std::chrono::duration<double>(end - start).count() << " seconds"
            << std::endl;
}
#endif // MASSCCS_V1_TIME_H
