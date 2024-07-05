#ifndef RMATRICES_HPP
#define RMATRICES_HPP

#include "RMatrix.hpp"

template<typename REAL = double>
class RMatricesT {
  std::vector<RMatrixT<REAL>> matrices{};
public:
  RMatricesT(int size) : matrices(size) {}

  void resize(int idx, int rows, int cols) { matrices[idx].resize(rows, cols); }

  void resize(int size) { matrices.resize(size); }

  int size() const { return matrices.size(); }

  RMatrixT<REAL> &operator[](int s) { return matrices[s]; }

  const RMatrixT<REAL> &operator[](int s) const { return matrices[s]; }

  int rows() const {
    int n = 0;
    for (int s = 0; s < matrices.size(); ++s) {
      if (matrices[s].rows() == 0) continue;
      if (n == 0) {
        n = matrices[s].rows();
      } else if (n != matrices[s].rows()) {
        THROW("wrong dim!") // TODO: this should be checked in advance
      }
    }
    return n;
  }

  int cols(int idx) const { return matrices[idx].cols(); }

  int cols() const {
    int n = 0;
    for (int s = 0; s < matrices.size(); ++s)
      n += cols(s);
    return n;
  }

  void set(int idx, const RMatrixT<REAL> &B) { matrices[idx] = B; }

  RMatrixT<REAL> Combine() const {
    RMatrixT<REAL> A(rows(), cols());
    int j = 0;
    for (int s = 0; s < size(); ++s)
      for (int k = 0; k < cols(s); ++j, ++k)
        for (int i = 0; i < rows(); ++i)
          A(i, j) = matrices[s][i][k];
    return A;
  }
};

using RMatrices = RMatricesT<double>;

#endif // RMATRICES_HPP
