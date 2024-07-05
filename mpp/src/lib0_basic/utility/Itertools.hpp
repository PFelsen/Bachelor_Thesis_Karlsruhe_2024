#ifndef ITERTOOLS_HPP
#define ITERTOOLS_HPP

#include <vector>
#include "Point.hpp"

namespace it {

template<std::size_t N>
class Ranges {
  int r[N];
  int state[N];
public:
  Ranges(int (&r)[N], int (&state)[N]) {
    std::copy(r, r + N, this->r);
    std::copy(state, state + N, this->state);
  }

  Ranges(int (&&r)[N]) {
    r[0]++;
    std::copy(r, r + N, this->r);
    std::fill(state, state + N, 0);
  }

  Ranges begin() { return *this; }

  Ranges end() {
    int rr[N];
    rr[0] = r[0] - 1;
    for (int i = 1; i < N; i++) {
      rr[i] = 0;
    }
    return Ranges(r, rr);
  }

  const auto &operator*() { return state; }

  Ranges operator++() {
    for (int i = N - 1; i >= 0; i--) {
      state[i]++;
      if (state[i] != r[i]) { break; }
      state[i] = 0;
    }
    return *this;
  }

  bool operator==(Ranges &other) {
    for (int i = 0; i < N; i++) {
      if (state[i] != other.state[i]) { return false; }
    }
    return true;
  }

  bool operator!=(Ranges &other) { return !this->operator==(other); }
};

class QuadratureIterator {

  int q;
  std::vector<Point> &qPoint;
  std::vector<double> &qWeight;
public:
  QuadratureIterator(std::vector<Point> &qPoint, std::vector<double> &qWeight, int q = 0) :
      qPoint(qPoint), qWeight(qWeight), q(q) {}

  QuadratureIterator begin() { return *this; }

  QuadratureIterator end() { return QuadratureIterator(qPoint, qWeight, qPoint.size()); }

  const std::tuple<int, Point, double> operator*() { return {q, qPoint[q], qWeight[q]}; }

  QuadratureIterator &operator++() {
    q++;
    return *this;
  }

  bool operator==(const QuadratureIterator qi) const { return qi.q == q; }

  bool operator!=(const QuadratureIterator qi) const { return qi.q != q; }
};
} // namespace it

#endif // ITERTOOLS_HPP
