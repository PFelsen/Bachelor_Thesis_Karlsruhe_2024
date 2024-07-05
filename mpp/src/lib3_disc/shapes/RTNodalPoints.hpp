#ifndef RTNODALPOINTS_HPP
#define RTNODALPOINTS_HPP

#include "ICell.hpp"

// TODO: Move to NodalPointProvider

template<typename T, int sDim, int tDim>
void RTNodalPointsTriangle(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int index, int N,
                           int currentN, int numberOfRekCalls) {
  if (currentN == 0) { z[index++] = c(); }
  for (int n = (numberOfRekCalls == 0) ? 1 : 0; n < currentN; n++) {
    z[index++] = c[0] + ((c[1] - PointT<T, sDim, tDim>(c[0])) * (T(n) / N))
                 + (c[1] + PointT<T, sDim, tDim>(c[2]) - 2 * PointT<T, sDim, tDim>(c[0]))
                       * (T(numberOfRekCalls) / N);
    z[index++] = c[1] + ((c[2] - PointT<T, sDim, tDim>(c[1])) * (T(n) / N))
                 + (c[2] + PointT<T, sDim, tDim>(c[0]) - 2 * PointT<T, sDim, tDim>(c[1]))
                       * (T(numberOfRekCalls) / N);
    z[index++] = c[2] + ((c[0] - PointT<T, sDim, tDim>(c[2])) * (T(n) / N))
                 + (c[1] + PointT<T, sDim, tDim>(c[0]) - 2 * PointT<T, sDim, tDim>(c[2]))
                       * (T(numberOfRekCalls) / N);
  }
  if (index != ((N + 1) * (N + 2)) / 2 - 3) {
    RTNodalPointsTriangle(c, z, index, N, currentN - 3, ++numberOfRekCalls);
  }
}

template<typename T, int sDim, int tDim>
void RTNodalPointsTriangle(const Cell &c, vector<PointT<T, sDim, tDim>> &z, int order) {
  RTNodalPointsTriangle(c, z, 0, order + 2, order + 2, 0);
}

template<class POINT>
vector<POINT> RTNodalPointsTriangle(const Cell &c, int order) {
  vector<POINT> z(3 * (order + 1) + (order * (order + 1)) / 2);
  RTNodalPointsTriangle(c, z, order);
  return z;
}


#endif // RTNODALPOINTS_HPP
