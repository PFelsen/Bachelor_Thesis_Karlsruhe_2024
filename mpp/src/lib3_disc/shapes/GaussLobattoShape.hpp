#ifndef GAUSSLOBATTOSHAPE_HPP
#define GAUSSLOBATTOSHAPE_HPP

#include "LagrangeShape.hpp"
#include "NodalPointProvider.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class GaussLobattoShape : public LagrangeShapeT<T, sDim, tDim> {
public:
  const int polynomialDegree;
  std::vector<T> localNodalPointsAsScalars;
  const int dim;

  GaussLobattoShape(int polynomialDegree, int dim) :
      LagrangeShapeT<T, sDim, tDim>(pow(polynomialDegree + 1, dim)),
      polynomialDegree(polynomialDegree), dim(dim) {
    localNodalPointsAsScalars = GAUSSLOBATTO_NODALPOINTS.at(polynomialDegree + 1);
  }

  const std::string Name() const {
    std::stringstream ss;
    ss << "GaussLobattoShape dim=" << dim << " degree " << polynomialDegree;
    return ss.str();
  }

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    if (dim == 1) {
      z = GaussLobattoNodalPointsInterval<T, sDim, tDim>(c, polynomialDegree + 1);
    } else if (dim == 2) {
      z = GaussLobattoNodalPointsQuadrilateral<T, sDim, tDim>(c, polynomialDegree + 1);
    } else if (dim == 3) {
      z = GaussLobattoNodalPointsHexaedron<T, sDim, tDim>(c, polynomialDegree + 1);
    } else {
      THROW("GL 1<=dim<=3")
    }
  }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override {
    auto &NP = localNodalPointsAsScalars;

    int i_idx[dim];
    for (int d = 0; d < dim; d++) {
      i_idx[d] = i % NP.size();
      i /= NP.size();
    }
    T value = 1;
    for (int d = 0; d < dim; d++) {
      for (int j = 0; j < NP.size(); j++) {
        if (i_idx[d] == j) continue;
        value *= (z[d] - NP[j]) / (NP[i_idx[d]] - NP[j]);
      }
    }
    return value;
  }

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override {
    auto &NP = localNodalPointsAsScalars;
    VectorFieldT<T, sDim> gradient;

    int i_idx[dim];
    for (int d = 0; d < dim; d++) {
      i_idx[d] = i % NP.size();
      i /= NP.size();
    }

    for (int d = 0; d < dim; d++) {
      for (int k = 0; k < NP.size(); k++) {
        if (k == i_idx[d]) continue;
        T toAdd = 1.0 / (NP[i_idx[d]] - NP[k]);
        for (int j = 0; j < NP.size(); j++) {
          if (j == i_idx[d] || j == k) continue;
          toAdd *= (z[d] - NP[j]) / (NP[i_idx[d]] - NP[j]);
        }
        gradient[d] += toAdd;
      }
      for (int dd = 0; dd < dim; dd++) {
        if (d == dd) continue;
        for (int j = 0; j < NP.size(); j++) {
          if (j == i_idx[dd]) continue;
          gradient[d] *= (z[dd] - NP[j]) / (NP[i_idx[dd]] - NP[j]);
        }
      }
    }
    return gradient;
  }
};


#endif // GAUSSLOBATTOSHAPE_HPP
