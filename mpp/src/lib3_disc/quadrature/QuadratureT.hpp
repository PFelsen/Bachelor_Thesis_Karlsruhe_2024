#ifndef QUADRATURET_HPP
#define QUADRATURET_HPP

#include "Celltype.hpp"
#include "Point.hpp"

#include <vector>

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QuadratureT {
protected:
  std::string name;
  std::vector<PointT<T, sDim, tDim>> z;
  std::vector<T> w;
public:
  QuadratureT(const std::string name = "", int numQPoints = 0) :
      name(name), z(numQPoints), w(numQPoints) {}

  const std::string &Name() const { return name; }

  int size() const { return z.size(); }

  const PointT<T, sDim, tDim> &QPoint(int i) const { return z[i]; }

  const std::vector<PointT<T, sDim, tDim>> &QPoints() const { return z; }

  const T &Weight(int i) const { return w[i]; }

  virtual ~QuadratureT() = default;
protected:
  void resize(int numQPoints) {
    z.resize(numQPoints);
    w.resize(numQPoints);
  }
};

template<typename T, int sDim, int tDim>
std::ostream &operator<<(std::ostream &s, const QuadratureT<T, sDim, tDim> &Q) {
  s << Q.Name() << ":" << endl;
  for (int q = 0; q < Q.size(); ++q)
    s << " " << Q.QPoint(q) << " <> " << Q.Weight(q) << endl;
  return s;
}

typedef QuadratureT<> Quadrature;

#ifdef BUILD_IA

typedef QuadratureT<IAInterval, SpaceDimension, TimeDimension> IAQuadrature;

#endif


template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class QTensorT : public QuadratureT<T, sDim, tDim> {
public:
  explicit QTensorT(std::vector<QuadratureT<T, sDim, tDim>> quads, bool special = false) {
    int new_size = 1;
    for (const auto &Q : quads) {
      new_size *= Q.size();
    }
    this->z.resize(new_size);
    this->w.resize(new_size);
    if (quads.size() == 2) {
      if (special) {
        int n = 0;
        for (int qY = 0; qY < quads[1].size(); ++qY) {
          for (int qX = 0; qX < quads[0].size(); ++qX) {
            this->z[n] = PointT<T, sDim, tDim>(quads[0].QPoint(qX)[0], quads[0].QPoint(qX)[1],
                                               quads[1].QPoint(qY)[0]);
            this->w[n++] = quads[0].Weight(qX) * quads[1].Weight(qY);
          }
        }
      } else {
        int n = 0;
        for (int qY = 0; qY < quads[1].size(); ++qY) {
          for (int qX = 0; qX < quads[0].size(); ++qX) {
            this->z[n] = PointT<T, sDim, tDim>(quads[0].QPoint(qX)[0], quads[1].QPoint(qY)[0]);
            this->w[n++] = quads[0].Weight(qX) * quads[1].Weight(qY);
          }
        }
      }
      this->name = "QTensor: [" + quads[0].Name() + ", " + quads[1].Name() + "]";
    } else if (quads.size() == 3) {
      int n = 0;
      for (int qZ = 0; qZ < quads[2].size(); ++qZ) {
        for (int qY = 0; qY < quads[1].size(); ++qY) {
          for (int qX = 0; qX < quads[0].size(); ++qX) {
            this->z[n] = PointT<T, sDim, tDim>(quads[0].QPoint(qX)[0], quads[1].QPoint(qY)[0],
                                               quads[2].QPoint(qZ)[0]);
            this->w[n++] = quads[0].Weight(qX) * quads[1].Weight(qY) * quads[2].Weight(qZ);
          }
        }
      }
      this->name =
          "QTensor: [" + quads[0].Name() + ", " + quads[1].Name() + ", " + quads[2].Name() + "]";
    } else {
      THROW("Only 2 and 3 dimensional TensorsQuads are implemented.")
    }
  }
};

using QTensor = QTensorT<>;

#endif // QUADRATURET_HPP
