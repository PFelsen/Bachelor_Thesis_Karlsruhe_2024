#ifndef SERENDIPITYSHAPES_HPP
#define SERENDIPITYSHAPES_HPP

#include "LagrangeShapes.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class SerendipityShapeT : public ShapeT<T, sDim, tDim> {
public:
  explicit SerendipityShapeT(int numNodalPoints) : ShapeT<T, sDim, tDim>(numNodalPoints) {}

  void fillValues(const vector<PointT<T, sDim, tDim>> &localQuadraturePoints) override {
    this->localValues.resize(localQuadraturePoints.size(), this->numNodalPoints);
    this->localGradient.resize(localQuadraturePoints.size(), this->numNodalPoints);
    for (int q = 0; q < localQuadraturePoints.size(); ++q)
      for (int i = 0; i < this->numNodalPoints; ++i) {
        this->localValues[q][i] = (*this)(localQuadraturePoints[q], i);
        this->localGradient[q][i] = this->LocalGradient(localQuadraturePoints[q], i);
      }
  }

  void
  fillFaceValues(const vector<vector<PointT<T, sDim, tDim>>> &localFaceQuadraturePoints) override {
    this->localFaceValues.resize(localFaceQuadraturePoints.size(),
                                 localFaceQuadraturePoints[0].size(), this->numNodalPoints);
    for (int face = 0; face < localFaceQuadraturePoints.size(); ++face)
      for (int q = 0; q < localFaceQuadraturePoints[face].size(); ++q)
        for (int i = 0; i < this->numNodalPoints; ++i) {
          this->localFaceValues[face][q][i] = (*this)(localFaceQuadraturePoints[face][q], i);
        }
  }
};

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P2QuadSerendipityT : public SerendipityShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P2QuadSerendipity"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override;

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const override;

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(c.Corners() + c.Faces());
    for (int i = 0; i < c.Corners(); ++i) {
      z[i] = c[i];
      z[i + c.Corners()] = c.Face(i);
    }
  }

  explicit P2QuadSerendipityT() : SerendipityShapeT<T, sDim, tDim>(8) {}
};

typedef P2QuadSerendipityT<> P2QuadSerendipity;

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class P2HexSerendipityT : public SerendipityShapeT<T, sDim, tDim> {
public:
  const std::string Name() const override { return "P2HexSerendipity"; }

  T operator()(const PointT<T, sDim, tDim> &z, int i) const override;

  VectorFieldT<T, sDim> LocalGradient(const PointT<T, sDim, tDim> &z, int i) const;

  void NodalPoints(const Cell &c, vector<PointT<T, sDim, tDim>> &z) const override {
    z.resize(c.Corners() + c.Edges());
    for (int i = 0; i < c.Corners(); ++i)
      z[i] = c[i];
    for (int i = 0; i < c.Edges(); ++i)
      z[i + c.Corners()] = c.Edge(i);
  }

  P2HexSerendipityT() : SerendipityShapeT<T, sDim, tDim>(20) {}
};

typedef P2HexSerendipityT<> P2HexSerendipity;

template<typename T, int sDim, int tDim>
T P2QuadSerendipityT<T, sDim, tDim>::operator()(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return ((1 - z[0]) * (1 - z[1]) * (1 - 2 * z[0] - 2 * z[1]));
  case 1:
    return (-z[0] * (1 - z[1]) * (1 - 2 * z[0] + 2 * z[1]));
  case 2:
    return (-z[0] * z[1] * (3 - 2 * z[0] - 2 * z[1]));
  case 3:
    return (-z[1] * (1 - z[0]) * (1 + 2 * z[0] - 2 * z[1]));
  case 4:
    return (4 * z[0] * (1 - z[0]) * (1 - z[1]));
  case 5:
    return (4 * z[0] * z[1] * (1 - z[1]));
  case 6:
    return (4 * z[0] * z[1] * (1 - z[0]));
  case 7:
    return (4 * z[1] * (1 - z[0]) * (1 - z[1]));
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
VectorFieldT<T, sDim>
P2QuadSerendipityT<T, sDim, tDim>::LocalGradient(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return VectorFieldT<T, sDim>((1 - z[1]) * (-3 + 4 * z[0] + 2 * z[1]),
                                 (1 - z[0]) * (-3 + 4 * z[1] + 2 * z[0]));
  case 1:
    return VectorFieldT<T, sDim>((1 - z[1]) * (-1 + 4 * z[0] - 2 * z[1]),
                                 -z[0] * (1 + 2 * z[0] - 4 * z[1]));
  case 2:
    return VectorFieldT<T, sDim>(-z[1] * (3 - 4 * z[0] - 2 * z[1]),
                                 -z[0] * (3 - 4 * z[1] - 2 * z[0]));
  case 3:
    return VectorFieldT<T, sDim>(-z[1] * (1 - 4 * z[0] + 2 * z[1]),
                                 (z[0] - 1) * (1 + 2 * z[0] - 4 * z[1]));
  case 4:
    return VectorFieldT<T, sDim>(4 * (1 - 2 * z[0]) * (1 - z[1]), -4 * z[0] * (1 - z[0]));
  case 5:
    return VectorFieldT<T, sDim>(4 * z[1] * (1 - z[1]), 4 * z[0] * (1 - 2 * z[1]));
  case 6:
    return VectorFieldT<T, sDim>(4 * z[1] * (1 - 2 * z[0]), 4 * z[0] * (1 - z[0]));
  case 7:
    return VectorFieldT<T, sDim>(-4 * z[1] * (1 - z[1]), 4 * (1 - z[0]) * (1 - 2 * z[1]));
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
T P2HexSerendipityT<T, sDim, tDim>::operator()(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return ((1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2])
            * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2]));
  case 1:
    return ((z[0]) * (1.0 - z[1]) * (1.0 - z[2])
            * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2]) - 5.0));
  case 2:
    return ((z[0]) * (z[1]) * (1.0 - z[2])
            * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0));
  case 3:
    return ((1.0 - z[0]) * (z[1]) * (1.0 - z[2])
            * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0));
  case 4:
    return ((1.0 - z[0]) * (1.0 - z[1]) * (z[2])
            * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (z[2]) - 5.0));
  case 5:
    return ((z[0]) * (1.0 - z[1]) * (z[2])
            * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (z[2]) - 5.0));
  case 6:
    return ((z[0]) * (z[1]) * (z[2]) * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (z[2]) - 5.0));
  case 7:
    return ((1.0 - z[0]) * (z[1]) * (z[2])
            * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (z[2]) - 5.0));
  case 8:
    return (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[0]);
  case 9:
    return (4.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[1]);
  case 10:
    return (4.0 * (1.0 - z[0]) * (z[1]) * (1.0 - z[2]) * z[0]);
  case 11:
    return (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[1]);
  case 12:
    return (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[2]);
  case 13:
    return (4.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]) * z[2]);
  case 14:
    return (4.0 * (z[0]) * (z[1]) * (1.0 - z[2]) * z[2]);
  case 15:
    return (4.0 * (1.0 - z[0]) * (z[1]) * (1.0 - z[2]) * z[2]);
  case 16:
    return (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (z[2]) * z[0]);
  case 17:
    return (4.0 * (z[0]) * (1.0 - z[1]) * (z[2]) * z[1]);
  case 18:
    return (4.0 * (1.0 - z[0]) * (z[1]) * (z[2]) * z[0]);
  case 19:
    return (4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (z[2]) * z[1]);
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
VectorFieldT<T, sDim>
P2HexSerendipityT<T, sDim, tDim>::LocalGradient(const PointT<T, sDim, tDim> &z, int i) const {
  switch (i) {
  case 0:
    return VectorFieldT<T, sDim>(-(1.0 - z[1]) * (1.0 - z[2])
                                         * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2])
                                     - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                                 -(1.0 - z[0]) * (1.0 - z[2])
                                         * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2])
                                     - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                                 -(1.0 - z[0]) * (1.0 - z[1])
                                         * (1.0 - 2.0 * z[0] - 2.0 * z[1] - 2.0 * z[2])
                                     - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - z[2]));
  case 1:
    return VectorFieldT<T, sDim>((1.0 - z[1]) * (1.0 - z[2])
                                         * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2])
                                            - 5.0)
                                     + 2.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                                 -z[0] * (1.0 - z[2])
                                         * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2])
                                            - 5.0)
                                     - 2.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                                 -z[0] * (1.0 - z[1])
                                         * (2.0 * (z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * (1.0 - z[2])
                                            - 5.0)
                                     - 2.0 * (z[0]) * (1.0 - z[1]) * (1.0 - z[2]));
  case 2:
    return VectorFieldT<T, sDim>(z[1] * (1.0 - z[2])
                                         * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                                     + 2.0 * z[0] * z[1] * (1.0 - z[2]),
                                 z[0] * (1.0 - z[2])
                                         * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                                     + 2.0 * z[0] * z[1] * (1.0 - z[2]),
                                 -z[0] * z[1]
                                         * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2]) - 5.0)
                                     - 2.0 * z[0] * z[1] * (1.0 - z[2]));
  case 3:
    return VectorFieldT<T, sDim>(-z[1] * (1.0 - z[2])
                                         * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2])
                                            - 5.0)
                                     - 2.0 * (1.0 - z[0]) * z[1] * (1.0 - z[2]),
                                 (1.0 - z[0]) * (1.0 - z[2])
                                         * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2])
                                            - 5.0)
                                     + 2.0 * (1.0 - z[0]) * z[1] * (1.0 - z[2]),
                                 -(1.0 - z[0]) * z[1]
                                         * (2.0 * (1.0 - z[0]) + 2.0 * (z[1]) + 2.0 * (1.0 - z[2])
                                            - 5.0)
                                     - 2.0 * (1.0 - z[0]) * z[1] * (1.0 - z[2]));
  case 4:
    return VectorFieldT<
        T, sDim>(-(1.0 - z[1]) * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                     - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[2],
                 -(1.0 - z[0]) * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                     - 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[2],
                 (1.0 - z[0]) * (1.0 - z[1])
                         * (2.0 * (1.0 - z[0]) + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                     + 2.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[2]);
  case 5:
    return VectorFieldT<T, sDim>((1.0 - z[1]) * z[2]
                                         * (2.0 * z[0] + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                                     + 2.0 * z[0] * (1.0 - z[1]) * z[2],
                                 -z[0] * z[2] * (2.0 * z[0] + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                                     - 2.0 * z[0] * (1.0 - z[1]) * z[2],
                                 z[0] * (1.0 - z[1])
                                         * (2.0 * z[0] + 2.0 * (1.0 - z[1]) + 2.0 * z[2] - 5.0)
                                     + 2.0 * z[0] * (1.0 - z[1]) * z[2]);
  case 6:
    return VectorFieldT<T, sDim>(z[1] * z[2] * (2.0 * z[0] + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                                     + 2.0 * z[0] * z[1] * z[2],
                                 z[0] * z[2] * (2.0 * z[0] + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                                     + 2.0 * z[0] * z[1] * z[2],
                                 z[0] * z[1] * (2.0 * (z[0]) + 2.0 * (z[1]) + 2.0 * z[2] - 5.0)
                                     + 2.0 * z[0] * z[1] * z[2]);
  case 7:
    return VectorFieldT<T, sDim>(-z[1] * z[2] * (2.0 * (1.0 - z[0]) + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                                     - 2.0 * (1.0 - z[0]) * z[1] * z[2],
                                 (1.0 - z[0]) * z[2]
                                         * (2.0 * (1.0 - z[0]) + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                                     + 2.0 * (1.0 - z[0]) * z[1] * z[2],
                                 (1.0 - z[0]) * z[1]
                                         * (2.0 * (1.0 - z[0]) + 2.0 * z[1] + 2.0 * z[2] - 5.0)
                                     + 2.0 * (1.0 - z[0]) * z[1] * z[2]);
  case 8:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - 2.0 * z[0]) * (1.0 - z[1]) * (1.0 - z[2]),
                                 -4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[0],
                                 -4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[0]);
  case 9:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[1],
                                 4.0 * (z[0]) * (1.0 - 2.0 * z[1]) * (1.0 - z[2]),
                                 -4.0 * (z[0]) * (1.0 - z[1]) * z[1]);
  case 10:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - 2.0 * z[0]) * z[1] * (1.0 - z[2]),
                                 4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[0],
                                 -4.0 * (1.0 - z[0]) * z[1] * z[0]);
  case 11:
    return VectorFieldT<T, sDim>(-4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[1],
                                 4.0 * (1.0 - z[0]) * (1.0 - 2.0 * z[1]) * (1.0 - z[2]),
                                 -4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[1]);
  case 12:
    return VectorFieldT<T, sDim>(-4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[2],
                                 -4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[2],
                                 4.0 * (1.0 - z[0]) * (1.0 - z[1]) * (1.0 - 2.0 * z[2]));
  case 13:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - z[1]) * (1.0 - z[2]) * z[2],
                                 -4.0 * (z[0]) * (1.0 - z[2]) * z[2],
                                 4.0 * z[0] * (1.0 - z[1]) * (1.0 - 2.0 * z[2]));
  case 14:
    return VectorFieldT<T, sDim>(4.0 * z[1] * (1.0 - z[2]) * z[2], 4.0 * z[0] * (1.0 - z[2]) * z[2],
                                 4.0 * z[0] * z[1] * (1.0 - 2.0 * z[2]));
  case 15:
    return VectorFieldT<T, sDim>(-4.0 * z[1] * (1.0 - z[2]) * z[2],
                                 4.0 * (1.0 - z[0]) * (1.0 - z[2]) * z[2],
                                 4.0 * (1.0 - z[0]) * z[1] * (1.0 - 2.0 * z[2]));
  case 16:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - 2.0 * z[0]) * (1.0 - z[1]) * z[2],
                                 -4.0 * (1.0 - z[0]) * z[2] * z[0],
                                 4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[0]);
  case 17:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - z[1]) * z[2] * z[1],
                                 4.0 * z[0] * (1.0 - 2.0 * z[1]) * z[2],
                                 4.0 * z[0] * (1.0 - z[1]) * z[1]);
  case 18:
    return VectorFieldT<T, sDim>(4.0 * (1.0 - 2.0 * z[0]) * z[1] * z[2],
                                 4.0 * (1.0 - z[0]) * z[2] * z[0],
                                 4.0 * (1.0 - z[0]) * z[1] * z[0]);
  case 19:
    return VectorFieldT<T, sDim>(-4.0 * (1.0 - z[1]) * z[2] * z[1],
                                 4.0 * (1.0 - z[0]) * (1.0 - 2.0 * z[1]) * z[2],
                                 4.0 * (1.0 - z[0]) * (1.0 - z[1]) * z[1]);
  default:
    THROW("Not implemented")
  }
}

template<typename T, int sDim, int tDim>
ShapeT<T, sDim, tDim> *createSerendipityShape(const CELLTYPE cellType) {
  switch (cellType) {
  case INTERVAL:
    return new IntervalShape<2, T, sDim, tDim>();
  case TRIANGLE:
    return new TriangularShape<2, T, sDim, tDim>();
  case QUADRILATERAL:
    return new P2QuadSerendipityT<T, sDim, tDim>();
  case TETRAHEDRON:
    return new P2TetT<T, sDim, tDim>();
  case HEXAHEDRON:
    return new P2HexSerendipityT<T, sDim, tDim>();
  default:
    THROW("Shape not implemented")
  }
}


#endif // SERENDIPITYSHAPES_HPP
