#ifndef CURLELEMENT_H
#define CURLELEMENT_H

#include "Element.hpp"
#include "RMatrix.hpp"

class CurlElement : public Element {
  const Shape &S;     // TODO: remove and use shape in Element!!
  double sign[216];   // TODO: should be a std::vector
  Point tangent[216]; // TODO: should be a std::vector
  ShapeValues<VectorField> vectorfield;
  ShapeValues<VectorField> curlfield;

  Point DofShape_int(double x, double y, double z, int i) const;
public:
  CurlElement(const VectorMatrixBase &base, const Cell &c);

  //  CurlElement(const VectorMatrixBase &base, const BasicElementT<> &baseElement);

  VectorField VectorValue(int q, int i) const { return vectorfield[q][i]; }

  VectorField CurlVector(int q, int i) const { return curlfield[q][i]; }

  VectorField VectorValue(int q, const Vector &u, int k = 0) const;

  VectorField VectorValue(const Point &, int k) const;

  VectorField VectorValue(const Point &, const Vector &, int k = 0) const;

  VectorField CurlVector(int q, const Vector &u, int k = 0) const;

  Point DofShape(const Point &z, int i) const { return DofShape_int(z[0], z[1], z[2], i); }

  Point TangentVector(int i) const { return tangent[i]; }

  Point Nodal_Point(int) const;

  Point Q_Point(int) const;

  Point LocalQ_Point(int) const;

  Point Edge_Corner(int, int) const;

  Point Face_Corner(int, int) const;

  //    short FaceNodalPoints (int i) const { return c.FaceEdges(i); }
  //    short FaceNodalPoint (int i, int j) const { return c.faceedge(i,j); }
  friend std::ostream &operator<<(std::ostream &s, const CurlElement &E);

  double GetSign(int i) const { return sign[i]; }

  Point ReferenceTangent(int) const;
};

#endif // CURLELEMENT_HPP
