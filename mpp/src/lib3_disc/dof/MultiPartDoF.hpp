#ifndef _MULTIPARTDOF_H_
#define _MULTIPARTDOF_H_

#include "Intval.hpp"
#include "LagrangeDoF.hpp"
#include "Quadrilateral.hpp"
#include "Triangle.hpp"

class MultiPartDoF : public LagrangeDoF {
  vector<Point> interfacePoints{};
public:
  MultiPartDoF(int degree, const Mesh &_M, int _m = 1, int n = 0, bool bnd = false) :
      LagrangeDoF(degree, _m, n, bnd) {}

  void AddInterfacePoint(Point p);

  void FinishInterfacePoints();

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  string Name() const override { return "MultiPartDoF"; }
};

#endif
