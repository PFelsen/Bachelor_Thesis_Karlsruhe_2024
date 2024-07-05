#ifndef EGDOF_HPP
#define EGDOF_HPP

// #include "DoFs.hpp"


// class EGDoF : public DoF {
//   short degree;
//   short dofSize;
// public:
//   explicit EGDoF(short CellStandardDegree, short dofSize = 1) :
//       DoF(0, false), degree(CellStandardDegree), dofSize(dofSize) {}
//
//   string Name() const override { return "EGDoF"; }
//
//   int TypeDoF(int) const override { return VERTEX; }
//
//   int GetNodalPoints(const Cell &c) const override;
//
//   void GetNodalPoints(const Cell &c, vector<Point> &z) const override;
//
//   int get_cell_deg(const Cell &c) const override {
//     return degree;
//   }
//
//   int CellDoFs(const Cell &c) const override;
//
//   void NodalDoFs(const Cell &c, vector<short> &z) const override;
//
//   int GetNodalPointsOnFace(const Cell &c, int face) const override;
//
//   int NodalPointOnFace(const Cell &c, int face, int k) const override;
// };


#endif // EGDOF_HPP
