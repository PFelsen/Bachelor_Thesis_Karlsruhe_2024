#ifndef WCDOF_HPP
#define WCDOF_HPP

#include "IDoF.hpp"

class WCDoF : public IDoF {
  int cell_std_deg;
  int face_std_dofs;
  int face_std_dual_dofs;
  int face_normal_dofs;
  std::unordered_map<Point, vector<short>> cell_deg;
  std::unordered_map<Point, vector<vector<short>>> face_dofs;
public:
  WCDoF(int cell_std_deg_, int face_std_dofs_, int face_std_dual_dofs_, int face_normal_dofs_);

  short NumberOfNodalPoints(const Cell &c) const override;

  std::vector<Point> GetNodalPoints(const Cell &c) const override;

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override;

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override;

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override;

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override;

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override;

  string Name() const override { return "WCDoF"; }

  int get_cell_deg(const Cell &c) const override;

  int get_cell_deg_WC(const Point &z, int ud) const;

  vector<short> get_face_dofs(const Point &z, int pd) const;

  void set_face_dofs(const Point &z, int pd, short dofs);
};

WCDoF *GetDoF();

#endif // WCDOF_HPP
