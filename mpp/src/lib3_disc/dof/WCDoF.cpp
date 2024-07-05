#include "WCDoF.hpp"

WCDoF::WCDoF(int cell_std_deg_, int face_std_dofs_, int face_std_dual_dofs_,
             int face_normal_dofs_) :
    IDoF(0, false), cell_std_deg(cell_std_deg_), face_std_dofs(face_std_dofs_),
    face_std_dual_dofs(face_std_dual_dofs_), face_normal_dofs(face_normal_dofs_) {}

short WCDoF::NumberOfNodalPoints(const Cell &c) const { return c.Faces(); }

std::vector<Point> WCDoF::GetNodalPoints(const Cell &c) const {
  std::vector<Point> z(NumberOfNodalPoints(c));
  for (int i = 0; i < c.Faces(); ++i)
    z[i] = c.Face(i);
  return z;
}

std::vector<short> WCDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  std::vector<Point> n = GetNodalPoints(c);
  std::vector<short> z(n.size());
  for (int i = 0; i < z.size(); ++i) {
    z[i] = 0;
    for (int pd = 0; pd < 2; ++pd)
      z[i] += get_face_dofs(n[i], pd).size();
    z[i] += get_face_dofs(n[i], 2).size();
  }
  return z;
}

short WCDoF::NumberOfNodalPointsOnFace(const Cell &c, int faceId) const { return 1; }

short WCDoF::IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const { return faceId; }

short WCDoF::NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const { THROW("TODO") }

short WCDoF::IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const { THROW("TODO") }

int WCDoF::get_cell_deg(const Cell &c) const {
  Point z = c();
  if (cell_deg.find(z) != cell_deg.end()) {
    const vector<short> &degs = cell_deg.find(z)->second;
    if (degs[0] > 0) return degs[0]; // Todo vorher [ud]
  }
  return cell_std_deg;
}

int WCDoF::get_cell_deg_WC(const Point &z, int ud) const {
  if (cell_deg.find(z) != cell_deg.end()) {
    const vector<short> &degs = cell_deg.find(z)->second;
    if (degs[ud] > 0) return degs[ud];
  }
  return cell_std_deg;
}

vector<short> WCDoF::get_face_dofs(const Point &z, int pd) const {
  if (face_dofs.find(z) != face_dofs.end()) {
    if (pd < 2) {
      const vector<vector<short>> &dofs = face_dofs.find(z)->second;
      if (dofs[pd].size() == 0 || dofs[pd][0] != -1) return dofs[pd];
    } else {
      const vector<vector<short>> &dofs = face_dofs.find(z)->second;
      if (dofs[2].size() == 0 || dofs[2][0] != -1) return dofs[2];
    }
  }
  int nr_std_dofs = 0;
  if (pd == 0) nr_std_dofs = face_std_dofs;
  else if (pd == 1) nr_std_dofs = face_std_dual_dofs;
  else nr_std_dofs = face_normal_dofs;
  vector<short> std_dofs(nr_std_dofs, 0.0);
  for (int i = 0; i < nr_std_dofs; i++)
    std_dofs[i] = i;
  return std_dofs;
}

void WCDoF::set_face_dofs(const Point &z, int pd, short dofs) {
  vector<vector<short>> *f_dofs = NULL;
  if (face_dofs.find(z) != face_dofs.end()) f_dofs = &face_dofs[z];
  else {
    f_dofs = &face_dofs[z];
    f_dofs->resize(3, vector<short>(1, -1));
  }
  vector<short> vec_dofs(dofs, 0);
  for (int i = 0; i < dofs; i++)
    vec_dofs[i] = i;
  if (pd < 2) (*f_dofs)[pd] = vec_dofs;
  else (*f_dofs)[2] = vec_dofs;
}

WCDoF *GetDoF() {
  int cell_std_deg = 3;
  int face_std_dofs = 3;
  int face_std_dual_dofs = 0;
  int face_std_normal_dofs = 0;
  Config::Get("degree", cell_std_deg);
  Config::Get("face_dofs", face_std_dofs);

  if (cell_std_deg == 2)
    if (face_std_dofs == 2) THROW("degree too small")

  Config::Get("ftestSN", face_std_dual_dofs);
  Config::Get("ftestN", face_std_normal_dofs);
  return new WCDoF(cell_std_deg, face_std_dofs, face_std_dual_dofs, face_std_normal_dofs);
}