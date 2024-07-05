#include "EGDoF.hpp"
#include "NodalPointProvider.hpp"

// int EGDoF::GetNodalPoints(const Cell &c) const {
//   return EquidistantNodalPointCount(c.ReferenceType(), degree) + 1;
// }
//
// void EGDoF::GetNodalPoints(const Cell &c, vector<Point> &z) const {
//   EquidistantNodalPoints(c, degree, z);
//   z.emplace_back(c());
// }
//
// int EGDoF::CellDoFs(const Cell &c) const {
//   return EquidistantNodalPointCount(c.ReferenceType(), degree) * dofSize + 1;
// }
//
// void EGDoF::NodalDoFs(const Cell &c, vector<short> &z) const {
//   z.resize(EquidistantNodalPointCount(c.ReferenceType(), degree) + 1);
//   for (int i = 0; i < z.size() - 1; ++i)
//     z[i] = dofSize;
//   z[c.Corners()] = 1;
// }
//
// int EGDoF::GetNodalPointsOnFace(const Cell &c, int face) const {
//   return EquidistantNodalPointCountOnFace(c.ReferenceType(), degree, face);
// }
//
// int EGDoF::NodalPointOnFace(const Cell &c, int face, int k) const {
//   return EquidistantNodalPointIdOnFace(c.ReferenceType(), degree, face, k);
// }
