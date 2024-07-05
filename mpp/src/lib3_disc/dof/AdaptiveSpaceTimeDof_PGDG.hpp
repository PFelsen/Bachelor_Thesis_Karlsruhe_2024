#ifndef ADAPTIVESPACETIMEDOF_PGDG_HPP
#define ADAPTIVESPACETIMEDOF_PGDG_HPP

#include "AdaptiveSpaceTimeDof.hpp"
#include "IDoF.hpp"
#include "LagrangeNodalPoints.hpp"
#include "NodalPointProvider.hpp"

class AdaptiveSpaceTimeDof_PGDG : public AdaptiveSpaceTimeDof {
  const int probDim;
  vector<vector<Point>> NP;
  vector<int> DOF;
public:
  AdaptiveSpaceTimeDof_PGDG(const Mesh &STM, const int &_dim, NodalPointProvider *npp) :
      AdaptiveSpaceTimeDof(1), probDim(_dim), DOF(7), NP(7) {
    cell c = STM.cells();
    if (c != STM.cells_end()) {
      CELLTYPE celltype = c.Type();
      for (int deg = 0; deg < 7; ++deg) {
        NP[deg] = npp->GetNodalPoints(SpaceCellType(celltype), deg);
        DOF[deg] = probDim * NP[deg].size();
      }
    }
  }

  string Name() const { return "AdaptiveSpaceTimeDof_PGDG"; }

  int get_dim() const { return probDim; }

  int get_m(const Point &p) const override {
    int deg = get_space_deg(p);
    return DOF[deg];
  }

  unsigned int get_m(const Cell &c) const { return get_m(c()); }

  short NumberOfNodalPoints(const int sDeg, const int tDeg) const override {
    return DOF[sDeg] * tDeg;
  }

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>(NumberOfNodalPoints(c), probDim);
  }

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override {
    THROW("NumberOfNodalPointsOnFace not implemented")
  }

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override {
    THROW("IdOfNodalPointOnFace not implemented")
  }

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override {
    THROW("NumberOfNodalPointsOnEdge not implemented")
  }

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId,
                             int k) const override{THROW("IdOfNodalPointOnEdge not implemented")}

  std::vector<short> AllocationSizesAtStoragePoints(const Cell &c) const override {
    return std::vector<short>{short(get_m(c()) * get_time_deg(c()))};
  }

  short NumberOfNodalPoints(const Cell &c) const override { return get_m(c()) * get_time_deg(c()); }

  std::vector<Point> GetNodalPoints(const Cell &c) const override {
    std::vector<Point> z(NumberOfNodalPoints(c));
    int spaceDeg = get_space_deg(c());
    int timeDeg = get_time_deg(c());
    double ti = (c.max() - c.min()) / (double)timeDeg;

    int cnt = 0;
    for (int l = 0; l < timeDeg; ++l) {
      for (int p = 0; p < probDim; ++p) {
        for (const Point &np : NP[spaceDeg]) {
          z[cnt] = c.SpaceCell().LocalToGlobal(np).CopyWithT(c.min() + (l + 1) * ti);
          cnt++;
        }
      }
    }
    return z;
  }

  void NodalPointsLocal(const int &deg, vector<Point> &z, const int p = 0) const {
    int nps = NP[deg].size();
    z.resize(nps);
    for (int i = 0; i < nps; ++i) {
      z[i] = NP[deg][i].CopyWithT(1.0);
    }
  }

  int NodalPointsLocal(const int &deg, const int p = 0) const { return NP[deg].size(); }
};

#endif // ADAPTIVESPACETIMEDOF_PGDG_HPP
