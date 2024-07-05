#ifndef ADAPTIVESPACETIMEDOF_DGDG_HPP
#define ADAPTIVESPACETIMEDOF_DGDG_HPP

#include "AdaptiveSpaceTimeDof.hpp"
#include "IDoF.hpp"
#include "LagrangeNodalPoints.hpp"
#include "NodalPointProvider.hpp"

class AdaptiveSpaceTimeDof_DGDG : public AdaptiveSpaceTimeDof {
  const int probDim;
  vector<vector<Point>> SpaceNPForDeg;
  vector<int> SpaceDoFCountForDeg;
public:
  AdaptiveSpaceTimeDof_DGDG(const Mesh &STM, const int &_dim, NodalPointProvider &npp) :
      AdaptiveSpaceTimeDof(0), probDim(_dim), SpaceDoFCountForDeg(7), SpaceNPForDeg(7) {
    cell c = STM.cells();
    if (c != STM.cells_end()) {
      CELLTYPE celltype = c.Type();
      for (int deg = 0; deg < 7; ++deg) {
        SpaceNPForDeg[deg] = npp.GetNodalPoints(SpaceCellType(celltype), deg);
        SpaceDoFCountForDeg[deg] = probDim * SpaceNPForDeg[deg].size();
      }
    }
  }

  string Name() const override { return "AdaptiveSpaceTimeDof_DGDG"; }

  int get_dim() const override { return probDim; }

  int get_m(const Point &p) const override {
    int deg = get_space_deg(p);
    return SpaceDoFCountForDeg[deg];
  }

  int get_m(const Cell &c) const { return get_m(c()); }

  short NumberOfNodalPoints(const int sDeg, const int tDeg) const override {
    return SpaceDoFCountForDeg[sDeg] * (1 + tDeg);
  }

  std::vector<short> AllocationSizesAtStoragePoints(const Cell &c) const override {
    return std::vector<short>{short(get_m(c()) * (1 + get_time_deg(c())))};
  }

  short NumberOfNodalPoints(const Cell &c) const override {
    return get_m(c()) * (1 + get_time_deg(c()));
  }

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>(NumberOfNodalPoints(c), probDim);
  }

  short NumberOfComponents() const override { return probDim; }

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

  vector<Point> GetNodalPoints(const Cell &c) const override {
    int spaceDeg = get_space_deg(c());
    int timeDeg = get_time_deg(c());
    vector<Point> z(NumberOfNodalPoints(c));
    double time_interval = (c.max() - c.min());
    double ti = time_interval * ((timeDeg == 0) ? 0.5 : 1.0 / timeDeg);

    int cnt = 0;
    double time = 0.0;
    for (int l = 0; l <= timeDeg; ++l) {
      for (int p = 0; p < probDim; ++p) {
        for (const Point &np : SpaceNPForDeg[spaceDeg]) {
          if (timeDeg == 0) {
            time = c.min() + ti;
          } else {
            time = c.min() + l * ti;
          }
          z[cnt] = c.SpaceCell().LocalToGlobal(np).CopyWithT(time);
          cnt++;
        }
      }
    }
    return z;
  }

  void NodalPointsLocal(const int &deg, vector<Point> &z, const int p = 0) const override {
    int nps = SpaceNPForDeg[deg].size();
    z.resize(nps);
    for (int i = 0; i < nps; ++i) {
      z[i] = SpaceNPForDeg[deg][i].CopyWithT(1.0);
    }
  }

  int NodalPointsLocal(const int &deg, const int p = 0) const override {
    return SpaceNPForDeg[deg].size();
  }
};

#endif // ADAPTIVESPACETIMEDOF_DGDG_HPP
