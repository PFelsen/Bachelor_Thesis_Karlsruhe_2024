#include "SpaceTimeDiscretization.hpp"
#include "Parallel.hpp"

template<typename T, int sDim, int tDim>
STDiscretizationT<T, sDim, tDim>::STDiscretizationT(const Meshes &meshes, DegreePair degree,
                                                    int pDim, bool singleEntries,
                                                    bool blockDiagonal) :
    IDiscretizationT<>(meshes, "SpaceTimeDiscretization"),
    cell_type(SpaceCellType(meshes.fine().cells().Type())), defaultDegree(degree), pDim(pDim),
    singleEntries(singleEntries), blockDiagonal(blockDiagonal) {}

template<typename T, int sDim, int tDim>
void STDiscretizationT<T, sDim, tDim>::CreateAdaptivityLevel(
    const std::unordered_map<Point, DegreePair> &map, LevelPair level) {
  if (IDiscretization::graphs.find(level) != IDiscretization::graphs.end()) {
    THROW("Adaptivitylevel already initialized: " + level.str())
  }
  IDiscretization::graphs[level] = std::move(createMatrixGraph(level, map));
}

template<typename T, int sDim, int tDim>
STDiscretizationT<T, sDim, tDim>::STDiscretizationT(const STDiscretizationT &otherDisc,
                                                    bool singleEntries, bool blockDiagonal) :
    IDiscretizationT<>(otherDisc.GetMeshes(), "SpaceTimeDiscretization"),
    cell_type(otherDisc.cell_type), otherDisc(&otherDisc), defaultDegree(otherDisc.defaultDegree),
    pDim(otherDisc.pDim), singleEntries(singleEntries), blockDiagonal(blockDiagonal) {
  if (cell_type != QUADRILATERAL) { THROW("Currently only QUADRILATERAL are supported.") }
}
