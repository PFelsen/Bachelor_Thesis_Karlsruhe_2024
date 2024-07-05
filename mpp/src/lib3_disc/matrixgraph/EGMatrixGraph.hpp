#ifndef EGMATRIXGRAPH_HPP
#define EGMATRIXGRAPH_HPP

#include "IMatrixGraph.hpp"

class EGMatrixGraph : public IMatrixGraph {
protected:
  void AddNBCells();

  void CellsOverlapEG();
public:
  EGMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof, int depth = 1);
};

#endif // EGMATRIXGRAPH_HPP
