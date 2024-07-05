#ifndef DGMATRIXGRAPH_HPP
#define DGMATRIXGRAPH_HPP

#include "IMatrixGraph.hpp"

class DGMatrixGraph : public IMatrixGraph {
protected:
  void AddNBCells();

  void AddCrossPoints();

  void SendProcSets();

  void CellsOverlap_dG1();
public:
  DGMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof);
};

#endif // DGMATRIXGRAPH_HPP
