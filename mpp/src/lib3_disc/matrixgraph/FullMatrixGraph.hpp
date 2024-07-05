#ifndef MPP_FULLMATRIXGRAPH_HPP
#define MPP_FULLMATRIXGRAPH_HPP

#include "IMatrixGraph.hpp"

// Todo finish this

class FullMatrixGraph : public IMatrixGraph {
public:
  FullMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof) :
      IMatrixGraph(mesh, std::move(_dof)) {}

  void FullCellsOverlap();
};


#endif // MPP_FULLMATRIXGRAPH_HPP
