#ifndef HYBRIDMATRIXGRAPH_HPP
#define HYBRIDMATRIXGRAPH_HPP

#include "MatrixGraph.hpp"

class HybridMatrixGraph : public MatrixGraph {
  std::unique_ptr<IDoF> shapeDoF = nullptr;
public:
  HybridMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _vectorDoF,
                    std::unique_ptr<IDoF> _shapeDoF);

  const IDoF &GetShapeDoF() const override { return *shapeDoF; }
};

#endif // HYBRIDMATRIXGRAPH_HPP
