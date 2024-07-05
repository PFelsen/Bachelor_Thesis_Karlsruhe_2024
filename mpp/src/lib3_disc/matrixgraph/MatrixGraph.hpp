#ifndef MATRIXGRAPH_HPP
#define MATRIXGRAPH_HPP

#include "IMatrixGraph.hpp"

class MatrixGraph : public IMatrixGraph {
public:
  MatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof, int depth = 1);
};

#endif // MATRIXGRAPH_HPP
