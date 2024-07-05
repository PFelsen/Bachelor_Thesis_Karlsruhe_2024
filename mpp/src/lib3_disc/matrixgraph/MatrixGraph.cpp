#include "MatrixGraph.hpp"

MatrixGraph::MatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _dof, int depth) :
    IMatrixGraph(mesh, std::move(_dof)) {
  AddCells(depth);
  AddOverlapCells(depth);
  Init();
}