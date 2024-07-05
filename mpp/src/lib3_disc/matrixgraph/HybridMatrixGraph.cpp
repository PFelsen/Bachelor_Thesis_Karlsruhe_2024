#include "HybridMatrixGraph.hpp"

HybridMatrixGraph::HybridMatrixGraph(const Mesh &mesh, std::unique_ptr<IDoF> _vectorDoF,
                                     std::unique_ptr<IDoF> _shapeDoF) :
    MatrixGraph(mesh, std::move(_vectorDoF)), shapeDoF(std::move(_shapeDoF)) {}