#include "LagrangeDisplacement.hpp"

#include <utility>
#include <vector>

#include "Mesh.hpp"
#include "Point.hpp"
#include "ProcSet.hpp"
#include "RowIterators.hpp"
#include "Vector.hpp"

void transformLagrangeDisplacement(MeshSynchronization &synchronization,
                                   const Vector &displacement) {
  const Mesh &mesh = displacement.GetMesh();
  const int dim = displacement.dim();
  const int rowCount = displacement.nR();
  const procset &procSetsEnd = displacement.procsets_end();

  std::vector<std::pair<Point, Point>> pointsToDisplacements;
  std::vector<std::pair<Point, ProcSet>> possibleDuplicates;
  pointsToDisplacements.reserve(rowCount);
  possibleDuplicates.reserve(mesh.ProcSetsCount());

  for (row row = displacement.rows(); row != displacement.rows_end(); ++row) {
    const Point &point = row();
    if (displacement.find_vertex(point) != displacement.vertices_end()) {
      Point pointDisplacement;
      for (int j = 0; j < dim; ++j) {
        pointDisplacement[j] = displacement(row, j);
      }
      pointsToDisplacements.emplace_back(point, pointDisplacement);

      // Find possible duplicates across procs. It shouldn't matter for Lagrange
      // but let's better be safe
      const procset &procSet = displacement.find_procset(point);
      if (procSet != procSetsEnd) { possibleDuplicates.emplace_back(point, *procSet); }
    }
  }

  synchronization.AddDeformation(pointsToDisplacements, possibleDuplicates);
}
