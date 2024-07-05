#include "DGDiscretization.hpp"

#include <utility>
#include <vector>

#include "Cell.hpp"
#include "DGVectorFieldElement.hpp"
#include "Mesh.hpp"
#include "Point.hpp"
#include "ProcSet.hpp"
#include "Vector.hpp"
#include "VectorField.hpp"

void transformDGDisplacement(MeshSynchronization &synchronization, const Vector &displacement) {
  static constexpr int minCornerCount = 2;
  const Mesh &mesh = displacement.GetMesh();
  const int dim = displacement.dim();
  const int rowCount = displacement.nR();
  const procset &procSetsEnd = displacement.procsets_end();

  std::vector<std::pair<Point, Point>> pointsToDisplacements;
  std::vector<std::pair<Point, ProcSet>> possibleDuplicates;
  pointsToDisplacements.reserve(rowCount * minCornerCount);
  possibleDuplicates.reserve(displacement.OverlapCount() * minCornerCount);

  for (row row = displacement.rows(); row != displacement.rows_end(); ++row) {
    cell cell = displacement.find_cell(row());
    if (cell != displacement.cells_end()) {
      DGVectorFieldElement element(displacement, *cell);
      const procset &procSet = displacement.find_procset(cell());

      for (int cornerIdx = 0; cornerIdx < cell.Corners(); ++cornerIdx) {
        // Collect displacements
        const Point &corner = cell.Corner(cornerIdx);
        Point cornerDisplacement;
        VectorField field = element.VectorValue(cell.LocalCorner(cornerIdx), displacement);
        for (int j = 0; j < dim; ++j) {
          cornerDisplacement[j] = field[j];
        }
        pointsToDisplacements.emplace_back(corner, cornerDisplacement);

        // Find possible duplicates across procs
        if (procSet != procSetsEnd) { possibleDuplicates.emplace_back(corner, *procSet); }
      }
    }
  }

  synchronization.AddDeformation(pointsToDisplacements, possibleDuplicates);
}
