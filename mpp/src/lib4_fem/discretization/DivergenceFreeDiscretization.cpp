#include "DivergenceFreeDiscretization.hpp"

#include <cstdint>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Cell.hpp"
#include "DivergenceFreeElement.hpp"
#include "Point.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"
#include "VectorField.hpp"

void transformDivergenceFreeData(DataSynchronization &synchronization,
                                 const std::vector<int> &indices, const Vector &data) {
  static constexpr int minCornerCount = 2;
  const int estimatedPointCount = data.GetMesh().CellCount() * minCornerCount;
  constexpr uint32_t dataSize = 2;

  std::unordered_set<const Point *> pointSet;
  std::vector<Point> points;
  std::vector<double> pointData;
  points.reserve(estimatedPointCount);
  pointSet.reserve(estimatedPointCount);
  pointData.reserve(estimatedPointCount * dataSize);

  for (cell cell = data.cells(); cell != data.cells_end(); ++cell) {
    const DivergenceFreeElement element(data, *cell);
    for (int i = 0; i < cell.Corners(); ++i) {
      const Point &corner = cell.Corner(i);
      if (!pointSet.contains(&corner)) {
        pointSet.insert(&corner);
        points.push_back(corner);

        const Velocity velocity = element.VelocityField(cell.LocalCorner(i), data);
        pointData.push_back(velocity[0]);
        pointData.push_back(velocity[1]);
      }
    }
  }

  BackendPlotData backendData = std::make_tuple(std::move(points), std::move(pointData), dataSize);
  synchronization.AddPointData(backendData);
}
