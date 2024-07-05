#include "RTDiscretization.hpp"

#include <cstdint>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Cell.hpp"
#include "Point.hpp"
#include "RTElement.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"
#include "VectorField.hpp"

void transformRTData(DataSynchronization &synchronization, const std::vector<int> &indices,
                     const Vector &data) {
  static constexpr int minCornerCount = 2;
  const int estimatedPointCount = data.GetMesh().CellCount() * minCornerCount;
  constexpr uint32_t dataSize = SpaceDimension;

  std::unordered_set<const Point *> pointSet;
  std::vector<Point> points;
  std::vector<double> pointData;
  points.reserve(estimatedPointCount);
  pointSet.reserve(estimatedPointCount);
  pointData.reserve(estimatedPointCount * dataSize);

  for (cell cell = data.cells(); cell != data.cells_end(); ++cell) {
    const RTElement element(data, *cell);
    for (int i = 0; i < cell.Corners(); ++i) {
      const Point &corner = cell.Corner(i);
      if (!pointSet.contains(&corner)) {
        pointSet.insert(&corner);
        points.push_back(corner);

        const Velocity velocity = element.VelocityField(cell.LocalCorner(i), data);
        for (int i = 0; i < dataSize; i++) {
          pointData.push_back(velocity[i]);
        }
      }
    }
  }

  BackendPlotData backendData = std::make_tuple(std::move(points), std::move(pointData), dataSize);
  synchronization.AddPointData(backendData);
}
