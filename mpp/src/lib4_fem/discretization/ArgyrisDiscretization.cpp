#include "ArgyrisDiscretization.hpp"

#include <cstdint>
#include <iterator>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "ArgyrisElement.hpp"
#include "Cell.hpp"
#include "Point.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"

void transformArgyrisData(DataSynchronization &synchronization, const std::vector<int> &indices,
                          const Vector &data) {
  static constexpr int minCornerCount = 2;
  const int estimatedPointCount = data.GetMesh().CellCount() * minCornerCount;
  const uint32_t dataSize = indices.size();
  std::unordered_set<const Point *> pointSet;
  std::vector<Point> points;
  std::vector<double> pointData;
  points.reserve(estimatedPointCount);
  pointSet.reserve(estimatedPointCount);
  pointData.reserve(estimatedPointCount * dataSize);

  for (cell cell = data.cells(); cell != data.cells_end(); ++cell) {
    const ArgyrisElement element(data, *cell);
    for (int i = 0; i < cell.Corners(); ++i) {
      const Point &corner = cell.Corner(i);
      if (!pointSet.contains(&corner)) {
        pointSet.insert(&corner);
        points.push_back(corner);

        std::transform(std::begin(indices), std::end(indices), std::back_inserter(pointData),
                       [&](const int value) {
                         return element.Value(cell.LocalCorner(i), data, value);
                       });
      }
    }
  }

  BackendPlotData backendData = std::make_tuple(std::move(points), std::move(pointData), dataSize);
  synchronization.AddPointData(backendData);
}
