#include "RTLagrangeDiscretization.hpp"

#include <cstdint>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Cell.hpp"
#include "Point.hpp"
#include "RTLagrangeElement.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"
#include "VectorField.hpp"

void transformRTLagrangeData(DataSynchronization &synchronization, const std::vector<int> &indices,
                             const Vector &data, const int degree) {
  static constexpr int minCornerCount = 2;
  const int cellCount = data.GetMesh().CellCount();
  const int estimatedPointCount = cellCount * minCornerCount;
  constexpr uint32_t pressureDataSize = 1;
  constexpr uint32_t velocityDataSize = SpaceDimension;
  const std::string &domainName = synchronization.GetDomainName();

  std::unordered_set<const Point *> pointSet;
  std::vector<Point> pressurePoints;
  std::vector<double> pressurePointData;
  std::vector<Point> velocityPoints;
  std::vector<double> velocityPointData;
  velocityPoints.reserve(estimatedPointCount);
  pointSet.reserve(estimatedPointCount);
  velocityPointData.reserve(estimatedPointCount * velocityDataSize);

  synchronization.SetDomainName(domainName + " - pressure");
  if (degree == 0) {
    for (cell cell = data.cells(); cell != data.cells_end(); ++cell) {
      // Velocity cell data
      const RTLagrangeElement element(data, *cell);
      for (int i = 0; i < cell.Corners(); ++i) {
        const Point &corner = cell.Corner(i);
        if (!pointSet.contains(&corner)) {
          pointSet.insert(&corner);
          velocityPoints.push_back(corner);

          const VectorField velocityField = element.VelocityField(cell.LocalCorner(i), data);
          for (int i = 0; i < velocityDataSize; i++) {
            velocityPointData.push_back(velocityField[i]);
          }
        }
      }

      // Pressure cell data
      std::unordered_set<const Point *> pressurePointSet;
      pressurePointSet.reserve(cellCount);
      pressurePoints.reserve(cellCount);
      pressurePointData.reserve(cellCount);
      const Point &center = cell.Center();
      if (!pressurePointSet.contains(&center)) {
        pressurePointSet.insert(&center);
        pressurePoints.push_back(center);
        pressurePointData.push_back(element.PressureValue(cell.LocalCenter(), data));
      }
    }

    BackendPlotData pressureBackendData =
        std::make_tuple(std::move(pressurePoints), std::move(pressurePointData), pressureDataSize);
    synchronization.AddCellData(pressureBackendData);
  } else {
    pressurePoints.reserve(estimatedPointCount);
    pressurePointData.reserve(estimatedPointCount);

    for (cell cell = data.cells(); cell != data.cells_end(); ++cell) {
      const RTLagrangeElement element(data, *cell);
      for (int i = 0; i < cell.Corners(); ++i) {
        const Point &corner = cell.Corner(i);
        if (!pointSet.contains(&corner)) {
          pointSet.insert(&corner);
          velocityPoints.push_back(corner);
          pressurePoints.push_back(corner);

          // Velocity cell data
          const VectorField velocityField = element.VelocityField(cell.LocalCorner(i), data);
          for (int i = 0; i < velocityDataSize; i++) {
            velocityPointData.push_back(velocityField[i]);
          }

          // Pressure point data
          pressurePointData.push_back(element.PressureValue(cell.LocalCorner(i), data));
        }
      }
    }

    BackendPlotData pressureBackendData =
        std::make_tuple(std::move(pressurePoints), std::move(pressurePointData), pressureDataSize);

    synchronization.AddPointData(pressureBackendData);
  }

  synchronization.SetDomainName(domainName + " - velocity");
  BackendPlotData velocityBackendData =
      std::make_tuple(std::move(velocityPoints), std::move(velocityPointData), velocityDataSize);
  synchronization.AddPointData(velocityBackendData);
  synchronization.SetDomainName(domainName);
}
