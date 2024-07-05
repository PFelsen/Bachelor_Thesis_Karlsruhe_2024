#include "IDiscretization.hpp"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <tuple>
#include <utility>
#include <vector>

#include "Mesh.hpp"
#include "Point.hpp"
#include "RowIterators.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"

void transformGenericCellData(DataSynchronization &synchronization, const std::vector<int> &indices,
                              const Vector &data) {
  const Mesh &mesh = data.GetMesh();
  const int rowCount = data.nR();
  const uint32_t dataSize = indices.size();

  std::vector<Point> points;
  std::vector<double> cellData;
  points.reserve(rowCount);
  cellData.reserve(rowCount * dataSize);

  for (row row = data.rows(); row != data.rows_end(); ++row) {
    const Point cell = row();
    if (data.find_cell(cell) != data.cells_end()) {
      points.push_back(cell);
      std::transform(std::begin(indices), std::end(indices), std::back_inserter(cellData),
                     [&](const int value) { return data(row, value); });
    }
  }

  BackendPlotData backendData = std::make_tuple(std::move(points), std::move(cellData), dataSize);
  synchronization.AddCellData(backendData);
}

void transformGenericPointData(DataSynchronization &synchronization,
                               const std::vector<int> &indices, const Vector &data) {
  const Mesh &mesh = data.GetMesh();
  const int rowCount = data.nR();
  const uint32_t dataSize = indices.size();

  std::vector<Point> points;
  std::vector<double> pointData;
  points.reserve(rowCount);
  pointData.reserve(rowCount * dataSize);

  for (row row = data.rows(); row != data.rows_end(); ++row) {
    const Point point = row();
    if (data.find_vertex(point) != data.vertices_end()) {
      points.push_back(point);

      std::transform(std::begin(indices), std::end(indices), std::back_inserter(pointData),
                     [&](const int value) { return data(row, value); });
    }
  }

  BackendPlotData backendData = std::make_tuple(std::move(points), std::move(pointData), dataSize);
  synchronization.AddPointData(backendData);
}
