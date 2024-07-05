#include "VtuPlot.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Assertion.hpp"
#include "Buffer.hpp"
#include "Cell.hpp"
#include "Celltype.hpp"
#include "ExchangeBuffer.hpp"
#include "IDiscretization.hpp"
#include "Mesh.hpp"
#include "MixedDoF.hpp"
#include "Parallel.hpp"
#include "PlotUtility.hpp"
#include "Point.hpp"
#include "ProcSet.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"
#include "Vertex.hpp"

Compression globalCompression = Compression::ASCII;

// Compression globalCompression = Compression::BASE64COMPRESSED;
bool globalParallelPlotting = true;
std::function<std::unique_ptr<IPlotBackend>()> globalPlotBackendCreator =
    std::make_unique<DefaultPlotBackend>;

class SerialMeshSynchronization : public MeshSynchronization {
private:
  int communicationSplit;
protected:
  void addCells(const Mesh &mesh) override {
    static constexpr int minCornerCount = 2;
    const int cellCount = mesh.CellCount();
    const int procCount = PPM->Size(communicationSplit);

    std::vector<Point> cells;
    std::vector<CELLTYPE> cellTypes;
    std::vector<size_t> cornerCounts;
    std::vector<Point> connectivities;
    cells.reserve(cellCount);
    cellTypes.reserve(cellCount);
    cornerCounts.reserve(cellCount);
    connectivities.reserve(cellCount * minCornerCount);

    // Collect cell information
    for (cell cell = mesh.cells(); cell != mesh.cells_end(); ++cell) {
      cells.push_back(cell());
      cellTypes.push_back(cell.ReferenceType());
      const int cornerCount = cell.Corners();
      cornerCounts.push_back(cornerCount);
      const auto &corners = cell->second->AsVector();
      connectivities.insert(std::end(connectivities), std::begin(corners), std::end(corners));
    }

    ExchangeBuffer exBuffer(communicationSplit);
    exBuffer.Send(0) << cells << cellTypes << cornerCounts << connectivities;

    // Send the data to proc 0 and merge it
    exBuffer.Communicate();
    if (!PPM->Master(communicationSplit)) { return; }

    CellInfo cellInfo;
    cellInfo.cells.reserve(cellCount * procCount);
    cellInfo.celltype.reserve(cellCount * procCount);
    cellInfo.cornerCount.reserve(cellCount * procCount);
    cellInfo.cellCorner.reserve(cellCount * procCount * minCornerCount);
    for (short q = 0; q < procCount; ++q) {
      std::vector<Point> receivedCells;
      std::vector<CELLTYPE> receivedCellTypes;
      std::vector<size_t> receivedCornerCount;
      std::vector<Point> receivedConnectivities;
      exBuffer.Receive(q) >> receivedCells;
      exBuffer.Receive(q) >> receivedCellTypes;
      exBuffer.Receive(q) >> receivedCornerCount;
      exBuffer.Receive(q) >> receivedConnectivities;

      cellInfo.cells.insert(std::end(cellInfo.cells), std::begin(receivedCells),
                            std::end(receivedCells));
      cellInfo.celltype.insert(std::end(cellInfo.celltype), std::begin(receivedCellTypes),
                               std::end(receivedCellTypes));
      cellInfo.cornerCount.insert(std::end(cellInfo.cornerCount), std::begin(receivedCornerCount),
                                  std::end(receivedCornerCount));
      cellInfo.cellCorner.insert(std::end(cellInfo.cellCorner), std::begin(receivedConnectivities),
                                 std::end(receivedConnectivities));
    }
    cellInfo.timeSteps = mesh.GetTimeMidpointValues();
    backend.AddCells(std::move(cellInfo));
  }
public:
  explicit SerialMeshSynchronization(IPlotBackend &backend, int communicationSplit) :
      MeshSynchronization(backend), communicationSplit(communicationSplit) {}

  void
  AddDeformation(const std::vector<std::pair<Point, Point>> &deformation,
                 const std::vector<std::pair<Point, ProcSet>> & /* possibleDuplicates */) override {
    const int procCount = PPM->Size(communicationSplit);
    ExchangeBuffer exBuffer(communicationSplit);
    exBuffer.Send(0) << deformation;

    exBuffer.Communicate();
    if (!PPM->Master(communicationSplit)) { return; }

    // Add duplicates. 'possibleDuplicates' are ignored as the
    // 'DisplacementBackend' takes care about duplicated points and it also is
    // aware of them as they are all located on the main proc
    std::unordered_map<Point, std::pair<Point, int>> deformationData;
    deformationData.reserve(deformation.size() * procCount);
    for (short q = 0; q < procCount; ++q) {
      std::vector<std::pair<Point, Point>> receivedDisplacements;
      exBuffer.Receive(q) >> receivedDisplacements;
      for (const auto &[point, receivedDisplacement] : receivedDisplacements) {
        auto &[pointDisplacement, count] = deformationData[point];
        pointDisplacement += receivedDisplacement;
        count++;
      }
    }

    backend.AddDisplacement(deformationData);

    BackendPlotData displacementPlotData;
    auto &[points, dataElements, dataSize] = displacementPlotData;
    for (auto &[point, deformed_data_entries] : deformationData) {
      points.push_back(point);
      const auto &[deformed_sum_of_point, deformed_count] = deformed_data_entries;
      const Point p = deformed_sum_of_point / deformed_count;
      for (size_t i = 0; i < SpaceDimension; i++) {
        dataElements.push_back(p[i]);
      }
      dataSize = SpaceDimension;
    }
    backend.AddPointData("DisplacementOnNodes", displacementPlotData);
  }
};

class ParallelMeshSynchronization : public MeshSynchronization {
private:
  int communicationSplit;
protected:
  void addCells(const Mesh &mesh) override {
    static constexpr int minCornerCount = 2;
    const int cellCount = mesh.CellCount();

    CellInfo cellInfo;

    std::vector<Point> &cells = cellInfo.cells;
    std::vector<CELLTYPE> &cellTypes = cellInfo.celltype;
    std::vector<size_t> &cornerCounts = cellInfo.cornerCount;
    std::vector<Point> &connectivity = cellInfo.cellCorner;
    cells.reserve(cellCount);
    cellTypes.reserve(cellCount);
    cornerCounts.reserve(cellCount);
    connectivity.reserve(cellCount * minCornerCount);
    cellInfo.timeSteps = mesh.GetTimeMidpointValues();

    for (cell cell = mesh.cells(); cell != mesh.cells_end(); ++cell) {
      cells.push_back(cell());
      cellTypes.push_back(cell.ReferenceType());
      cornerCounts.push_back(cell.Corners());
      const auto &corners = cell->second->AsVector();
      connectivity.insert(std::end(connectivity), std::begin(corners), std::end(corners));
    }

    connectivity.shrink_to_fit();

    backend.AddCells(std::move(cellInfo));
  }
public:
  explicit ParallelMeshSynchronization(IPlotBackend &backend, int communicationSplit) :
      MeshSynchronization(backend), communicationSplit(communicationSplit) {}

  void AddDeformation(const std::vector<std::pair<Point, Point>> &deformation,
                      const std::vector<std::pair<Point, ProcSet>> &possibleDuplicates) override {
    static constexpr int minCornerCount = 2;
    const int procCount = PPM->Size(communicationSplit);
    const int currentProc = PPM->Proc(communicationSplit);

    // Build initial displacement map for the current proc
    std::unordered_map<Point, std::pair<Point, int>> deformationData;
    deformationData.reserve(deformation.size());
    for (const auto &[point, displacement] : deformation) {
      // Add duplicates on current process local dataset
      auto &[pointDisplacement, count] = deformationData[point];
      pointDisplacement += displacement;
      count++;
    }

    // Collect and communicate possible duplicates in order to create the median
    // of them across processes
    std::vector<std::unordered_map<Point, std::pair<Point, int>>> deformationDataToSend(procCount);
    for (const auto &[point, procSet] : possibleDuplicates) {
      const auto it = deformationData.find(point);
      if (it != std::end(deformationData)) {
        for (const intproc proc : procSet) {
          if (currentProc == proc) {
            // Skip sending to self
            continue;
          }
          deformationDataToSend[proc].insert(*it);
        }
      }
    }

    ExchangeBuffer exBuffer(communicationSplit);
    for (short proc = 0; proc < procCount; proc++) {
      exBuffer.Send(proc) << deformationDataToSend[proc];
    }

    exBuffer.Communicate();
    for (short proc = 0; proc < procCount; ++proc) {
      while (exBuffer.Receive(proc).size() < exBuffer.ReceiveSize(proc)) {
        std::unordered_map<Point, std::pair<Point, int>> receivedDeformationData;
        exBuffer.Receive(proc) >> receivedDeformationData;

        // Add incoming duplicates
        for (const auto &[point, displacementDataElement] : receivedDeformationData) {
          // Skip points which do not belong to this proc
          auto it = deformationData.find(point);
          if (it != std::end(deformationData)) {
            const auto &[receivedPointDisplacement, receivedCount] = displacementDataElement;
            auto &[pointDisplacement, count] = it->second;
            pointDisplacement += receivedPointDisplacement;
            count += receivedCount;
          }
        }
      }
    }

    backend.AddDisplacement(deformationData);

    BackendPlotData displacementPlotData;
    auto &[points, dataElements, dataSize] = displacementPlotData;
    for (auto &[point, deformed_data_entries] : deformationData) {
      points.push_back(point);
      const auto &[deformed_sum_of_point, deformed_count] = deformed_data_entries;
      const Point p = deformed_sum_of_point / deformed_count;
      for (size_t i = 0; i < SpaceDimension; i++) {
        dataElements.push_back(p[i]);
      }
      dataSize = SpaceDimension;
    }
    backend.AddPointData("DisplacementOnNodes", displacementPlotData);
  }
};

class SerialDataSynchronization : public DataSynchronization {
private:
  int communicationSplit;
  const Mesh &mesh;

  void sendAndReceiveData(BackendPlotData &data) const {
    const int procCount = PPM->Size(communicationSplit);

    ExchangeBuffer exBuffer(communicationSplit);
    exBuffer.Send(0) << data;

    exBuffer.Communicate();
    if (!PPM->Master(communicationSplit)) { return; }

    auto &[points, dataElements, _dataSize] = data;
    const size_t pointCount = points.size();
    const size_t dataElementCount = dataElements.size();

    points.clear();
    dataElements.clear();
    points.reserve(pointCount * procCount);
    dataElements.reserve(dataElementCount * procCount);

    for (short q = 0; q < procCount; ++q) {
      BackendPlotData receivedData;
      exBuffer.Receive(q) >> receivedData;
      // Skip data size and data dim as they are expected to be the same on
      // every proc, therefore "data" already contains the correct values
      const auto &[receivedPoints, receivedDataElements, _receivedDataSize] = receivedData;
      points.insert(std::end(points), std::begin(receivedPoints), std::end(receivedPoints));
      dataElements.insert(std::end(dataElements), std::begin(receivedDataElements),
                          std::end(receivedDataElements));
    }
  }
public:
  explicit SerialDataSynchronization(IPlotBackend &backend, const Mesh &mesh,
                                     int communicationSplit) :
      DataSynchronization(backend), communicationSplit(communicationSplit), mesh(mesh) {}

  void AddCellData(BackendPlotData &data) override {
    sendAndReceiveData(data);

    if (!PPM->Master(communicationSplit)) { return; }
    backend.AddCellData(domainName, data);
  }

  void AddPointData(BackendPlotData &data) override {
    auto &[points, dataElements, dataSize] = data;

    // Remove points which do not belong to this process
    std::vector<Point> procPoints;
    std::vector<double> procDataElements;
    procPoints.reserve(points.size());
    procDataElements.reserve(dataElements.size());
    for (int i = 0; i < points.size(); ++i) {
      if (!mesh.GetProcSets().master(points[i])) { continue; }
      procPoints.push_back(points[i]);
      const int off = i * dataSize;
      for (int j = 0; j < dataSize; ++j) {
        procDataElements.push_back(dataElements[off + j]);
      }
    }
    points = procPoints;
    dataElements = procDataElements;

    sendAndReceiveData(data);

    if (!PPM->Master(communicationSplit)) { return; }
    backend.AddPointData(domainName, data);
  }
};

class ParallelDataSynchronization : public DataSynchronization {
public:
  explicit ParallelDataSynchronization(IPlotBackend &backend) : DataSynchronization(backend) {}

  void AddCellData(BackendPlotData &data) override { backend.AddCellData(domainName, data); }

  void AddPointData(BackendPlotData &data) override { backend.AddPointData(domainName, data); }
};

void VtuPlot::initialize(const Mesh &init) {
  // Check if already initialized
  if (meshPtr == &init) { return; }
  if (meshPtr != nullptr) { THROW("Plotting: Meshes don't match!") }

  meshPtr = &init;
  communicationSplit = init.CommSplit();

  // Select synchronization
  IPlotBackend &ownedBackend = *backend;
  if (plotConfig.parallelPlotting) {
    meshSynchronization =
        std::make_unique<ParallelMeshSynchronization>(ownedBackend, init.CommSplit());
    dataSynchronization = std::make_unique<ParallelDataSynchronization>(ownedBackend);
  } else {
    meshSynchronization =
        std::make_unique<SerialMeshSynchronization>(ownedBackend, init.CommSplit());
    dataSynchronization =
        std::make_unique<SerialDataSynchronization>(ownedBackend, init, init.CommSplit());
  }
  // Add mesh data
  meshSynchronization->addCells(init);

  // Add subdomain and procload to plot
  dataSynchronization->SetDomainName("Subdomain");
  const int cellCount = init.CellCount();
  std::vector<Point> points;
  std::vector<double> sdData;
  std::vector<double> procloadData;
  points.reserve(cellCount);
  sdData.reserve(cellCount);
  procloadData.reserve(cellCount);

  for (cell cell = init.cells(); cell != init.cells_end(); ++cell) {
    points.push_back(cell());
    sdData.push_back(cell.Subdomain());
    procloadData.push_back(PPM->Proc(0)); // Is 0 correct?
  }
  std::vector<Point> procloadPoints(points);

  BackendPlotData subdomainBackendData = std::make_tuple(std::move(points), std::move(sdData), 1);
  dataSynchronization->AddCellData(subdomainBackendData);

  dataSynchronization->SetDomainName("ProcLoad");
  BackendPlotData procloadBackendData =
      std::make_tuple(std::move(procloadPoints), std::move(procloadData), 1);
  dataSynchronization->AddCellData(procloadBackendData);
}

void VtuPlot::addData(const string &name, const Vector &data, std::vector<int> indices) {
  initialize(data.GetMesh());

  dataCount++;

  const IDiscretization &discToPlot = data.GetDisc();
  if (indices.empty()) {
    int dim;
    if (typeid(data.GetDoF()) == typeid(MixedDoF)) {
      dim = std::min(3, int(MixedDoF::Cast(data.GetDoF()).NumberOfComponents(0)));
    } else {
      dim = std::min(3, int(data.GetDoF().NumberOfComponents()));
    }

    for (int i = 0; i < dim; ++i) {
      indices.emplace_back(i);
    }
  }

  dataSynchronization->SetDomainName(name);
  discToPlot.transformData(*dataSynchronization, indices, data);
}

// TODO MOVE THIS SOME PLACE MORE NICE
Buffer &operator<<(Buffer &buffer, const PointData &pointData) {
  buffer << pointData.data;
  return buffer;
}

Buffer &operator>>(Buffer &buffer, PointData &pointData) {
  std::vector<double> localVector;
  buffer >> localVector;
  pointData = PointData(std::move(localVector));
  return buffer;
}

// END

void VtuPlot::DeformGrid(const Vector &displacement) {
  initialize(displacement.GetMesh());

  const IDiscretization &discToPlot = displacement.GetDisc();
  discToPlot.transformDisplacement(*meshSynchronization, displacement);
}

VtuPlot &operator<<(VtuPlot &plot, const Vector &v) {
  plot.AddData("value" + std::to_string(plot.DataCount()), v);
  return plot;
}

VtuPlot &operator<<(VtuPlot &plot, PlotStatus status) {
  plot.PlotFile();
  return plot;
}

VtuPlot &operator<<(VtuPlot &plot, const string &name) {
  plot.FileName(name);
  return plot;
}

VtuPlot &operator<<(VtuPlot &plot, const std::pair<std::string, PlotStatus> &saveoptions) {
  return plot << saveoptions.first << saveoptions.second;
}

void VtuPlot::AddData(const string &name, const Vector &data) {
  addData(name, data, std::vector<int>{});
}

void VtuPlot::AddData(const string &name, const Vector &data, int index) {
  addData(name, data, std::vector<int>{index});
}

void VtuPlot::AddDataMultipleEntries(const string &name, const Vector &data,
                                     const vector<int> &indices) {
  for (int i : indices) {
    std::string new_name = name + "_" + std::to_string(i);
    AddData(new_name, data, i);
  }
}

void VtuPlot::AddData(const string &name, const Vector &data, const vector<int> &indices) {
  addData(name, data, indices);
}
