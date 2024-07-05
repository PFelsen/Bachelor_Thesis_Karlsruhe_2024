#ifndef VTUPLOT_H
#define VTUPLOT_H

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Vtu11.hpp>

#include "Assertion.hpp"
#include "Celltype.hpp"
#include "Mesh.hpp"
#include "Parallel.hpp"
#include "PlotUtility.hpp"
#include "Point.hpp"
#include "Synchronization.hpp"
#include "Vector.hpp"

typedef std::unordered_map<Point, PointData> PointDataMap;

struct VtuUnstructuredMesh {
  std::uint_fast8_t pointDim_ = 3;
  std::vector<double> points_;
  std::vector<vtu11::VtkIndexType> connectivity_;
  std::vector<vtu11::VtkIndexType> offsets_;
  std::vector<vtu11::VtkCellType> types_;

  const std::vector<double> &points() const { return points_; }

  const std::vector<vtu11::VtkIndexType> &connectivity() const { return connectivity_; }

  const std::vector<vtu11::VtkIndexType> &offsets() const { return offsets_; }

  const std::vector<vtu11::VtkCellType> &types() const { return types_; }

  std::uint_fast8_t pointDimension() const { return pointDim_; }

  size_t numberOfPoints() const { return points_.size() / pointDim_; }

  size_t numberOfCells() const { return types_.size(); }
};

struct VtuInfo {
  VtuUnstructuredMesh mesh;
  std::vector<vtu11::DataSetInfo> dataSetInfo;
  std::vector<vtu11::DataSetData> dataSetData;
  std::string postFix = "";
};

struct PlotMeshProvider {
  std::vector<CELLTYPE> cellTypes;
  std::vector<size_t> connectivity;
  std::vector<size_t> offset;
  std::vector<Point> points;
  std::vector<Point> cells;
};

struct DataEntry {
  vtu11::DataSetData values;
  size_t componentSize;
};

struct PlotDataEntry {
  std::unordered_map<std::string, DataEntry> pointData;
  std::unordered_map<std::string, DataEntry> cellData;
};

enum VtkCellTypes {
  VTK_UNKNOWN = 0,
  VTK_VERTEX = 1,
  VTK_POLY_VERTEX = 2,
  VTK_LINE = 3,
  VTK_POLY_LINE = 4,
  VTK_TRIANGLE = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON = 7,
  VTK_PIXEL = 8,
  VTK_QUAD = 9,
  VTK_TETRA = 10,
  VTK_VOXEL = 11,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE = 13,
  VTK_PYRAMID = 14,
  VTK_QUADRATIC_EDGE = 21,
  VTK_QUADRATIC_TRIANGLE = 22,
  VTK_QUADRATIC_QUAD = 23,
  VTK_QUADRATIC_TETRA = 24,
  VTK_QUADRATIC_HEXAHEDRON = 25
};

struct CellInfo {
  std::vector<Point> cells;
  std::vector<Point> cellCorner;
  std::vector<CELLTYPE> celltype;
  std::vector<size_t> cornerCount;
  std::vector<double> timeSteps;
};

class IPlotBackend {
private:
  // Keep this disabled to force the use of 3D points as long as
  // https://gitlab.kitware.com/paraview/paraview/-/issues/19058
  // remains unresolved
  static constexpr bool kEnableLowerPointDim = false;
public:
  virtual ~IPlotBackend() = default;

  virtual std::vector<VtuInfo> GetVtuData() = 0;

  virtual void AddCells(CellInfo &&cellInfo) = 0;

  virtual void AddCellData(const std::string &name, const BackendPlotData &cellData) = 0;

  virtual void AddPointData(const std::string &name, const BackendPlotData &cellData) = 0;

  virtual void AddDisplacement(const std::unordered_map<Point, std::pair<Point, int>> &vector) {}

  /**
   * @brief Get the ID of the data piece
   *
   * @param communicationSplit
   * @return int ID or -1 if this piece should be ignored
   */
  virtual int GetPieceId(const int communicationSplit, const VtuInfo &info) {
    return PPM->Proc(communicationSplit);
  }

  virtual int GetPieceCount(const int communicationSplit, const VtuInfo &info) {
    return PPM->Size(communicationSplit);
  }

  inline static void buildDataSets(const vtu11::DataSetType type,
                                   const std::unordered_map<std::string, DataEntry> &dataEntrys,
                                   std::vector<vtu11::DataSetInfo> &dataSetInfos,
                                   std::vector<vtu11::DataSetData> &dataSetData) {
    const auto oldSizes = dataSetInfos.size();
    dataSetInfos.resize(oldSizes + dataEntrys.size());
    dataSetData.resize(oldSizes + dataEntrys.size());
    auto infoIterator = dataSetInfos.begin() + oldSizes;
    auto dataIterator = dataSetData.begin() + oldSizes;
    for (const auto &[name, entry] : dataEntrys) {
      *(infoIterator++) = std::make_tuple(name, type, entry.componentSize);
      *(dataIterator++) = entry.values;
    }
  }

  inline static void buildDataSets(const PlotDataEntry &dataEntry,
                                   std::vector<vtu11::DataSetInfo> &dataSetInfos,
                                   std::vector<vtu11::DataSetData> &dataSetData) {
    buildDataSets(vtu11::DataSetType::CellData, dataEntry.cellData, dataSetInfos, dataSetData);
    buildDataSets(vtu11::DataSetType::PointData, dataEntry.pointData, dataSetInfos, dataSetData);
  }

  static constexpr VtkCellTypes CellTypeToVTK(const CELLTYPE type) {
    switch (type) {
    case POINT:
      return VTK_VERTEX;
    case SPACETIME_INTERVAL:
    case INTERVAL:
    case TINTERVAL:
      return VTK_LINE;
    case SPACETIME_TRIANGLE:
    case TRIANGLE:
      return VTK_TRIANGLE;
    case SPACETIME_QUADRILATERAL:
    case QUADRILATERAL:
    case QUADRILATERAL2:
      return VTK_QUAD;
    case TETRAHEDRON:
      return VTK_TETRA;
    case PYRAMID:
      return VTK_PYRAMID;
    case SPACETIME_HEXAHEDRON:
    case HEXAHEDRON:
    case HEXAHEDRON20:
    case HEXAHEDRON27:
      return VTK_HEXAHEDRON;
    default:
      return VTK_UNKNOWN;
    }
  }

  inline static std::vector<double> PointToDouble(const std::span<const Point> points,
                                                  const size_t pointSize = SpaceDimension) {
    Assert(pointSize <= 3);
    std::vector<double> pointData;
    if constexpr (kEnableLowerPointDim) {
      pointData.resize(points.size() * pointSize);
    } else {
      pointData.resize(points.size() * 3);
    }
    size_t currentIndex = 0;
    for (const auto &point : points) {
      for (size_t i = 0; i < pointSize; i++) {
        pointData[currentIndex + i] = point[i];
      }
      if constexpr (kEnableLowerPointDim) {
        currentIndex += pointSize;
      } else {
        for (size_t i = pointSize; i < 3; i++) {
          pointData[currentIndex + i] = 0;
        }
        currentIndex += 3;
      }
    }
    return pointData;
  }

  inline static DataEntry &FindOrAdd(std::unordered_map<std::string, DataEntry> &dataEntries,
                                     const std::string &name, const size_t componentSize) {
    Assert(componentSize > 0);
    auto iteratorFind = dataEntries.find(name);
    if (iteratorFind == std::end(dataEntries)) {
      iteratorFind = dataEntries.emplace(name, DataEntry{{}, componentSize}).first;
    }
    DataEntry &dataEntry = iteratorFind->second;
    Assert(dataEntry.componentSize == componentSize);
    return dataEntry;
  }

  template<std::invocable<const Point &> Invocable>
  inline static void InsertIntoData(DataEntry &outputEntry, const BackendPlotData &plotData,
                                    const Invocable &&invocable) {
    const auto &[pointsToMatch, dataInsert, componentSize] = plotData;
    Assert(outputEntry.componentSize == componentSize);

    const auto beginRawData = outputEntry.values.begin();
    auto beginIncomingPointData = dataInsert.begin();
    for (const Point &point : pointsToMatch) {
      const auto destinationIterator = beginRawData + (invocable(point) * componentSize);
      std::copy(beginIncomingPointData, beginIncomingPointData + componentSize,
                destinationIterator);
      beginIncomingPointData += componentSize;
    }
  }

  inline static std::unordered_map<Point, size_t> createCache(const std::vector<Point> &points) {
    std::unordered_map<Point, size_t> map;
    for (size_t i = 0; i < points.size(); i++) {
      map[points[i]] = i;
    }
    return map;
  }

  inline static void InsertIntoDataLinear(DataEntry &outputEntry, const BackendPlotData &plotData,
                                          const std::vector<Point> &pointTo) {
    outputEntry.values.resize(pointTo.size() * outputEntry.componentSize);
    const auto beginPointIterator = std::begin(pointTo);
    const auto endPointIterator = std::end(pointTo);
    InsertIntoData(outputEntry, plotData, [&](auto &point) {
      const auto found = std::find(beginPointIterator, endPointIterator, point);
      Assert(found != endPointIterator);
      return std::distance(beginPointIterator, found);
    });
  }

  inline static void InsertIntoDataConstant(DataEntry &outputEntry, const BackendPlotData &plotData,
                                            const std::unordered_map<Point, size_t> &pointTo,
                                            const size_t count) {
    outputEntry.values.resize(count * outputEntry.componentSize);
    const auto endIterator = std::end(pointTo);
    InsertIntoData(outputEntry, plotData, [&](auto &point) {
      const auto found = pointTo.find(point);
      Assert(found != endIterator);
      return found->second;
    });
  }

  inline static void buildMesh(const PlotMeshProvider &provider, VtuUnstructuredMesh &mesh,
                               const size_t pointSize = SpaceDimension) {
    mesh.types_.resize(provider.cellTypes.size());
    std::transform(std::begin(provider.cellTypes), std::end(provider.cellTypes),
                   std::begin(mesh.types_), [](auto cellType) { return CellTypeToVTK(cellType); });
    if constexpr (kEnableLowerPointDim) { mesh.pointDim_ = pointSize; }
    mesh.points_ = PointToDouble(provider.points, pointSize);
    mesh.connectivity_ = std::vector<vtu11::VtkIndexType>(std::begin(provider.connectivity),
                                                          std::end(provider.connectivity));
    mesh.offsets_ =
        std::vector<vtu11::VtkIndexType>(std::begin(provider.offset), std::end(provider.offset));
  }
};

class DefaultPlotBackend : public IPlotBackend {
public:
  PlotMeshProvider provider;
  PlotDataEntry entrys;
  std::unordered_map<Point, size_t> cellCache;
  std::unordered_map<Point, size_t> pointCache;

  virtual void AddCells(CellInfo &&cellInfo) override {
    provider.cells = std::move(cellInfo.cells);
    provider.cellTypes = std::move(cellInfo.celltype);
    provider.connectivity.resize(cellInfo.cellCorner.size());
    std::transform(std::begin(cellInfo.cellCorner), std::end(cellInfo.cellCorner),
                   std::begin(provider.connectivity), [&](const Point &point) mutable {
                     if (!pointCache.contains(point)) {
                       const auto index = provider.points.size();
                       provider.points.push_back(point);
                       pointCache[point] = index;
                       return index;
                     }
                     return pointCache.at(point);
                   });
    provider.offset.resize(cellInfo.cornerCount.size());
    std::transform_inclusive_scan(std::begin(cellInfo.cornerCount), std::end(cellInfo.cornerCount),
                                  std::begin(provider.offset), std::plus<>(), std::identity());
    cellCache = createCache(provider.cells);
  }

  virtual std::vector<VtuInfo> GetVtuData() override {
    std::vector<VtuInfo> returnVector(1);
    VtuInfo &info = returnVector[0];
    buildMesh(provider, info.mesh);
    buildDataSets(entrys, info.dataSetInfo, info.dataSetData);
    return returnVector;
  }

  virtual void AddCellData(const std::string &name, const BackendPlotData &cellData) override {
    auto &entry = FindOrAdd(entrys.cellData, name, std::get<2>(cellData));
    InsertIntoDataConstant(entry, cellData, cellCache, provider.cells.size());
  }

  virtual void AddPointData(const std::string &name, const BackendPlotData &cellData) override {
    auto &entry = FindOrAdd(entrys.pointData, name, std::get<2>(cellData));
    InsertIntoDataConstant(entry, cellData, pointCache, provider.points.size());
  }
};

class DefaultSpaceTimeBackend : public DefaultPlotBackend {
public:
  virtual std::vector<VtuInfo> GetVtuData() override {
    std::vector<VtuInfo> returnVector(1);
    VtuInfo &info = returnVector[0];
    buildMesh(provider, info.mesh, SpaceDimension + TimeDimension);
    buildDataSets(entrys, info.dataSetInfo, info.dataSetData);
    for (auto &celltype : info.mesh.types_) {
      if (celltype == VTK_TRIANGLE) {
        celltype = VTK_WEDGE;
      } else if (celltype == VTK_QUAD) {
        celltype = VTK_HEXAHEDRON;
      } else if (celltype == VTK_LINE) {
        celltype = VTK_PIXEL;
      } else {
        Exit("Could not promote the given cell type! Are you using the wrong "
             "backend?");
      }
    }
    return returnVector;
  }
};

class DisplacementBackend : public DefaultPlotBackend {
  std::vector<Point> oldData;
public:
  virtual void AddCells(CellInfo &&cellInfo) override {
    DefaultPlotBackend::AddCells(std::move(cellInfo));
    oldData = provider.points;
  }

  virtual void AddDisplacement(
      const std::unordered_map<Point, std::pair<Point, int>> &deformationData) override {
    std::transform(std::begin(oldData), std::end(oldData), std::begin(provider.points),
                   [&](const Point &point) {
                     auto it = deformationData.find(point);
                     if (it != std::end(deformationData)) {
                       const auto &[d, count] = it->second;
                       return point + (d / count);
                     }
                     return point;
                   });
  }
};

class SpaceTimeBackend : public IPlotBackend {
public:
  std::vector<double> keys;
  std::vector<PlotDataEntry> entrys;
  std::vector<PlotMeshProvider> meshes;
  std::vector<std::unordered_map<Point, size_t>> cellCaches;
  std::vector<std::unordered_map<Point, size_t>> pointCaches;
  std::unordered_map<Point, std::vector<double>> cellToTimeCache;

  static inline std::ptrdiff_t FindOrInsert(const double key, std::vector<double> keys,
                                            std::invocable<double, double> auto &&near,
                                            std::invocable<std::ptrdiff_t> auto &&insertion) {
    const auto keyBegin = keys.begin();
    auto foundIterator = std::lower_bound(keyBegin, keys.end(), key);
    if (foundIterator == keys.end() || !near(*foundIterator, key)) {
      foundIterator = std::upper_bound(foundIterator, keys.end(), key);
      if (foundIterator == keys.end() || !near(*foundIterator, key)) {
        keys.insert(foundIterator, key);
        insertion(std::distance(keyBegin, foundIterator));
      }
    }
    return std::distance(keyBegin, foundIterator);
  }

  std::ptrdiff_t Find(const double key) {
    return FindOrInsert(
        key, keys, [](double a, double b) { return mpp_geometry::isNear(a, b); },
        [&](std::ptrdiff_t distance) {
          throw std::runtime_error("Could not find " + to_string(key));
        });
  }

  virtual void AddCells(CellInfo &&cellInfo) override {
    const auto cellCount = cellInfo.cells.size();
    Assert(cellCount == cellInfo.celltype.size());
    Assert(cellCount == cellInfo.cornerCount.size());
    auto corners = cellInfo.cellCorner.begin();
    const auto timeBegin = cellInfo.timeSteps.begin();
    const auto timeEnd = cellInfo.timeSteps.end();
    std::stable_sort(timeBegin, timeEnd);

    const auto reserveSize = cellInfo.timeSteps.size();
    keys = cellInfo.timeSteps;
    cellCaches.resize(reserveSize);
    pointCaches.resize(reserveSize);
    meshes.resize(reserveSize);
    entrys.resize(reserveSize);
    for (auto [cell, type, corner] =
             std::make_tuple(cellInfo.cells.begin(), cellInfo.celltype.begin(),
                             cellInfo.cornerCount.begin());
         cell < cellInfo.cells.end(); cell++, type++, corner++) {
      const size_t cornerCount = *corner;
      const auto [minIter, maxIter] =
          std::minmax_element(corners, corners + cornerCount,
                              [](auto &point1, auto &point2) { return point1.t() < point2.t(); });
      const auto minTime = minIter->t();
      const auto maxTime = maxIter->t();
      const auto currentTime = cell->t();
      const auto lowerBound =
          std::max(std::lower_bound(timeBegin, timeEnd, minTime) - 1, timeBegin);
      const auto upperBound = std::min(std::upper_bound(lowerBound, timeEnd, maxTime) + 1, timeEnd);
      for (auto iter = lowerBound; iter != upperBound; iter++) {
        const auto currenTime = *iter;
        if (currenTime > maxTime || currenTime < minTime) continue;
        // We found a plottable time index which should contain this cell
        // So we insert
        const auto index = Find(currenTime);
        PlotMeshProvider &provider = meshes[index];
        const auto &cellPoint = *cell;
        auto &values = cellToTimeCache[cellPoint];
        values.push_back(currenTime);
        auto &cellCache = cellCaches[index];
        cellCache[cellPoint] = provider.cells.size();
        provider.cells.push_back(cellPoint);
        provider.cellTypes.push_back(*type);
        const auto previous = provider.offset.empty() ? (size_t)0 : provider.offset.back();
        const auto actualCornerCount = cornerCount / 2;
        provider.offset.push_back(previous + actualCornerCount);
        auto &pointCache = pointCaches[index];
        std::transform(corners, corners + actualCornerCount,
                       std::back_inserter(provider.connectivity), [&](const auto &point) mutable {
                         if (pointCache.contains(point)) {
                           return pointCache[point];
                         } else {
                           const auto nextIndex = provider.points.size();
                           provider.points.push_back(point);
                           pointCache[point] = nextIndex;
                           return nextIndex;
                         }
                       });
      }

      corners += cornerCount;
    }
    auto cacheOutput = cellCaches.begin();
    for (const auto &mesh : meshes) {
      *cacheOutput++ = createCache(mesh.cells);
    }
  }

  virtual std::vector<VtuInfo> GetVtuData() override {
    Assert(keys.size() == entrys.size());
    Assert(entrys.size() == meshes.size());
    std::vector<VtuInfo> infos(keys.size());
    for (size_t i = 0; i < keys.size(); i++) {
      VtuInfo &info = infos[i];
      buildMesh(meshes[i], info.mesh);
      buildDataSets(entrys[i], info.dataSetInfo, info.dataSetData);
      info.postFix = "_" + std::to_string(i);
    }
    return infos;
  }

  virtual void AddCellData(const std::string &name, const BackendPlotData &cellData) override {
    const auto &[points, data, size] = cellData;
    auto dataIter = data.begin();
    const auto keyBegin = keys.begin();
    const auto keyEnd = keys.end();
    for (const auto &point : points) {
      for (const auto timePoint : cellToTimeCache[point]) {
        const auto index = Find(timePoint);
        const auto &cellCache = cellCaches[index];
        const auto cellSize = meshes[index].cells.size();
        auto &entryPair = entrys[index].cellData;
        auto &entry = entryPair[name];
        entry.componentSize = size;
        entry.values.resize(size * cellSize);
        std::copy(dataIter, dataIter + size, entry.values.begin() + size * cellCache.at(point));
      }
      dataIter += size;
    }
  }

  virtual void AddPointData(const std::string &name, const BackendPlotData &pointData) override {
    /**
     * This is currently not possible but could be implemented to support
     * point values in discontinuous cases
     */
    ERROR("Point values not implementet for space time case!");
  }

  int GetPieceId(const int communicationSplit, const VtuInfo &info) override {
    const bool hasData = info.mesh.numberOfCells() > 0;
    const int sum = PPM->ExclusiveScanSum(static_cast<int>(hasData), communicationSplit);
    if (!hasData) { return -1; }

    // ExclusiveScanSum returns undefined values for proc with lowest rank
    if (PPM->Master(communicationSplit)) { return 0; }

    return sum;
  }

  int GetPieceCount(const int communicationSplit, const VtuInfo &info) override {
    const bool hasData = info.mesh.numberOfCells() > 0;
    return PPM->SumOnCommSplit(static_cast<int>(hasData), communicationSplit);
  }
};

enum Format { VTU };

enum Compression { ASCII, BASE64, BASE64COMPRESSED };

using namespace std::string_literals;
const std::array compressionString = {"ascii"s, "base64inline"s, "rawbinarycompressed"s};

template<>
inline Compression parse<Compression>(const std::string &value) {
  const auto found = std::ranges::find(compressionString, value);
  Assert(found != compressionString.end());
  return (Compression)std::distance(compressionString.begin(), found);
}

extern Compression globalCompression;
extern bool globalParallelPlotting;
extern std::function<std::unique_ptr<IPlotBackend>()> globalPlotBackendCreator;

struct PlotConfig {
  Compression compression = globalCompression;
  bool parallelPlotting = globalParallelPlotting;
  std::function<std::unique_ptr<IPlotBackend>()> plotBackendCreator = globalPlotBackendCreator;
  Format format = VTU;

  PlotConfig WithParallelPlotting() {
    parallelPlotting = true;
    return *this;
  }

  PlotConfig WithoutParallelPlotting() {
    parallelPlotting = false;
    return *this;
  }
};

class VtuPlot {
  std::string plotName;
  const PlotConfig plotConfig;
  const std::unique_ptr<IPlotBackend> backend;

  Mesh const *meshPtr = nullptr;
  int dataCount = 0;
  int communicationSplit = 0;

  std::unique_ptr<MeshSynchronization> meshSynchronization;
  std::unique_ptr<DataSynchronization> dataSynchronization;
protected:
  std::string dataPath = Config::GetPlotPath();

  void initialize(const Mesh &init);

  void addData(const string &name, const Vector &data, std::vector<int> indices);

  inline std::map<std::string, std::string>
  getActiveData(const std::vector<vtu11::DataSetInfo> &dataSetInfo, const vtu11::DataSetType type) {
    const auto dataSetEnd = dataSetInfo.end();
    static const std::set<std::string> lookup = {"ProcLoad", "Subdomain"};
    const auto foundCell = std::ranges::find_if(dataSetInfo, [&](const auto &value) {
      // TODO fix this in python script
      if (lookup.contains(std::get<0>(value))) return false;
      return std::get<1>(value) == type;
    });
    std::map<std::string, std::string> cellAttribute;
    if (foundCell != dataSetEnd) {
      const auto &value = *foundCell;
      const auto size = std::get<2>(value);
      cellAttribute[size > 1 ? "Vectors" : "Scalars"] = std::get<0>(value);
    }
    return cellAttribute;
  }
public:
  VtuPlot(const std::string &plotName = "default", const PlotConfig &plotConfig = PlotConfig()) :
      plotName(plotName), plotConfig(plotConfig), backend(plotConfig.plotBackendCreator()) {}

  explicit VtuPlot(const Mesh &init, const std::string &plotName = "default",
                   const PlotConfig &plotConfig = PlotConfig()) : VtuPlot(plotName, plotConfig) {
    initialize(init);
  }

  explicit VtuPlot(const Vector &init, const std::string &plotName = "default",
                   const PlotConfig &plotConfig = PlotConfig()) :
      VtuPlot(init.GetMesh(), plotName, plotConfig) {}

  void PlotFile(const std::string &fileNameIn) {
    const std::vector<VtuInfo> vtuInfos = backend->GetVtuData();
    const std::string &writeMode = compressionString[plotConfig.compression];

    std::vector<vtu11::DataSetInfo> const *headerInfo = nullptr;
    std::map<std::string, std::string> headerPointAttribute;
    std::map<std::string, std::string> headerCellAttribute;

    if (plotConfig.parallelPlotting) {
      // Find first info set which contains all the data we need for the pvtu
      // file and prioritize it
      for (const auto &vtuInfo : vtuInfos) {
        if (vtuInfo.mesh.numberOfCells() > 0) {
          headerInfo = &vtuInfo.dataSetInfo;
          headerPointAttribute = getActiveData(vtuInfo.dataSetInfo, vtu11::DataSetType::PointData);
          headerCellAttribute = getActiveData(vtuInfo.dataSetInfo, vtu11::DataSetType::CellData);
          break;
        }
      }
    }

    // Fallback to first element
    if (headerInfo == nullptr && !vtuInfos.empty()) {
      const auto &vtuInfo = vtuInfos[0];

      headerInfo = &vtuInfo.dataSetInfo;
      const std::map<std::string, std::string> pointAttribute =
          getActiveData(vtuInfo.dataSetInfo, vtu11::DataSetType::PointData);
      const std::map<std::string, std::string> cellAttribute =
          getActiveData(vtuInfo.dataSetInfo, vtu11::DataSetType::CellData);
    }

    for (const auto &vtuInfo : vtuInfos) {
      std::string fileName = fileNameIn + vtuInfo.postFix;

      if (plotConfig.parallelPlotting) {
        const int pieceCount = backend->GetPieceCount(communicationSplit, vtuInfo);
        const int pieceId = backend->GetPieceId(communicationSplit, vtuInfo);

        // Write main pvtu
        if (pieceId == 0) {
          vtu11::writePVtu(dataPath, fileName, vtuInfo.mesh, *headerInfo, pieceCount,
                           headerPointAttribute, headerCellAttribute);
        }

        // Wait until the output folder is created
        PPM->Barrier(communicationSplit);

        if (pieceId >= 0) {
          while (true) {
            try {
              vtu11::writePartition(dataPath, fileName, vtuInfo.mesh, vtuInfo.dataSetInfo,
                                    vtuInfo.dataSetData, pieceId, writeMode);
              break;
            } catch (const std::runtime_error &e) {
              using namespace std::chrono_literals;
              if (errno != ENFILE) { throw; }
              std::this_thread::sleep_for(50ms);
            }
          }
        }
      } else {
        if (!PPM->Master(communicationSplit)) { return; }
        const std::string filepath = dataPath + fileName + ".vtu";
        std::filesystem::create_directories(dataPath);
        vtu11::writeVtu(filepath, vtuInfo.mesh, vtuInfo.dataSetInfo, vtuInfo.dataSetData, writeMode,
                        headerPointAttribute, headerCellAttribute);
      }
    }
  }

  int DataCount() const { return dataCount; }

  void FileName(const string &name) { plotName = name; }

  const std::string &FileName() const { return plotName; }

  void ChangeDataPath(const string &dataPath) { this->dataPath = dataPath; }

  void DeformGrid(const Vector &displacement);

  void AddData(const string &name, const Vector &data);

  void AddData(const string &name, const Vector &data, int index);

  void AddData(const string &name, const Vector &data, const std::vector<int> &indices);

  void AddDataMultipleEntries(const string &name, const Vector &data, const vector<int> &indices);

  void PlotFile() { PlotFile(plotName); }
};

VtuPlot &operator<<(VtuPlot &plot, const Vector &v);

VtuPlot &operator<<(VtuPlot &plot, PlotStatus status);

VtuPlot &operator<<(VtuPlot &plot, const string &name);

VtuPlot &operator<<(VtuPlot &plot, const std::pair<std::string, PlotStatus> &saveoptions);

#endif // VTUPLOT_H
