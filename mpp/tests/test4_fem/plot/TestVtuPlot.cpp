#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Cell.hpp>
#include <Config.hpp>
#include <DGDiscretization.hpp>
#include <IDiscretization.hpp>
#include <LagrangeDiscretization.hpp>
#include <Mesh.hpp>
#include <Meshes.hpp>
#include <MeshesCreator.hpp>
#include <Plotting.hpp>
#include <Point.hpp>
#include <ProcSet.hpp>
#include <RowIterators.hpp>
#include <Vector.hpp>
#include <VtuPlot.hpp>
#include <gtest/gtest-param-test.h>
#include <pugiconfig.hpp>
#include <pugixml.hpp>

#if USE_SPACETIME
#include <SpaceTimeDiscretization_DGDG.hpp>
#include <SpaceTimeTools.hpp>
#endif

#include "TestEnvironment.hpp"

namespace fs = std::filesystem;

namespace {

/** @see IPlotBackend#kEnableLowerPointDim */
constexpr bool kEnableLowerPointDim = false;

void checkXmlAttributes(
    const std::vector<std::pair<std::string, std::string>> &attributeExpectedPairs,
    const pugi::xml_node &node) {
  for (const auto &[name, expected] : attributeExpectedPairs) {
    EXPECT_EQ(node.attribute(name.c_str()).as_string(), expected)
        << "Node: " << node << ", for attribute=" << name;
  }
}

std::vector<double> ParseStrDataArray(const std::string &values) {
  auto itrBeforeNum = values.begin();
  std::vector<double> output;

  while (itrBeforeNum != values.end()) {
    const auto itrAfterNum = std::find(itrBeforeNum + 1, values.end(), ' ');
    std::string number;
    number.reserve(std::distance(itrBeforeNum, itrAfterNum));
    for (auto iter = itrBeforeNum; iter != itrAfterNum; iter++) {
      const auto character = *iter;
      if (!std::isspace(character)) number.push_back(character);
    }
    if (!number.empty()) { output.push_back(std::stod(number)); }
    itrBeforeNum = itrAfterNum;
  }
  return output;
}

std::vector<Point> DataArrayToPoints(const std::vector<double> &data, const size_t dim) {
  const size_t dataSize = data.size();
  EXPECT_EQ(0, dataSize % dim);

  const size_t numPoints = dataSize / dim;
  std::vector<Point> points(numPoints);

  auto dataIt = std::begin(data);
  auto pointIt = std::begin(points);

  for (size_t i = 0; i < numPoints; ++i) {
    Point point;
    for (size_t j = 0; j < dim; ++j) {
      point[j] = *dataIt++;
    }
    *pointIt++ = point;
  }

  return points;
}

} // namespace

class TestVtuPlot : public TestWithParam<std::string> {
protected:
  std::unique_ptr<Meshes> meshes;
  std::unique_ptr<VtuPlot> plot;
  const Mesh &mesh;
  std::string file;
  int numVertices;
  int numCorners;
  int numCells;
  pugi::xml_document doc;
  pugi::xml_node fileNode;
  pugi::xml_node pieceNode;
  pugi::xml_node cellsNode;
  pugi::xml_node pointsNode;
  pugi::xml_node cellData;
  pugi::xml_node pointData;

  std::string format = "ascii";

  TestVtuPlot(
      MeshesCreator &&meshesCreator = MeshesCreator(GetParam()).WithPLevel(4).WithLevel(4)) :
      meshes(meshesCreator.CreateUnique()), mesh(meshes->fine()) {
    const auto folder = Config::GetPlotPath();
    if (PPM->Master()) {
      if (fs::exists(folder)) fs::remove_all(folder);
    }
    PPM->Barrier(0);

    file = folder + GetParam() + ".vtu";

    numCells = mesh.SpaceCellCountGeometry(0);
    numCorners = mesh.cells().SpaceCell().Corners();
    numVertices = mesh.SpaceVertexCountGeometry();

    // meshes->PrintInfo();
    std::pair<double, double> xBorders = {0, 2000};
    std::pair<double, double> yBorders = {0, 2000};
    int dimX = 11;
    int dimY = 11;
    std::vector<int> boundary = {2, 2, 2, 2};
    CoarseGeometry Test = RectangularGeometry2D(xBorders, yBorders, dimX, dimY, boundary);
  }

  pugi::xml_node checkDataArrayAttributes(
      const std::string &name, const pugi::xml_node &parent, const std::string &type = "Float32",
      const std::vector<std::pair<std::string, std::string>> &more = {}) const {
    const pugi::xml_node currentNode =
        parent.find_child_by_attribute("DataArray", "Name", name.c_str());
    EXPECT_NE(currentNode.type(), pugi::node_null)
        << "Node not found with Name=" << name << " in Parent=" << parent.name();
    std::vector vector = more;
    vector.emplace_back("Name", name);
    vector.emplace_back("format", format);
    vector.emplace_back("type", type);
    checkXmlAttributes(vector, currentNode);
    return currentNode;
  }

  std::vector<double>
  checkDataArray(const std::string &name, const pugi::xml_node &parent, const size_t size,
                 const std::string &type = "Float32",
                 const std::vector<std::pair<std::string, std::string>> &more = {}) const {
    const auto currentNode = checkDataArrayAttributes(name, parent, type, more);
    const auto dataRead = ParseStrDataArray(currentNode.child_value());
    EXPECT_EQ(dataRead.size(), size) << "Node: " << name;
    return dataRead;
  }

  void checkDataArray(const std::string &name, const pugi::xml_node &parent, const size_t size,
                      const double data, const std::string &type = "Float32",
                      const std::vector<std::pair<std::string, std::string>> &more = {}) const {
    const auto actualData = checkDataArray(name, parent, size, type, more);
    for (const double actual : actualData) {
      ASSERT_DOUBLE_EQ(actual, data) << "Node " << name;
    }
  }

  void checkDataArray(const std::string &name, const pugi::xml_node &parent,
                      const std::vector<double> &data, const std::string &type = "Float32") const {
    const auto actualData = checkDataArray(name, parent, data.size(), type);
    for (size_t i = 0; i < data.size(); i++) {
      ASSERT_DOUBLE_EQ(data[i], actualData[i]) << "Node: " << name << std::endl << "Index: " << i;
    }
  }

  void ReadFile() {
    if (!PPM->Master(0)) return;
    ASSERT_TRUE(fs::exists(file)) << file;
    const auto result = doc.load_file(file.c_str());
    ASSERT_TRUE(result) << file;

    fileNode = doc.child("VTKFile");
    pieceNode = fileNode.child("UnstructuredGrid").child("Piece");
    cellsNode = pieceNode.child("Cells");
    pointsNode = pieceNode.child("Points");
    cellData = pieceNode.child("CellData");
    pointData = pieceNode.child("PointData");
  }

  void TestMesh() {
    if (!PPM->Master(0)) return;

    checkXmlAttributes({{"byte_order", "LittleEndian"},
                        {"type", "UnstructuredGrid"},
                        {"version", "0.1"}},
                       fileNode);

    checkXmlAttributes({{"NumberOfCells", std::to_string(numCells)},
                        {"NumberOfPoints", std::to_string(numVertices)}},
                       pieceNode);

    checkDataArray("connectivity", cellsNode, numCells * numCorners, "Int64");

    std::vector<double> offsetExpected;
    offsetExpected.resize(numCells);
    for (size_t i = 0; i < numCells; ++i) {
      offsetExpected[i] = (i + 1) * numCorners;
    }
    checkDataArray("offsets", cellsNode, offsetExpected, "Int64");

    const double cellTypeExpected = VtkCellType(SpaceCellType(meshes->fine().cells().Type()));
    checkDataArray("types", cellsNode, numCells, cellTypeExpected, "Int8");

    checkDataArray("Subdomain", cellData, numCells, 1);

    // checkDataArray("ProcLoad", cellData, numCells, 0);
  }
};

class TestVtuPlotWithMesh : public TestVtuPlot {
public:
  TestVtuPlotWithMesh() : TestVtuPlot() {
    plot = std::make_unique<VtuPlot>(mesh);
    PPM->Barrier();
    plot->PlotFile(GetParam());
    PPM->Barrier();

    ReadFile();
  }
};

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithMesh,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                // #ifdef USE_SPACETIME
                                //   , "SpaceTimeSquare"
                                // #endif
                                ));

TEST_P(TestVtuPlotWithMesh, TestVtuFile) { TestMesh(); }

class TestVtuPlotWithLagrangeDisc : public TestVtuPlot {
protected:
  std::shared_ptr<const LagrangeDiscretization> disc;
  Vector v;
  int dataSize;
  int degree;
public:
  TestVtuPlotWithLagrangeDisc(int degree, int size) :
      TestVtuPlot(), disc(std::make_shared<LagrangeDiscretization>(*meshes, degree, size)),
      v(Vector(1.0, disc)), degree(degree), dataSize(size) {
    plot = std::make_unique<VtuPlot>(v);
    plot->AddData("testvalue", v);
    PPM->Barrier();
    plot->PlotFile(GetParam());
    PPM->Barrier();

    ReadFile();
  }

  TestVtuPlotWithLagrangeDisc(int degree, int size, int index) :
      TestVtuPlot(), disc(std::make_shared<LagrangeDiscretization>(*meshes, degree, size)),
      v(Vector(1.0, disc)), degree(degree), dataSize(1) {
    plot = std::make_unique<VtuPlot>(v);
    plot->AddData("testvalue", v, index);
    PPM->Barrier();
    plot->PlotFile(GetParam());
    PPM->Barrier();

    ReadFile();
  }

  TestVtuPlotWithLagrangeDisc(int degree, int size, const std::vector<int> &index) :
      TestVtuPlot(), disc(std::make_shared<LagrangeDiscretization>(*meshes, degree, size)),
      v(Vector(1.0, disc)), degree(degree), dataSize(index.size()) {
    plot = std::make_unique<VtuPlot>(v);
    plot->AddDataMultipleEntries("testvalue", v, index);
    PPM->Barrier();
    plot->PlotFile(GetParam());
    PPM->Barrier();

    ReadFile();
  }

  void TestVertexScalarData() {
    if (!PPM->Master(0)) return;
    if constexpr (kEnableLowerPointDim) {
      checkXmlAttributes({{"NumberOfComponents", std::to_string(SpaceDimension)},
                          {"format", format},
                          {"type", "Float64"}},
                         pointsNode.child("DataArray"));
    } else {
      checkXmlAttributes({{"NumberOfComponents", "3"}, {"format", format}, {"type", "Float64"}},
                         pointsNode.child("DataArray"));
    }

    checkDataArray("testvalue", pointData, numVertices * dataSize, 1);
  }

  void TestCellData() {
    if (!PPM->Master(0)) return;
    if constexpr (kEnableLowerPointDim) {
      checkXmlAttributes({{"NumberOfComponents", std::to_string(SpaceDimension)},
                          {"format", format},
                          {"type", "Float64"}},
                         pointsNode.child("DataArray"));
    } else {
      checkXmlAttributes({{"NumberOfComponents", "3"}, {"format", format}, {"type", "Float64"}},
                         pointsNode.child("DataArray"));
    }
    checkDataArray("testvalue", cellData, numCells * dataSize, 1);
  }
};

class TestVtuPlotWithLagrangeScalarData : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeScalarData() : TestVtuPlotWithLagrangeDisc(1, 1) {}
};

TEST_P(TestVtuPlotWithLagrangeScalarData, TestVtuFile) { TestVertexScalarData(); }

class TestVtuPlotWithLagrangeVectorData : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeVectorData() : TestVtuPlotWithLagrangeDisc(1, 3) {}
};

TEST_P(TestVtuPlotWithLagrangeVectorData, TestVtuFile) { TestVertexScalarData(); }

class TestVtuPlotWithLagrangeDoFVertexDataSingle : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeDoFVertexDataSingle() : TestVtuPlotWithLagrangeDisc(1, 5, 2) {}
};

TEST_P(TestVtuPlotWithLagrangeDoFVertexDataSingle, TestVtuFile) { TestVertexScalarData(); }

class TestVtuPlotWithLagrangeDoFVertexDataMultiple : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeDoFVertexDataMultiple() : TestVtuPlotWithLagrangeDisc(1, 7) {}
};

TEST_P(TestVtuPlotWithLagrangeDoFVertexDataMultiple, TestVtuFile) {
  // TODO
}

class TestVtuPlotWithLagrangeDoFCellScalar : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeDoFCellScalar() : TestVtuPlotWithLagrangeDisc(0, 1) {}
};

TEST_P(TestVtuPlotWithLagrangeDoFCellScalar, TestVtuFile) { TestCellData(); }

class TestVtuPlotWithLagrangeDoFCellVector : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeDoFCellVector() : TestVtuPlotWithLagrangeDisc(0, 3) {}
};

TEST_P(TestVtuPlotWithLagrangeDoFCellVector, TestVtuFile) { TestCellData(); }

class TestVtuPlotWithLagrangeDoFCellScalarSingle : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeDoFCellScalarSingle() : TestVtuPlotWithLagrangeDisc(0, 5, 2) {}
};

TEST_P(TestVtuPlotWithLagrangeDoFCellScalarSingle, TestVtuFile) { TestCellData(); }

class TestVtuPlotWithLagrangeDoFCellMultiple : public TestVtuPlotWithLagrangeDisc {
public:
  TestVtuPlotWithLagrangeDoFCellMultiple() : TestVtuPlotWithLagrangeDisc(0, 7, {0, 2, 3, 5}) {}
};

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeScalarData,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeVectorData,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeDoFVertexDataSingle,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeDoFVertexDataMultiple,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeDoFCellScalar,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeDoFCellVector,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithLagrangeDoFCellScalarSingle,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

std::vector<double> createTestVector(const double offset) {
  std::vector<double> vector;
  vector.resize(SpaceDimension + TimeDimension);
  for (size_t i = 0; i < SpaceDimension + TimeDimension; i++) {
    vector[i] = i * offset;
  }
  return vector;
}

TEST(TestVtuBackend, pointToDouble) {
  std::vector<Point> points;
  std::vector<double> emptyData = IPlotBackend::PointToDouble(points);
  EXPECT_TRUE(emptyData.empty());

  EXPECT_TRUE(IPlotBackend::PointToDouble({}).empty());

  std::vector<std::vector<double>> pointData = {createTestVector(0), createTestVector(0.5),
                                                createTestVector(0.8)};
  points.push_back(Point(pointData[0]));
  points.push_back(Point(pointData[1]));
  points.push_back(Point(pointData[2]));

  std::vector<double> rawData = IPlotBackend::PointToDouble(points);
  if constexpr (kEnableLowerPointDim) {
    EXPECT_EQ(rawData.size(), points.size() * SpaceDimension);
    for (size_t i = 0; i < pointData.size(); i++) {
      const auto &data = pointData[i];
      for (size_t j = 0; j < SpaceDimension; j++) {
        EXPECT_EQ(rawData[j + i * SpaceDimension], data[j]);
      }
    }
  } else {
    EXPECT_EQ(rawData.size(), points.size() * 3);
    for (size_t i = 0; i < pointData.size(); i++) {
      const auto &data = pointData[i];
      for (size_t j = 0; j < SpaceDimension; j++) {
        EXPECT_EQ(rawData[j + i * 3], data[j]);
      }
      for (size_t j = SpaceDimension; j < 3; j++) {
        EXPECT_EQ(rawData[j + i * 3], 0);
      }
    }
  }
}

std::vector<Point> makePointList(const std::vector<double> &allDouble) {
  std::vector<Point> points;
  points.reserve(allDouble.size());
  for (double d : allDouble)
    points.push_back(Point(createTestVector(d)));
  return points;
}

TEST(TestVtuBackend, IPlotBackendDefaultTests) {
  std::vector<Point> points = makePointList({0, 2, 0.444});

  BackendPlotData data;
  auto &[pointsToAdd, values, component] = data;
  component = 3;
  pointsToAdd.resize(points.size());
  values.resize(component * points.size());
  for (size_t i = 0; i < points.size(); i++) {
    pointsToAdd[i] = points[points.size() - i - 1];
    values[component * i] = -((double)i);
  }

  auto cache = IPlotBackend::createCache(points);
  EXPECT_EQ(cache.size(), points.size());

  for (size_t i = 0; i < points.size(); i++) {
    EXPECT_EQ(cache[points[i]], i);
  }

  DataEntry entry;
  entry.componentSize = 3;
  IPlotBackend::InsertIntoDataConstant(entry, data, cache, points.size());
  EXPECT_EQ(entry.values.size(), points.size() * component);
  for (size_t i = 0; i < points.size(); i++) {
    const double index = points.size() - i - 1;
    EXPECT_DOUBLE_EQ(entry.values[i * component], -index);
    EXPECT_DOUBLE_EQ(entry.values[i * component + 1], 0);
    EXPECT_DOUBLE_EQ(entry.values[i * component + 2], 0);
  }

  entry.values.clear();
  IPlotBackend::InsertIntoDataLinear(entry, data, points);
  EXPECT_EQ(entry.values.size(), points.size() * component);
  for (size_t i = 0; i < points.size(); i++) {
    const double index = points.size() - i - 1;
    EXPECT_DOUBLE_EQ(entry.values[i * component], -index);
    EXPECT_DOUBLE_EQ(entry.values[i * component + 1], 0);
    EXPECT_DOUBLE_EQ(entry.values[i * component + 2], 0);

    std::unordered_map<std::string, DataEntry> dataEntryForNames;
    auto &foundEntry = IPlotBackend::FindOrAdd(dataEntryForNames, "test", 3);
    EXPECT_EQ(foundEntry.componentSize, 3);
    EXPECT_TRUE(foundEntry.values.empty());

    foundEntry.values.resize(6);
    std::fill(foundEntry.values.begin(), foundEntry.values.end(), 2);

    auto &foundEntryNext = IPlotBackend::FindOrAdd(dataEntryForNames, "test", 3);
    EXPECT_EQ(foundEntryNext.componentSize, foundEntry.componentSize);
    EXPECT_EQ(foundEntryNext.values, foundEntry.values);
    EXPECT_EQ(foundEntryNext.values.size(), 6);
  }
}

void testDataFind(std::vector<Point> &points,
                  std::function<void(std::string &, BackendPlotData &)> &&function,
                  std::unordered_map<std::string, DataEntry> &dataEntry) {
  for (size_t i = 0; i < points.size(); i++) {
    BackendPlotData internalBackendData;
    auto &[currentPoints, data, component] = internalBackendData;
    component = 3;
    const auto reference = points.size() - i - 1;
    currentPoints.push_back(points[reference]);
    data.reserve(component);
    for (size_t j = 0; j < component; j++) {
      data.push_back((double)(j + reference * component));
    }
    std::string test = "test";
    function(test, internalBackendData);
  }

  const auto foundIterator = dataEntry.find("test");
  EXPECT_NE(foundIterator, std::end(dataEntry));
  const auto &cellEntry = foundIterator->second;
  EXPECT_EQ(cellEntry.componentSize, 3);
  EXPECT_EQ(cellEntry.values.size(), points.size() * cellEntry.componentSize);
  for (size_t i = 0; i < cellEntry.values.size(); i++) {
    EXPECT_DOUBLE_EQ(cellEntry.values[i], (double)i);
  }
}

TEST(TestVtuBackend, DefaultBackendTesting) {
  DefaultPlotBackend plotBackend;

  std::vector<Point> points = makePointList({0, 0.22222, 0.444});
  std::vector<Point> pointsCorner = makePointList({0, 0.1, 0.3, 0.5});

  CellInfo cellInfo;
  cellInfo.cornerCount = {2, 1, 1};
  cellInfo.cells = points;
  cellInfo.cellCorner = pointsCorner;
  cellInfo.celltype = {CELLTYPE::QUADRILATERAL, CELLTYPE::QUADRILATERAL, CELLTYPE::QUADRILATERAL};
  plotBackend.AddCells(std::move(cellInfo));
  EXPECT_EQ(plotBackend.cellCache.size(), points.size());
  EXPECT_EQ(plotBackend.pointCache.size(), pointsCorner.size());

  testDataFind(
      plotBackend.provider.cells,
      [&](auto &string, auto &backend) { plotBackend.AddCellData(string, backend); },
      plotBackend.entrys.cellData);

  testDataFind(
      plotBackend.provider.points,
      [&](auto &string, auto &backend) { plotBackend.AddPointData(string, backend); },
      plotBackend.entrys.pointData);

  const VtuInfo plotData = plotBackend.GetVtuData()[0];
  EXPECT_EQ(plotData.dataSetInfo.size(), 2);
  EXPECT_EQ(plotData.dataSetData.size(), 2);
  EXPECT_EQ(plotData.mesh.connectivity_.size(), 4);
  EXPECT_EQ(plotData.mesh.offsets_.size(), 3);
  if constexpr (kEnableLowerPointDim) {
    EXPECT_EQ(plotData.mesh.pointDim_, SpaceDimension);
  } else {
    EXPECT_EQ(plotData.mesh.pointDim_, 3);
  }
  EXPECT_EQ(plotData.mesh.points_.size(), 4 * plotData.mesh.pointDim_);
  EXPECT_EQ(plotData.mesh.types_.size(), 3);

  for (size_t i = 0; i < 3; i++) {
    EXPECT_EQ(plotData.mesh.types_[i], VTK_QUAD);
  }
  EXPECT_EQ(plotData.mesh.offsets_[0], 2);
  for (size_t i = 0; i < 4; i++) {
    EXPECT_EQ(plotData.mesh.connectivity_[i], i);
  }
  for (size_t i = 0; i < pointsCorner.size(); i++) {
    auto &currentpoint = pointsCorner[i];
    for (size_t j = 0; j < SpaceDimension; j++) {
      const double value = plotData.mesh.points_[i * plotData.mesh.pointDim_ + j];
      EXPECT_DOUBLE_EQ(currentpoint[j], value);
    }
  }
}

#if USE_SPACETIME && SpaceDimension < 3
void testMemoryTransitivity(const PlotMeshProvider &plotMesh, const size_t count) {
  EXPECT_EQ(plotMesh.cells.size(), count);
  EXPECT_EQ(plotMesh.cells.size(), plotMesh.cellTypes.size());
  EXPECT_EQ(plotMesh.cellTypes.size(), plotMesh.offset.size());
}

void testSpaceTimeMemoryTransitivity(SpaceTimeBackend &backend) {
  EXPECT_EQ(backend.keys.size(), backend.entrys.size());
  EXPECT_EQ(backend.entrys.size(), backend.meshes.size());
  EXPECT_EQ(backend.meshes.size(), backend.pointCaches.size());
  EXPECT_EQ(backend.pointCaches.size(), backend.cellCaches.size());
}

TEST(TestSpaceTimeBackend, FunctionTest) {
  SpaceTimeBackend backendSpaceTime;

  CellInfo cellInfo;
  cellInfo.cornerCount = {2, 2, 2};
  cellInfo.cells = {{3.5, 3, 3, 3}, {1.5, 1, 1, 1}, {2.5, 2, 2, 2}};
  cellInfo.celltype = {CELLTYPE::INTERVAL, CELLTYPE::INTERVAL, CELLTYPE::INTERVAL};
  cellInfo.cellCorner = {{3, 3, 3, 3}, {4, 3, 3, 3}, {1, 1, 1, 1},
                         {2, 1, 1, 1}, {2, 2, 2, 2}, {3, 2, 2, 2}};
  cellInfo.timeSteps = {1, 2, 3};
  CellInfo cellCopy = cellInfo;
  backendSpaceTime.AddCells(std::move(cellCopy));

  EXPECT_EQ(backendSpaceTime.keys.size(), 3);
  testSpaceTimeMemoryTransitivity(backendSpaceTime);

  for (const auto &point : cellInfo.cells) {
    const size_t value = (size_t)point.t() - 1;
    const auto index = backendSpaceTime.Find(point.t());
    EXPECT_EQ(index, value);
    const auto &cache = backendSpaceTime.cellCaches.at(index);
    EXPECT_TRUE(cache.contains(point));
    auto &mesh = backendSpaceTime.meshes[index];
    for (const auto cellType : mesh.cellTypes)
      EXPECT_EQ(cellType, CELLTYPE::INTERVAL);
  }

  const auto &oneMesh = backendSpaceTime.meshes[0];
  testMemoryTransitivity(oneMesh, 1);
  EXPECT_EQ(oneMesh.cells[0], cellInfo.cells[1]);
  EXPECT_EQ(oneMesh.offset[0], 1);
  EXPECT_EQ(oneMesh.connectivity[0], 0);

  const auto &twoMesh = backendSpaceTime.meshes[1];
  testMemoryTransitivity(twoMesh, 1);
  EXPECT_EQ(twoMesh.cells[0], cellInfo.cells[2]);
  EXPECT_EQ(twoMesh.offset[0], 1);
  EXPECT_EQ(twoMesh.connectivity[0], 0);

  const auto &threeMesh = backendSpaceTime.meshes[2];
  testMemoryTransitivity(threeMesh, 1);
  EXPECT_EQ(threeMesh.cells[0], cellInfo.cells[0]);
  EXPECT_EQ(threeMesh.offset[0], 1);
  EXPECT_EQ(threeMesh.connectivity[0], 0);
}

TEST(TestSpaceTime, SpaceTimeSquareTriangles) {
  MeshesCreator meshCreator("SpaceTimeSquareTriangles");
  meshCreator.WithLevel(4);
  meshCreator.WithPLevel(4);
  const auto meshes = meshCreator.CreateUnique();
  meshes->PrintInfo();
  auto discretization = std::make_shared<STDiscretization_DGDG>(*meshes, DegreePair{1, 1}, 2);

  Vector usedVector(discretization);
  fillSpaceTimeVector(usedVector, 1, [](const Point &x, const cell &c, int component) {
    if (component == 0) { return x[0]; }
    if (component == 1) { return x.t(); }
    return 0.0;
  });

  PlotConfig config;
  config.plotBackendCreator = std::make_unique<SpaceTimeBackend>;
  VtuPlot vtuPlot("SpaceTimeSquareTriangles", config);
  vtuPlot.AddData("test", usedVector);
  vtuPlot.PlotFile();

  config.parallelPlotting = true;
  VtuPlot vtuPlotParallel("SpaceTimeSquareTrianglesParallel", config);
  vtuPlotParallel.AddData("test", usedVector);
  vtuPlotParallel.PlotFile();

  PlotConfig config2;
  config2.plotBackendCreator = std::make_unique<DefaultSpaceTimeBackend>;
  VtuPlot vtuPlot2("SpaceTimeSquareTrianglesDefault", config2);
  vtuPlot2.AddData("test", usedVector);
  vtuPlot2.PlotFile();

  config2.parallelPlotting = true;
  VtuPlot vtuPlotParallel2("SpaceTimeSquareTrianglesDefaultParallel", config2);
  vtuPlotParallel2.AddData("test", usedVector);
  vtuPlotParallel2.PlotFile();
}

TEST(TestSpaceTime, SpaceTimeSquareQuads) {
  MeshesCreator meshCreator("SpaceTimeSquare");
  meshCreator.WithLevel(4);
  meshCreator.WithPLevel(4);
  const auto meshes = meshCreator.CreateUnique();
  meshes->PrintInfo();
  auto discretization = std::make_shared<STDiscretization_DGDG>(*meshes, DegreePair{1, 1}, 2);

  Vector usedVector(discretization);
  fillSpaceTimeVector(usedVector, 3, [](const Point &x, const cell &c, int component) {
    if (component == 0) { return 1; }
    if (component == 1) { return 2; }
    if (component == 2) { return 3; }
    return 4;
  });

  PlotConfig config;
  config.plotBackendCreator = std::make_unique<SpaceTimeBackend>;
  VtuPlot vtuPlot("SpaceTimeSquare", config);
  vtuPlot.AddData("test", usedVector, 0);
  vtuPlot.PlotFile();

  config.parallelPlotting = true;
  VtuPlot vtuPlotParallel("SpaceTimeSquareParallel", config);
  vtuPlotParallel.AddData("test", usedVector, 0);
  vtuPlotParallel.PlotFile();

  PlotConfig config2;
  config2.plotBackendCreator = std::make_unique<DefaultSpaceTimeBackend>;
  VtuPlot vtuPlot2("SpaceTimeSquareDefault", config2);
  vtuPlot2.AddData("test", usedVector, 0);
  vtuPlot2.PlotFile();

  config2.parallelPlotting = true;
  VtuPlot vtuPlotParallel2("SpaceTimeSquareDefaultParallel", config2);
  vtuPlotParallel2.AddData("test", usedVector, 0);
  vtuPlotParallel2.PlotFile();
}

#endif

TEST(TestPlotting, multiplelPlotting) {
  VtuPlot &plot1 = mpp::plot("test1");
  ASSERT_EQ(plot1.FileName(), "test1");
  VtuPlot &plot2 = mpp::plot("test2");
  ASSERT_EQ(plot2.FileName(), "test2");
  ASSERT_NE(plot1.FileName(), plot2.FileName());

  auto creator = MeshesCreator("Triangle").WithPLevel(4).WithLevel(4);
  auto meshes = creator.CreateUnique();
  auto discretization = std::make_shared<LagrangeDiscretization>(*meshes, 1, 1);
  Vector vectorData(1.0, discretization);
  mpp::PlotData("testPlotData", "testPlotData", vectorData);
  std::string plotPath = Config::GetPlotPath();
  if (!PPM->Master(0)) return;
  ASSERT_TRUE(fs::exists(plotPath + "testPlotData.vtu"));
}

TEST(TestPlotVtuFile, plotFile) {
  std::string plotPath = Config::GetPlotPath();
  if (PPM->Master(0)) {
    if (!fs::exists(plotPath)) fs::create_directories(plotPath);

    for (size_t i = 0; i <= 3; i++) {
      if (fs::exists(plotPath + "test" + std::to_string(i))) {
        fs::remove(plotPath + "test" + std::to_string(i));
      }
    }
  }

  VtuPlot vtuPlot1("test1");
  EXPECT_EQ(vtuPlot1.FileName(), "test1");
  vtuPlot1.PlotFile();
  PPM->Barrier();

  EXPECT_TRUE(std::filesystem::exists(plotPath + "test1.vtu"));

  VtuPlot vtuPlot2("test2");
  EXPECT_EQ(vtuPlot2.FileName(), "test2");
  vtuPlot2.PlotFile("test3");
  PPM->Barrier();

  EXPECT_TRUE(std::filesystem::exists(plotPath + "test3.vtu"));

  if (!PPM->Master(0)) return;
  pugi::xml_document doc;
  std::string filePath = plotPath + "test1.vtu";
  pugi::xml_parse_result result = doc.load_file(filePath.c_str());
  EXPECT_TRUE(result);
  pugi::xml_node parent = doc.child("VTKFile").child("UnstructuredGrid").child("Piece");
  auto pointsDataArrayNode = parent.child("Points").child("DataArray");
  auto cellsDataArrayNode1 = parent.child("Cells").child("DataArray");
  auto cellsDataArrayNode2 = cellsDataArrayNode1.next_sibling();
  auto cellsDataArrayNode3 = cellsDataArrayNode2.next_sibling();

  // Checking file1
  checkXmlAttributes({{"byte_order", "LittleEndian"},
                      {"type", "UnstructuredGrid"},
                      {"version", "0.1"}},
                     doc.child("VTKFile"));

  checkXmlAttributes({{"NumberOfCells", "0"}, {"NumberOfPoints", "0"}}, parent);

  if constexpr (kEnableLowerPointDim) {
    checkXmlAttributes({{"NumberOfComponents", std::to_string(SpaceDimension)},
                        {"format", "ascii"},
                        {"type", "Float64"}},
                       pointsDataArrayNode);
  } else {
    checkXmlAttributes({{"NumberOfComponents", "3"}, {"format", "ascii"}, {"type", "Float64"}},
                       pointsDataArrayNode);
  }

  checkXmlAttributes({{"Name", "connectivity"}, {"format", "ascii"}, {"type", "Int64"}},
                     cellsDataArrayNode1);
  checkXmlAttributes({{"Name", "offsets"}, {"format", "ascii"}, {"type", "Int64"}},
                     cellsDataArrayNode2);
  checkXmlAttributes({{"Name", "types"}, {"format", "ascii"}, {"type", "Int8"}},
                     cellsDataArrayNode3);
}

class TestVtuPlotWithDisplacement : public TestWithParam<std::string> {
protected:
  const std::unique_ptr<Meshes> meshes;
  const Mesh &mesh;
  const std::shared_ptr<const IDiscretization> disc;

  Vector displacement;
  std::unordered_map<Point, std::vector<double>> displacementData;
public:
  explicit TestVtuPlotWithDisplacement(
      const std::function<std::shared_ptr<IDiscretization>(Meshes &)> &discCreator,
      MeshesCreator &&meshCreator = MeshesCreator()) :
      meshes(meshCreator.WithMeshName(GetParam()).WithPLevel(4).WithLevel(4).CreateUnique()),
      mesh(meshes->fine()), disc(discCreator(*meshes)), displacement(0.0, disc) {}

  void Plot(PlotConfig &config, const std::string &plotPath) {
    if (PPM->Master()) {
      if (fs::exists(plotPath)) fs::remove_all(plotPath);
    }
    PPM->Barrier();

    VtuPlot plot = VtuPlot(mesh, GetParam() + "_deformed", config);
    plot.DeformGrid(displacement);
    plot.ChangeDataPath(plotPath);
    plot.PlotFile(GetParam() + "_deformed");

    // Sync fs cache as the result will be read on multiple procs

    PPM->Barrier();
  }

  void TestDisplacement(const std::string &deformedPath) {
    pugi::xml_document document;
    const auto result = document.load_file(deformedPath.c_str());
    EXPECT_TRUE(result);

    const pugi::xml_node pieceNode =
        document.child("VTKFile").child("UnstructuredGrid").child("Piece");
    const pugi::xml_node pointsNode = pieceNode.child("Points").child("DataArray");

    const int pointCount = pieceNode.attribute("NumberOfPoints").as_int();
    const int pointDim = pointsNode.attribute("NumberOfComponents").as_int();

    const std::vector<Point> deformedPoints =
        DataArrayToPoints(ParseStrDataArray(pointsNode.child_value()), pointDim);
    EXPECT_EQ(pointCount, deformedPoints.size());

    const std::unordered_set<Point> deformedPointsSet(std::begin(deformedPoints),
                                                      std::end(deformedPoints));

    for (const auto &[point, displacement] : displacementData) {
      Point displacedPoint(point);
      for (size_t i = 0; i < displacement.size(); ++i) {
        displacedPoint[i] += displacement[i];
      }
      bool it = deformedPointsSet.contains(displacedPoint);
      EXPECT_TRUE(deformedPointsSet.contains(displacedPoint))
          << "Plot does not contain displaced point " << displacedPoint << " which originally was "
          << point;
    }
  }
};

class TestVtuPlotWithDisplacementLagrange : public TestVtuPlotWithDisplacement {
private:
  inline void fillDisplacement() {
    static constexpr int minCornerCount = 2;
    const int proc = PPM->Proc();
    const int dim = displacement.dim();
    const procset procSetEnd = displacement.procsets_end();
    displacementData.reserve(displacement.nR());

    for (row row = displacement.rows(); row != displacement.rows_end(); ++row) {
      const Point point = row();
      if (displacement.find_vertex(point) != displacement.vertices_end()) {
        std::vector<double> displacementDataEntry(dim);
        const bool existsOnOneProcOnly = displacement.find_procset(point) == procSetEnd;

        for (int j = 0; j < dim; ++j) {
          const double jEntry = point[j];
          displacement(row, j) = jEntry;
          displacementDataEntry[j] = jEntry;
        }
        displacement(row, 0) *= 2;
        displacementDataEntry[0] *= 2;

        displacementData[point] = displacementDataEntry;
      }
    }
  }
public:
  explicit TestVtuPlotWithDisplacementLagrange(int degree = 1, int size = SpaceDimension) :
      TestVtuPlotWithDisplacement([=](Meshes &meshes) {
        return std::make_unique<LagrangeDiscretization>(meshes, degree, size);
      }) {
    fillDisplacement();
  }
};

class TestVtuPlotWithDisplacementDG : public TestVtuPlotWithDisplacement {
private:
  inline void fillDisplacement() {
    static constexpr int minCornerCount = 2;
    const int proc = PPM->Proc();
    const int displacementOfAxis = proc + 1;
    const int dim = displacement.dim();
    const procset procSetEnd = displacement.procsets_end();
    displacementData.reserve(displacement.nR() * minCornerCount);

    const std::vector<double> displacementDataEntry(dim, displacementOfAxis);

    for (row row = displacement.rows(); row != displacement.rows_end(); ++row) {
      const cell cell = displacement.find_cell(row());
      if (cell != displacement.cells_end()) {
        const bool existsOnOneProcOnly = displacement.find_procset(cell()) == procSetEnd;

        for (int cornerIdx = 0; cornerIdx < cell.Corners(); ++cornerIdx) {
          const int offset = cornerIdx * dim;
          for (int j = 0; j < dim; ++j) {
            displacement(row, offset + j) = displacementOfAxis;
          }

          if (existsOnOneProcOnly) {
            // Only test displaced points which do not exist on multiple
            // processes
            displacementData[cell.Corner(cornerIdx)] = displacementDataEntry;
          }
        }
      }
    }
  }
public:
  explicit TestVtuPlotWithDisplacementDG(int degree = 1, int size = SpaceDimension) :
      TestVtuPlotWithDisplacement(
          [=](Meshes &meshes) { return std::make_unique<DGDiscretization>(meshes, degree, size); },
          MeshesCreator()) {
    fillDisplacement();
  }
};

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithDisplacementLagrange,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Triangle", "Square", "Square2Triangles"
#endif
#if SpaceDimension >= 3
                                ,
                                "Tetrahedron", "Hexahedron"
#endif
                                ));

TEST_P(TestVtuPlotWithDisplacementLagrange, TestDisplacementSerial) {
  const std::string plotPath = Config::GetPlotPath() + "displacement_lagrange/";
  const std::string plotName = GetParam() + "_deformed";
  const std::string serialDeformedPath = plotPath + plotName + ".vtu";

  PlotConfig config{
      .compression = Compression::ASCII,
      .parallelPlotting = false,
      .plotBackendCreator = std::make_unique<DisplacementBackend>,
  };

  Plot(config, plotPath);
  TestDisplacement(serialDeformedPath);
}

TEST_P(TestVtuPlotWithDisplacementLagrange, TestDisplacementParallel) {
  const std::string plotPath = Config::GetPlotPath() + "displacement_lagrange/";
  const std::string plotName = GetParam() + "_deformed";
  const std::string parallelDeformedPath =
      plotPath + plotName + "/" + plotName + "_" + std::to_string(PPM->Proc()) + ".vtu";

  PlotConfig config{
      .compression = Compression::ASCII,
      .parallelPlotting = true,
      .plotBackendCreator = std::make_unique<DisplacementBackend>,
  };

  Plot(config, plotPath);
  TestDisplacement(parallelDeformedPath);
}

INSTANTIATE_TEST_SUITE_P(TestVtuPlot, TestVtuPlotWithDisplacementDG,
                         Values("Interval"
#if SpaceDimension >= 2
                                ,
                                "Square"
#endif
#if SpaceDimension >= 3
                                ,
                                "Hexahedron"
#endif
                                ));

TEST_P(TestVtuPlotWithDisplacementDG, TestDisplacementSerial) {
  const std::string plotPath = Config::GetPlotPath() + "displacement_dg/";
  const std::string plotName = GetParam() + "_deformed";
  const std::string serialDeformedPath = plotPath + plotName + ".vtu";

  PlotConfig config{
      .compression = Compression::ASCII,
      .parallelPlotting = false,
      .plotBackendCreator = std::make_unique<DisplacementBackend>,
  };

  Plot(config, plotPath);
  TestDisplacement(serialDeformedPath);
}

TEST_P(TestVtuPlotWithDisplacementDG, TestDisplacementParallel) {
  const std::string plotPath = Config::GetPlotPath() + "displacement_dg/";
  const std::string plotName = GetParam() + "_deformed";
  const std::string parallelDeformedPath =
      plotPath + plotName + "/" + plotName + "_" + std::to_string(PPM->Proc()) + ".vtu";

  PlotConfig config{
      .compression = Compression::ASCII,
      .parallelPlotting = true,
      .plotBackendCreator = std::make_unique<DisplacementBackend>,
  };

  Plot(config, plotPath);
  TestDisplacement(parallelDeformedPath);
}

int main(int argc, char **argv) {
  MppTest mppTest = MppTest(MppTestBuilder(argc, argv)
                                .WithScreenLogging()
                                .WithParallelListeners()
                                .WithConfigEntry("Overlap", "STCellsWithFaces")
                                .WithConfigEntry("MeshVerbose", 2)
                                .WithConfigEntry("MeshesVerbose", 1)
                                .WithConfigEntry("Compression", "ascii")
                                .WithConfigEntry("ParallelPlotting", false)
                                .WithPPM());

  return mppTest.RUN_ALL_MPP_TESTS();
}