#ifndef COARSEGEOMETRY_HPP
#define COARSEGEOMETRY_HPP

#include <memory>
#include "Celltype.hpp"
#include "CheckOrientation.hpp"
#include "Config.hpp"
#include "DataSet.hpp"
#include "Functools.hpp"
#include "Parallel.hpp"

// TODO: Move to Config!
static std::string geometrySavePath() {
  if (ReplaceGeometries) { return Config::GetGeoPath(); }
  return std::string(ProjectBuildDir) + "/data/geo/";
}

class Buffer;

class DataContainer;

class M_ifstream;

using Coordinates = std::vector<Point>;

class CellId {
  CELLTYPE type;
  int subdomain;
  std::vector<int> cornerIndices{};
public:
  CellId(){};

  CellId(CELLTYPE type, int subdomain, std::vector<int> &&cornerIndices) :
      type(type), subdomain(subdomain), cornerIndices(std::move(cornerIndices)) {}

  CELLTYPE Type() const { return type; }

  int Subdomain() const { return subdomain; }

  const std::vector<int> &CornerIndices() const { return cornerIndices; }

  bool operator==(const CellId &other) {
    std::vector<int> diff;
    std::set_difference(cornerIndices.begin(), cornerIndices.end(), other.cornerIndices.begin(),
                        other.cornerIndices.end(), std::inserter(diff, diff.begin()));

    return (type == other.type) && (subdomain == other.subdomain) && diff.empty();
  }

  friend Buffer &operator<<(Buffer &b, const CellId &id);

  friend Buffer &operator>>(Buffer &b, CellId &id);

  friend std::ostream &operator<<(std::ostream &s, const CellId &i);

  friend class CoarseGeometry;
};

using CellIds = std::vector<CellId>;

class FaceId {
  int bndType{};
  std::vector<int> faceIndices{};
public:
  FaceId() = default;

  FaceId(int bndType, std::vector<int> &&faceIndices) :
      bndType(bndType), faceIndices(std::move(faceIndices)) {}

  int BndType() const { return bndType; }

  const std::vector<int> FaceIndices() const { return faceIndices; }

  bool operator==(const FaceId &other) {
    std::vector<int> diff;
    std::set_difference(faceIndices.begin(), faceIndices.end(), other.faceIndices.begin(),
                        other.faceIndices.end(), std::inserter(diff, diff.begin()));

    return (bndType == other.bndType) && diff.empty();
  }

  friend Buffer &operator<<(Buffer &b, const FaceId &id);

  friend Buffer &operator>>(Buffer &b, FaceId &ids);

  friend std::ostream &operator<<(std::ostream &s, const FaceId &i);
};

using FaceIds = std::vector<FaceId>;

using VertexDataList = std::vector<DataContainer>;

using CellDataList = std::vector<DataContainer>;

using TimeSteps = std::vector<double>;

using CoarseGeometryData = std::tuple<Coordinates, CellIds, FaceIds>;

class CoarseGeometry {
  // CoarseGeometry <- geoName, subdomains, timeSteps, scaling, transformation

  friend class MeshesCreator;
protected:
  std::string geoName{};

  std::unordered_map<std::string, std::string> geoInformation{};

  std::string getInformation(const std::string &key) {
    if (geoInformation.find(key) == geoInformation.end()) return "";
    return geoInformation[key];
  }

  /// Maximum length of a line in the .geo file
  const int len = 256;

  double scalingFactor = 1.0;

  int refineTimeStepCount = 0;

  int dimension = 0;

  Coordinates coordinates{};

  CellIds cellIds{};

  FaceIds faceIds{};

  VertexDataList vdataList{};

  CellDataList cdataList{};
public:
  TimeSteps timeSteps{};
protected:
  M_ifstream loadGeometryFile(const std::string &name);

  /// Functions to read old geo file format
  bool readGeometryFile(std::istream &geoFileStream);

  static int calculateSpaceDim(std::vector<Point> &coordinates);

  static CELLTYPE getCellType(int calculatedSpaceDim, int tp, int numberOfCorners);

  void readHeader(std::istream &geoFileStream, char *L);

  bool insertInformation(string &&line);

  void readVerteces(std::istream &geoFileStream, char *L);

  bool insertCoordinate(string &&line);

  bool insertVtuCellId(string &&line, string celltype, string globalSubdomain);

  void readCells(std::istream &geoFileStream, char *L);

  bool insertLegacyCellId(string &&line, int calculatedSpaceDim);

  bool insertVtuFaceId(string &&line, string globalBoundary);

  void readFaces(std::istream &geoFileStream, char *L);

  bool insertLegacyFaceId(string &&line);

  void readData(std::istream &geoFileStream, char *L);

  bool insertVtuData(std::string &&line, std::vector<DataContainer> &dataList);

  bool insertLegacyData(std::string &&line, std::vector<DataContainer> &dataList);

  void readTimes(std::istream &geoFileStream, char *L);

  void insertTimeSteps(std::string &&line);

  TimeSteps refineTimeSteps(const TimeSteps &timeSteps);

  void check();
public:
  explicit CoarseGeometry(const std::string &geometryFilename);

  explicit CoarseGeometry(CoarseGeometryData &&data, std::string name = "Initialized by data") :
      CoarseGeometry(std::move(std::get<0>(data)), std::move(std::get<1>(data)),
                     std::move(std::get<2>(data)), {}, name) {}

  explicit CoarseGeometry(Coordinates &&coordinates, CellIds &&cellIds, FaceIds &&faceIds,
                          std::string name = "Initialized by data") :
      CoarseGeometry(std::move(coordinates), std::move(cellIds), std::move(faceIds), {}, name) {}

  explicit CoarseGeometry(Coordinates &&coordinates, CellIds &&cellIds, FaceIds &&faceIds,
                          TimeSteps &&timeSteps, std::string name = "Initialized by data") :
      CoarseGeometry(std::move(coordinates), std::move(cellIds), std::move(faceIds),
                     std::move(timeSteps), {}, {}, name) {}

  explicit CoarseGeometry(Coordinates &&coordinates, CellIds &&cellIds, FaceIds &&faceIds,
                          VertexDataList &&vertexData, CellDataList &&cellData,
                          std::string name = "Initialized by data") :
      CoarseGeometry(std::move(coordinates), std::move(cellIds), std::move(faceIds), {},
                     std::move(vertexData), std::move(cellData), name) {}

  explicit CoarseGeometry(Coordinates &&coordinates, CellIds &&cellIds, FaceIds &&faceIds,
                          TimeSteps &&timeSteps, VertexDataList &&vertexData,
                          CellDataList &&cellData, std::string name = "Initialized by data");

  virtual bool HasMapGeometry() const { return false; }

  virtual Point operator()(const Point &x) const { return x; }

  // Todo: maybe change this to linear affine or even nonlinear transformation (this should be done
  // via operator()!)
  void Scale(double scalingFactor);

  virtual ~CoarseGeometry() = default;

  std::string Name() const { return geoName; }

  const Coordinates &GetCoordinates() const { return coordinates; }

  const CellIds &GetCellIds() const { return cellIds; }

  const FaceIds &GetFaceIds() const { return faceIds; }

  const VertexDataList &GetVertexDataList() const { return vdataList; }

  const CellDataList &GetCellDataList() const { return cdataList; }

  std::vector<Point> Corners(const CellId &cellId) const;

  Point FaceCenter(const FaceId &faceId) const;

  int Dim() const { return dimension; }

  void CommunicateGeometry(); // TODO: make protected

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const CoarseGeometry &agg);


  bool SaveGeometryFile(const string &geometryFilename);

  bool withData();
};

class RectangularGeometry2D : public CoarseGeometry {
  Coordinates createPoints(const std::pair<double, double> &xBorders,
                           const std::pair<double, double> &yBorders, int dimX, int dimY);

  CellIds createCells(int dimX, int dimY);

  FaceIds createFaces(int dimX, int dimY, std::vector<int> boundaries);
public:
  RectangularGeometry2D(const std::pair<double, double> &xBorders,
                        const std::pair<double, double> &yBorders, int dimX, int dimY,
                        std::vector<int> boundaries, std::string name = "") :
      RectangularGeometry2D(xBorders, yBorders, dimX, dimY, std::move(boundaries), {}, {},
                            std::move(name)) {}

  RectangularGeometry2D(const std::pair<double, double> &xBorders,
                        const std::pair<double, double> &yBorders, int dimX, int dimY,
                        std::vector<int> boundaries, TimeSteps &&timeSteps, std::string name = "") :
      RectangularGeometry2D(xBorders, yBorders, dimX, dimY, std::move(boundaries),
                            std::move(timeSteps), {}, {}, std::move(name)) {}

  RectangularGeometry2D(const std::pair<double, double> &xBorders,
                        const std::pair<double, double> &yBorders, int dimX, int dimY,
                        std::vector<int> boundaries, VertexDataList &&vertexData,
                        CellDataList &&cellData, std::string name = "") :
      RectangularGeometry2D(xBorders, yBorders, dimX, dimY, std::move(boundaries), {},
                            std::move(vertexData), std::move(cellData), std::move(name)) {}

  RectangularGeometry2D(const std::pair<double, double> &xBorders,
                        const std::pair<double, double> &yBorders, int dimX, int dimY,
                        std::vector<int> boundaries, TimeSteps &&timeSteps,
                        VertexDataList &&vertexData, CellDataList &&cellData,
                        std::string name = "") :
      CoarseGeometry(createPoints(xBorders, yBorders, dimX, dimY), createCells(dimX, dimY),
                     createFaces(dimX, dimY, std::move(boundaries)), std::move(timeSteps),
                     std::move(vertexData), std::move(cellData), std::move(name)) {}
};

CoarseGeometry *CreateCoarseGeometry(const std::string &meshName, VertexDataList &&vertexData = {},
                                     CellDataList &&cellData = {});

std::unique_ptr<CoarseGeometry> CreateCoarseGeometryUnique(const std::string &meshName,
                                                           VertexDataList &&vertexData = {},
                                                           CellDataList &&cellData = {});

std::shared_ptr<CoarseGeometry> CreateCoarseGeometryShared(const std::string &meshName,
                                                           VertexDataList &&vertexData = {},
                                                           CellDataList &&cellData = {});


#endif // of #ifndef COARSEGEOMETRY_HPP
