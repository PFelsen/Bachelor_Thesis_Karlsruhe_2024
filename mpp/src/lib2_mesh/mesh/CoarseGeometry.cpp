#include "CoarseGeometry.hpp"

#include "Assertion.hpp"
#include "Buffer.hpp"
#include "Celltype.hpp"
#include "M_IOFiles.hpp"
#include "Parallel.hpp"

Buffer &operator<<(Buffer &b, const CellId &id) {
  return b << id.type << id.subdomain << id.cornerIndices;
}

Buffer &operator>>(Buffer &b, CellId &id) {
  return b >> id.type >> id.subdomain >> id.cornerIndices;
}

std::ostream &operator<<(std::ostream &s, const CellId &id) {
  s << "type " << id.type << " sd " << id.subdomain << " id";
  for (int cornerIndex : id.cornerIndices)
    s << " " << cornerIndex;
  return s;
}

Buffer &operator<<(Buffer &b, const FaceId &id) { return b << id.bndType << id.faceIndices; }

Buffer &operator>>(Buffer &b, FaceId &id) { return b >> id.bndType >> id.faceIndices; }

std::ostream &operator<<(std::ostream &s, const FaceId &id) {
  s << "bnd " << id.bndType << " id";
  for (int faceIndice : id.faceIndices)
    s << " " << faceIndice;
  return s;
}

M_ifstream CoarseGeometry::loadGeometryFile(const std::string &name) {
  std::string filePath = Config::GetGeoPath() + name + ".geo";
  if (!FileExists(filePath))
    Exit("No mesh with name: " + name + " at " + Config::GetGeoPath()) return {filePath.c_str()};
}

bool CoarseGeometry::readGeometryFile(std::istream &geoFileStream) {
  char L[len];
  geoFileStream.getline(L, len);

  readHeader(geoFileStream, L);
  readVerteces(geoFileStream, L);
  readCells(geoFileStream, L);

  readFaces(geoFileStream, L);
  readData(geoFileStream, L);
  readTimes(geoFileStream, L);
  return true;
}

void CoarseGeometry::readTimes(std::istream &geoFileStream, char *L) {
  if (strncasecmp("times", L, 5) == 0) {
    geoFileStream.getline(L, len);
    insertTimeSteps(L);
  }
}

void CoarseGeometry::readData(std::istream &geoFileStream, char *L) {
  const auto &defaultVData = getInformation("vertexdata");
  if (!defaultVData.empty()) {
    DataContainer data(splitToDouble(defaultVData, " "));
    for (int i = 0; i < coordinates.size(); ++i) {
      vdataList.emplace_back(data);
    }
  } else {
    if (strncasecmp("vdata", L, 5) == 0) geoFileStream.getline(L, len);
    if (geoInformation.find("withdata") != geoInformation.end()) {
      while (insertVtuData(L, vdataList))
        geoFileStream.getline(L, len);
      if (withData() && vdataList.empty()) {
        for (int i = 0; i < coordinates.size(); ++i) {
          vdataList.emplace_back(0);
        }
      }
    } else {
      while (insertLegacyData(L, vdataList))
        geoFileStream.getline(L, len);
    }
  }

  const auto &defaultCData = getInformation("celldata");
  if (!defaultCData.empty()) {
    DataContainer data(splitToDouble(defaultCData, " "));
    for (int i = 0; i < cellIds.size(); ++i) {
      cdataList.emplace_back(data);
    }
  } else {
    if (strncasecmp("cdata", L, 5) == 0) geoFileStream.getline(L, len);
    if (geoInformation.find("withdata") != geoInformation.end()) {
      while (insertVtuData(L, cdataList))
        geoFileStream.getline(L, len);
      if (withData() && cdataList.empty()) {
        for (int i = 0; i < cellIds.size(); ++i) {
          cdataList.emplace_back(0);
        }
      }
    } else {
      while (insertLegacyData(L, cdataList))
        geoFileStream.getline(L, len);
    }
  }
}

void CoarseGeometry::readHeader(std::istream &geoFileStream, char *L) {
  /*
   HEADER
  Cellformat=Vtu
  Celltype=Hexahedron
  Facetype=Quadrilateral
  VertexData = [0 0 0]
  CellData=[1 0 0 0 1 0 0 0 1]
  Subdomain = 1
  Boundary=1
   */
  if (strncasecmp("header", L, 6) == 0) geoFileStream.getline(L, len);
  while (insertInformation(L))
    geoFileStream.getline(L, len);
}

void CoarseGeometry::readVerteces(std::istream &geoFileStream, char *L) {
  if (strncasecmp("points", L, 6) == 0) geoFileStream.getline(L, len);
  while (insertCoordinate(std::string(L)))
    geoFileStream.getline(L, len);
}

void CoarseGeometry::readCells(std::istream &geoFileStream, char *L) {
  if (strncasecmp("cells", L, 5) == 0) geoFileStream.getline(L, len);
  if (geoInformation.find("cellformat") != geoInformation.end()) {
    std::string format = lowerCase(geoInformation["cellformat"]);
    std::string celltype = getInformation("celltype");
    std::string globalSubdomain = getInformation("subdomain");

    if (format == "vtu") {
      while (insertVtuCellId(std::string(L), celltype, globalSubdomain)) {
        geoFileStream.getline(L, len);
      }
    }
  } else {
    int dim = calculateSpaceDim(coordinates);
    while (insertLegacyCellId(std::string(L), dim)) {
      geoFileStream.getline(L, len);
    }
  }
}

void CoarseGeometry::readFaces(std::istream &geoFileStream, char *L) {
  if (strncasecmp("faces", L, 5) == 0) geoFileStream.getline(L, len);
  if (geoInformation.find("cellformat") != geoInformation.end()) {
    std::string format = lowerCase(geoInformation["cellformat"]);
    std::string globalboundary = getInformation("boundary");

    if (format == "vtu") {
      while (insertVtuFaceId(std::string(L), globalboundary)) {
        geoFileStream.getline(L, len);
      }
    }
  } else {
    while (insertLegacyFaceId(L))
      geoFileStream.getline(L, len);
  }
}

int CoarseGeometry::calculateSpaceDim(std::vector<Point> &coordinates) {
  int calculatedSpaceDim = 1;
  if constexpr (SpaceDimension > 1) {
    for (auto &coordinate : coordinates)
      if (coordinate[1] != 0) calculatedSpaceDim = 2;
  }
  if constexpr (SpaceDimension > 2) {
    for (auto &coordinate : coordinates)
      if (coordinate[2] != 0) calculatedSpaceDim = 3;
  }
  return calculatedSpaceDim;
}

CELLTYPE CoarseGeometry::getCellType(int calculatedSpaceDim, int tp, int numberOfCorners) {
  if (calculatedSpaceDim == 1) return INTERVAL;
  else if (calculatedSpaceDim == 2) {
    if (tp == 3 && numberOfCorners == 3) return TRIANGLE;
    if (tp == 4 && numberOfCorners == 4) return QUADRILATERAL;
    if (tp == 4 && numberOfCorners == 8) return QUADRILATERAL2;
  } else {
    if (tp == 4 && numberOfCorners == 4) return TETRAHEDRON;
    if (tp == 8 && numberOfCorners == 8) return HEXAHEDRON;
    if (tp == 8 && numberOfCorners == 20) return HEXAHEDRON20;
    if (tp == 8 && numberOfCorners == 27) return HEXAHEDRON27;
  }
  THROW("No corresponding CellType")
}

bool CoarseGeometry::insertCoordinate(string &&line) {
  try {
    auto z = splitToDouble(line, " ");
    if (z.empty()) return false;
    coordinates.emplace_back(z);
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::insertVtuCellId(string &&line, string celltype, string globalSubdomain) {
  try {
    auto lineContent = splitToInt(line, " ");
    if (lineContent.empty()) return false;

    int cType = lineContent[0];
    if (!celltype.empty()) {
      try {
        cType = std::stoi(celltype);
      } catch (std::exception &e) { cType = CellnameToVtu(celltype); }
    }
    int subdomain =
        globalSubdomain.empty() ? lineContent[(int)celltype.empty()] : std::stoi(globalSubdomain);

    int nodes = lineContent.size() - (int)celltype.empty() - (int)globalSubdomain.empty();

    std::vector<int> indices(nodes);
    std::uninitialized_move(lineContent.begin() + lineContent.size() - nodes, lineContent.end(),
                            indices.begin());

    cellIds.emplace_back(vtuToCelltype(cType), subdomain, std::move(indices));
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::insertLegacyCellId(string &&line, int calculatedSpaceDim) {
  try {
    auto lineContent = splitToInt(line, " ");
    if (lineContent.empty()) return false;

    int nodes = lineContent[0];
    int subdomain = lineContent[1];

    if (nodes > lineContent.size() - 2) { return false; }

    std::vector<int> indices(nodes);
    std::uninitialized_move(lineContent.begin() + 2, lineContent.end(), indices.begin());

    cellIds.emplace_back(getCellType(calculatedSpaceDim, nodes, indices.size()), subdomain,
                         std::move(indices));
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::insertVtuFaceId(string &&line, string globalBoundary) {
  try {
    auto lineContent = splitToInt(line, " ");
    if (lineContent.empty()) return false;

    int boundary = globalBoundary.empty() ? lineContent[0] : std::stoi(globalBoundary);

    int nodes = lineContent.size() - (int)globalBoundary.empty();

    std::vector<int> indices(nodes);
    std::uninitialized_move(lineContent.begin() + lineContent.size() - nodes, lineContent.end(),
                            indices.begin());

    faceIds.emplace_back(FaceId(boundary, std::move(indices)));
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::insertLegacyFaceId(string &&line) {
  try {
    auto lineContent = splitToInt(line, " ");
    if (lineContent.empty()) return false;

    int nodes = lineContent[0];
    int boundary = lineContent[1];

    if (nodes > lineContent.size() - 2) { return false; }

    std::vector<int> indices(nodes);
    std::uninitialized_move(lineContent.begin() + 2, lineContent.end(), indices.begin());

    faceIds.emplace_back(FaceId(boundary, std::move(indices)));
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::insertVtuData(std::string &&line, std::vector<DataContainer> &dataList) {
  try {
    auto data = splitToDouble(line, " ");
    if (data.empty()) return false;

    dataList.emplace_back(DataContainer(std::move(data)));
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::insertLegacyData(std::string &&line, std::vector<DataContainer> &dataList) {
  try {
    auto lineContent = split(line, " ");

    int size = std::stoi(lineContent[0]);

    if (size == 0) {
      dataList.emplace_back(DataContainer(0));
      return true;
    }
    if (size != lineContent.size() - 1) { return false; }

    std::vector<double> data;
    std::transform(lineContent.cbegin() + 1, lineContent.cend(), std::back_inserter(data),
                   [](const std::string &s) { return std::stod(s); });

    dataList.emplace_back(DataContainer(std::move(data)));
    return true;
  } catch (std::exception &e) { return false; }
}

void CoarseGeometry::insertTimeSteps(std::string &&line) {
  std::vector<string> ts_strings = split(line, " ");
  ts_strings = ft::filter(ts_strings, [](auto &s) { return !s.empty(); });
  timeSteps.resize(ts_strings.size(), 0.0);
  std::transform(ts_strings.begin(), ts_strings.end(), timeSteps.begin(), ft::stod);
}

// TODO is this function necessary?
std::vector<double> CoarseGeometry::refineTimeSteps(const std::vector<double> &timeSteps) {
  //  refining scale times the time steps
  Exit("Not tested!");
  vector<vector<double>> init_timesteps(refineTimeStepCount + 1);
  init_timesteps[0] = timeSteps;
  int SF_inv = 2.0;
  for (unsigned int S = 1; S < refineTimeStepCount + 1; ++S) {
    unsigned int m = init_timesteps[S - 1].size();
    init_timesteps[S].resize(2 * m - 1);
    init_timesteps[S][0] = init_timesteps[S - 1][0];
    init_timesteps[S][m + (m - 1) * (SF_inv - 1) - 1] = init_timesteps[S - 1][m - 1];
    for (unsigned int k = 0; k < m - 1; ++k) {
      init_timesteps[S][k * SF_inv + 1] = init_timesteps[S - 1][k + 1];
      double dt = init_timesteps[S - 1][k + 2] - init_timesteps[S - 1][k + 1];
      for (unsigned int l = 1; l < SF_inv; ++l)
        init_timesteps[S][k * SF_inv + l + 1] =
            init_timesteps[S - 1][k + 1] + dt * 1 / double(SF_inv) * l;
    }
  }
  return init_timesteps[refineTimeStepCount];
}

void CoarseGeometry::check() {
  // check if all indices of cells are contained in coordinates
  for (const auto &cell_id : cellIds) {
    const std::vector<int> &indices = cell_id.CornerIndices();
    for (auto idx : indices) {
      if (idx >= coordinates.size()) THROW("At least one cell index not found in coordinates")
    }
  }

  // check if all indices of faces are contained in coordinates
  for (const auto &face_id : faceIds) {
    const std::vector<int> &indices = face_id.FaceIndices();
    for (auto idx : indices) {
      if (idx >= coordinates.size()) THROW("At least one face index not found in coordinates")
    }
  }

#ifdef USE_DATAMESH
  // check data mesh sizes
  if (coordinates.size() != vdataList.size())
    THROW("Length of VData does not match length of coordinates")
  if (cellIds.size() != cdataList.size()) THROW("Length of CData does not match length of cells")
#endif

  // check orientation
  for (auto &cell_id : cellIds) {
    std::vector<Point> corners = Corners(cell_id);
    CheckOrientation::Check(cell_id.Type(), corners, cell_id.cornerIndices);
  }
}

CoarseGeometry::CoarseGeometry(const string &geometryFilename) {
  if (PPM->Master(0)) {
    M_ifstream geoFileStream = loadGeometryFile(geometryFilename);
    readGeometryFile(geoFileStream);
    if (scalingFactor != 1.0) {
      mout << "  scaling geometry by " << scalingFactor << endl;
      Scale(scalingFactor);
    }
    geoFileStream.close();

    if (geoInformation.empty() || getInformation("checked").empty()
        || getInformation("checked") == "false") {
      check();
      SaveGeometryFile(geometryFilename);
    }
  }
  CommunicateGeometry();
  geoName = geometryFilename;

  if (cellIds.empty()) { THROW("cellIds are empty, can't init dimension"); }
  dimension = CellDim(cellIds[0].Type());
}

CoarseGeometry::CoarseGeometry(Coordinates &&coordinates, CellIds &&cellIds, FaceIds &&faceIds,
                               TimeSteps &&timeSteps, VertexDataList &&vertexData,
                               CellDataList &&cellData, std::string name) :
    coordinates(std::move(coordinates)), cellIds(std::move(cellIds)), faceIds(std::move(faceIds)),
    timeSteps(std::move(timeSteps)) {
  if (vertexData.size() == 1) vdataList = VertexDataList(this->coordinates.size(), vertexData[0]);
  else vdataList = std::move(vertexData);
  if (cellData.size() == 1) cdataList = CellDataList(this->cellIds.size(), cellData[0]);
  else cdataList = std::move(cellData);
  geoName = name;
  check();

  if (this->cellIds.empty()) { THROW("cellIds are empty, can't init dimension"); }
  dimension = CellDim(this->cellIds[0].Type());
}

void CoarseGeometry::Scale(double scalingFactor) {
  for (auto &z : coordinates)
    z *= scalingFactor;
}

std::vector<Point> CoarseGeometry::Corners(const CellId &cellId) const {
  const std::vector<int> &ids = cellId.CornerIndices();
  std::vector<Point> corners(ids.size());
  for (int i = 0; i < ids.size(); ++i)
    corners[i] = coordinates[ids[i]];
  return corners;
}

Point CoarseGeometry::FaceCenter(const FaceId &faceId) const {
  const std::vector<int> &ids = faceId.FaceIndices();
  Point faceCenter;
  for (int id : ids)
    faceCenter += coordinates[id];
  faceCenter *= (1.0 / ids.size());
  return faceCenter;
}

void CoarseGeometry::CommunicateGeometry() {
  if (PPM->Size(0) != 1) {
    ExchangeBuffer exBuffer(0);
    if (PPM->Master(0)) {
      for (int q = 1; q < PPM->Size(0); ++q) {
        exBuffer.Send(q) << cellIds;
        exBuffer.Send(q) << faceIds;
        exBuffer.Send(q) << timeSteps;
        exBuffer.Send(q) << coordinates;
      }
    }
    exBuffer.Communicate();
    if (!PPM->Master(0)) {
      exBuffer.Receive(0) >> cellIds;
      exBuffer.Receive(0) >> faceIds;
      exBuffer.Receive(0) >> timeSteps;
      exBuffer.Receive(0) >> coordinates;
    }
    exBuffer.ClearBuffers();
  }
}

template<typename S>
LogTextStream<S> &operator<<(LogTextStream<S> &s, const CoarseGeometry &M) {
  s << "POINTS: " << M.coordinates.size() << endl;
  s << M.coordinates << "CELLS: " << M.cellIds.size();
  s << endl << M.cellIds << "FACES: " << M.faceIds.size();
  s << endl << M.faceIds;
#ifdef USE_DATAMESH
  s << "VDATA: " << M.vdataList.size();
  s << endl << M.vdataList << "CDATA: " << M.cdataList.size();
  s << endl;
  s << M.cdataList;
#endif
  s << endl;
  return s;
}

bool CoarseGeometry::insertInformation(string &&line) {
  trim(line);
  if (line.empty()) return true;

  std::vector<std::string> lineType = split(line, "=");
  if (lineType.size() != 2) { return false; }

  trim(lineType[0]);
  trim(lineType[1]);
  geoInformation.try_emplace(lowerCase(lineType[0]), lineType[1]);

  return true;
}

bool CoarseGeometry::SaveGeometryFile(const string &geometryFilename) {
  if (!geoInformation.empty() && geoInformation.find("checked") != geoInformation.end()
      && geoInformation["checked"] == "true")
    return true;
  try {
    std::ofstream file(geometrySavePath() + geometryFilename + ".geo");
    file.precision(14);

    file << "HEADER:" << std::endl;
    if (geoInformation.find("checked") == geoInformation.end()) file << "Checked=true" << std::endl;

    file << "WithData=" << (withData() ? "true" : "false") << std::endl;

    if (geoInformation.find("cellformat") == geoInformation.end())
      file << "Cellformat=Vtu" << std::endl;
    for (const auto &info : geoInformation) {
      if (info.first != "checked" || info.first != "withdata" || info.first != "cellformat")
        file << capitalize(info.first) << "=" << info.second << std::endl;
    }

    file << "POINTS:" << std::endl;
    int currentIndex = 0;
    for (auto v : coordinates) {
      file << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    file << "CELLS:" << std::endl;
    bool globaCelllType = geoInformation.find("celltype") != geoInformation.end();
    bool globalSubdomain = geoInformation.find("subdomain") != geoInformation.end();
    for (const auto &c : cellIds) {
      file << (globaCelllType ? "" : std::to_string(CelltypeToVtu(c.Type())) + " ");
      file << (globalSubdomain ? "" : std::to_string(c.Subdomain()) + " ");

      for (auto i : c.CornerIndices()) {
        file << i << " ";
      }
      file << std::endl;
    }

    file << "FACES:" << std::endl;
    bool globalFaceType = geoInformation.find("facetype") != geoInformation.end();
    bool globalBoundary = geoInformation.find("boundary") != geoInformation.end();
    // If Face type and boundary are set globally, we don't have to add the to .geo file
    if (!globalFaceType || !globalBoundary) {
      for (auto f : faceIds) {
        file << (globalBoundary ? "" : std::to_string(f.BndType()) + " ");

        for (auto i : f.FaceIndices()) {
          file << i << " ";
        }
        file << std::endl;
      }
    }

    if (withData()) {
      file << "VDATA:" << std::endl;
      for (auto data : vdataList) {
        for (auto d : data) {
          file << d << " ";
        }
        file << std::endl;
      }


      file << "CDATA:" << std::endl;
      for (auto data : cdataList) {
        for (auto d : data) {
          file << d << " ";
        }
        file << std::endl;
      }
    }
    file.close();
    return true;
  } catch (std::exception &e) { return false; }
}

bool CoarseGeometry::withData() {
  return (geoInformation.find("withdata") == geoInformation.end()) ?
             !(vdataList.empty() && cdataList.empty()) :
             (geoInformation["withdata"] == "true");
}

Coordinates RectangularGeometry2D::createPoints(const std::pair<double, double> &xBorders,
                                                const std::pair<double, double> &yBorders, int dimX,
                                                int dimY) {
  Coordinates coord;
  std::vector<double> xCoord = ft::linspace(xBorders.first, xBorders.second, dimX);
  std::vector<double> yCoord = ft::linspace(yBorders.first, yBorders.second, dimY);
  for (size_t j = 0; j < dimY; j++) {
    for (size_t i = 0; i < dimX; i++) {
      coord.emplace_back(xCoord[i], yCoord[j]);
    }
  }
  return coord;
}

CellIds RectangularGeometry2D::createCells(int dimX, int dimY) {
  CellIds ids;
  for (int j = 0; j < dimY - 1; j++) {
    for (int i = 0; i < dimX - 1; i++) {
      ids.emplace_back(
          CellId(QUADRILATERAL, 1,
                 {i + j * dimX, i + 1 + j * dimX, i + dimX + 1 + j * dimX, i + dimX + j * dimX}));
    }
  }
  return ids;
}

FaceIds RectangularGeometry2D::createFaces(int dimX, int dimY, std::vector<int> boundaries) {
  FaceIds ids;
  for (int i = 0; i < dimX - 1; ++i) {
    ids.emplace_back(FaceId(boundaries[0], {i, i + 1})); // left boundary
  }
  for (int j = 0; j < dimY - 1; ++j) {
    ids.emplace_back(
        FaceId(boundaries[2], {dimX - 1 + j * dimX, dimX - 1 + (j + 1) * dimX})); // bottom boundary
  }
  for (int i = 0; i < dimX - 1; ++i) {
    ids.emplace_back(
        FaceId(boundaries[1], {dimX * dimY - 1 - i, dimX * dimY - 1 - i - 1})); // right boundary
  }
  for (int j = 0; j < dimY - 1; ++j) {
    ids.emplace_back(FaceId(boundaries[3], {dimX * (dimY - 1) - j * dimX,
                                            dimX * (dimY - 1) - (j + 1) * dimX})); // top boundary
  }
  return ids;
}

class IntervalCoarseGeometry : public CoarseGeometry {
public:
  explicit IntervalCoarseGeometry(VertexDataList &&vertexData = {}, CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}}, CellIds{{INTERVAL, 1, {0, 1}}},
                     FaceIds{{1, {0}}, {1, {1}}}, std::move(vertexData), std::move(cellData),
                     "Interval") {}
};

class IntervalCoarseGeometryRefined : public CoarseGeometry {
public:
  explicit IntervalCoarseGeometryRefined(VertexDataList &&vertexData = {},
                                         CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}},
                     CellIds{{INTERVAL, 1, {0, 1}}, {INTERVAL, 1, {1, 2}}},
                     FaceIds{{1, {0}}, {1, {2}}}, std::move(vertexData), std::move(cellData),
                     "IntervalRefined") {}
};

class IntervalCoarseGeometry2Refined : public CoarseGeometry {
public:
  explicit IntervalCoarseGeometry2Refined(VertexDataList &&vertexData = {},
                                          CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {0.25, 0.0}, {0.5, 0.0}, {0.75, 0.0}, {1.0, 0.0}},
                     CellIds{{INTERVAL, 1, {0, 1}},
                             {INTERVAL, 1, {1, 2}},
                             {INTERVAL, 1, {2, 3}},
                             {INTERVAL, 1, {3, 4}}},
                     FaceIds{{1, {0}}, {1, {4}}}, std::move(vertexData), std::move(cellData),
                     "Interval2Refined") {}
};

class IntervalCoarseGeometry3Refined : public CoarseGeometry {
public:
  explicit IntervalCoarseGeometry3Refined(VertexDataList &&vertexData = {},
                                          CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0},
                                 {0.125, 0.0},
                                 {0.25, 0.0},
                                 {0.375, 0.0},
                                 {0.5, 0.0},
                                 {0.625, 0.0},
                                 {0.75, 0.0},
                                 {0.875, 0.0},
                                 {1.0, 0.0}},
                     CellIds{{INTERVAL, 1, {0, 1}},
                             {INTERVAL, 1, {1, 2}},
                             {INTERVAL, 1, {2, 3}},
                             {INTERVAL, 1, {3, 4}},
                             {INTERVAL, 1, {4, 5}},
                             {INTERVAL, 1, {5, 6}},
                             {INTERVAL, 1, {6, 7}},
                             {INTERVAL, 1, {7, 8}}},
                     FaceIds{{1, {0}}, {1, {8}}}, std::move(vertexData), std::move(cellData),
                     "Interval3Refined") {}
};

class TriangleCoarseGeometry : public CoarseGeometry {
public:
  explicit TriangleCoarseGeometry(VertexDataList &&vertexData = {}, CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
                     CellIds{{TRIANGLE, 1, {0, 1, 2}}},
                     FaceIds{{1, {0, 1}}, {1, {1, 2}}, {1, {2, 0}}}, std::move(vertexData),
                     std::move(cellData), "Triangle") {}
};

class SquareTurnedCoarseGeometry : public CoarseGeometry {
public:
  explicit SquareTurnedCoarseGeometry(VertexDataList &&vertexData = {},
                                      CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.70710678118, 0.0}, // sqrt(2)/2 to obtain unit length
                                 {0.0, 0.70710678118},
                                 {-0.70710678118, 0.0},
                                 {0.0, -0.70710678118}},
                     CellIds{{QUADRILATERAL, 1, {0, 1, 2, 3}}},
                     FaceIds{{1, {0, 1}}, {1, {1, 2}}, {1, {2, 3}}, {1, {3, 0}}},
                     std::move(vertexData), std::move(cellData), "SquareTurned") {}
};

class UnitSquareCoarseGeometry : public CoarseGeometry {
public:
  explicit UnitSquareCoarseGeometry(VertexDataList &&vertexData = {},
                                    CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}},
                     CellIds{{QUADRILATERAL, 1, {0, 1, 2, 3}}},
                     FaceIds{{1, {0, 1}}, {0, {1, 2}}, {2, {2, 3}}, {0, {3, 0}}},
                     std::move(vertexData), std::move(cellData), "UnitSquare") {}
};

class Square2TrianglesCoarseGeometry : public CoarseGeometry {
public:
  explicit Square2TrianglesCoarseGeometry(VertexDataList &&vertexData = {},
                                          CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}},
                     CellIds{{TRIANGLE, 1, {0, 1, 3}}, {TRIANGLE, 1, {2, 3, 1}}},
                     FaceIds{{1, {0, 1}}, {1, {1, 2}}, {1, {2, 3}}, {1, {3, 0}}},
                     std::move(vertexData), std::move(cellData), "Square2Triangles") {}
};

class Square4TrianglesCoarseGeometry : public CoarseGeometry {
public:
  explicit Square4TrianglesCoarseGeometry(VertexDataList &&vertexData = {},
                                          CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.5, 0.5}},
                     CellIds{{TRIANGLE, 1, {0, 1, 4}},
                             {TRIANGLE, 1, {1, 2, 4}},
                             {TRIANGLE, 1, {2, 3, 4}},
                             {TRIANGLE, 1, {3, 0, 4}}},
                     FaceIds{{1, {0, 1}}, {1, {1, 2}}, {1, {2, 3}}, {1, {3, 0}}},
                     std::move(vertexData), std::move(cellData), "Square4Triangles") {}
};

class Square4TrianglesCoarseGeometryTurned : public CoarseGeometry {
public:
  explicit Square4TrianglesCoarseGeometryTurned(VertexDataList &&vertexData = {},
                                                CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.70710678118, 0.0}, // sqrt(2)/2 to obtain unit length
                                 {0.0, 0.70710678118},
                                 {-0.70710678118, 0.0},
                                 {0.0, -0.70710678118},
                                 {0.0, 0.0}},
                     CellIds{{TRIANGLE, 1, {0, 1, 4}},
                             {TRIANGLE, 1, {1, 2, 4}},
                             {TRIANGLE, 1, {2, 3, 4}},
                             {TRIANGLE, 1, {3, 0, 4}}},
                     FaceIds{{1, {0, 1}}, {1, {1, 2}}, {1, {2, 3}}, {1, {3, 0}}},
                     std::move(vertexData), std::move(cellData), "Square4TrianglesTurned") {}
};

class SquareTriangleCoarseGeometry : public CoarseGeometry {
public:
  explicit SquareTriangleCoarseGeometry(VertexDataList &&vertexData = {},
                                        CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {2.0, 0.0}},
                     CellIds{{QUADRILATERAL, 1, {0, 1, 2, 3}}, {TRIANGLE, 1, {1, 4, 2}}},
                     FaceIds{{1, {0, 1}}, {1, {1, 4}}, {1, {4, 2}}, {1, {2, 3}}, {1, {3, 0}}},
                     std::move(vertexData), std::move(cellData), "SquareTriangle") {}
};

class TetrahedronCoarseGeometry : public CoarseGeometry {
public:
  explicit TetrahedronCoarseGeometry(VertexDataList &&vertexData = {},
                                     CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0, 0.0},
                                 {1.0, 0.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0}},
                     CellIds{{TETRAHEDRON, 1, {0, 1, 2, 3}}},
                     FaceIds{{1, {0, 1, 2}}, {1, {0, 1, 3}}, {1, {1, 2, 3}}, {1, {2, 0, 3}}},
                     std::move(vertexData), std::move(cellData), "Tetrahedron") {}
};

class HexahedronCoarseGeometry : public CoarseGeometry {
public:
  explicit HexahedronCoarseGeometry(VertexDataList &&vertexData = {},
                                    CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0, 0.0},
                                 {1.0, 0.0, 0.0},
                                 {1.0, 1.0, 0.0},
                                 {0.0, 1.0, 0.0},
                                 {0.0, 0.0, 1.0},
                                 {1.0, 0.0, 1.0},
                                 {1.0, 1.0, 1.0},
                                 {0.0, 1.0, 1.0}},
                     CellIds{{HEXAHEDRON, 1, {0, 1, 2, 3, 4, 5, 6, 7}}},
                     FaceIds{{1, {0, 1, 2, 3}},
                             {1, {0, 1, 5, 4}},
                             {1, {1, 2, 6, 5}},
                             {1, {2, 3, 7, 6}},
                             {1, {3, 0, 4, 7}},
                             {1, {4, 5, 6, 7}}},
                     std::move(vertexData), std::move(cellData), "Hexahedron") {}
};

class SpaceTimeUnitInterval : public CoarseGeometry {
public:
  explicit SpaceTimeUnitInterval(VertexDataList &&vertexData = {}, CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0}, {1.0}}, CellIds{{INTERVAL, 1, {0, 1}}},
                     FaceIds{{1, {0}}, {1, {1}}}, TimeSteps{0.0, 1.0}, std::move(vertexData),
                     std::move(cellData), "SpaceTimeUnitInterval") {}
};

class STSquareCoarseGeometry : public RectangularGeometry2D {
public:
  STSquareCoarseGeometry() :
      RectangularGeometry2D(std::pair<double, double>({0, 1}), std::pair<double, double>({0, 1}), 2,
                            2, std::vector<int>({1, 1, 1, 1}), TimeSteps{0.0, 1.0}, "STSquare") {}
};

class STSquareTrianglesCoarseGeometry : public CoarseGeometry {
public:
  explicit STSquareTrianglesCoarseGeometry(VertexDataList &&vertexData = {},
                                           CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}},
                     CellIds{{TRIANGLE, 1, {0, 1, 2}}, {TRIANGLE, 1, {0, 2, 3}}}, FaceIds{},
                     TimeSteps{0.0, 1.0}, std::move(vertexData), std::move(cellData),
                     "STSquareTrianglesCoarseGeometry") {}
};

class STSquareCenteredCoarseGeometry : public RectangularGeometry2D {
public:
  STSquareCenteredCoarseGeometry() :
      RectangularGeometry2D(std::pair<double, double>({-0.5, 0.5}),
                            std::pair<double, double>({0, 1}), 2, 2, std::vector<int>({1, 1, 1, 1}),
                            TimeSteps{0.0, 0.5, 1.0}, "STSquareCentered") {}
};

class STTubeSquaresCoarseGeometry : public RectangularGeometry2D {
public:
  STTubeSquaresCoarseGeometry() :
      RectangularGeometry2D(std::pair<double, double>({-2, 4}), std::pair<double, double>({0, 1}),
                            7, 2, std::vector<int>({1, 1, 1, 1}), TimeSteps{0, 1, 2, 3, 4},
                            "STTubeSquares") {}
};

class FWIExampleCoarseGeometry : public RectangularGeometry2D {
public:
  FWIExampleCoarseGeometry() :
      RectangularGeometry2D(std::pair<double, double>({-1.125, 1.5}),
                            std::pair<double, double>({0, 1.25}), 21, 10,
                            std::vector<int>({1, 1, 1, 0}),
                            TimeSteps{0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125,
                                      1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0},
                            "FWIExampleCoarseGeometry") {}
};

class SpaceTimeUnitCube : public CoarseGeometry {
public:
  explicit SpaceTimeUnitCube(VertexDataList &&vertexData = {}, CellDataList &&cellData = {}) :
      CoarseGeometry(Coordinates{{0, 0, 0},
                                 {1, 0, 0},
                                 {1, 1, 0},
                                 {0, 1, 0},
                                 {0, 0, 1},
                                 {1, 0, 1},
                                 {1, 1, 1},
                                 {0, 1, 1}},
                     CellIds{{HEXAHEDRON, 0, {0, 1, 2, 3, 4, 5, 6, 7}}},
                     FaceIds{{1, {0}}, {1, {1}}}, TimeSteps{0.0, 1.0}, std::move(vertexData),
                     std::move(cellData), "SpaceTimeUnitCube") {}
};

CoarseGeometry *CreateCoarseGeometry(const std::string &meshName, VertexDataList &&vertexData,
                                     CellDataList &&cellData) {
  if (meshName == "Interval" || meshName == "Interval_DirichletBC") {
    return new IntervalCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "IntervalRefined" || meshName == "IntervalRefined_DirichletBC") {
    return new IntervalCoarseGeometryRefined(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Interval2Refined" || meshName == "Interval2Refined_DirichletBC") {
    return new IntervalCoarseGeometry2Refined(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Interval3Refined" || meshName == "Interval3Refined_DirichletBC") {
    return new IntervalCoarseGeometry3Refined(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Triangle" || meshName == "Triangle_DirichletBC") {
    return new TriangleCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Square" || meshName == "Square_DirichletBC") {
    return new RectangularGeometry2D(std::pair<double, double>({0, 1}),
                                     std::pair<double, double>({0, 1}), 2, 2,
                                     std::vector<int>({1, 1, 1, 1}), std::move(vertexData),
                                     std::move(cellData), "Square");
  }

  if (meshName == "Square10x1") {
    return new RectangularGeometry2D(std::pair<double, double>({0, 10}),
                                     std::pair<double, double>({0, 1}), 11, 2,
                                     std::vector<int>({1, 1, 1, 1}), std::move(vertexData),
                                     std::move(cellData), "Square10x1");
  }

  if (meshName == "SquareTurned" || meshName == "SquareTurned_DirichletBC") {
    return new SquareTurnedCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "UnitSquare") {
    return new UnitSquareCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "SquareRefined" || meshName == "SquareRefined_DirichletBC") {
    return new RectangularGeometry2D(std::pair<double, double>({0, 1}),
                                     std::pair<double, double>({0, 1}), 3, 3,
                                     std::vector<int>({1, 1, 1, 1}), std::move(vertexData),
                                     std::move(cellData), "SquareRefined");
  }

  if (meshName == "Square2Triangles" || meshName == "Square2Triangles_DirichletBC") {
    return new Square2TrianglesCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "SquareTriangle" || meshName == "SquareTriangle_DirichletBC") {
    return new SquareTriangleCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Square4Triangles" || meshName == "Square4Triangles_DirichletBC") {
    return new Square4TrianglesCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Square4TrianglesTurned" || meshName == "Square4TrianglesTurned_DirichletBC") {
    return new Square4TrianglesCoarseGeometryTurned(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Tetrahedron" || meshName == "Tetrahedron_DirichletBC") {
    return new TetrahedronCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "Hexahedron" || meshName == "Hexahedron_DirichletBC") {
    return new HexahedronCoarseGeometry(std::move(vertexData), std::move(cellData));
  }

  if (meshName == "SpaceTimeSquare" || meshName == "STSquare") {
    return new STSquareCoarseGeometry();
  }

  if (meshName == "SpaceTimeSquareTriangles") { return new STSquareTrianglesCoarseGeometry(); }

  if (meshName == "STTubeSquares") { return new STTubeSquaresCoarseGeometry(); }

  if (meshName == "SpaceTimeSquareCentered" || meshName == "STSquareCentered") {
    return new STSquareCenteredCoarseGeometry();
  }

  if (meshName == "SpaceTimeUnitInterval") { return new SpaceTimeUnitInterval(); }

  if (meshName == "SpaceTimeUnitCube") { return new SpaceTimeUnitCube(); }

  if (meshName == "FWIExampleGeometry") { return new FWIExampleCoarseGeometry(); }

  return new CoarseGeometry(meshName);
}

std::unique_ptr<CoarseGeometry> CreateCoarseGeometryUnique(const std::string &meshName,
                                                           VertexDataList &&vertexData,
                                                           CellDataList &&cellData) {
  return std::unique_ptr<CoarseGeometry>(
      CreateCoarseGeometry(meshName, std::move(vertexData), std::move(cellData)));
}

std::shared_ptr<CoarseGeometry> CreateCoarseGeometryShared(const std::string &meshName,
                                                           VertexDataList &&vertexData,
                                                           CellDataList &&cellData) {
  return std::shared_ptr<CoarseGeometry>(
      CreateCoarseGeometry(meshName, std::move(vertexData), std::move(cellData)));
}
