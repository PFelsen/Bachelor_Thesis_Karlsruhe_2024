#include "VtuConverter.hpp"
#include "MeshesCreator.hpp"

#include <vtkCellData.h>
#include <vtkPointData.h>

#include <utility>

std::string loadFile(std::string name) {
  std::string fileName = std::string(ProjectSourceDir) + "/" + name;
  std::cout << "Loading: " << fileName << std::endl;

  if (!file_exists(fileName.c_str())) { Warning("File not found") return 0; }
  return fileName;
}

std::vector<Point> CornersFromIds(const CellId &cellId, const Coordinates &coordinates) {
  const std::vector<int> &ids = cellId.CornerIndices();
  std::vector<Point> corners(ids.size());
  for (int i = 0; i < ids.size(); ++i)
    corners[i] = coordinates[ids[i]];
  return corners;
}

CellId VtuConverter::readVtkCell(vtkCell &cell, int sd) {
  auto &cellCornerIds = *(cell.GetPointIds());
  std::vector<int> ids(cellCornerIds.GetNumberOfIds());
  for (int i = 0; i < ids.size(); ++i) {
    ids[i] = cellCornerIds.GetId(i);
  }

  return CellId(vtuToCelltype(cell.GetCellType()), sd, std::move(ids));
}

void VtuConverter::readFaces(const CellId &cellID, const std::vector<Point> &meshPoints,
                             std::vector<FaceId> &faceIds) {
  auto ids = cellID.CornerIndices();
  auto corners = CornersFromIds(cellID, meshPoints);

  std::map<Point, int> PointID{};
  for (int i = 0; i < corners.size(); ++i) {
    PointID[corners[i]] = ids[i];
  }

  auto c = std::unique_ptr<Cell>(CreateCell(cellID.Type(), cellID.Subdomain(), corners));

  for (int f = 0; f < c->Faces(); ++f) {
    for (const auto &bndData : boundaryData) {
      const auto &bndCondition = bndData.first;
      const auto &bndPoints = bndData.second;
      bool onBoundary = true;
      for (int fc = 0; fc < c->FaceCorners(f); ++fc) {
        auto faceCorner = c->FaceCorner(f, fc);
        if (std::find_if(bndPoints.begin(), bndPoints.end(),
                         [&faceCorner](const Point &p) { return p.isNear(faceCorner, 1e-4); })
            == bndPoints.end()) {
          onBoundary = false;
          break;
        }
      }
      if (!onBoundary) continue;

      std::vector<int> fIndices(c->FaceCorners(f));
      for (int fc = 0; fc < c->FaceCorners(f); ++fc) {
        fIndices[fc] = PointID[c->FaceCorner(f, fc)];
      }

      faceIds.emplace_back(cellID.Subdomain() < 4000 ? bndCondition.first : bndCondition.second,
                           std::move(fIndices));
    }
  }
}

void VtuConverter::smoothBoundaryData(FaceIds &faceIds) {
  auto isApplicable = [this](const FaceId &faceId) {
    return std::find(smoothingBndCnd.begin(), smoothingBndCnd.end(), faceId.BndType())
           != smoothingBndCnd.end();
  };
  FaceIds applicableIds;
  std::copy_if(faceIds.begin(), faceIds.end(), std::back_inserter(applicableIds), isApplicable);

  faceIds.erase(std::remove_if(faceIds.begin(), faceIds.end(),
                               [faceIds, applicableIds, isApplicable](const FaceId &faceId1) {
                                 if (!isApplicable(faceId1)) return false;
                                 for (auto &faceId2 : applicableIds) {
                                   if (faceId1.BndType() == faceId2.BndType()) { continue; }
                                   std::vector<int> fIntersection;
                                   auto fIds1 = faceId1.FaceIndices();
                                   std::sort(fIds1.begin(), fIds1.end());
                                   auto fIds2 = faceId2.FaceIndices();
                                   std::sort(fIds2.begin(), fIds2.end());
                                   std::set_intersection(fIds1.begin(), fIds1.end(), fIds2.begin(),
                                                         fIds2.end(),
                                                         std::back_inserter(fIntersection));
                                   if (!fIntersection.empty()
                                       && fIntersection.size() != fIds1.size())
                                     return true;
                                 }
                                 return false;
                               }),
                faceIds.end());
}

// M++ Mesh creation
CoarseGeometry VtuConverter::extractGeometry(const vtkSmartPointer<vtkUnstructuredGrid> &grid,
                                             const std::string &name) {
  auto &pointData = *(grid->GetPoints()->GetData());


  std::vector<Point> meshPoints(pointData.GetNumberOfTuples());
  Coordinates trimmedPoints(meshPoints.size());
  std::vector<DataContainer> trimmedVertexData(meshPoints.size());
  if (vertexData.empty()) {
    vertexData = std::vector<std::vector<double>>(meshPoints.size());
    for (int i = 0; i < vertexData.size(); ++i) {
      vertexData[i] = {-10.0, 0.0, 0.0};
    }
  }

  for (int pointID = 0; pointID < pointData.GetNumberOfTuples(); ++pointID) {
    // auto usedPoint = usedVerteces.find(pointID);

    // if (usedPoint != usedVerteces.end()) {
    auto point = pointData.GetTuple(pointID);
    meshPoints[pointID] = Point(point[0], point[1], point[2]);
    // usedVerteces.erase(usedPoint);
    // } else {
    //   meshPoints[pointID] = Infty;
    // }
  }

  // Trim Point Array
  // std::unordered_map<int, int> pointIndexConversion{};
  // int conversionIndex = 0;
  for (int i = 0; i < meshPoints.size(); ++i) {
    trimmedPoints[i] = Point(meshPoints[i][0], meshPoints[i][1], meshPoints[i][2]);
    trimmedVertexData[i] = DataContainer{std::move(vertexData[i])};
    /*if (meshPoints[i] == Infty)
      continue;

    vertexIndex[meshPoints[i]] = conversionIndex;
    trimmedPoints[conversionIndex] = std::vector<double>{meshPoints[i][0], meshPoints[i][1],
                                                         meshPoints[i][2]};
    trimmedVertexData[conversionIndex] = DataContainer{std::move(vertexData[i])};
    pointIndexConversion[i] = conversionIndex++;*/
  }

  // Get Maximum Cell Dimension to exclude Faces
  /*int maxDim = 0;
  for (int cellID = 0; cellID < grid->GetNumberOfCells(); ++cellID) {
    maxDim = max(maxDim, grid->GetCell(cellID)->GetCellDimension());
  }
  int realCells = 0;
  for (int cellID = 0; cellID < grid->GetNumberOfCells(); ++cellID) {
    realCells += (maxDim == grid->GetCell(cellID)->GetCellDimension());
  }*/

  // Read only real cells
  CellIds meshCells{};
  std::vector<DataContainer> sortedCellData{};

  FaceIds meshFaces{};

  // std::unordered_set<int> usedVerteces{};
  for (int cellID = 0; cellID < grid->GetNumberOfCells(); ++cellID) {
    auto cell = grid->GetCell(cellID);
    // if (cell->GetCellDimension() < maxDim)
    //   continue;
    auto convertedCell =
        readVtkCell(*cell, subdomainData.empty() ? 0 : convertSD(subdomainData[cellID]));
    /* for (int i = 2; i < convertedCell.size(); ++i) {
       usedVerteces.insert(convertedCell[i]);
     }*/
    readFaces(convertedCell, trimmedPoints, meshFaces);

    meshCells.emplace_back(std::move(convertedCell));
    sortedCellData.emplace_back(std::move(cellData[cellID]));
  }

  // Convert Indices in Cell Array
  /*for (auto &cell: meshCells) {
    for (int i = 2; i < cell.size(); ++i) {
      cell[i] = pointIndexConversion[cell[i]];
    }
  }

  std::cout << meshCells.size() << " - " << sortedCellData.size() << std::endl;*/

  smoothBoundaryData(meshFaces);

  CoarseGeometry geo(std::move(trimmedPoints), std::move(meshCells), std::move(meshFaces),
                     std::move(trimmedVertexData), std::move(sortedCellData), name);
  return geo;
}

bool VtuConverter::WriteGeoFile(const Mesh &M, const string &path) {
  try {
    ofstream file(path + M.Name() + ".geo");


    file << "POINTS:" << std::endl;
    int currentIndex = 0;
    for (vertex v = M.vertices(); v != M.vertices_end(); ++v, ++currentIndex) {
      vertexIndex[v()] = currentIndex;
      indexVertex[currentIndex] = v();
      file << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    file << "CELLS:" << std::endl;
    currentIndex = 0;
    for (cell c = M.cells(); c != M.cells_end(); ++c, ++currentIndex) {
      indexCell[currentIndex] = c();
      file << c.Corners() << " " << c.Subdomain();
      for (int i = 0; i < c.Corners(); ++i) {
        file << " " << vertexIndex[c.Corner(i)];
      }
      file << std::endl;
    }

    file << "FACES:" << std::endl;
    for (cell c = M.cells(); c != M.cells_end(); ++c) {
      for (int f = 0; f < c.Faces(); ++f) {
        auto bndFace = M.find_bnd_face(c.Face(f));
        if (bndFace != M.bnd_faces_end()) {
          if (bndFace.Part() == 0) continue;
          file << c.FaceCorners(f) << " " << bndFace.Part();
          for (int j = 0; j < c.FaceCorners(f); ++j) {
            file << " " << vertexIndex[c.FaceCorner(f, j)];
          }
          file << std::endl;
        }
      }
    }

    file << "VDATA:" << std::endl;
    for (int i = 0; i < indexVertex.size(); ++i) {
      DataContainer data = M.find_vertex(indexVertex[i]).GetData();
      file << data.size();
      for (auto d : data) {
        file << " " << d;
      }
      file << std::endl;
    }

    file << "CDATA:" << std::endl;
    for (int i = 0; i < indexCell.size(); ++i) {
      DataContainer data = M.find_cell(indexCell[i]).GetData();
      file << data.size();
      for (auto d : data) {
        file << " " << d;
      }
      file << std::endl;
    }

    return true;
  } catch (std::exception &e) { return false; }
}

VtuConverter::VtuConverter(const string &filePath) {
  extractFileName(filePath);

  std::string extension = "";
  if (filePath.find_last_of(".") != std::string::npos) {
    extension = filePath.substr(filePath.find_last_of("."));
  }
  // Drop the case of the extension
  std::transform(extension.begin(), extension.end(), extension.begin(), tolower);

  if (extension == ".vtu") {
    vtkNew<vtkXMLUnstructuredGridReader> reader;
    reader->SetFileName(filePath.c_str());
    reader->Update();
    unstructuredGrid = reader->GetOutput();
  }

  unstructuredGrid->Print(std::cout);
}

const Mesh &VtuConverter::Create() {
  geometry = std::make_shared<CoarseGeometry>(extractGeometry(unstructuredGrid, fileName));

  meshes = MeshesCreator("This does nothing").WithCoarseGeometry(geometry).CreateUnique();

  return meshes->coarse();
}

void VtuConverter::extractFileName(std::string filePath) {
  filePath = filePath.substr(0, filePath.find(".vtu"));
  size_t found = filePath.find_last_of("/");
  if (found != filePath.npos) { filePath = filePath.substr(found + 1); }
  fileName = filePath;
}

VtuConverter &VtuConverter::UseCellArrayAsSubdomain(int cellArrayIndex) {
  auto vtuData = unstructuredGrid->GetCellData()->GetArray(cellArrayIndex);
  subdomainData.resize(vtuData->GetNumberOfTuples());
  for (int i = 0; i < subdomainData.size(); ++i) {
    auto tuple = vtuData->GetTuple(i);
    subdomainData[i] = (int)(vtuData->GetTuple(i))[0];
  }
  return *this;
}

VtuConverter &VtuConverter::UseArrayForCellData(int cellArrayIndex) {
  auto vtuData = unstructuredGrid->GetCellData()->GetArray(cellArrayIndex);
  if (cellData.empty()) { cellData.resize(vtuData->GetNumberOfTuples()); }
  for (int i = 0; i < cellData.size(); ++i) {
    auto tuple = vtuData->GetTuple(i);
    for (int j = 0; j < vtuData->GetNumberOfComponents(); ++j) {
      cellData[i].emplace_back(tuple[j]);
    }
  }
  return *this;
}

VtuConverter &VtuConverter::UseArrayForVertexData(int pointArrayIndex) {
  auto vtuData = unstructuredGrid->GetPointData()->GetArray(pointArrayIndex);
  if (vertexData.empty()) { vertexData.resize(vtuData->GetNumberOfTuples()); }
  for (int i = 0; i < vertexData.size(); ++i) {
    auto tuple = vtuData->GetTuple(i);
    for (int j = 0; j < vtuData->GetNumberOfComponents(); ++j) {
      vertexData[i].emplace_back(tuple[j]);
    }
  }
  return *this;
}

VtuConverter &VtuConverter::WithExcitationFromFile(int excitationArrayIndex, double duration,
                                                   double amplitude) {
  auto vtuData = unstructuredGrid->GetPointData()->GetArray(excitationArrayIndex);
  if (vertexData.empty()) { vertexData.resize(vtuData->GetNumberOfTuples()); }
  for (int i = 0; i < vertexData.size(); ++i) {
    auto excitation = vtuData->GetTuple(i)[0];
    vertexData[i] = excitation < 0 ? std::vector<double>{-10.0, 0.0, 0.0} :
                                     std::vector<double>{excitation, duration, amplitude};
  }
  return *this;
}

VtuConverter &VtuConverter::WithBoundaryFromData(std::pair<int, int> bc, vector<Point> &&points) {
  boundaryData.try_emplace(bc, std::move(points));
  return *this;
}

VtuConverter &VtuConverter::WithBoundarySmoothing(std::vector<int> bndConditions) {
  smoothingBndCnd = std::move(bndConditions);
  return *this;
}

VtuConverter &VtuConverter::WithSubdomainConversion(std::function<int(int)> conversionFct) {
  convertSD = std::move(conversionFct);
  return *this;
}
