#ifndef VTUCONVERTER_HPP
#define VTUCONVERTER_HPP

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <utility>

#include "Meshes.hpp"

class VtuConverter {
  std::string fileName{};
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  std::shared_ptr<CoarseGeometry> geometry;
  std::unique_ptr<Meshes> meshes;

  void extractFileName(std::string filePath);

  std::unordered_map<Point, int> vertexIndex{};
  std::unordered_map<int, Point> indexVertex{};
  std::unordered_map<int, Point> indexCell{};

  std::vector<int> subdomainData{};
  std::vector<int> smoothingBndCnd{};

  std::function<int(int)> convertSD = [](int sd) { return sd; };

  struct pair_hash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &pair) const {
      return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
  };

  std::unordered_map<std::pair<int, int>, std::vector<Point>, pair_hash> boundaryData;

  CellId readVtkCell(vtkCell &cell, int sd);

  std::vector<std::vector<double>> cellData{};
  std::vector<std::vector<double>> vertexData{};

  CoarseGeometry extractGeometry(const vtkSmartPointer<vtkUnstructuredGrid> &grid,
                                 const string &name);

  void smoothBoundaryData(FaceIds &faceIds);
public:
  explicit VtuConverter(const std::string &filePath);

  VtuConverter &UseCellArrayAsSubdomain(int cellArrayIndex);

  VtuConverter &UseArrayForCellData(int cellArrayIndex);

  VtuConverter &UseArrayForVertexData(int pointArrayIndex);

  const Mesh &Create();

  bool WriteGeoFile(const Mesh &M, const string &path);

  VtuConverter &WithExcitationFromFile(int excitationArrayIndex, double duration = 0.003,
                                       double amplitude = 30.0);

  VtuConverter &WithBoundaryFromData(std::pair<int, int> bc, vector<Point> &&points);

  void readFaces(const CellId &cellID, const std::vector<Point> &meshPoints,
                 std::vector<FaceId> &faceIds);


  VtuConverter &WithBoundarySmoothing(std::vector<int> bndConditions);
  VtuConverter &WithSubdomainConversion(std::function<int(int)> conversionFct);
};

std::string loadFile(std::string name);


#endif // VTUCONVERTER_HPP
