#ifndef CHECKORIENTATION_HPP
#define CHECKORIENTATION_HPP

#include "Celltype.hpp"
#include "Point.hpp"

class CheckOrientation {
public:
  static void CheckInterval(std::vector<Point> &corners, std::vector<int> &cornerIndices);

  static void CheckTriangle(std::vector<Point> &corners, std::vector<int> &cornerIndices);

  static void CheckQuadrilateral(std::vector<Point> &corners, std::vector<int> &cornerIndices);

  static void CheckTetrahedron(std::vector<Point> &corners, std::vector<int> &cornerIndices);

  static void CheckHexahedron(std::vector<Point> &corners, std::vector<int> &cornerIndices);
  static void Check(CELLTYPE type, std::vector<Point> &corners, std::vector<int> &cornerIndices);
};

#endif // CHECKORIENTATION_HPP
