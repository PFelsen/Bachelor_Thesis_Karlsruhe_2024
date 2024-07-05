#ifndef CELLTYPE_HPP
#define CELLTYPE_HPP

#include "Point.hpp"

enum CELLTYPE {
  NONE = 0,
  POINT = 1,
  INTERVAL = 122,
  TINTERVAL = 144,
  TRIANGLE = 233,
  QUADRILATERAL = 244,
  QUADRILATERAL2 = 248,
  TETRAHEDRON = 344,
  PYRAMID = 355,
  PRISM = 366,
  HEXAHEDRON = 388,
  HEXAHEDRON20 = 3820,
  HEXAHEDRON27 = 3827,

  SPACETIME_INTERVAL = 1122,
  SPACETIME_QUADRILATERAL = 1244,
  SPACETIME_HEXAHEDRON = 1388,

  SPACETIME_TRIANGLE = 1233
};

std::string to_string(CELLTYPE type);

constexpr bool isSpaceTimeCellType(const CELLTYPE &type) {
  return type == SPACETIME_INTERVAL || type == SPACETIME_QUADRILATERAL
         || type == SPACETIME_HEXAHEDRON || type == SPACETIME_TRIANGLE;
}

constexpr CELLTYPE SpaceTimeCellType(const CELLTYPE &type) {
  if (type == QUADRILATERAL) {
    return SPACETIME_QUADRILATERAL;
  } else if (type == INTERVAL) {
    return SPACETIME_INTERVAL;
  } else if (type == HEXAHEDRON) {
    return SPACETIME_HEXAHEDRON;
  } else if (type == TRIANGLE) {
    return SPACETIME_TRIANGLE;
  }
  return type;
}

constexpr CELLTYPE SpaceCellType(const CELLTYPE &type) {
  if (type == SPACETIME_QUADRILATERAL) {
    return QUADRILATERAL;
  } else if (type == SPACETIME_INTERVAL) {
    return INTERVAL;
  } else if (type == SPACETIME_HEXAHEDRON) {
    return HEXAHEDRON;
  } else if (type == SPACETIME_TRIANGLE) {
    return TRIANGLE;
  }
  return type;
}

constexpr CELLTYPE FaceCellType(const CELLTYPE &type) {
  switch (type) {
  case INTERVAL:
    return POINT;
  case TRIANGLE:
    return INTERVAL;
  case QUADRILATERAL:
  case SPACETIME_INTERVAL:
    return INTERVAL;
  case TETRAHEDRON:
    return TRIANGLE;
  case HEXAHEDRON:
  case SPACETIME_QUADRILATERAL:
    return QUADRILATERAL;
  case SPACETIME_HEXAHEDRON:
    return HEXAHEDRON;
  default:
    THROW("Celltype not implemented!")
  }
}

constexpr std::array<std::pair<CELLTYPE, int>, 5>
    celltype_vtu{std::pair<CELLTYPE, int>{INTERVAL, 3}, std::pair<CELLTYPE, int>{TRIANGLE, 5},
                 std::pair<CELLTYPE, int>{QUADRILATERAL, 9},
                 std::pair<CELLTYPE, int>{TETRAHEDRON, 10},
                 std::pair<CELLTYPE, int>{HEXAHEDRON, 12}};

static const std::array<std::pair<std::string, int>, 5>
    cellname_vtu{std::pair<std::string, int>{"Interval", 3},
                 std::pair<std::string, int>{"Triangle", 5},
                 std::pair<std::string, int>{"Quadrilateral", 9},
                 std::pair<std::string, int>{"Tetrahedron", 10},
                 std::pair<std::string, int>{"Hexahedron", 12}};

constexpr CELLTYPE vtuToCelltype(int vtuCellType) {
  for (auto pair : celltype_vtu) {
    if (pair.second == vtuCellType) return pair.first;
  }
  return NONE;
}

constexpr int CelltypeToVtu(CELLTYPE type) {
  for (const auto &pair : celltype_vtu) {
    if (pair.first == type) return pair.second;
  }
  return -1;
}

static int CellnameToVtu(const std::string &name) {
  for (const auto &pair : cellname_vtu) {
    if (pair.first == name) return pair.second;
  }
  return -1;
}

constexpr int VtkCellType(CELLTYPE cType) { return CelltypeToVtu(cType); }

constexpr int CellDim(CELLTYPE cType) {
  switch (cType) {
  case NONE:
  case POINT:
    return 0;
  case INTERVAL:
  case TINTERVAL:
  case SPACETIME_INTERVAL:
    return 1;
  case TRIANGLE:
  case QUADRILATERAL:
  case QUADRILATERAL2:
  case SPACETIME_QUADRILATERAL:
  case SPACETIME_TRIANGLE:
    return 2;
  case TETRAHEDRON:
  case PRISM:
  case PYRAMID:
  case HEXAHEDRON:
  case HEXAHEDRON20:
  case HEXAHEDRON27:
  case SPACETIME_HEXAHEDRON:
    return 3;
  default:
    THROW("Celltype not implemented!")
  }
}

#endif // CELLTYPE_HPP
