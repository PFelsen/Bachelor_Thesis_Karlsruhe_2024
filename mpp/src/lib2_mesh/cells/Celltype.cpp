#include "Celltype.hpp"
#include "Assertion.hpp"

std::string to_string(CELLTYPE type) {
  switch (type) {
  case INTERVAL:
    return "Interval";
  case TRIANGLE:
    return "Triangle";
  case QUADRILATERAL:
    return "Square";
  case TETRAHEDRON:
    return "Tetrahedron";
  case HEXAHEDRON:
    return "Hexahedron";
  default:
    THROW("Celltype not implemented!")
  }
}
