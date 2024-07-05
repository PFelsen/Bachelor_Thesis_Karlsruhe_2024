#include "Quadrature.hpp"


#include <unordered_map>

std::unordered_map<std::string, Quadrature *> Quads;

const QemptyT<> quadEmpty;

const QpointT<> quadPoint;

template<>
const QuadratureT<> &GetQuadratureT(const CELLTYPE cellType, int exactUpTo, int exactUpTo1,
                                    int exactUpTo2) {
  if (cellType == NONE) return quadEmpty;
  if (cellType == POINT) return quadPoint;

  int type = cellType;
  std::string key = std::to_string(type)
                        .append("_")
                        .append(std::to_string(exactUpTo))
                        .append("_")
                        .append(std::to_string(exactUpTo1))
                        .append("_")
                        .append(std::to_string(exactUpTo2));

  if (Quads.empty()) Quads.reserve(1);
  auto quad = Quads.find(key);

  if (quad != Quads.end()) return *(quad->second);

  switch (cellType) {
  case INTERVAL:
    Quads[key] = new Qint(exactUpTo);
    break;
  case TRIANGLE:
    Quads[key] = new QtriSym(exactUpTo);
    break;
  case QUADRILATERAL:
  case SPACETIME_INTERVAL:
    Quads[key] = new Qquad(exactUpTo, exactUpTo1);
    break;
  case TETRAHEDRON:
    Quads[key] = new QtetSym(exactUpTo);
    break;
  case PRISM:
    Quads[key] = new QpriSym(exactUpTo, exactUpTo1);
    break;
  case HEXAHEDRON:
  case SPACETIME_QUADRILATERAL:
    Quads[key] = new Qhex(exactUpTo, exactUpTo1, exactUpTo2);
    break;
  default:
    THROW("No Quadrature implemented for given cell type!")
  }

  quad = Quads.find(key);
  if (quad != Quads.end()) return *(quad->second);

  THROW("Quadrature not implemented!")
}

#ifdef BUILD_IA

std::unordered_map<std::string, QuadratureT<IAInterval, SpaceDimension, TimeDimension> *> IAQuads;

const QemptyT<IAInterval, SpaceDimension, TimeDimension> IAquadEmpty;

const QpointT<IAInterval, SpaceDimension, TimeDimension> IAquadPoint;

template<>
const QuadratureT<IAInterval, SpaceDimension, TimeDimension> &
GetQuadratureT(const CELLTYPE cellType, int exactUpTo, int exactUpTo1, int exactUpTo2) {
  if (cellType == NONE) return IAquadEmpty;
  if (cellType == POINT) return IAquadPoint;

  int type = cellType;
  std::string key = std::to_string(type)
                        .append("_")
                        .append(std::to_string(exactUpTo))
                        .append("_")
                        .append(std::to_string(exactUpTo1))
                        .append("_")
                        .append(std::to_string(exactUpTo2));

  if (IAQuads.empty()) IAQuads.reserve(1);
  auto quad = IAQuads.find(key);

  if (quad != IAQuads.end()) return *(quad->second);

  switch (cellType) {
  case INTERVAL:
    IAQuads[key] = new IAQint(exactUpTo);
    break;
  case TRIANGLE:
    IAQuads[key] = new IAQtriSym(exactUpTo);
    break;
  case QUADRILATERAL:
    IAQuads[key] = new IAQquad(exactUpTo, exactUpTo1);
    break;
  case TETRAHEDRON:
    IAQuads[key] = new IAQtetSym(exactUpTo);
    break;
  case PRISM:
    IAQuads[key] = new IAQpriSym(exactUpTo, exactUpTo1);
    break;
  case HEXAHEDRON:
    IAQuads[key] = new IAQhex(exactUpTo, exactUpTo1, exactUpTo2);
    break;
  default:
    THROW("No interval Quadrature implemented for given cell type!")
  }

  quad = IAQuads.find(key);
  if (quad != IAQuads.end()) return *(quad->second);

  THROW("Quadrature not implemented!")
}

#endif

const Quadrature &GetQuadrature(const CELLTYPE cellType, int exactUpTo, int exactUpTo1,
                                int exactUpTo2) {
  return GetQuadratureT<double, SpaceDimension, TimeDimension>(cellType, exactUpTo, exactUpTo1,
                                                               exactUpTo2);
}

template<typename T, int sDim, int tDim>
const QuadratureT<T, sDim, tDim> &GetGLQuadratureT(const CELLTYPE cellType, int npcount) {
  if (cellType == NONE) return quadEmpty;
  if (cellType == POINT) return quadPoint;

  int type = cellType;
  std::string key =
      std::string("GL").append(std::to_string(type)).append("_").append(std::to_string(npcount));

  if (Quads.empty()) Quads.reserve(1);
  auto quad = Quads.find(key);

  if (quad != Quads.end()) return *(quad->second);

  switch (cellType) {
  case INTERVAL:
    Quads[key] = new QintGaussLobatto(npcount);
    break;
  case QUADRILATERAL: {
    const QuadratureT<> &oneD = GetGLQuadratureT<T, sDim, tDim>(INTERVAL, npcount);
    Quads[key] = new QTensor({oneD, oneD});
    break;
  }
  case HEXAHEDRON: {
    const QuadratureT<> &oneD = GetGLQuadratureT<T, sDim, tDim>(INTERVAL, npcount);
    Quads[key] = new QTensor({oneD, oneD, oneD});
    break;
  }
  default:
    THROW("No Quadrature implemented for given cell type!")
  }

  quad = Quads.find(key);
  if (quad != Quads.end()) return *(quad->second);

  THROW("Quadrature not implemented!")
}

const Quadrature &GetGLQuadrature(const CELLTYPE cellType, int npcount) {
  return GetGLQuadratureT<double, SpaceDimension, TimeDimension>(cellType, npcount);
}

template<class T, int sDim, int tDim>
const QuadratureT<T, sDim, tDim> &GetQuadrature(const std::string &name) {
  if (name == "Qempty") return GetQuadratureT<T, sDim, tDim>(NONE);
  if (name == "Qpoint1") return quadPoint;
  if (name == "Qint1") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 1);
  if (name == "Qint2") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 3);
  if (name == "Qint3") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 5);
  if (name == "Qint4") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 7);
  if (name == "Qint5") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 9);
  if (name == "Qint6") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 11);
  if (name == "Qint7") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 13);
  if (name == "Qint8") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 15);
  if (name == "Qint9") return GetQuadratureT<T, sDim, tDim>(INTERVAL, 17);
  if (name == "Qtri1") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 1);
  if (name == "Qtri3") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 2);
  if (name == "Qtri4") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 3);
  if (name == "Qtri6") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 4);
  if (name == "Qtri7") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 5);
  if (name == "Qtri12") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 6);
  if (name == "Qtri13") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 7);
  if (name == "Qtri16") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 8);
  if (name == "Qtri25") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 10);
  if (name == "Qtri36") return GetQuadratureT<T, sDim, tDim>(TRIANGLE, 12);
  if (name == "Qquad1") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 1);
  if (name == "Qquad4") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 3);
  if (name == "Qquad9") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 5);
  if (name == "Qquad16") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 7);
  if (name == "Qquad25") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 9);
  if (name == "Qquad36") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 11);
  if (name == "Qquad49") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 13);
  if (name == "Qquad64") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 15);
  if (name == "Qquad81") return GetQuadratureT<T, sDim, tDim>(QUADRILATERAL, 17);
  if (name == "Qtet1") return GetQuadratureT<T, sDim, tDim>(TETRAHEDRON, 1);
  if (name == "Qtet4") return GetQuadratureT<T, sDim, tDim>(TETRAHEDRON, 2);
  if (name == "Qtet11") return GetQuadratureT<T, sDim, tDim>(TETRAHEDRON, 4);
  if (name == "Qtet15") return GetQuadratureT<T, sDim, tDim>(TETRAHEDRON, 5);
  if (name == "Qpri6") return GetQuadratureT<T, sDim, tDim>(PRISM, 2);
  if (name == "Qpri8") return GetQuadratureT<T, sDim, tDim>(PRISM, 3);
  if (name == "Qhex1") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 1);
  if (name == "Qhex8") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 3);
  if (name == "Qhex27") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 5);
  if (name == "Qhex64") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 7);
  if (name == "Qhex125") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 9);
  if (name == "Qhex216") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 11);
  if (name == "Qhex343") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 13);
  if (name == "Qhex512") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 15);
  if (name == "Qhex729") return GetQuadratureT<T, sDim, tDim>(HEXAHEDRON, 17);
  THROW("Quadrature not implemented")
}

const Quadrature &GetQuadrature(const std::string &name) {
  return GetQuadrature<double, SpaceDimension, TimeDimension>(name);
}

void clearQuadrature() {
  for (auto quad = Quads.begin(); quad != Quads.end(); ++quad) {
    delete quad->second;
  }
  Quads.clear();

#ifdef BUILD_IA
  for (auto quad = IAQuads.begin(); quad != IAQuads.end(); ++quad)
    delete quad->second;
  IAQuads.clear();
#endif
}
