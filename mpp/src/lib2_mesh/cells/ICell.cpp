#include "ICell.hpp"

template<>
double Cell::LocalFaceAreaT<double>(int faceId) const {
  return LocalFaceArea(faceId);
}

template<>
const Point &Cell::LocalFaceNormalT(int faceId) const {
  return LocalFaceNormal(faceId);
}

#ifdef BUILD_IA

template<>
IAInterval Cell::LocalFaceAreaT<IAInterval>(int faceId) const {
  return LocalFaceAreaIA(faceId);
}

template<>
const IAPoint &Cell::LocalFaceNormalT(int faceId) const {
  return LocalFaceNormalIA(faceId);
}

#endif