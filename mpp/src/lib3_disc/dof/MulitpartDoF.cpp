#include "MultiPartDoF.hpp"

void MultiPartDoF::AddInterfacePoint(Point p) { interfacePoints.emplace_back(p); }

void MultiPartDoF::FinishInterfacePoints() {
  // Removing duplicate points by first sorting, ...
  std::sort(interfacePoints.begin(), interfacePoints.end());
  // ... placing duplicate elements at the end of the vector ...
  auto last = std::unique(interfacePoints.begin(), interfacePoints.end());
  // ... and then removing the end items.
  interfacePoints.erase(last, interfacePoints.end());
}

std::vector<short> MultiPartDoF::DoFSizesAtNodalPoints(const Cell &c) const {
  std::vector<Point> sp = GetNodalPoints(c);
  std::vector<short> allocSizes(sp.size());
  for (unsigned int i = 0; i < sp.size(); ++i) {
    if (std::find(interfacePoints.begin(), interfacePoints.end(), sp[i]) != interfacePoints.end())
      allocSizes[i] = 2 * short(m);
    else allocSizes[i] = short(m);
  }
  return allocSizes;
}