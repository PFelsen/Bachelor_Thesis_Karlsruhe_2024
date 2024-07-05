#include "IDoF.hpp"

SubVectorMask::SubVectorMask() {
  for (int i = 0; i < MaxPointTypes; ++i) {
    tp[i] = 0;
    m[i] = 0;
  }
}

SubVectorMask::SubVectorMask(const int *n, const char *c) {
  int k = 0;
  for (int i = 0; i < MaxPointTypes; ++i) {
    if (n[i] > 0) {
      m[i] = new bool[n[i]];
      for (int j = 0; j < n[i]; ++j)
        m[i][j] = (c[k++] == '1');
    } else m[i] = 0;
    tp[i] = n[i];
  }
}

SubVectorMask::SubVectorMask(const SubVectorMask &s) {
  for (int i = 0; i < MaxPointTypes; ++i) {
    int n = tp[i] = s.tp[i];
    if (n) {
      m[i] = new bool[n];
      for (int j = 0; j < n; ++j)
        m[i][j] = s.m[i][j];
    } else m[i] = 0;
  }
}

SubVectorMask &SubVectorMask::operator=(const SubVectorMask &s) {
  for (int i = 0; i < MaxPointTypes; ++i) {
    int n = tp[i] = s.tp[i];
    if (n) {
      m[i] = new bool[n];
      for (int j = 0; j < n; ++j)
        m[i][j] = s.m[i][j];
    } else m[i] = 0;
  }
  return *this;
}

SubVectorMask::~SubVectorMask() {
  for (int i = 0; i < MaxPointTypes; ++i)
    if (m[i]) delete[] m[i];
}

std::ostream &operator<<(std::ostream &s, const SubVectorMask &S) {
  for (int i = 0; i < MaxPointTypes; ++i)
    for (int j = 0; j < S.tp[i]; ++j)
      s << S.m[i][j];
  return s;
}

short IDoF::NumberOfStoragePoints(const Cell &c) const { return NumberOfNodalPoints(c); }

std::vector<Point> IDoF::GetStoragePoints(const Cell &c) const { return GetNodalPoints(c); }

std::vector<short> IDoF::AllocationSizesAtStoragePoints(const Cell &c) const {
  return DoFSizesAtNodalPoints(c);
}

std::vector<short> IDoF::AccumulatedAllocationSizes(const Cell &c) const {
  std::vector<short> accuSizes = AllocationSizesAtStoragePoints(c);
  for (int i = 1; i < accuSizes.size(); ++i) {
    accuSizes[i] = accuSizes[i - 1] + accuSizes[i];
  }
  return accuSizes;
}

std::vector<Point> IDoF::GetNodalPointsOnFace(const Cell &c, int faceId) const {
  std::vector<Point> np = GetNodalPoints(c);
  std::vector<Point> npOnFace(NumberOfNodalPointsOnFace(c, faceId));
  for (int k = 0; k < npOnFace.size(); ++k)
    npOnFace[k] = np[IdOfNodalPointOnFace(c, faceId, k)];
  return npOnFace;
}

std::vector<Point> IDoF::GetNodalPointsOnEdge(const Cell &c, int edgeId) const {
  std::vector<Point> np = GetStoragePoints(c);
  std::vector<Point> npOnEdge(NumberOfNodalPointsOnEdge(c, edgeId));
  for (int k = 0; k < npOnEdge.size(); ++k)
    npOnEdge[k] = np[IdOfNodalPointOnEdge(c, edgeId, k)];
  return npOnEdge;
}

short IDoF::NumberOfStoragePointsOnFace(const Cell &c, int faceId) const {
  return NumberOfNodalPointsOnFace(c, faceId);
}

short IDoF::IdOfStoragePointOnFace(const Cell &c, int faceId, int k) const {
  return IdOfNodalPointOnFace(c, faceId, k);
}

std::vector<Point> IDoF::GetStoragePointsOnFace(const Cell &c, int faceId) const {
  std::vector<Point> sp = GetStoragePoints(c);
  std::vector<Point> spOnFace(NumberOfStoragePointsOnFace(c, faceId));
  for (int k = 0; k < spOnFace.size(); ++k)
    spOnFace[k] = sp[IdOfStoragePointOnFace(c, faceId, k)];
  return spOnFace;
}

short IDoF::NumberOfStoragePointsOnEdge(const Cell &c, int edgeId) const {
  return NumberOfNodalPointsOnEdge(c, edgeId);
}

short IDoF::IdOfStoragePointOnEdge(const Cell &c, int edgeId, int k) const {
  return IdOfNodalPointOnEdge(c, edgeId, k);
}

std::vector<Point> IDoF::GetStoragePointsOnEdge(const Cell &c, int edgeId) const {
  std::vector<Point> sp = GetStoragePoints(c);
  std::vector<Point> spOnEdge(NumberOfStoragePointsOnEdge(c, edgeId));
  for (int k = 0; k < spOnEdge.size(); ++k)
    spOnEdge[k] = sp[IdOfStoragePointOnEdge(c, edgeId, k)];
  return spOnEdge;
}

// void IDoF::AccumulatedDoFs(const Cell &c, vector<short> &accumulatedNodalDoFs) const {
//   std::vector<short> nodalDoFs;
//   NodalDoFs(c, nodalDoFs);
//   accumulatedNodalDoFs.resize(nodalDoFs.size());
//   accumulatedNodalDoFs[0] = nodalDoFs[0];
//   for (int i = 1; i < nodalDoFs.size(); ++i) {
//     accumulatedNodalDoFs[i] = accumulatedNodalDoFs[i - 1] + nodalDoFs[i];
//   }
// }

// void IDoF::GetNodalPointsOnFace(const Cell &c, int face, vector<Point> &z) const {
//   vector<Point> w;
//   GetNodalPoints(c, w);
//   z.resize(GetNodalPointsOnFace(c, face));
//   for (int j = 0; j < z.size(); ++j)
//     z[j] = w[NodalPointOnFace(c, face, j)];
// }

// void IDoF::NodalDoFsOnFace(const Cell &c, int face, vector<short> &z) const {
//   vector<short> w;
//   NodalDoFs(c, w);
//   z.resize(GetNodalPointsOnFace(c, face));
//   for (int j = 0; j < z.size(); ++j)
//     z[j] = w[NodalPointOnFace(c, face, j)];
// }

void IDoF::TypeDoFs(int *n) const {
  for (int i = 0; i < MaxPointTypes; ++i)
    n[i] = TypeDoF(i);
}

const SubVectorMask &IDoF::GetSubVector(const char *name) const {
  std::map<string, SubVectorMask>::const_iterator s = sub.find(name);
  if (s == sub.end()) Exit(string("wrong subvector name ") + string(name));
  return s->second;
}

void IDoF::AddSubVector(const char *name, const char *c, const int *n) {
  sub[name] = SubVectorMask(n, c);
}

void IDoF::AddSubVector(const char *name, const char *c) {
  int n[MaxPointTypes];
  TypeDoFs(n);
  sub[name] = SubVectorMask(n, c);
}

template<typename S>
LogTextStream<S> &operator<<(LogTextStream<S> &s, const IDoF &D) {
  s << D.Name() << ": ";
  int n[MaxPointTypes];
  D.TypeDoFs(n);
  for (int i = 0; i < MaxPointTypes; ++i)
    s << n[i];
  return s << endl << D.Sub();
}