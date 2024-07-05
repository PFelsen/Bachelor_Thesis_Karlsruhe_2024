#include "VectorAccess.hpp"
#include "MultiPartDoF.hpp"

void ShapeVectorIteratorGeneral::Incr() {
  id++;
  if (id == data->accumulatedDoF) {
    data++;
    entry = data->vectorEntry;
  } else {
    entry += numberOfComponents;
  }
}

ShapeVectorIteratorGeneral::ShapeVectorIteratorGeneral(int id) : ShapeId(id) {}

ShapeVectorIteratorGeneral::ShapeVectorIteratorGeneral(int numberOfComponents,
                                                       VectorAccessDataGeneral *data) :
    ShapeId(0), numberOfComponents(numberOfComponents), entry(data[0].vectorEntry),
    data(std::move(data)) {}

ShapeVectorIteratorEquallyDistributed::ShapeVectorIteratorEquallyDistributed(int id) :
    ShapeId(id) {}

ShapeVectorIteratorEquallyDistributed::ShapeVectorIteratorEquallyDistributed(
    int numberOfComponents, VectorAccessDataEquallyDistributed *data) :
    ShapeId(0), numberOfComponents(numberOfComponents), data(std::move(data)) {}

RowBndValues::RowBndValues(Vector &u, const Cell &c) : BF(u.GetMesh(), c), a(0), b(0) {
  rows R(u.GetMatrixGraph(), c);
  a.resize(R.size());
  b.resize(R.size());
  for (int i = 0; i < R.size(); ++i) {
    a[i] = u(R[i]);
    b[i] = u.D(R[i]);
  }
}

MixedRowValues::MixedRowValues(Vector &u, const Cell &c, const rows &R) : a(u.NumberOfDoFs()) {
  const MixedDoF &mixedDoF = MixedDoF::Cast(u.GetDoF());
  for (int n = 0; n < u.NumberOfDoFs(); ++n) {
    const vector<SPData> &data_n = mixedDoF.StoragePointData(n, c);
    a[n].resize(data_n.size());
    for (int i = 0; i < a[n].size(); ++i)
      a[n][i] = u(R[data_n[i].index]) + data_n[i].shift;
  }
}

MixedRowBndValues::MixedRowBndValues(Vector &u, const Cell &c, const rows &R) :
    BF(u.GetMesh(), c), a(u.NumberOfDoFs()), b(u.NumberOfDoFs()) {
  const MixedDoF &mixedDoF = MixedDoF::Cast(u.GetDoF());
  for (int n = 0; n < u.NumberOfDoFs(); ++n) {
    const vector<SPData> &data_n = mixedDoF.StoragePointData(n, c);
    a[n].resize(data_n.size());
    b[n].resize(data_n.size());

    for (int i = 0; i < a[n].size(); ++i) {
      a[n][i] = u(R[data_n[i].index]) + data_n[i].shift;
      b[n][i] = u.D(R[data_n[i].index]) + data_n[i].shift;
    }
  }
}

DGSizedRowBndValues::DGSizedRowBndValues(Vector &u, const Cell &c, int elemN) :
    BF(u.GetMesh(), c), elemSize(elemN) {
  row r = u.find_row(c());
  a = u(r);
  b = u.D(r);
}

Scalar &DGSizedRowBndValues::operator()(int i, int l) { return (a[i + l * elemSize]); }

bool &DGSizedRowBndValues::D(int i, int l) { return b[i + l * elemSize]; }

bool DGSizedRowBndValues::onBnd() const { return BF.onBnd(); }

int DGSizedRowBndValues::bc(int i) const { return BF[i]; }

DGRowValues::DGRowValues(Vector &u, const row &r, int elemN) : a(u(r)), elemSize(elemN) {}

DGRowValues::DGRowValues(Vector &u, const Cell &c, int elemN) :
    DGRowValues(u, u.find_row(c()), elemN) {}

Scalar &DGRowValues::operator()(int i, int l) { return (a[i + l * elemSize]); }

// RowFaceValues::RowFaceValues(int face, Vector &u, const Cell &c, const rows &R, int n) {
//   const vector<NPData> &data_face_n = u.NodalPointFaceData(n, c, face);
//   a.resize(data_face_n.size());
//   for (int i = 0; i < a.size(); ++i)
//     a[i] = u(R[data_face_n[i].index]) + data_face_n[i].shift;
// }
//
// RowBndFaceValues::RowBndFaceValues(int face, Vector &u, const Cell &c, const rows &R, int n) {
//   const vector<NPData> &data_face_n = u.NodalPointFaceData(n, c, face);
//   a.resize(data_face_n.size());
//   b.resize(data_face_n.size());
//   for (int i = 0; i < a.size(); ++i) {
//     a[i] = u(R[data_face_n[i].index]) + data_face_n[i].shift;
//     b[i] = u.D(R[data_face_n[i].index]) + data_face_n[i].shift;
//   }
// }

PartRowBndValues::PartRowBndValues(Vector &u, const rows &R, const Cell &c, int mpIndex) :
    BF(u.GetMesh(), c), a(R.size()), b(R.size()) {
  for (unsigned int i = 0; i < R.size(); ++i) {
    int q = 0;
    if (R[i].n() > 1) q = mpIndex;
    a[i] = u(R[i]) + q;
    b[i] = u.D(R[i]) + q;
  }
}
