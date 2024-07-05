#include "MatrixAccess.hpp"

void ShapeMatrixEntryIteratorGeneral::Incr() {
  id++;
  if (id == *accuDoFs) {
    accuDoFs++;
    entries++;
    entry = *entries;
    submatrixLength++;
  } else {
    entry += numberOfComponents;
  }
}

ShapeMatrixEntryIteratorGeneral::ShapeMatrixEntryIteratorGeneral(int id) : ShapeId(id) {}

ShapeMatrixEntryIteratorGeneral::ShapeMatrixEntryIteratorGeneral(short numberOfComponents,
                                                                 const short *accumulatedDoFs,
                                                                 double **entries,
                                                                 const int *submatrixLength) :
    ShapeId(0), numberOfComponents(numberOfComponents), entry(entries[0]),
    accuDoFs(std::move(accumulatedDoFs)), entries(std::move(entries)),
    submatrixLength(std::move(submatrixLength)) {}

void ShapeMatrixRowIteratorGeneral::Incr() {
  id++;
  if (id == *accuDoFs) {
    accuDoFs++;
    rows++;
    row = *rows;
  } else {
    for (int n = 0; n < row.size(); ++n)
      row[n] += numberOfComponents * submatrixLength->operator[](n);
  }
}

ShapeMatrixRowIteratorGeneral::ShapeMatrixRowIteratorGeneral(int id) : ShapeId(id) {}

ShapeMatrixRowIteratorGeneral::ShapeMatrixRowIteratorGeneral(
    short numberOfComponents, const std::vector<int> &submatrixLength,
    const std::vector<short> &accumulatedDoFs, std::vector<double *> *rows) :
    ShapeId(0), numberOfComponents(numberOfComponents), row(rows[0]),
    submatrixLength(&submatrixLength), accumulatedDoFs(&accumulatedDoFs),
    accuDoFs(&accumulatedDoFs[0]), rows(std::move(rows)) {}

RowEntries::RowEntries(Matrix &A, const rows &R) :
    a(R.size(), vector<Scalar *>(R.size())), n(R.size()) {
  for (int i = 0; i < R.size(); ++i) {
    n[i] = R[i].n();
    for (int j = 0; j < R.size(); ++j) {
      a[i][j] = A(R[i], R[j]);
    }
  }
}

MixedRowEntries::MixedRowEntries(Matrix &A, const Cell &c, const rows &R) :
    a(A.NumberOfDoFs(), vector<vector<vector<Scalar *>>>(A.NumberOfDoFs())),
    length(A.NumberOfDoFs()) {

  const MixedDoF &mixedDoF = MixedDoF::Cast(A.GetDoF());
  for (int n = 0; n < A.NumberOfDoFs(); ++n) {
    const vector<SPData> &data_n = mixedDoF.StoragePointData(n, c);
    int nps_n = data_n.size();
    length[n].resize(nps_n);
    for (int i = 0; i < nps_n; ++i)
      length[n][i] = R[data_n[i].index].n();
    for (int m = 0; m < A.NumberOfDoFs(); ++m) {
      const vector<SPData> &data_m = mixedDoF.StoragePointData(m, c);
      int nps_m = data_m.size();
      a[n][m].resize(nps_n);
      for (int i = 0; i < nps_n; ++i) {
        a[n][m][i].resize(nps_m);
        for (int j = 0; j < nps_m; ++j) {
          a[n][m][i][j] = A(R[data_n[i].index], R[data_m[j].index])
                          + data_n[i].shift * R[data_m[j].index].n() + data_m[j].shift;
        }
      }
    }
  }
}

MixedRowEntries::MixedRowEntries(Matrix &A, const Cell &c, const Cell &cf) :
    a(A.NumberOfDoFs(), vector<vector<vector<Scalar *>>>(A.NumberOfDoFs())),
    length(A.NumberOfDoFs()) {
  rows R(A.GetMatrixGraph(), c);
  rows R2(A.GetMatrixGraph(), cf);

  const MixedDoF &mixedDoF = MixedDoF::Cast(A.GetDoF());
  for (int n = 0; n < A.NumberOfDoFs(); ++n) {
    const vector<SPData> &data_n = mixedDoF.StoragePointData(n, c);
    int nps_n = data_n.size();
    length[n].resize(nps_n);
    for (int i = 0; i < nps_n; ++i)
      length[n][i] = R[data_n[i].index].n();
    for (int m = 0; m < A.NumberOfDoFs(); ++m) {
      const vector<SPData> &data_m = mixedDoF.StoragePointData(m, cf);
      int nps_m = data_m.size();
      a[n][m].resize(nps_n);
      for (int i = 0; i < nps_n; ++i) {
        a[n][m][i].resize(nps_m);
        for (int j = 0; j < nps_m; ++j) {
          a[n][m][i][j] = A(R[data_n[i].index], R2[data_m[j].index])
                          + data_n[i].shift * R2[data_m[j].index].n() + data_m[j].shift;
        }
      }
    }
  }
}

void DGRowEntries::initialize(Matrix &A, const row &r, const row &rf) {
  n = r.n();
  nf = rf.n();
  a = A(r, rf);
}

DGRowEntries::DGRowEntries(Matrix &A, const row &r, const row &rf) { initialize(A, r, rf); }

DGRowEntries::DGRowEntries(Matrix &A, const Cell &c, const Cell &cf, bool prev) : prev(prev) {
  time_deg = prev ? A.GetDoF().get_time_deg(cf) : 0;
  initialize(A, A.find_row(c()), A.find_row(cf()));
}

Scalar &DGRowEntries::operator()(int k, int l) {
  if (k * nf + l >= n * nf) {
    mout << OUT(n) << OUT(nf) << OUT(k) << OUT(l) << endl;
    Exit("too large")
  }
  if (prev) return a[k * nf + l + (time_deg - 1) * (nf / time_deg)];
  return a[k * nf + l];
}

DGSizedRowEntries::DGSizedRowEntries(Matrix &A, const row &R, int elemN, const row &Rf, int elemNf,
                                     bool prev) :
    n(R.n()), nf(Rf.n()), elemsize(elemN), elemsize_nghbr(elemNf), prev(prev),
    time_deg(prev ? A.GetDoF().get_time_deg(Rf()) : 0) {
  a = A(R, Rf);
}

DGSizedRowEntries::DGSizedRowEntries(Matrix &A, const Cell &c, int elemN, const Cell &cf,
                                     int elemNf, bool prev) :
    DGSizedRowEntries(A, A.find_row(c()), elemN, A.find_row(cf()), elemNf, prev) {}

DGSizedRowEntries::DGSizedRowEntries(Matrix &A, const Cell &c, int elemN) :
    DGSizedRowEntries(A, A.find_row(c()), elemN, A.find_row(c()), elemN) {}

Scalar &DGSizedRowEntries::operator()(int i, int j, int k, int l) {
  int k_1 = k * elemsize + i;
  int l_1 = l * elemsize_nghbr + j;
  /*if (k_1 * nf + l_1 >= n * nf) {
    mout << OUT (n) << OUT (nf) << OUT (k_1) << OUT (l_1) << endl;
    Exit ("too large")
  }
  if (prev) return a[k_1 * nf + l_1 + (time_deg - 1) * (nf / time_deg)];*/
  return a[k_1 * nf + l_1];
}

PartRowEntries::PartRowEntries(Matrix &A, const rows &R, int p) :
    a(R.size(), vector<Scalar *>(R.size())), n(R.size()), q(R.size()) {
  for (unsigned int i = 0; i < R.size(); ++i) {
    q[i] = 0;
    n[i] = R[i].n();
    if (n[i] > 1) q[i] = p;
  }
  for (unsigned int i = 0; i < R.size(); ++i) {
    for (unsigned int j = 0; j < R.size(); ++j) {
      a[i][j] = A(R[i], R[j]) + q[i] * n[j] + q[j];
    }
  }
}
