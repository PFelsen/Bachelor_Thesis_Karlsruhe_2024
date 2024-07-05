#ifndef MATRIXACCESS_HPP
#define MATRIXACCESS_HPP

#include "Matrix.hpp"
#include "ShapeIterator.hpp"

///-------------------------------------------------------------------------------------------------
/// Several basic classes for Matrix iterators are defined below
///-------------------------------------------------------------------------------------------------
template<typename ROW_ITER, typename ENTRY_ITER>
class MatrixAccessBaseT {
  double &operator()(ROW_ITER &i, ENTRY_ITER &j) { return j.MatrixEntry(); }

  double &operator()(ROW_ITER &i, ENTRY_ITER &j, int k, int l) { return j.MatrixEntry(k, l); }

  double &MatrixAccess(ROW_ITER &i, ENTRY_ITER &j) { return j.MatrixEntry(); }

  double &MatrixAccess(ROW_ITER &i, ENTRY_ITER &j, int k, int l) { return j.MatrixEntry(k, l); }
};

template<typename ELEMENT, typename ROW_ITER, typename ENTRY_ITER>
class MatrixAccessT : public MatrixAccessBaseT<ROW_ITER, ENTRY_ITER> {
protected:
  std::vector<std::vector<double *>> matrixEntries{};
  int numberOfComponents = 1;
public:
  MatrixAccessT(const ELEMENT &baseElement, Matrix &matrix) {
    numberOfComponents = baseElement.GetVectorMatrixBase().NumberOfComponents();
    matrixEntries = std::vector<std::vector<double *>>{};
    for (int i = 0; i < baseElement.size(); ++i) {
      std::vector<double *> row_i{};
      for (int j = 0; j < baseElement.size(); ++j)
        row_i.push_back(matrix(baseElement[i], baseElement[j]));
      matrixEntries.push_back(row_i);
    }
  }
};

template<typename ELEMENT, typename ROW_ITER, typename ENTRY_ITER>
class MixedMatrixAccessT : public MatrixAccessBaseT<ROW_ITER, ENTRY_ITER> {
protected:
  std::vector<std::vector<std::vector<std::vector<double *>>>> matrixEntries{};
  std::vector<std::vector<int>> submatrixLength{};
  std::vector<int> numberOfComponents{};
public:
  MixedMatrixAccessT(const ELEMENT &mixedElement, Matrix &matrix) {
    const std::vector<std::vector<SPData>> &nodalPointData = mixedElement.NodalPointData();
    for (int n = 0; n < mixedElement.GetVectorMatrixBase().NumberOfDoFs(); ++n) {
      numberOfComponents.push_back(mixedElement.GetDoF().NumberOfComponents(n));
      std::vector<int> submatrixLength_n{};
      std::vector<std::vector<std::vector<double *>>> matrixEntries_n{};
      for (int i = 0; i < nodalPointData[n].size(); ++i) {
        submatrixLength_n.push_back(mixedElement[nodalPointData[n][i].index].n());
        std::vector<std::vector<double *>> matrixEntries_n_i{};
        for (int m = 0; m < mixedElement.GetVectorMatrixBase().NumberOfDoFs(); ++m) {
          std::vector<double *> matrixEntries_n_i_m{};
          for (int j = 0; j < nodalPointData[m].size(); ++j) {
            matrixEntries_n_i_m.push_back(matrix(mixedElement[nodalPointData[n][i].index],
                                                 mixedElement[nodalPointData[m][j].index])
                                          + nodalPointData[n][i].shift
                                                * mixedElement[nodalPointData[m][j].index].n()
                                          + nodalPointData[m][j].shift);
          }
          matrixEntries_n_i.push_back(matrixEntries_n_i_m);
        }
        matrixEntries_n.push_back(matrixEntries_n_i);
      }
      submatrixLength.push_back(submatrixLength_n);
      matrixEntries.push_back(matrixEntries_n);
    }
  }
};

///-------------------------------------------------------------------------------------------------
/// Shape-Matrix iterators for different types of DoFs are defined below
///-------------------------------------------------------------------------------------------------

///-------------------------------------------------------------------------------------------------
/// Classes for providing Shape-Matrix iterators for different types of DoFs are
/// defined below
///-------------------------------------------------------------------------------------------------

/**
 * Provides the iteration over the shape id's and the access to a matrix of a
 * general discretization which is based on a single DoF, i.e., DoFs contains
 * only one DoF. This class computes the the storage id's using the matrix
 * graph. Thus, in principle it can be used for every discretization (with a
 * single DoF), however in some cases it more efficient to use a different
 * version. For instance, if the number of components is equally distributed
 * over the nodal points (i.e., the storage points) like for the
 * LagrangeDiscretization one better uses the class
 * MatrixAccessEquallyDistributed.
 */
struct ShapeMatrixEntryIteratorGeneral : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;
private:
  short numberOfComponents = -1;
  const short *accuDoFs = nullptr;
  double **entries = nullptr;
  double *entry = nullptr;
  const int *submatrixLength = nullptr;

  void Incr();
public:
  ShapeMatrixEntryIteratorGeneral(int id);

  ShapeMatrixEntryIteratorGeneral(short numberOfComponents, const short *accumulatedDoFs,
                                  double **entries, const int *submatrixLength);

  double &MatrixEntry() { return *entry; }

  double &MatrixEntry(int k, int l) { return entry[k * (*submatrixLength) + l]; }

  int NumberOfComponents() const { return numberOfComponents; }
};

struct ShapeMatrixRowIteratorGeneral : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;

  using Iter = ShapeMatrixEntryIteratorGeneral;
private:
  short numberOfComponents = 1;
  const std::vector<int> *submatrixLength = nullptr;   // as reference
  const std::vector<short> *accumulatedDoFs = nullptr; // as reference
  const short *accuDoFs = nullptr;
  std::vector<double *> *rows = nullptr;
  std::vector<double *> row{};

  void Incr();
public:
  ShapeMatrixRowIteratorGeneral(int id);

  ShapeMatrixRowIteratorGeneral(short numberOfComponents, const std::vector<int> &submatrixLength,
                                const std::vector<short> &accumulatedDoFs,
                                std::vector<double *> *rows);

  ShapeMatrixRowIteratorGeneral &IterateOverShapeAndMatrixRow() { return *this; }

  auto begin() {
    return ShapeIterator(
        Iter(numberOfComponents, accumulatedDoFs->data(), row.data(), submatrixLength->data()));
  }

  auto end() { return ShapeIterator(Iter(accumulatedDoFs->back())); }

  int NumberOfComponents() const { return numberOfComponents; }
};

template<typename ELEMENT>
class MatrixAccessGeneral :
    public MatrixAccessT<ELEMENT, ShapeMatrixRowIteratorGeneral, ShapeMatrixEntryIteratorGeneral> {
  using Iter = ShapeMatrixRowIteratorGeneral;

  std::vector<short> accumulatedDoFs{};
  std::vector<int> submatrixLength{};
public:
  MatrixAccessGeneral(const ELEMENT &baseElement, Matrix &matrix) :
      MatrixAccessT<ELEMENT, ShapeMatrixRowIteratorGeneral,
                    ShapeMatrixEntryIteratorGeneral>(baseElement, matrix) {
    accumulatedDoFs = baseElement.GetVectorMatrixBase().GetDoF().AccumulatedAllocationSizes(
        baseElement.GetCell());
    for (int i = 0; i < baseElement.size(); ++i)
      submatrixLength.push_back(baseElement[i].n());
    // Add dummy vector to avoid additional if statement in increment operator
    this->matrixEntries.push_back(std::vector<double *>{});
  }

  MatrixAccessGeneral &IterateOverShapeAndMatrix() { return *this; }

  auto begin() {
    return ShapeIterator(Iter(this->numberOfComponents, submatrixLength, accumulatedDoFs,
                              this->matrixEntries.data()));
  }

  auto end() { return ShapeIterator(Iter(accumulatedDoFs.back())); }
};

/**
 * Use this version if the number of components is equally distributed over the
 * nodal points (i.e., the storage points) like the LagrangeDiscretization. Note
 * that in this case the shape functions (in some way) are connected to the
 * nodal points.
 */
struct ShapeMatrixEntryIteratorEquallyDistributed : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;
private:
  short numberOfComponents = -1;
  double **entries = nullptr;

  void Incr() {
    id++;
    entries++;
  }
public:
  ShapeMatrixEntryIteratorEquallyDistributed(int id) : ShapeId(id) {}

  ShapeMatrixEntryIteratorEquallyDistributed(short numberOfComponents, double **entries) :
      ShapeId(0), numberOfComponents(numberOfComponents), entries(std::move(entries)) {}

  double &MatrixEntry() { return **entries; }

  double &MatrixEntry(int k, int l) { return (*entries)[k * numberOfComponents + l]; }

  int NumberOfComponents() const { return numberOfComponents; }
};

struct ShapeMatrixRowIteratorEquallyDistributed : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;

  using Iter = ShapeMatrixEntryIteratorEquallyDistributed;
private:
  short numberOfComponents = 1;
  std::vector<double *> *rows = nullptr;

  void Incr() {
    id++;
    rows++;
  }
public:
  ShapeMatrixRowIteratorEquallyDistributed(int id) : ShapeId(id) {}

  ShapeMatrixRowIteratorEquallyDistributed(short numberOfComponents, std::vector<double *> *rows) :
      ShapeId(0), numberOfComponents(numberOfComponents), rows(std::move(rows)) {}

  ShapeMatrixRowIteratorEquallyDistributed &IterateOverShapeAndMatrixRow() { return *this; }

  auto begin() { return ShapeIterator(Iter(numberOfComponents, rows->data())); }

  auto end() { return ShapeIterator(Iter(rows->size())); }

  int NumberOfComponents() const { return numberOfComponents; }
};

template<typename ELEMENT>
class MatrixAccessEquallyDistributed :
    public MatrixAccessT<ELEMENT, ShapeMatrixRowIteratorEquallyDistributed,
                         ShapeMatrixEntryIteratorEquallyDistributed> {
  using Iter = ShapeMatrixRowIteratorEquallyDistributed;
public:
  MatrixAccessEquallyDistributed(const ELEMENT &baseElement, Matrix &matrix) :
      MatrixAccessT<ELEMENT, ShapeMatrixRowIteratorEquallyDistributed,
                    ShapeMatrixEntryIteratorEquallyDistributed>(baseElement, matrix) {}

  MatrixAccessEquallyDistributed &IterateOverShapeAndMatrix() { return *this; }

  auto begin() { return ShapeIterator(Iter(this->numberOfComponents, this->matrixEntries.data())); }

  auto end() { return ShapeIterator(Iter(this->matrixEntries.size())); }
};

/**
 * Provides the iteration over the shape id and the access to a vector of a
 * general discretization which possibly is based on more than a single DoF. If
 * the discretization contains only one single DoF it is more efficient to use
 * the class VectorAccessGeneral. This class computes the the storage id's using
 * the matrix graph. Thus, in principle it can be used for every discretization
 * (with more than one DoF), however in some cases it more efficient to use a
 * different version. For instance, if the number of components is equally
 * distributed over the nodal points (i.e., the storage points) for ALL DoF
 * contained in the discretization (e.g. TaylorHoodDiscretization) one better
 * uses the class MixedVectorAccessEquallyDistributed.
 */
// struct ShapeMatrixEntryIteratorGeneral : public ShapeId {
//   template<typename DATA>
//   friend struct ShapeIterator;
// private:
//   short numberOfComponents = -1;
//   const short *accuDoFs = nullptr;
//   double **entries = nullptr;
//   double *entry = nullptr;
//   const int *submatrixLength = nullptr;
//
//   void Incr();
// public:
//   ShapeMatrixEntryIteratorGeneral(int id);
//
//   ShapeMatrixEntryIteratorGeneral(short numberOfComponents, const short
//   *accumulatedDoFs,
//                                   double **entries, const int
//                                   *submatrixLength);
//
//   double &MatrixEntry() { return *entry; }
//
//   double &MatrixEntry(int k, int l) { return entry[k * (*submatrixLength) +
//   l]; }
// };
//
// struct ShapeMatrixRowIteratorGeneral : public ShapeId {
//   template<typename DATA>
//   friend struct ShapeIterator;
//
//   using Iter = ShapeMatrixEntryIteratorGeneral;
// private:
//   short numberOfComponents = 1;
//   const int *submatrixLength = nullptr;
//   const std::vector<short> *accumulatedDoFs = nullptr; // as reference
//   const short *accuDoFs = nullptr;
//   std::vector<double *> *rows = nullptr;
//   std::vector<double *> row{};
//
//   void Incr();
//
// public:
//   ShapeMatrixRowIteratorGeneral(int id);
//
//   ShapeMatrixRowIteratorGeneral(short numberOfComponents, const int
//   *submatrixLength,
//                                 const std::vector<short> &accumulatedDoFs,
//                                 std::vector<double *> *rows);
//
//   ShapeMatrixRowIteratorGeneral &IterateOverShapeAndMatrixRow() { return
//   *this; }
//
//   auto begin() {
//     return ShapeIterator(Iter(numberOfComponents, accumulatedDoFs->data(),
//     row.data(), submatrixLength));
//   }
//
//   auto end() { return ShapeIterator(Iter(accumulatedDoFs->back())); }
// };

// template<typename ELEMENT>
// class MixedMatrixAccessGeneral : public MixedMatrixAccessT<ELEMENT,
//                                                            ShapeMatrixRowIteratorGeneral,
//                                                            ShapeMatrixEntryIteratorGeneral>
//                                                            {
//   using Iter = ShapeMatrixRowIteratorGeneral;
//
//   std::vector<std::vector<short>> accumulatedDoFs{};
//   std::vector<std::vector<int>> submatrixLength{};
// public:
//   MixedMatrixAccessGeneral(const ELEMENT &mixedElement, Matrix &matrix)
//       : MixedMatrixAccessT<ELEMENT, ShapeMatrixRowIteratorGeneral,
//       ShapeMatrixEntryIteratorGeneral> (mixedElement, matrix) {
//     for (int n = 0; n < mixedElement.GetVectorMatrixBase().NumberOfDoFs();
//     ++n) {
//       std::vector<short> accumulatedDoFs_n{};
//       mixedElement.GetVectorMatrixBase().AccumulatedDoFs(n,
//       mixedElement.GetCell(), accumulatedDoFs_n);
//       accumulatedDoFs.push_back(accumulatedDoFs_n);
//       std::vector<int> submatrixLength_n{};
//       for (int i = 0; i < mixedElement.size(); ++i)
//         submatrixLength_n.push_back(mixedElement[i].n());
//       submatrixLength.push_back(submatrixLength_n);
//       //Add dummy vector to avoid additional if statement in increment
//       operator this->matrixEntries[n].push_back(std::vector<double *>{});
//     }
//   }
//
//   auto IterateOverShapeAndMatrix(int n) { return IteratorProvider{*this, n};
//   }
//
//   auto begin(int n) {
//     return ShapeIterator(Iter(this->numberOfComponents[n],
//     this->numberOfComponents, submatrixLength.data(), accumulatedDoFs,
//     this->matrixEntries.data()));
//   }
//
//   auto end(int n) { return ShapeIterator(Iter(accumulatedDoFs[n].back())); }
// };

/**
 * Use this version if the number of components is equally distributed over the
 * nodal points (i.e., the storage points) for ALL DoF contained in the
 * discretization (e.g. TaylorHoodDiscretization). Note that in this case the
 * shape functions (in some way) are connected to the nodal points. However, a
 * DoF contained in the discretization not necessarily needs to have a shape
 * function corresponding to each nodal point (cf. TaylorHood P2-P1 where only
 * at the vertices all shape functions have a nodal point)
 */
struct MixedShapeMatrixEntryIteratorEquallyDistributed : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;
private:
  int numberOfComponents = -1;
  double **entries = nullptr;
  int *submatrixLength = nullptr;

  void Incr() {
    id++;
    entries++;
    submatrixLength++;
  }
public:
  MixedShapeMatrixEntryIteratorEquallyDistributed(int id) : ShapeId(id) {}

  MixedShapeMatrixEntryIteratorEquallyDistributed(int numberOfComponents, int *submatrixLength,
                                                  double **entries) :
      ShapeId(0), numberOfComponents(numberOfComponents),
      submatrixLength(std::move(submatrixLength)), entries(std::move(entries)) {}

  double &MatrixEntry() { return **entries; }

  double &MatrixEntry(int k, int l) { return (*entries)[k * (*submatrixLength) + l]; }

  int NumberOfComponents() const { return numberOfComponents; }
};

struct MixedShapeMatrixRowIteratorEquallyDistributed : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;

  template<typename T>
  friend struct IteratorProvider;
private:
  int numberOfComponents_n = -1;
  std::vector<int> *numberOfComponents = nullptr;           // as reference
  std::vector<std::vector<int>> *submatrixLength = nullptr; // as reference
  std::vector<std::vector<double *>> *rows = nullptr;

  void Incr() {
    id++;
    rows++;
  }
public:
  MixedShapeMatrixRowIteratorEquallyDistributed(int id) : ShapeId(id) {}

  MixedShapeMatrixRowIteratorEquallyDistributed(int numberOfComponents_n,
                                                std::vector<int> &numberOfComponents,
                                                std::vector<std::vector<int>> &submatrixLength,
                                                std::vector<std::vector<double *>> *rows) :
      ShapeId(0), numberOfComponents_n(numberOfComponents_n),
      numberOfComponents(&numberOfComponents), submatrixLength(&submatrixLength),
      rows(std::move(rows)) {}

  auto IterateOverShapeAndMatrixRow(int m) { return IteratorProvider{*this, m}; }

  auto begin(int m) {
    return ShapeIterator(
        MixedShapeMatrixEntryIteratorEquallyDistributed(numberOfComponents->operator[](m),
                                                        submatrixLength->operator[](m).data(),
                                                        rows->operator[](m).data()));
  }

  auto end(int m) {
    return ShapeIterator(
        MixedShapeMatrixEntryIteratorEquallyDistributed(rows->operator[](m).size()));
  }

  int NumberOfComponents() const { return numberOfComponents_n; }
};

template<typename ELEMENT>
class MixedMatrixAccessEquallyDistributed :
    public MixedMatrixAccessT<ELEMENT, MixedShapeMatrixRowIteratorEquallyDistributed,
                              MixedShapeMatrixEntryIteratorEquallyDistributed> {
  template<typename T>
  friend struct IteratorProvider;

  using Iter = MixedShapeMatrixRowIteratorEquallyDistributed;
public:
  MixedMatrixAccessEquallyDistributed(const ELEMENT &mixedElement, Matrix &matrix) :
      MixedMatrixAccessT<ELEMENT, MixedShapeMatrixRowIteratorEquallyDistributed,
                         MixedShapeMatrixEntryIteratorEquallyDistributed>(mixedElement, matrix) {}

  auto IterateOverShapeAndMatrix(int n) { return IteratorProvider{*this, n}; }

  auto begin(int n) {
    return ShapeIterator(Iter(this->numberOfComponents[n], this->numberOfComponents,
                              this->submatrixLength, this->matrixEntries[n].data()));
  }

  auto end(int n) { return ShapeIterator(Iter(this->matrixEntries[n].size())); }
};

class RowEntries {
  vector<vector<Scalar *>> a;
  vector<int> n;
public:
  RowEntries(Matrix &A, const rows &R);

  RowEntries(Matrix &A, const Cell &c) : RowEntries(A, rows(A.GetMatrixGraph(), c)) {}

  Scalar &operator()(int i, int j, int k = 0, int l = 0) { return a[i][j][k * n[j] + l]; }
};

class MixedRowEntries {
  vector<vector<vector<vector<Scalar *>>>> a;
  vector<vector<int>> length;
public:
  MixedRowEntries(Matrix &A, const Cell &c, const rows &R);

  MixedRowEntries(Matrix &A, const Cell &c) : MixedRowEntries(A, c, rows(A.GetMatrixGraph(), c)) {}

  MixedRowEntries(Matrix &A, const Cell &c, const Cell &cf);

  Scalar &operator()(int n, int m, int i, int j, int k = 0, int l = 0) {
    return a[n][m][i][j][k * length[m][j] + l];
  }
};

class DGRowEntries {
public:
  Scalar *a;
  int n;
  int nf;
  bool prev;
  int time_deg;

  void initialize(Matrix &A, const row &r, const row &rf);

  DGRowEntries(Matrix &A, const row &r, const row &rf);

  DGRowEntries(Matrix &A, const row &r) : DGRowEntries(A, r, r) {}

  DGRowEntries(Matrix &A, const Cell &c, const Cell &cf, bool prev = false);

  DGRowEntries(Matrix &A, const Cell &c) : DGRowEntries(A, A.find_row(c())) {}

  Scalar &operator()(int k, int l);
};

class DGSizedRowEntries {
  Scalar *a;
  int n;
  int nf;
  int elemsize;
  int elemsize_nghbr;
  bool prev;
  int time_deg;
public:
  DGSizedRowEntries(Matrix &A, const row &R, int elemN, const row &Rf, int elemNf,
                    bool prev = false);

  DGSizedRowEntries(Matrix &A, const Cell &c, int elemN, const Cell &cf, int elemNf,
                    bool prev = false);

  DGSizedRowEntries(Matrix &A, const Cell &c, int elemN);

  Scalar &operator()(int i, int j, int k = 0, int l = 0);
};

class PartRowEntries {
  vector<vector<Scalar *>> a;
  vector<int> n;
  vector<int> q;
public:
  PartRowEntries(Matrix &A, const rows &R, int p);

  Scalar &operator()(int i, int j, int k = 0, int l = 0) { return a[i][j][k * n[j] + l]; }
};

#endif // MATRIXACCESS_HPP
