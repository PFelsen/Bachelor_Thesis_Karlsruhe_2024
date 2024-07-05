#ifndef VECTORACCESS_HPP
#define VECTORACCESS_HPP

#include "MixedDoF.hpp"
#include "ShapeIterator.hpp"
#include "Vector.hpp"

///-------------------------------------------------------------------------------------------------
/// Several basic classes for Vector iterators are defined below
///-------------------------------------------------------------------------------------------------
template<typename ITERATOR>
class VectorAccessBaseT {
  double &operator()(ITERATOR &i) { return i.VectorEntry(); }

  double &operator()(ITERATOR &i, int k) { return i.VectorEntry(k); }

  double &VectorAccess(ITERATOR &i) { return i.VectorEntry(); }

  double &VectorAccess(ITERATOR &i, int k) { return i.VectorEntry(k); }
};

struct VectorAccessDataGeneral {
  double *vectorEntry{};
  short accumulatedDoF{};
};

struct VectorAccessDataEquallyDistributed {
  double *vectorEntry{};
};

///-------------------------------------------------------------------------------------------------
/// Shape-Vector iterators for different types of DoFs are defined below
///-------------------------------------------------------------------------------------------------
struct ShapeVectorIteratorGeneral : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;
private:
  double *entry = nullptr;

  int numberOfComponents{};
  VectorAccessDataGeneral *data = nullptr;

  void Incr();
public:
  ShapeVectorIteratorGeneral(int id);

  ShapeVectorIteratorGeneral(int numberOfComponents, VectorAccessDataGeneral *data);

  double &VectorEntry() { return *entry; }

  double &VectorEntry(int k) { return entry[k]; }

  int NumberOfComponents() const { return numberOfComponents; }
};

struct ShapeVectorIteratorEquallyDistributed : public ShapeId {
  template<typename DATA>
  friend struct ShapeIterator;
private:
  int numberOfComponents{};
  VectorAccessDataEquallyDistributed *data = nullptr;

  void Incr() {
    id++;
    data++;
  }
public:
  ShapeVectorIteratorEquallyDistributed(int id);

  ShapeVectorIteratorEquallyDistributed(int numberOfComponents,
                                        VectorAccessDataEquallyDistributed *data);

  double &VectorEntry() { return *data->vectorEntry; }

  double &VectorEntry(int k) { return (data->vectorEntry)[k]; }

  int NumberOfComponents() const { return numberOfComponents; }
};

///-------------------------------------------------------------------------------------------------
/// Classes for providing Shape-Vector iterators for different types of DoFs are defined below
///-------------------------------------------------------------------------------------------------

/**
 * Provides the iteration over the shape id and the access to a vector of a general discretization
 * which is based on a single DoF, i.e., DoFs contains only one DoF. This class computes the
 * the storage id's using the matrix graph. Thus, in principle it can be used for every
 * discretization (with a single DoF), however in some cases it more efficient to use a different
 * version. For instance, if the number of components is equally distributed over the nodal points
 * (i.e., the storage points) like for the LagrangeDiscretization one better uses the
 * class VectorAccessEquallyDistributed.
 */
template<typename ELEMENT>
class VectorAccessGeneral : public VectorAccessBaseT<ShapeVectorIteratorGeneral> {
  using Iter = ShapeVectorIteratorGeneral;

  int numberOfComponents{};
  std::vector<VectorAccessDataGeneral> data{};
public:
  VectorAccessGeneral(const ELEMENT &baseElement, Vector &vector) {
    numberOfComponents = baseElement.GetVectorMatrixBase().NumberOfComponents();
    short accuDofs = 0;
    for (int i = 0; i < baseElement.size(); ++i) {
      accuDofs += baseElement.n(i);
      data.push_back(VectorAccessDataGeneral{vector(baseElement[i]), accuDofs});
    }
  }

  VectorAccessGeneral &IterateOverShapeAndVector() { return *this; }

  auto begin() { return ShapeIterator(Iter(numberOfComponents, data.data())); }

  auto end() { return ShapeIterator(Iter(data.back().accumulatedDoF)); }
};

/**
 * Use this version if the number of components is equally distributed over the nodal points
 * (i.e., the storage points) like the LagrangeDiscretization. Note that in this case the
 * shape functions (in some way) are connected to the nodal points.
 */
template<typename ELEMENT>
class VectorAccessEquallyDistributed :
    public VectorAccessBaseT<ShapeVectorIteratorEquallyDistributed> {
  using Iter = ShapeVectorIteratorEquallyDistributed;

  int numberOfComponents{};
  std::vector<VectorAccessDataEquallyDistributed> data{};
public:
  VectorAccessEquallyDistributed(const ELEMENT &baseElement, Vector &vector) {
    numberOfComponents = baseElement.GetVectorMatrixBase().NumberOfComponents(0);
    for (int i = 0; i < baseElement.size(); ++i)
      data.push_back(VectorAccessDataEquallyDistributed{vector(baseElement[i])});
  }

  VectorAccessEquallyDistributed &IterateOverShapeAndVector() { return *this; }

  auto begin() { return ShapeIterator(Iter(numberOfComponents, data.data())); }

  auto end() { return ShapeIterator(Iter(data.size())); }
};

/**
 * Provides the iteration over the shape id and the access to a vector of a general discretization
 * which possibly is based on more than a single DoF. If the discretization contains only one
 * single DoF it is more efficient to use the class VectorAccessGeneral. This class computes the
 * the storage id's using the matrix graph. Thus, in principle it can be used for every
 * discretization (with more than one DoF), however in some cases it more efficient to use a
 * different version. For instance, if the number of components is equally distributed over the
 * nodal points (i.e., the storage points) for ALL DoF contained in the discretization (e.g.
 * TaylorHoodDiscretization) one better uses the class MixedVectorAccessEquallyDistributed.
 */
template<typename ELEMENT>
class MixedVectorAccessGeneral : public VectorAccessBaseT<ShapeVectorIteratorGeneral> {
  template<typename T>
  friend struct IteratorProvider;

  using Iter = ShapeVectorIteratorGeneral;

  std::vector<int> numberOfComponents{};
  std::vector<std::vector<VectorAccessDataGeneral>> data{};
public:
  MixedVectorAccessGeneral(const ELEMENT &mixedElement, Vector &vector) {
    const std::vector<std::vector<SPData>> &nodalPointData = mixedElement.NodalPointData();
    for (int n = 0; n < mixedElement.GetVectorMatrixBase().NumberOfDoFs(); ++n) {
      numberOfComponents.push_back(mixedElement.GetVectorMatrixBase().NumberOfComponents(n));
      std::vector<short> accumulatedDoFs_n =
          mixedElement.GetDoF().AccumulatedDoFs(n, mixedElement.GetCell());
      std::vector<VectorAccessDataGeneral> data_n{};
      for (int i = 0; i < nodalPointData[n].size(); ++i)
        data_n.push_back(VectorAccessDataGeneral{vector(mixedElement[nodalPointData[n][i].index])
                                                     + nodalPointData[n][i].shift,
                                                 accumulatedDoFs_n[i]});
      data.push_back(data_n);
    }
  }

  auto IterateOverShapeAndVector(int n) { return IteratorProvider{*this, n}; }

  auto begin(int n) { return ShapeIterator(Iter(numberOfComponents[n], data[n].data())); }

  auto end(int n) { return ShapeIterator(Iter(data[n].back().accumulatedDoF)); }
};

/**
 * Use this version if the number of components is equally distributed over the nodal points
 * (i.e., the storage points) for ALL DoF contained in the discretization (e.g.
 * TaylorHoodDiscretization). Note that in this case the shape functions (in some way) are
 * connected to the nodal points. However, a DoF contained in the discretization not
 * necessarily needs to have a shape function corresponding to each nodal point (cf. TaylorHood
 * P2-P1 where only at the vertices all shape functions have a nodal point)
 */
template<typename ELEMENT>
class MixedVectorAccessEquallyDistributed :
    public VectorAccessBaseT<ShapeVectorIteratorEquallyDistributed> {
  template<typename T>
  friend struct IteratorProvider;

  using Iter = ShapeVectorIteratorEquallyDistributed;

  std::vector<int> numberOfComponents{};
  std::vector<std::vector<VectorAccessDataEquallyDistributed>> data{};
public:
  MixedVectorAccessEquallyDistributed(const ELEMENT &mixedElement, Vector &vector) {
    const std::vector<std::vector<SPData>> &nodalPointData = mixedElement.NodalPointData();
    for (int n = 0; n < mixedElement.GetVectorMatrixBase().NumberOfDoFs(); ++n) {
      numberOfComponents.push_back(mixedElement.GetDoF().NumberOfComponents(n));
      std::vector<VectorAccessDataEquallyDistributed> data_n{};
      for (int i = 0; i < nodalPointData[n].size(); ++i)
        data_n.push_back(VectorAccessDataEquallyDistributed{
            vector(mixedElement[nodalPointData[n][i].index]) + nodalPointData[n][i].shift});
      data.push_back(data_n);
    }
  }

  auto IterateOverShapeAndVector(int n) { return IteratorProvider{*this, n}; }

  auto begin(int n) { return ShapeIterator(Iter(numberOfComponents[n], data[n].data())); }

  auto end(int n) { return ShapeIterator(Iter(data[n].size())); }
};

///-------------------------------------------------------------------------------------------------
/// Old Vector access classes are defined below
///-------------------------------------------------------------------------------------------------
class RowValues {
  vector<Scalar *> a;
public:
  RowValues() : a(0) {}

  RowValues(Vector &u, const rows &R) : a(R.size()) {
    for (int i = 0; i < R.size(); ++i)
      a[i] = u(R[i]);
  }

  RowValues(Vector &u, const Cell &c) : RowValues(u, rows(u.GetMatrixGraph(), c)) {}

  Scalar &operator()(int i, int k = 0) { return a[i][k]; }
};

class RowValuesVector {
  std::vector<RowValues> rowvalues{};
public:
  RowValuesVector(Vectors &u, const rows &R) {
    for (int n = 0; n < u.size(); ++n)
      rowvalues.emplace_back(RowValues(u[n], R));
  }

  RowValuesVector(Vectors &u, const Cell &c) : RowValuesVector(u, rows(u.GetMatrixGraph(), c)) {}

  RowValues &operator[](int n) { return rowvalues[n]; }
};

class RowBndValues {
  BFParts BF;
  vector<Scalar *> a;
  vector<bool *> b;
public:
  RowBndValues(Vector &u, const Cell &c);

  Scalar &operator()(int i, int k = 0) { return a[i][k]; }

  bool &D(int i, int k = 0) { return b[i][k]; }

  bool onBnd() const { return BF.onBnd(); }

  int bc(int i) const { return BF[i]; }
};

/**
 * Sets the RowValues corresponding to the nodal points
 * n: index of the shapes (ordered as in the discretization/DoFs)
 * i: index of the shape function (corresponding to n)
 * k: index to number the dimension (if vector valued solutions are needed)
 */
class MixedRowValues {
  vector<vector<Scalar *>> a;
public:
  MixedRowValues(Vector &u, const Cell &c, const rows &R);

  MixedRowValues(Vector &u, const Cell &c) : MixedRowValues(u, c, rows(u.GetMatrixGraph(), c)) {}

  Scalar &operator()(int n, int i, int k = 0) { return a[n][i][k]; }
};

class MixedRowValuesVector {
  std::vector<MixedRowValues> rowvalues{};
public:
  MixedRowValuesVector(Vectors &u, const Cell &c, const rows &R) {
    for (int n = 0; n < u.size(); ++n)
      rowvalues.emplace_back(MixedRowValues(u[n], c, R));
  }

  MixedRowValuesVector(Vectors &u, const Cell &c) :
      MixedRowValuesVector(u, c, rows(u.GetMatrixGraph(), c)) {}

  MixedRowValues &operator[](int n) { return rowvalues[n]; }
};

/**
 * Sets the RowValues corresponding to the nodal points
 * n: index of the shapes (ordered as in the discretization/DoFs)
 * i: index of the shape function (corresponding to n)
 * k: index to number the dimension (if vector valued solutions are needed)
 */
class MixedRowBndValues {
  BFParts BF;
  vector<vector<Scalar *>> a;
  vector<vector<bool *>> b;
public:
  MixedRowBndValues(Vector &u, const Cell &c, const rows &R);

  MixedRowBndValues(Vector &u, const Cell &c) :
      MixedRowBndValues(u, c, rows(u.GetMatrixGraph(), c)) {}

  Scalar &operator()(int n, int i, int k = 0) { return a[n][i][k]; }

  bool &D(int n, int i, int k = 0) { return b[n][i][k]; }

  bool onBnd() const { return BF.onBnd(); }

  int bc(int i) const { return BF[i]; }
};

class DGRowValues {
  Scalar *a{nullptr};
  int elemSize{1};
public:
  DGRowValues() {}

  DGRowValues(Vector &u, const row &r, int elemN);

  DGRowValues(Vector &u, const Cell &c, int elemN);

  Scalar &operator()(int i, int l);
};

class DGSizedRowBndValues {
  BFParts BF;
  Scalar *a;
  bool *b;
  int elemSize;
public:
  DGSizedRowBndValues(Vector &u, const Cell &c, int elemN);

  Scalar &operator()(int i, int l);

  bool &D(int i, int l = 0);

  bool onBnd() const;

  int bc(int i) const;
};

class PartRowBndValues {
  BFParts BF;
  vector<Scalar *> a;
  vector<bool *> b;
public:
  PartRowBndValues(Vector &u, const rows &R, const Cell &c, int mpIndex = 0);

  Scalar &operator()(int i, int k = 0) { return a[i][k]; }

  bool &D(int i, int k = 0) { return b[i][k]; }

  bool onBnd() const { return BF.onBnd(); }

  int bc(int i) const { return BF[i]; }
};

#endif // VECTORACCESS_HPP
