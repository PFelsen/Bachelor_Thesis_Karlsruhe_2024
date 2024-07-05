#ifndef ELEMENTDATA_HPP
#define ELEMENTDATA_HPP

#include "Algebra.hpp"
#include "IDiscretization.hpp"

#ifdef BUILD_IA

#include "IAInterval.hpp"

#endif


template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
struct ElementDataT {
  BFParts bfParts;
  const QuadratureT<TT, sDim, tDim> &quad;
  vector<TT> quadWeights{};
  vector<PointT<TT, sDim, tDim>> quadPoints{};
#ifdef AFFINE_TRANSFORMATION
  TransformationT<TT, sDim> trafo;
#else
  vector<TransformationT<TT, sDim>> trafos{};
#endif

  ElementDataT(const Cell &c, const QuadratureT<TT, sDim, tDim> &quad, const Mesh &mesh) :
      ElementDataT(quad, mesh, c) {

#ifdef AFFINE_TRANSFORMATION
    trafo = c.GetTransformation(quad.QPoint(0));
    for (int q = 0; q < quad.size(); ++q) {
      quadWeights.push_back(trafo.Det() * quad.Weight(q));
      quadPoints.push_back(c.LocalToGlobal(quad.QPoint(q)));
    }
#else
    trafos.reserve(quad.size());
    for (int q = 0; q < quad.size(); ++q) {
      trafos.push_back(c.GetTransformation(quad.QPoint(q)));
      quadWeights.push_back(trafos[q].Det() * quad.Weight(q));
      quadPoints.push_back(c.LocalToGlobal(quad.QPoint(q)));
    }
#endif
  }

  const TransformationT<TT, sDim> &GetTransformation(int q) const {
#ifdef AFFINE_TRANSFORMATION
    return trafo;
#else
    return trafos[q];
#endif
  }

  virtual TT Area() const {
    TT area{0.0};
    for (auto weight : this->quadWeights)
      area += weight;
    return area;
  }

  virtual const PointT<TT, sDim, tDim> &QLocal(int q) const {
    THROW("QLocal is not defined for volume elements");
  }

  virtual const PointT<TT, sDim, tDim> &QNormal(int q) const {
    THROW("QNormal is not defined for volume elements");
  }

  virtual short Sign() const { THROW("Sign is not defined for volume elements"); }

  virtual int FaceID() const { THROW("FaceID is not defined for volume elements"); }
protected:
  ElementDataT(const QuadratureT<TT, sDim, tDim> &quad, const Mesh &mesh, const Cell &c) :
      quad(quad), bfParts(mesh, c) {
    quadWeights.reserve(quad.size());
    quadPoints.reserve(quad.size());
  }
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
struct FaceElementDataT : public ElementDataT<TT, sDim, tDim> {
  vector<PointT<TT, sDim, tDim>> quadNormals{};
  vector<PointT<TT, sDim, tDim>> quadPointCellLocals{};
  short sign{1};
  int faceID;

  TT area{0};

  FaceElementDataT(const Cell &c, int face, const QuadratureT<TT, sDim, tDim> &quad,
                   const Mesh &mesh) : ElementDataT<TT, sDim, tDim>(quad, mesh, c), faceID(face) {
    this->quadNormals.reserve(quad.size());
    this->quadPointCellLocals.reserve(quad.size());

    TT refArea{};
    for (int q = 0; q < quad.size(); ++q)
      refArea += quad.Weight(q);
    refArea = c.LocalFaceAreaT<TT>(faceID) / refArea;

#ifdef AFFINE_TRANSFORMATION
    auto qCellLocal = c.FaceLocalToLocalOriented(faceID, quad.QPoint(0));
    quadPointCellLocals.push_back(qCellLocal);
    this->trafo = c.GetTransformation(quadPointCellLocals[0]);
    computeQPointData(refArea, c, 0);
    for (int q = 1; q < quad.size(); ++q) {
      quadPointCellLocals.push_back(c.FaceLocalToLocalOriented(faceID, quad.QPoint(q)));
      computeQPointData(refArea, c, q);
    }
#else
    this->trafos.resize(quad.size());
    for (int q = 0; q < quad.size(); ++q) {
      auto qCellLocal = c.FaceLocalToLocalOriented(faceID, quad.QPoint(q));
      this->trafos[q] = c.GetTransformation(qCellLocal);
      quadPointCellLocals.push_back(qCellLocal);
      computeQPointData(refArea, c, q);
    }
#endif
    if (mid(quadNormals[0]) < Origin) sign = -1;
    for (auto weight : this->quadWeights)
      area += weight;
  }

  TT Area() const override { return area; }

  const PointT<TT, sDim, tDim> &QLocal(int q) const override { return quadPointCellLocals[q]; }

  const PointT<TT, sDim, tDim> &QNormal(int q) const override { return quadNormals[q]; }

  short Sign() const override { return sign; }

  int FaceID() const override { return faceID; }
private:
  void computeQPointData(TT refArea, const Cell &c, int q) {
    this->quadPoints.push_back(c.LocalToGlobal(quadPointCellLocals[q]));
    auto qNormal = this->GetTransformation(q)
                   * PointT<TT, sDim, tDim>(c.LocalFaceNormalT<TT, sDim, tDim>(faceID));
    TT w = norm(qNormal);
    this->quadWeights.push_back(this->GetTransformation(q).Det() * w * refArea
                                * this->quad.Weight(q));
    quadNormals.push_back(qNormal / w);
  }
};

template<typename DATA>
struct QuadIterator {
protected:
  DATA data;
public:
  QuadIterator(DATA &&data) : data(std::move(data)) {}

  QuadIterator &operator++() {
    data.Incr();
    return *this;
  }

  bool operator==(const QuadIterator &other) const { return data.id == other.data.id; }

  bool operator!=(const QuadIterator &other) const { return data.id != other.data.id; }

  typename DATA::Refs operator*() const { return data.References(); }
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
struct QuadPairT {
  int id;
  const TT *weight = nullptr;

  void Incr() {
    id++;
    weight++;
  }

  struct Refs {
    int id;
    const TT &weight;
  };

  Refs References() const { return Refs{id, *weight}; }

  struct Iter {
    const ElementDataT<TT, sDim, tDim> &elementData;

    QuadIterator<QuadPairT<TT, sDim, tDim>> begin() const {
      return QuadIterator<QuadPairT<TT, sDim, tDim>>(
          QuadPairT<TT, sDim, tDim>{0, &elementData.quadWeights[0]});
    }

    QuadIterator<QuadPairT<TT, sDim, tDim>> end() const {
      return QuadIterator<QuadPairT<TT, sDim, tDim>>(
          QuadPairT<TT, sDim, tDim>{elementData.quad.size(), nullptr});
    }
  };
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
struct QuadTripleT {
  int id;
  const TT *weight = nullptr;
  const PointT<TT, sDim, tDim> *point = nullptr;

  void Incr() {
    id++;
    weight++;
    point++;
  }

  struct Refs {
    int id;
    const TT &weight;
    const PointT<TT, sDim, tDim> &point;
  };

  Refs References() const { return Refs{id, *weight, *point}; }

  struct Iter {
    const ElementDataT<TT, sDim, tDim> &elementData;

    QuadIterator<QuadTripleT<TT, sDim, tDim>> begin() const {
      return QuadIterator<QuadTripleT<TT, sDim, tDim>>(
          QuadTripleT<TT, sDim, tDim>{0, &elementData.quadWeights[0], &elementData.quadPoints[0]});
    }

    QuadIterator<QuadTripleT<TT, sDim, tDim>> end() const {
      return QuadIterator<QuadTripleT<TT, sDim, tDim>>(
          QuadTripleT<TT, sDim, tDim>{elementData.quad.size(), nullptr, nullptr});
    }
  };
};

template<typename TT = double, int sDim = SpaceDimension, int tDim = TimeDimension>
struct QuadQuadrupleT {
  int id;
  const TT *weight = nullptr;
  const PointT<TT, sDim, tDim> *point = nullptr;
  const PointT<TT, sDim, tDim> *localPoint = nullptr;

  void Incr() {
    id++;
    weight++;
    point++;
    localPoint++;
  }

  struct Refs {
    int id;
    const TT &weight;
    const PointT<TT, sDim, tDim> &point;
    const PointT<TT, sDim, tDim> &localPoint;
  };

  Refs References() const { return Refs{id, *weight, *point, *localPoint}; }

  struct Iter {
    const ElementDataT<TT, sDim, tDim> &elementData;

    QuadIterator<QuadQuadrupleT<TT, sDim, tDim>> begin() const {
      return QuadIterator<QuadQuadrupleT<TT, sDim, tDim>>(
          QuadQuadrupleT<TT, sDim, tDim>{0, &elementData.quadWeights[0], &elementData.quadPoints[0],
                                         &elementData.quad.QPoints()[0]});
    }

    QuadIterator<QuadQuadrupleT<TT, sDim, tDim>> end() const {
      return QuadIterator<QuadQuadrupleT<TT, sDim, tDim>>(
          QuadQuadrupleT<TT, sDim, tDim>{elementData.quad.size(), nullptr, nullptr, nullptr});
    }
  };
};

#endif // ELEMENTDATA_HPP
