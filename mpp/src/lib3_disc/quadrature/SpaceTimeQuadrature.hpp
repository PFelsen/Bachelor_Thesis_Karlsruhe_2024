#ifndef SPACETIME_SPACETIMEQUADRATURE_HPP
#define SPACETIME_SPACETIMEQUADRATURE_HPP

#include "Quadrature.hpp"

template<typename T = double, int sDim = SpaceDimension, int tDim = TimeDimension>
class SpaceTimeQuadratureT : public QuadratureT<T, sDim, tDim> {
protected:
  const QuadratureT<T, sDim, tDim> &SQ;
  const QuadratureT<T, sDim, tDim> &TQ;
public:
  explicit SpaceTimeQuadratureT(const QuadratureT<T, sDim, tDim> &SQ,
                                const QuadratureT<T, sDim, tDim> &TQ) :
      QuadratureT<T, sDim, tDim>("SpaceTimeQuadrature: [" + SQ.Name() + ", " + TQ.Name() + "]",
                                 SQ.size() * TQ.size()),
      SQ(SQ), TQ(TQ) {
    int n = 0;
    for (int tq = 0; tq < TQ.size(); ++tq) {
      double t = TQ.QPoint(tq)[0];
      for (int sq = 0; sq < SQ.size(); ++sq) {
        this->z[n] = SQ.QPoint(sq).CopyWithT(t);
        this->w[n++] = SQ.Weight(sq) * TQ.Weight(tq);
      }
    }
  }

  const QuadratureT<T, sDim, tDim> &GetSpaceQuad() const { return SQ; }

  const QuadratureT<T, sDim, tDim> &GetTimeQuad() const { return TQ; }

  virtual ~SpaceTimeQuadratureT() = default;
};

using SpaceTimeQuadrature = SpaceTimeQuadratureT<>;

#endif