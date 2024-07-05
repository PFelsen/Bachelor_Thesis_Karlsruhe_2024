#include "Results.hpp"
#include "ViscoAcousticElements.hpp"

SpaceTimeViscoAcousticElement::SpaceTimeViscoAcousticElement
    (const STDiscretization &disc, const VectorMatrixBase &g, const cell &tc, int nL)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(disc.GetCellShape(deg.space)),
    TS(disc.GetTimeShape(deg.time)),
    TestSpace(disc.GetTimeShape(deg.time - 1, true)),
    space_value(SS.values()),
    time_value(TS.values()),
    testspace_value(TestSpace.values()),
    dim(tc.dim()),
    SpaceQ(disc.GetCellQuad(deg.space)),
    TimeQ(disc.GetTimeQuad(deg.time)),
    numL(nL) {

  qWeight.resize(nQ());
  space_qWeight.resize(nSpaceQ());
  time_qWeight.resize(nTimeQ());
  qPoint.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  dtVelocityST.resize(nQ() * j_dimension());
  divVelocityST.resize(nQ() * j_dimension());
  velocityST_TestSpace.resize(nQ() * i_dimension());

  pressureST.resize(nQ() * j_dimension());
  dtPressureST.resize(nQ() * j_dimension());
  gradPressureST.resize(nQ() * j_dimension());
  pressureST_TestSpace.resize(nQ() * i_dimension());

  time_gradient.resize(TimeQ.size());
  std::vector<std::vector<VectorField> > gradient(nQ());

  T = tc.GetTransformation(SpaceQ.QPoint(0));
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      gradient[tq * SpaceQ.size() + sq].resize(SS.size());
      for (int i = 0; i < SS.size(); ++i) {
        gradient[tq * SpaceQ.size() + sq][i] = T * SS.LocalGradient(sq, i);
      }
    }

    time_gradient[tq].resize(TS.size());
    for (int i = 0; i < TS.size(); ++i) {
      time_gradient[tq][i] = (1. / T.DetTime()) * TS.LocalGradient(tq, i);
    }
  }

  //Quadraturgewichte
  for (int tq = 0; tq < TimeQ.size(); ++tq)
    time_qWeight[tq] = T.DetTime() * TimeQ.Weight(tq);
  for (int sq = 0; sq < SpaceQ.size(); ++sq)
    space_qWeight[sq] = T.DetSpace() * SpaceQ.Weight(sq);
  //Quadraturpunkte
  for (int tq = 0; tq < TimeQ.size(); ++tq)
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      int index = tq * SpaceQ.size() + sq;
      qWeight[index] = space_qWeight[sq] * time_qWeight[tq];
      qPoint[index] = tc.LocalToGlobal(SpaceQ.QPoint(sq).CopyWithT(TimeQ.QPoint(tq)[0]));
    }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      int idQ = (tq * nSpaceQ() + sq);
      for (int si = 0; si < SS.size(); ++si) {
        //velocity
        for (int k = 0; k < dim; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index =
                idQ * j_dimension() + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST[index][k] = space_value[sq][si] * double(time_value[tq][ti]);
            dtVelocityST[index][k] = space_value[sq][si] * double(time_gradient[tq][ti][0]);
            divVelocityST[index] = gradient[sq][si][k] * double(time_value[tq][ti]);
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index =
                idQ * i_dimension() + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST_TestSpace[index][k] = space_value[sq][si] * double(testspace_value[tq][ti]);
          }
        }

        //pressure
        for (int k = 0; k <= numL; ++k) {
          int id_k = (dim + k) * SS.size();
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = idQ * j_dimension() + ti * SS.size() * (1 + numL + dim) + id_k + si;
            pressureST[index] = space_value[sq][si] * double(time_value[tq][ti]);
            dtPressureST[index] = space_value[sq][si] * double(time_gradient[tq][ti][0]);
            gradPressureST[index] = gradient[sq][si] * double(time_value[tq][ti]);
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = idQ * i_dimension() + ti * SS.size() * (1 + numL + dim) + id_k + si;
            pressureST_TestSpace[index] = space_value[sq][si] * double(testspace_value[tq][ti]);
          }
        }
      }
    }
  }
}

