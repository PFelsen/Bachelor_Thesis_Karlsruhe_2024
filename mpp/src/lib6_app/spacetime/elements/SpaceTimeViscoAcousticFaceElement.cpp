#include "SpaceTimeViscoAcousticFaceElement.hpp"

SpaceTimeViscoAcousticFaceElement::SpaceTimeViscoAcousticFaceElement
  (const STDiscretization &D,
   const VectorMatrixBase &g,
   const cell &tc,
   int f_id,
   int nL)
  :
  deg(g.GetDoF().GetDegree(*tc)),
  SS(D.GetCellShape(deg.space)),
  TS(D.GetTimeShape(deg.time)),
  TestSpace(D.GetTimeShape(deg.time - 1, true)),
  faceqpoints(D.FaceQPoints()),
  TC(tc), f(g.find_face(tc.Face(f_id))), r(g.find_row(tc())),
  fid(f_id), dim(tc.dim()),
  SpaceQ(D.GetCellFaceQuad(f_id, deg.space)),
  TimeQ(D.GetTimeFaceQuad(f_id, deg.time)),
  numL(nL) {

  ww = TC.SpaceCell().LocalFaceArea(fid);

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());


  velocityST.resize(nQ() * j_dimension());
  velocityST_TestSpace.resize(nQ() * i_dimension());
  pressureST.resize(nQ() * j_dimension());
  pressureST_TestSpace.resize(nQ() * i_dimension());
  T = TC.GetTransformation(faceqpoints[fid][0]);
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      int idx = tq * SpaceQ.size() + sq;
      qPoint[idx] = TC.LocalToGlobal(faceqpoints[fid][sq].CopyWithT(TimeQ.QPoint(tq)[0]));
      qNormal[idx] = T * TC.SpaceCell().LocalFaceNormal(fid);
      double w = norm(qNormal[idx]);
      qWeight[idx] = T.Det() * w * ww * SpaceQ.Weight(sq) * TimeQ.Weight(tq);
      qNormal[idx] /= w;
    }
  }

  ShapeValues<Scalar> spaceValue(SpaceQ.size(), SS.size());
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int i = 0; i < SS.size(); ++i) spaceValue[sq][i] = SS(fid, sq, i);
  }
  ShapeValues<Scalar> timeValue(TimeQ.size(), TS.size());
  ShapeValues<Scalar> testspaceValue(TimeQ.size(), TestSpace.size());
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int i = 0; i < TS.size(); ++i)
      timeValue[tq][i] = TS(TimeQ.QPoint(tq), i);
    for (int i = 0; i < TestSpace.size(); ++i)
      testspaceValue[tq][i] = TestSpace(TimeQ.QPoint(tq), i);
  }
  int offset = SS.size() * (1 + numL + dim);
  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      int offset1 = tq * nSpaceQ() + sq;
      for (int si = 0; si < SS.size(); ++si) {
        for (int k = 0; k < dim; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = offset1 * j_dimension() + ti * offset + k * SS.size() + si;
            velocityST[index][k] = spaceValue[sq][si] * double(timeValue[tq][ti]);
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = offset1 * i_dimension() + ti * offset + k * SS.size() + si;
            velocityST_TestSpace[index][k] = spaceValue[sq][si] * double(testspaceValue[tq][ti]);
          }
        }
        for (int k = 0; k <= numL; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = offset1 * j_dimension() + ti * offset + (dim + k) * SS.size() + si;
            pressureST[index] = spaceValue[sq][si] * double(timeValue[tq][ti]);
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = offset1 * i_dimension() + ti * offset + (dim + k) * SS.size() + si;
            pressureST_TestSpace[index] = spaceValue[sq][si] * double(testspaceValue[tq][ti]);
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    for (int j = 0; j < j_dimension(); ++j) {
      j_zero[q * j_dimension() + j] = (std::abs(Pressure(q, j)) == 0)
        && (norm(Velocity(q, j)) == 0);
    }
    for (int i = 0; i < i_dimension(); ++i) {
      i_zero[q * i_dimension() + i] = ((std::abs(PressureTestSpace(q, i)) == 0)
        && (norm(VelocityTestSpace(q, i)) == 0));
    }
  }
}
