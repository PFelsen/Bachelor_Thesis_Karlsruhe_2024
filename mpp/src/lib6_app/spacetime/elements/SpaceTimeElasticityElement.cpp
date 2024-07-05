#include "SpaceTimeElasticityElement.hpp"

SpaceTimeElasticityElement::SpaceTimeElasticityElement(const STDiscretization &disc,
                                                       const VectorMatrixBase &g,
                                                       const cell &tc)
  :
  deg(g.GetDoF().GetDegree(*tc)),
  SS(disc.GetCellShape(deg.space)),
  TS(disc.GetTimeShape(deg.time)),
  TestSpace(disc.GetTimeShape(deg.time - 1)),
  spaceValue(SS.values()),
  timeValue(TS.values()),
  testspaceValue(TestSpace.values()),
  TC(tc), dim(tc.dim()),
  SpaceQ(disc.GetCellQuad(deg.space)),
  TimeQ(disc.GetTimeQuad(deg.time)) {

  qWeight.resize(nQ());
  space_qWeight.resize(nSpaceQ());
  time_qWeight.resize(nTimeQ());
  qPoint.resize(nQ());
  T.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  dtvelocityST.resize(nQ() * j_dimension());
  strainST.resize(nQ() * j_dimension());
  velocityST_TestSpace.resize(nQ() * i_dimension());
  stressST.resize(nQ() * j_dimension());
  dtstressST.resize(nQ() * j_dimension());
  divstressST.resize(nQ() * j_dimension());
  stressST_TestSpace.resize(nQ() * i_dimension());

  time_gradient.resize(TimeQ.size());
  vector<vector<VectorField> > gradient;
  gradient.resize(nQ());

  double a = tc.min();
  double b = tc.max();
  timedet = b - a;

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      T[tq * SpaceQ.size() + sq] = tc.GetTransformation(SpaceQ.QPoint(sq));
      gradient[tq * SpaceQ.size() + sq].resize(SS.size());
      for (int i = 0; i < SS.size(); ++i) {
        gradient[tq * SpaceQ.size() + sq][i] = T[sq] * SS.LocalGradient(sq, i);
      }
    }

    time_gradient[tq].resize(TS.size());
    for (int i = 0; i < TS.size(); ++i) {
      time_gradient[tq][i] = (1. / timedet) * TS.LocalGradient(tq, i);
    }
  }

  //Quadraturgewichte
  for (int tq = 0; tq < TimeQ.size(); ++tq)
    time_qWeight[tq] = timedet * TimeQ.Weight(tq);
  for (int sq = 0; sq < SpaceQ.size(); ++sq)
    space_qWeight[sq] = T[sq].Det() * SpaceQ.Weight(sq);
  //Quadraturpunkte
  for (int tq = 0; tq < TimeQ.size(); ++tq)
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      qWeight[tq * SpaceQ.size() + sq] =
        T[sq].Det() * SpaceQ.Weight(sq) * timedet * TimeQ.Weight(tq);
      qPoint[tq * SpaceQ.size() + sq] = tc[SpaceQ.QPoint(sq)].
        CopyWithT(a + timedet * TimeQ.QPoint(tq)[0]);
    }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int k = 0; k < mod(); ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index =
              (tq * nSpaceQ() + sq) * j_dimension() + ti * SS.size() * mod()
                + SS.size() * k + si;
            velocityST[index] = zero;
            dtvelocityST[index] = zero;
            strainST[index] = Zero;
            stressST[index] = Zero;
            dtstressST[index] = Zero;
            divstressST[index] = zero;
            if (k < dim) {
              stressST[index][k][k] =
                spaceValue[sq][si] * double(timeValue[tq][ti]);
              dtstressST[index][k][k] =
                spaceValue[sq][si] * double(time_gradient[tq][ti][0]);
              divstressST[index][k] =
                gradient[sq][si][k] * double(timeValue[tq][ti]);
            } else if (k < mod() - dim) {
              if (dim == 2) {
                stressST[index][0][1] =
                  spaceValue[sq][si] * double(timeValue[tq][ti]);
                stressST[index][1][0] = stressST[index][0][1];

                dtstressST[index][0][1] = spaceValue[sq][si]
                  * double(time_gradient[tq][ti][0]);
                dtstressST[index][1][0] = dtstressST[index][0][1];

                divstressST[index][0] =
                  gradient[sq][si][1] * double(timeValue[tq][ti]);
                divstressST[index][1] =
                  gradient[sq][si][0] * double(timeValue[tq][ti]);
              }
              if (dim == 3) {
                Exit("not implemented");//TODO
              }
            } else {
              int k_ind = k + dim - mod();
              velocityST[index][k_ind] =
                spaceValue[sq][si] * double(timeValue[tq][ti]);
              dtvelocityST[index][k_ind] =
                spaceValue[sq][si] * double(time_gradient[tq][ti][0]);
              Tensor Dv = Zero;
              Dv[k_ind] = gradient[sq][si] * double(timeValue[tq][ti]);
              strainST[index] = sym(Dv);
            }
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index =
              (tq * nSpaceQ() + sq) * i_dimension() + ti * SS.size() * mod()
                + SS.size() * k + si;
            velocityST_TestSpace[index] = zero;
            stressST_TestSpace[index] = Zero;
            if (k < dim) {
              stressST_TestSpace[index][k][k] =
                spaceValue[sq][si] * double(testspaceValue[tq][ti]);
            } else if (k < mod() - dim) {
              if (dim == 2) {
                stressST_TestSpace[index][0][1] =
                  spaceValue[sq][si] * double(testspaceValue[tq][ti]);
                stressST_TestSpace[index][1][0] =
                  stressST_TestSpace[index][0][1];
              }
              /*if (dim==3) {
                 TODO
              }*/

            } else {
              int k_ind = k + dim - mod();
              velocityST_TestSpace[index][k_ind] =
                spaceValue[sq][si] * double(testspaceValue[tq][ti]);
            }
          }
        }
      }
    }
  }
}

SpaceTimeElasticityFaceElement::SpaceTimeElasticityFaceElement
  (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc, int f_id)
  :
  deg(g.GetDoF().GetDegree(*tc)),
  SS(D.GetCellShape(deg.space)),
  TS(D.GetTimeShape(deg.time)),
  TestSpace(D.GetTimeShape(deg.time - 1)),
  timeValue(TS.values()),
  testspaceValue(TestSpace.values()),
  TC(tc), f(g.find_face(tc.Face(f_id))),
  fid(f_id), dim(tc.dim()),
  SpaceQ(D.GetCellFaceQuad(f_id, deg.space)),
  TimeQ(D.GetTimeFaceQuad(f_id, deg.time)) {
  ww = TC.LocalFaceArea(fid);

  double a = TC.min();
  double b = TC.max();
  timedet = b - a;

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  qTangent.resize(nQ());
  T.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  velocityST_TestSpace.resize(nQ() * i_dimension());
  stressST.resize(nQ() * j_dimension());
  stressST_TestSpace.resize(nQ() * i_dimension());

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      const Point &x = SpaceQ.QPoint(sq);
      Point y = Local(x);
      qLocal[tq * SpaceQ.size() + sq] = y;

      qPoint[tq * SpaceQ.size() + sq] = TC[y].
        CopyWithT(a + timedet * TimeQ.QPoint(tq)[0]);

      T[tq * SpaceQ.size() + sq] = TC.GetTransformation(y);
      qNormal[tq * SpaceQ.size() + sq] =
        T[tq * SpaceQ.size() + sq] * TC.LocalFaceNormal(fid);
      double w = norm(qNormal[tq * SpaceQ.size() + sq]);
      qWeight[tq * SpaceQ.size() + sq] = T[tq * SpaceQ.size() + sq].Det()
        * w * ww * SpaceQ.Weight(sq) * timedet * TimeQ.Weight(tq);
      qNormal[tq * SpaceQ.size() + sq] /= w;
      qTangent[tq * SpaceQ.size() + sq] =
        TC.FaceCorner(f_id, 1) - TC.FaceCorner(f_id, 0);
      w = norm(qTangent[tq * SpaceQ.size() + sq]);
      qTangent[tq * SpaceQ.size() + sq] /= w;
    }
  }

  double value[SpaceQ.size()][SS.size()];
  for (int q = 0; q < SpaceQ.size(); ++q) {
    for (int i = 0; i < SS.size(); ++i) {
      value[q][i] = SS(QLocal(q), i);
    }
  }

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      T[tq * SpaceQ.size() + sq] = tc.GetTransformation(SpaceQ.QPoint(sq));
    }
  }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int k = 0; k < mod(); ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index =
              (tq * nSpaceQ() + sq) * j_dimension() + ti * SS.size() * mod()
                + SS.size() * k + si;
            velocityST[index] = zero;
            stressST[index] = Zero;
            if (k < dim) {
              stressST[index][k][k] =
                value[sq][si] * double(timeValue[tq][ti]);
            } else if (k < mod() - dim) {
              if (dim == 2) {
                stressST[index][0][1] =
                  value[sq][si] * double(timeValue[tq][ti]);
                stressST[index][1][0] = stressST[index][0][1];
              }

            } else {
              int k_ind = k + dim - mod();
              velocityST[index][k_ind] =
                value[sq][si] * double(timeValue[tq][ti]);
            }
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index =
              (tq * nSpaceQ() + sq) * i_dimension() + ti * SS.size() * mod()
                + SS.size() * k + si;
            velocityST_TestSpace[index] = zero;
            stressST_TestSpace[index] = Zero;
            if (k < dim) {
              stressST_TestSpace[index][k][k] =
                value[sq][si] * double(testspaceValue[tq][ti]);
            } else if (k < mod() - dim) {
              if (dim == 2) {
                stressST_TestSpace[index][0][1] =
                  value[sq][si] * double(testspaceValue[tq][ti]);
                stressST_TestSpace[index][1][0] =
                  stressST_TestSpace[index][0][1];
              }
            } else {
              int k_ind = k + dim - mod();
              velocityST_TestSpace[index][k_ind] =
                value[sq][si] * double(testspaceValue[tq][ti]);
            }
          }
        }
      }
    }
  }
}


Point SpaceTimeElasticityFaceElement::Local(const Point &local) const {
  switch (TC.Type()) {
    case TRIANGLE:
      return (1 - local[0]) * LocalFaceCorner(0) + local[0] * LocalFaceCorner(1);
    case TETRAHEDRON: Exit("do it");
    case QUADRILATERAL:
      return (1 - local[0]) * LocalFaceCorner(0) + local[0] * LocalFaceCorner(1);
    case HEXAHEDRON: Exit("do it");
  }
  Exit("Not implemented");
}
