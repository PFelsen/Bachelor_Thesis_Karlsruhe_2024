#include "SpaceTimeViscoElasticElement.hpp"

SpaceTimeViscoElasticElement::SpaceTimeViscoElasticElement
    (const STDiscretization &disc, const VectorMatrixBase &g, const cell &tc, int nL)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(disc.GetCellShape(deg.space)),
    SS2(SS),
    TS(disc.GetTimeShape(deg.time)),
    TestSpace(disc.GetTimeShape(deg.time - 1, true)),
    spaceValue(SS.values()),
    space_value2(spaceValue),
    timeValue(TS.values()),
    testspaceValue(TestSpace.values()),
    TC(tc), dim(tc.dim()),
    SpaceQ(disc.GetCellQuad(deg.space)),
    TimeQ(disc.GetTimeQuad(deg.time)),
    numL(nL) {

  qWeight.resize(nQ());
  space_qWeight.resize(nSpaceQ());
  time_qWeight.resize(nTimeQ());
  qPoint.resize(nQ());
  T.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  dtVelocityST.resize(nQ() * j_dimension());
  strainST.resize(nQ() * j_dimension());
  velocityST_TestSpace.resize(nQ() * i_dimension());

  stressST.resize(nQ() * j_dimension());
  dtStressST.resize(nQ() * j_dimension());
  divStressST.resize(nQ() * j_dimension());
  stressST_TestSpace.resize(nQ() * i_dimension());
  /*stressST.resize( 1+numL );
  dtStressST.resize( 1+numL );
  divStressST.resize( 1+numL );
  stressST_TestSpace.resize( 1+numL );
  for(k=0;k<1+numL){
      stressST[k].resize( nQ()*j_dimension() );
      dtStressST[k].resize( nQ()*j_dimension() );
      divStressST[k].resize( nQ()*j_dimension() );
      stressST_TestSpace[k].resize( nQ()*i_dimension() );
  }*/

  time_gradient.resize(TimeQ.size());
  vector<vector<VectorField>> gradient;
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
      qWeight[tq * SpaceQ.size() + sq]
          = T[sq].Det() * SpaceQ.Weight(sq) * timedet * TimeQ.Weight(tq);
      qPoint[tq * SpaceQ.size() + sq] = tc[SpaceQ.QPoint(sq)].
          CopyWithT(a + timedet * TimeQ.QPoint(tq)[0]);
    }

  int S1S2 = SS.size() * (dim + td()) + SS2.size() * td() * numL;
  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        //velocity
        for (int k = 0; k < dim; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * S1S2 + k * SS.size() + si;
            velocityST[index][k]
                = spaceValue[sq][si] * double(timeValue[tq][ti]);
            dtVelocityST[index][k]
                = spaceValue[sq][si] * double(time_gradient[tq][ti][0]);
            Tensor Dv = Zero;
            Dv[k] = gradient[sq][si] * double(timeValue[tq][ti]);
            strainST[index] = sym(Dv);
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * i_dimension()
                        + ti * S1S2 + k * SS.size() + si;
            velocityST_TestSpace[index][k]
                = spaceValue[sq][si] * double(testspaceValue[tq][ti]);
          }
        }

        //stress
        for (int k = 0; k < td(); ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * S1S2 + (dim + k) * SS.size() + si;
            if (k < dim) {
              stressST[index][k][k] = spaceValue[sq][si] * double(timeValue[tq][ti]);
              dtStressST[index][k][k] = spaceValue[sq][si] * double(time_gradient[tq][ti][0]);
              divStressST[index][k] = gradient[sq][si][k] * double(timeValue[tq][ti]);
            } else if (dim == 2) {
              stressST[index][0][1] = spaceValue[sq][si] * double(timeValue[tq][ti]);
              stressST[index][1][0] = stressST[index][0][1];
              dtStressST[index][0][1] = spaceValue[sq][si] * double(time_gradient[tq][ti][0]);
              dtStressST[index][1][0] = dtStressST[index][0][1];
              divStressST[index][0] = gradient[sq][si][1] * double(timeValue[tq][ti]);
              divStressST[index][1] = gradient[sq][si][0] * double(timeValue[tq][ti]);
            } else Exit("ToDo dim=3")
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * i_dimension()
                        + ti * S1S2 + (dim + k) * SS.size() + si;
            if (k < dim) {
              stressST_TestSpace[index][k][k] =
                  spaceValue[sq][si] * double(testspaceValue[tq][ti]);
            } else if (dim == 2) {
              stressST_TestSpace[index][0][1] =
                  spaceValue[sq][si] * double(testspaceValue[tq][ti]);
              stressST_TestSpace[index][1][0] = stressST_TestSpace[index][0][1];
            } else Exit("ToDo dim=3")
          }
        }
      }

      for (int nL = 1; nL <= numL; nL++) {
        for (int si = 0; si < SS2.size(); ++si) {
          //visco stress
          for (int k = 0; k < td(); ++k) {
            for (int ti = 0; ti < TS.size(); ++ti) {
              int index = (tq * nSpaceQ() + sq) * j_dimension()
                          + ti * S1S2 + (dim + td()) * SS.size()
                          + k * SS2.size() + si;
              if (k < dim) {
                stressST[index][k][k] = space_value2[sq][si] * double(timeValue[tq][ti]);
                dtStressST[index][k][k] = space_value2[sq][si] * double(time_gradient[tq][ti][0]);

                divStressST[index][k] = gradient[sq][si][k] * double(timeValue[tq][ti]);
              } else if (dim == 2) {
                stressST[index][0][1] = space_value2[sq][si] * double(timeValue[tq][ti]);
                stressST[index][1][0] = stressST[index][0][1];
                dtStressST[index][0][1] = space_value2[sq][si] * double(time_gradient[tq][ti][0]);
                dtStressST[index][1][0] = dtStressST[index][0][1];

                divStressST[index][0] = gradient[sq][si][1] * double(timeValue[tq][ti]);
                divStressST[index][1] = gradient[sq][si][0] * double(timeValue[tq][ti]);
              } else Exit("ToDo dim=3")
            }
            for (int ti = 0; ti < TestSpace.size(); ++ti) {
              int index = (tq * nSpaceQ() + sq) * i_dimension()
                          + ti * S1S2 + (dim + td()) * SS.size()
                          + k * SS2.size() + si;
              if (k < dim) {
                stressST_TestSpace[index][k][k] =
                    space_value2[sq][si] * double(testspaceValue[tq][ti]);
              } else if (dim == 2) {
                stressST_TestSpace[index][0][1] =
                    space_value2[sq][si] * double(testspaceValue[tq][ti]);
                stressST_TestSpace[index][1][0] = stressST_TestSpace[index][0][1];
              } else Exit("ToDo dim=3")
            }
          }
        }
      }
    }
  }
}

SpaceTimeViscoElasticFaceElement::SpaceTimeViscoElasticFaceElement
    (const STDiscretization &D,
     const VectorMatrixBase &g,
     const cell &tc,
     int f_id,
     int nL)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(D.GetCellShape(deg.space)),
    SS2(SS),
    TS(D.GetTimeShape(deg.time)),
    TestSpace(D.GetTimeShape(deg.time - 1, true)),
    faceqpoints(D.FaceQPoints()),
    TC(tc), f(g.find_face(tc.Face(f_id))),
    fid(f_id), dim(tc.dim()),
    SpaceQ(D.GetCellFaceQuad(f_id, deg.space)),
    TimeQ(D.GetTimeFaceQuad(f_id, deg.time)),
    numL(nL) {

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
      /** Transformation **/
      qPoint[tq * SpaceQ.size() + sq] = TC[faceqpoints[fid][sq]].
          CopyWithT(a + timedet * TimeQ.QPoint(tq)[0]);

      T[tq * SpaceQ.size() + sq] = TC.GetTransformation(faceqpoints[fid][sq]);
      qNormal[tq * SpaceQ.size() + sq] = T[tq * SpaceQ.size() + sq] * TC.LocalFaceNormal(fid);
      double w = norm(qNormal[tq * SpaceQ.size() + sq]);
      qWeight[tq * SpaceQ.size() + sq] = T[tq * SpaceQ.size() + sq].Det()
                                         * w * ww * SpaceQ.Weight(sq) * timedet * TimeQ.Weight(tq);
      qNormal[tq * SpaceQ.size() + sq] /= w;
      qTangent[tq * SpaceQ.size() + sq] = TC.FaceCorner(f_id, 1) - TC.FaceCorner(f_id, 0);
      w = norm(qTangent[tq * SpaceQ.size() + sq]);
      qTangent[tq * SpaceQ.size() + sq] /= w;
    }
  }

  double spaceValue[SpaceQ.size()][SS.size()];
  double space_value2[SpaceQ.size()][SS.size()];
  //ShapeValue space_value2;
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int i = 0; i < SS.size(); ++i) spaceValue[sq][i] = SS(fid, sq, i);
    //for (int i=0; i<SS2.size(); ++i)
    //    space_value2[sq][i] = SS2(fid,sq,i);
  }
  double timeValue[TimeQ.size()][TS.size()];
  double testspaceValue[TimeQ.size()][TS.size()];
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int i = 0; i < TS.size(); ++i)
      timeValue[tq][i] = TS(TimeQ.QPoint(tq), i);
    for (int i = 0; i < TestSpace.size(); ++i)
      testspaceValue[tq][i] = TestSpace(TimeQ.QPoint(tq), i);
  }

  int S1S2 = SS.size() * (dim + td()) + SS2.size() * td() * numL;
  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        //velocity
        for (int k = 0; k < dim; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * j_dimension() + ti * S1S2 + k * SS.size() + si;
            velocityST[index][k] = spaceValue[sq][si] * double(timeValue[tq][ti]);
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * i_dimension()
                        + ti * S1S2 + k * SS.size() + si;
            velocityST_TestSpace[index][k] = spaceValue[sq][si] * double(testspaceValue[tq][ti]);
          }
        }

        //stress
        for (int k = 0; k < td(); ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * S1S2 + (dim + k) * SS.size() + si;
            if (k < dim) {
              stressST[index][k][k] = spaceValue[sq][si] * double(timeValue[tq][ti]);
            } else if (dim == 2) {
              stressST[index][0][1] = spaceValue[sq][si] * double(timeValue[tq][ti]);
              stressST[index][1][0] = stressST[index][0][1];
            } else Exit("ToDo dim=3")
          }
          for (int ti = 0; ti < TestSpace.size(); ++ti) {
            int index = (tq * nSpaceQ() + sq) * i_dimension()
                        + ti * S1S2 + (dim + k) * SS.size() + si;
            if (k < dim) {
              stressST_TestSpace[index][k][k] =
                  spaceValue[sq][si] * double(testspaceValue[tq][ti]);
            } else if (dim == 2) {
              stressST_TestSpace[index][0][1] =
                  spaceValue[sq][si] * double(testspaceValue[tq][ti]);
              stressST_TestSpace[index][1][0] = stressST_TestSpace[index][0][1];
            } else Exit("ToDo dim=3")
          }
        }
      }
      for (int nL = 1; nL <= numL; nL++) {
        for (int si = 0; si < SS2.size(); ++si) {
          //visco stress
          for (int k = 0; k < td(); ++k) {
            for (int ti = 0; ti < TS.size(); ++ti) {
              int index = (tq * nSpaceQ() + sq) * j_dimension()
                          + ti * S1S2 + (dim + td()) * SS.size()
                          + k * SS2.size() + si;
              if (k < dim) {
                stressST[index][k][k] = space_value2[sq][si] * double(timeValue[tq][ti]);
              } else if (dim == 2) {
                stressST[index][0][1] = space_value2[sq][si] * double(timeValue[tq][ti]);
                stressST[index][1][0] = stressST[index][0][1];
              } else Exit("ToDo dim=3")
            }
            for (int ti = 0; ti < TestSpace.size(); ++ti) {
              int index = (tq * nSpaceQ() + sq) * i_dimension()
                          + ti * S1S2 + (dim + td()) * SS.size()
                          + k * SS2.size() + si;
              if (k < dim) {
                stressST_TestSpace[index][k][k] =
                    space_value2[sq][si] * double(testspaceValue[tq][ti]);
              } else if (dim == 2) {
                stressST_TestSpace[index][0][1] =
                    space_value2[sq][si] * double(testspaceValue[tq][ti]);
                stressST_TestSpace[index][1][0] = stressST_TestSpace[index][0][1];
              } else Exit("ToDo dim=3")
            }
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    for (int j = 0; j < i_dimension(); ++j)
      j_zero[q * j_dimension() + j] =
          ((norm(Stress(q, j)) == 0) && (norm(Velocity(q, j)) == 0));
    for (int i = 0; i < i_dimension(); ++i)
      i_zero[q * i_dimension() + i] =
          ((norm(StressTestSpace(q, i)) == 0)
           && (norm(VelocityTestSpace(q, i)) == 0));
  }
}