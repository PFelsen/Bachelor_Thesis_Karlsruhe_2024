#include "DGTViscoAcousticElementsGLGL.hpp"
// Todo remove this header again and the function
// transformPointLocaltoGlobal() in project
#include "SpaceTimeTools.hpp"

SpaceTimeViscoAcousticDGTElementGLGL::SpaceTimeViscoAcousticDGTElementGLGL
    (const STDiscretization &disc, const VectorMatrixBase &g, const cell &tc, int nL)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(disc.GetCellShape(deg.space)),
    TS(disc.GetTimeShape(deg.time)),
    spaceValue(SS.values()),
    timeValue(TS.values()),
    TC(tc), dim(tc.dim()),
    SpaceQ(disc.GetCellQuad(deg.space)),
    TimeQ(disc.GetTimeQuad(deg.time)),
    numL(nL) {

  qWeight.resize(nQ());
  space_qWeight.resize(nSpaceQ());
  time_qWeight.resize(nTimeQ());
  qPoint.resize(nQ());
  std::vector<Transformation> T(nQ());

  velocityST.resize(nQ() * j_dimension());
  dtVelocityST.resize(nQ() * j_dimension());
  divVelocityST.resize(nQ() * j_dimension());

  pressureST.resize(nQ() * j_dimension());
  dtPressureST.resize(nQ() * j_dimension());
  gradPressureST.resize(nQ() * j_dimension());

  time_gradient.resize(TimeQ.size());
  std::vector<std::vector<VectorField> > gradient(nQ());

  double a = tc.min();
  double b = tc.max();
  timedet = b - a;

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      T[tq * SpaceQ.size() + sq] = tc.SpaceCell().GetTransformation(SpaceQ.QPoint(sq));
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


  for (int tq = 0; tq < TimeQ.size(); ++tq)
    time_qWeight[tq] = timedet * TimeQ.Weight(tq);
  for (int sq = 0; sq < SpaceQ.size(); ++sq)
    space_qWeight[sq] = T[sq].DetSpace() * SpaceQ.Weight(sq);

  for (int tq = 0; tq < TimeQ.size(); ++tq)
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      int index = tq * SpaceQ.size() + sq;
      qWeight[index] = space_qWeight[sq] * time_qWeight[tq];
      qPoint[index] = tc.LocalToGlobal(SpaceQ.QPoint(sq)).CopyWithT(a + timedet * TimeQ.QPoint(tq)[0]);
    }


  if (nTimeQ() != TS.size()) {
    Exit("Time: GLQuad should have same amount of qpoints as GLShape Nodalpoints.")
  }
  if (nSpaceQ() != SS.size()) {
    Exit("Space: GLQuad should have same amount of qpoints as GLShape Nodalpoints.")
  }
  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      int offset1 = (tq * nSpaceQ() + sq) * j_dimension();
      int offset2 = offset1 + tq * SS.size() * (1 + numL + dim);
      int offset = offset2 + sq;
      for (int k = 0; k <= numL; ++k) {
        int index = offset + (dim + k) * SS.size();
        pressureST[index] = 1;
      }
      for (int k = 0; k < dim; ++k) {
        int index = offset + k * SS.size();
        velocityST[index][k] = 1;
      }
      for (int si = 0; si < SS.size(); ++si) {
        int offset3 = offset2 + si;
        for (int k = 0; k < dim; ++k) {
          int index = offset3 + k * SS.size();
          divVelocityST[index] = gradient[sq][si][k];
        }
        for (int k = 0; k <= numL; ++k) {
          int index = offset3 + (dim + k) * SS.size();
          gradPressureST[index] = gradient[sq][si];
        }
      }
      for (int ti = 0; ti < TS.size(); ++ti) {
        int offset4 =
            (tq * nSpaceQ() + sq) * j_dimension() + ti * SS.size() * (1 + numL + dim) + sq;
        for (int k = 0; k < dim; ++k) {
          int index = offset4 + k * SS.size();
          dtVelocityST[index][k] = time_gradient[tq][ti][0];
        }
        for (int k = 0; k <= numL; ++k) {
          int index = offset4 + (dim + k) * SS.size();
          dtPressureST[index] = time_gradient[tq][ti][0];
        }
      }
    }
  }
}


void SpaceTimeViscoAcousticDGTFaceElementGLGL::updateCell(cell TC) {

  double a = TC.min();
  double b = TC.max();
  double timedet = b - a;
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int tq = 0; tq < TimeQ.size(); ++tq) {

      /** Transformation **/
      qPoint[tq * SpaceQ.size() + sq] = TC[faceqpoints[fid][sq]].
        CopyWithT(a + timedet * TimeQ.QPoint(tq)[0]);

      T[tq * SpaceQ.size() + sq] = TC.GetTransformation(faceqpoints[fid][sq]);
      qNormal[tq * SpaceQ.size() + sq] = T[tq * SpaceQ.size() + sq] * TC.LocalFaceNormal(fid);
      double w = norm(qNormal[tq * SpaceQ.size() + sq]);
      qWeight[tq * SpaceQ.size() + sq] = T[tq * SpaceQ.size() + sq].Det()
                                         * w * ww * SpaceQ.Weight(sq) * timedet * TimeQ.Weight(tq);
      qNormal[tq * SpaceQ.size() + sq] /= w;
    }
  }
}

SpaceTimeViscoAcousticDGTFaceElementGLGL::SpaceTimeViscoAcousticDGTFaceElementGLGL
    (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc, int fid, int nL)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(D.GetCellShape(deg.space)),
    TS(D.GetTimeShape(deg.time)),
    faceqpoints(D.FaceQPoints()),
    ww(tc.LocalFaceArea(fid)),
    fid(fid), dim(tc.dim()),
    SpaceQ(D.GetCellFaceQuad(fid, deg.space)),
    TimeQ(D.GetTimeFaceQuad(fid, deg.time)),
    numL(nL) {

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  T.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  pressureST.resize(nQ() * j_dimension());

  updateCell(tc);


  double spaceValue[SpaceQ.size()][SS.size()];
  static bool space_once[4] = {true, true, true, true};
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int si = 0; si < SS.size(); ++si) {
      spaceValue[sq][si] = SS(fid, sq, si);
    }
  }

  double timeValue[TimeQ.size()][TS.size()];
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int i = 0; i < TS.size(); ++i)
      timeValue[tq][i] = TS(TimeQ.QPoint(tq), i);
  }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          double value = spaceValue[sq][si] * timeValue[tq][ti];
          for (int k = 0; k < dim; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;

            velocityST[index][k] = value;

          }
          for (int k = 0; k <= numL; ++k) {

            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * SS.size() * (dim + 1 + numL) + (dim + k) * SS.size()
                        + si;
            pressureST[index] = value;
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    non_zero_indices_for_q.emplace_back(std::vector<int>{});
    for (int j = 0; j < j_dimension(); ++j) {
      j_zero[q * j_dimension() + j] =
          ((std::abs(Pressure(q, j)) == 0) && (norm(Velocity(q, j)) == 0));


    }
    for (int i = 0; i < i_dimension(); ++i) {
      i_zero[q * i_dimension() + i] =
          ((std::abs(PressureTestSpace(q, i)) == 0)
           && (norm(VelocityTestSpace(q, i)) == 0));
      if (!i_zero[q * i_dimension() + i]) {
        non_zero_indices_for_q[q].push_back(i);
      }
    }
  }
}


SpaceTimeViscoAcousticDGTFaceElementGLGL::SpaceTimeViscoAcousticDGTFaceElementGLGL
    (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc,
     int fid, int nL, string dgdg)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(D.GetCellShape(deg.space)),
    TS(D.GetTimeShape(deg.time)),
    faceqpoints(D.FaceQPoints()),
    ww(tc.LocalFaceArea(fid)),
    fid(fid), dim(tc.dim()),
    SpaceQ(D.GetCellQuad(deg.space)),
    TimeQ(GetQuadrature("Qint1")),
    numL(nL) {
  face f = g.GetMesh().find_face(tc.Face(fid));

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  T.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  pressureST.resize(nQ() * j_dimension());

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      T[tq * SpaceQ.size() + sq] = tc.GetTransformation(SpaceQ.QPoint(sq));
      /** Transformation **/
      qPoint[tq * SpaceQ.size() + sq] = tc[SpaceQ.QPoint(sq)].CopyWithT(f().t());
      qWeight[tq * SpaceQ.size() + sq] = T[sq].Det() * SpaceQ.Weight(sq);
    }
  }

  double spaceValue[SpaceQ.size()][SS.size()];
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int i = 0; i < SS.size(); ++i) spaceValue[sq][i] = SS(sq, i);
  }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {

        //velocity
        for (int k = 0; k < dim; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            double timeValue = 0.0;
            if (fid == tc.Faces() - 2) timeValue = TS(Point(0.0), ti);
            if (fid == tc.Faces() - 1) timeValue = TS(Point(1.0), ti);
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST[index][k] = spaceValue[sq][si] * timeValue;
          }
        }
        //pressure
        for (int k = 0; k <= numL; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            double timeValue = 0.0;
            if (fid == tc.Faces() - 2) timeValue = TS(Point(0.0), ti);
            if (fid == tc.Faces() - 1) timeValue = TS(Point(1.0), ti);
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * SS.size() * (dim + 1 + numL) + (dim + k) * SS.size()
                        + si;
            pressureST[index] = spaceValue[sq][si] * timeValue;
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    for (int j = 0; j < j_dimension(); ++j)
      j_zero[q * j_dimension() + j] = ((std::abs(Pressure(q, j)) == 0) &&
                                       (norm(Velocity(q, j)) == 0));
    for (int i = 0; i < i_dimension(); ++i) {
      i_zero[q * i_dimension() + i] = ((std::abs(PressureTestSpace(q, i)) == 0)
                                       && (norm(VelocityTestSpace(q, i)) == 0));
    }
  }
}

SpaceTimeViscoAcousticDGTFaceElementGLGL::SpaceTimeViscoAcousticDGTFaceElementGLGL
    (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc,
     int fid, int nL, string dgdg, const cell &tcf)
    :
    deg(g.GetDoF().GetDegree(*tc)),
    SS(D.GetCellShape(deg.space)),
    TS(D.GetTimeShape(deg.time)),
    faceqpoints(D.FaceQPoints()),
    ww(tc.LocalFaceArea(fid)),
    fid(fid), dim(tc.dim()),
    SpaceQ(D.GetCellQuad(g.GetDoF().GetDegree(*tcf).space)),
    TimeQ(GetQuadrature("Qint1")),
    numL(nL) {
  face f = g.GetMesh().find_face(tc.Face(fid));

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  T.resize(nQ());

  velocityST.resize(nQ() * j_dimension());
  pressureST.resize(nQ() * j_dimension());

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      T[tq * SpaceQ.size() + sq] = tc.GetTransformation(SpaceQ.QPoint(sq));
      /** Transformation **/
      qPoint[tq * SpaceQ.size() + sq] = tc[SpaceQ.QPoint(sq)].CopyWithT(f().t());
      qWeight[tq * SpaceQ.size() + sq] = T[sq].Det() * SpaceQ.Weight(sq);
    }
  }

  double spaceValue[SpaceQ.size()][SS.size()];
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    Point qpLoc;
    transformPointGlobalToLocal(qPoint[sq], qpLoc, tc);
    for (int i = 0; i < SS.size(); ++i) spaceValue[sq][i] = SS(qpLoc, i);
  }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {

        //velocity
        for (int k = 0; k < dim; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            double timeValue = 0.0;
            if (fid == tc.Faces() - 2) timeValue = TS(Point(0.0), ti);
            if (fid == tc.Faces() - 1) timeValue = TS(Point(1.0), ti);
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST[index][k] = spaceValue[sq][si] * timeValue;
          }
        }
        //pressure
        for (int k = 0; k <= numL; ++k) {
          for (int ti = 0; ti < TS.size(); ++ti) {
            double timeValue = 0.0;
            if (fid == tc.Faces() - 2) timeValue = TS(Point(0.0), ti);
            if (fid == tc.Faces() - 1) timeValue = TS(Point(1.0), ti);
            int index = (tq * nSpaceQ() + sq) * j_dimension()
                        + ti * SS.size() * (dim + 1 + numL) + (dim + k) * SS.size()
                        + si;
            pressureST[index] = spaceValue[sq][si] * timeValue;
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    for (int j = 0; j < j_dimension(); ++j)
      j_zero[q * j_dimension() + j] = ((std::abs(Pressure(q, j)) == 0) &&
                                       (norm(Velocity(q, j)) == 0));
    for (int i = 0; i < i_dimension(); ++i) {
      i_zero[q * i_dimension() + i] = ((std::abs(PressureTestSpace(q, i)) == 0)
                                       && (norm(VelocityTestSpace(q, i)) == 0));
    }
  }
}
