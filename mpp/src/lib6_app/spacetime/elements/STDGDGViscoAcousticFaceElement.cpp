#include "STDGDGViscoAcousticFaceElement.hpp"

SpaceTimeViscoAcousticDGTFaceElement::SpaceTimeViscoAcousticDGTFaceElement
  (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc, int f_id, int nL)
  :
  deg(g.GetDoF().GetDegree(tc())),
  shape(GetSpaceTimeShape(g, *tc)),
  faceqpoints(D.FaceQPoints()),
  TC(tc), f(g.find_face(tc.Face(f_id))),  r(g.find_row(tc())),
  fid(f_id), dim(tc.dim()),
  SpaceQ(D.GetCellFaceQuad(f_id, deg.space)),
  TimeQ(D.GetTimeFaceQuad(f_id, deg.time)),
  numL(nL) {

  double ww = TC.SpaceCell().LocalFaceArea(fid);
  double a = TC.min();
  double b = TC.max();
  det_time = b - a;

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  velocityST.resize(nQ() * j_dimension());
  pressureST.resize(nQ() * j_dimension());

  const Transformation &T = tc.SpaceCell().GetTransformation(faceqpoints[fid][0]);

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      Point localQPoint = faceqpoints[fid][sq].CopyWithT(TimeQ.QPoint(tq)[0]);
      qPoint[tq * SpaceQ.size() + sq] = TC.LocalToGlobal(localQPoint);

      qNormal[tq * SpaceQ.size() + sq] = T * TC.SpaceCell().LocalFaceNormal(fid);
      double w = norm(qNormal[tq * SpaceQ.size() + sq]);
      qWeight[tq * SpaceQ.size() + sq] = T.Det()
        * w * ww * SpaceQ.Weight(sq) * det_time * TimeQ.Weight(tq);
      qNormal[tq * SpaceQ.size() + sq] /= w;
    }
  }
  const Shape &SS = shape.GetSpaceShape();
  ShapeValues<double> space_values(SpaceQ.size(), SS.size());
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int i = 0; i < SS.size(); ++i) {
      space_values[sq][i] = SS(fid, sq, i);
    }
  }
  const Shape &TS = shape.GetTimeShape();
  ShapeValues<double> time_values(TimeQ.size(), TS.size());
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int i = 0; i < TS.size(); ++i) {
      time_values[tq][i] = TS(TimeQ.QPoint(tq), i);
    }
  }

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          for (int k = 0; k < dim; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
              + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST[index][k] = space_values[sq][si] * time_values[tq][ti];
          }
          for (int k = 0; k <= numL; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
              + ti * SS.size() * (dim + 1 + numL) + (dim + k) * SS.size()
              + si;
            pressureST[index] = space_values[sq][si] * time_values[tq][ti];
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    for (int j = 0; j < j_dimension(); ++j)
      j_zero[q * j_dimension() + j] =
        ((std::abs(Pressure(q, j)) == 0) && (norm(Velocity(q, j)) == 0));
    for (int i = 0; i < i_dimension(); ++i)
      i_zero[q * i_dimension() + i] =
        ((std::abs(Pressure(q, i)) == 0) && (norm(Velocity(q, i)) == 0));
  }
}

SpaceTimeViscoAcousticDGTFaceElement::SpaceTimeViscoAcousticDGTFaceElement
  (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc,
   int f_id, int nL, string dgdg)
  :
  deg(g.GetDoF().GetDegree(*tc)),
  shape(GetSpaceTimeShape(g, *tc)),
  faceqpoints(D.FaceQPoints()),
  TC(tc), f(g.find_face(tc.Face(f_id))),  r(g.find_row(tc())),
  fid(f_id), dim(tc.dim()),
  SpaceQ(D.GetCellQuad(deg.space)),
  TimeQ(GetQuadrature("Qint1")),
  numL(nL) {

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  velocityST.resize(nQ() * j_dimension());
  pressureST.resize(nQ() * j_dimension());

  const Transformation &T = tc.SpaceCell().GetTransformation(SpaceQ.QPoint(0));

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      qPoint[tq * SpaceQ.size() + sq] = TC.SpaceCell().
        LocalToGlobal(SpaceQ.QPoint(sq)).
        CopyWithT(f().t());
      qWeight[tq * SpaceQ.size() + sq] = T.Det() * SpaceQ.Weight(sq);
    }
  }

  const Shape &SS = shape.GetSpaceShape();
  ShapeValues<double> space_values(SpaceQ.size(), SS.size());
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    for (int i = 0; i < SS.size(); ++i) space_values[sq][i] = SS(sq, i);
  }
  const Shape &TS = shape.GetTimeShape();

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          double timeValue = 0.0;
          if (fid == TC.Faces() - 2) timeValue = TS(Point(0.0), ti);
          if (fid == TC.Faces() - 1) timeValue = TS(Point(1.0), ti);
          for (int k = 0; k < dim; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
              + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST[index][k] = space_values[sq][si] * timeValue;
          }
          for (int k = 0; k <= numL; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
              + ti * SS.size() * (dim + 1 + numL) + (dim + k) * SS.size()
              + si;
            pressureST[index] = space_values[sq][si] * timeValue;
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
    for (int i = 0; i < i_dimension(); ++i)
      i_zero[q * i_dimension() + i] = ((std::abs(Pressure(q, i)) == 0)
        && (norm(Velocity(q, i)) == 0));
  }
}

SpaceTimeViscoAcousticDGTFaceElement::SpaceTimeViscoAcousticDGTFaceElement
  (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc,
   int f_id, int nL, string dgdg, const cell &tcf)
  :
  deg(g.GetDoF().GetDegree(*tcf)),
  shape(GetSpaceTimeShape(g, *tc)),
  faceqpoints(D.FaceQPoints()),
  TC(tc), f(g.find_face(tc.Face(f_id))), r(g.find_row(tc())),
  fid(f_id), dim(tc.dim()),
  SpaceQ(D.GetCellQuad(deg.space)),
  TimeQ(GetQuadrature("Qint1")),
  numL(nL) {

  //Assert(tc.Faces() - 1 == f_id);
  //Assert(tc.Face(tc.Faces() - 1) == tcf.Face(tcf.Faces() - 2));

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  velocityST.resize(nQ() * j_dimension());
  pressureST.resize(nQ() * j_dimension());

  const Transformation &T = tc.SpaceCell().GetTransformation(SpaceQ.QPoint(0));

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      qPoint[tq * SpaceQ.size() + sq] = TC.SpaceCell().
        LocalToGlobal(SpaceQ.QPoint(sq)).
        CopyWithT(f().t());
      qWeight[tq * SpaceQ.size() + sq] = T.Det() * SpaceQ.Weight(sq);
    }
  }

  const Shape &SS = shape.GetSpaceShape();
  ShapeValues<double> space_values(SpaceQ.size(), SS.size());
  for (int sq = 0; sq < SpaceQ.size(); ++sq) {
    Point qpLoc = TC.GlobalToLocal(qPoint[sq]);
    for (int i = 0; i < SS.size(); ++i) space_values[sq][i] = SS(qpLoc, i);
  }
  const Shape &TS = shape.GetTimeShape();

  for (int tq = 0; tq < nTimeQ(); ++tq) {
    for (int sq = 0; sq < nSpaceQ(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          double timeValue = 0.0;
          if (fid == TC.Faces() - 2) timeValue = TS(Point(0.0), ti);
          if (fid == TC.Faces() - 1) timeValue = TS(Point(1.0), ti);
          for (int k = 0; k < dim; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
              + ti * SS.size() * (1 + numL + dim) + k * SS.size() + si;
            velocityST[index][k] = space_values[sq][si] * timeValue;
          }

          for (int k = 0; k <= numL; ++k) {
            int index = (tq * nSpaceQ() + sq) * j_dimension()
              + ti * SS.size() * (dim + 1 + numL) + (dim + k) * SS.size()
              + si;
            pressureST[index] = space_values[sq][si] * timeValue;
          }
        }
      }
    }
  }

  j_zero.resize(nQ() * j_dimension());
  i_zero.resize(nQ() * i_dimension());
  for (int q = 0; q < nQ(); ++q) {
    for (int j = 0; j < j_dimension(); ++j) {
      j_zero[q * j_dimension() + j] = ((std::abs(Pressure(q, j)) == 0) &&
        (norm(Velocity(q, j)) == 0));
    }
    for (int i = 0; i < i_dimension(); ++i) {
      i_zero[q * i_dimension() + i] = ((std::abs(Pressure(q, i)) == 0)
        && (norm(Velocity(q, i)) == 0));
    }
  }
}

const Point &SpaceTimeViscoAcousticDGTFaceElement::QLocal(int q) const { return qLocal[q]; }

const Point &SpaceTimeViscoAcousticDGTFaceElement::QPoint(int q) const { return qPoint[q]; }

const Point &SpaceTimeViscoAcousticDGTFaceElement::QNormal(int q) const { return qNormal[q]; }

double SpaceTimeViscoAcousticDGTFaceElement::QWeight(int q) const { return qWeight[q]; }

int SpaceTimeViscoAcousticDGTFaceElement::nSpaceQ() const { return SpaceQ.size(); }

int SpaceTimeViscoAcousticDGTFaceElement::nTimeQ() const { return TimeQ.size(); }

int SpaceTimeViscoAcousticDGTFaceElement::nQ() const { return SpaceQ.size() * TimeQ.size(); }

bool SpaceTimeViscoAcousticDGTFaceElement::is_zero_test(int nq, int i) const { return i_zero[nq * i_dimension() + i]; }

bool SpaceTimeViscoAcousticDGTFaceElement::is_zero(int nq, int j) const { return j_zero[nq * j_dimension() + j]; }

int SpaceTimeViscoAcousticDGTFaceElement::i_dimension() const { return shape.size() * (1 + dim + numL); }

int SpaceTimeViscoAcousticDGTFaceElement::j_dimension() const { return shape.size() * (1 + dim + numL); }

int SpaceTimeViscoAcousticDGTFaceElement::variable(int j) const {
  const Shape &SS = shape.GetSpaceShape();
  int tmp = j % (SS.size() * (dim + 1 + numL));
  if (tmp < dim * SS.size()) return 0;
  for (int nL = 0; nL <= numL; nL++)
    if (tmp - SS.size() * dim < (nL + 1) * SS.size()) return 1 + nL;
  Exit("BLUPP");
}

VectorField SpaceTimeViscoAcousticDGTFaceElement::Velocity(int nq, int j) {
  return velocityST[nq * j_dimension() + j];
}

double SpaceTimeViscoAcousticDGTFaceElement::Pressure(int nq, int j) {
  return pressureST[nq * j_dimension() + j];
}

double SpaceTimeViscoAcousticDGTFaceElement::Pressure(int nq, const double *u) {
  double P = 0;
  for (int j = 0; j < j_dimension(); ++j)
    P += u[j] * Pressure(nq, j);
  return P;
}

VectorField SpaceTimeViscoAcousticDGTFaceElement::Velocity(Point QP, int i) {
  int var = variable(i);
  if (var != 0) return zero;
  const Shape &SS = shape.GetSpaceShape();
  const Shape &TS = shape.GetTimeShape();
  int ti = i / (SS.size() * (1 + dim + numL));
  int si = i % (SS.size() * (1 + dim + numL));
  int ssi = si / SS.size();
  int ssj = si % SS.size();

  VectorField V = zero;
  V[ssi] = 1.0;
  V *= SS(QP, ssj) * TS(Point(QP.t()), ti);
  return V;
}

VectorField SpaceTimeViscoAcousticDGTFaceElement::Velocity(int nq, const double *u) {
  VectorField V = zero;
  for (int j = 0; j < j_dimension(); ++j)
    V += u[j] * Velocity(nq, j);
  return V;
}

double SpaceTimeViscoAcousticDGTFaceElement::Pressure(Point QP, int i) {
  int var = variable(i);
  if (var == 0) return 0.0;
  const Shape &SS = shape.GetSpaceShape();
  const Shape &TS = shape.GetTimeShape();
  int ti = i / (SS.size() * (1 + dim + numL));
  int si = i % SS.size();

  Scalar P = SS(QP, si) * TS(Point(QP.t()), ti);
  return P;
}



