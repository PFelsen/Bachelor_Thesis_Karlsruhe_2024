#include "STDGDGTransportFaceElement.hpp"

const SpaceTimeShape &GetSpaceTimeShape(const VectorMatrixBase &base, const Cell &c);

SpaceTimeTransportDGTFaceElement::SpaceTimeTransportDGTFaceElement
    (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc, int f_id) :
    deg(g.GetDoF().GetDegree(tc())),
    shape(GetSpaceTimeShape(g, *tc)),
    faceqpoints(D.FaceQPoints()),
    TC(tc), f(g.find_face(tc.Face(f_id))), r(g.find_row(tc())),
    fid(f_id), dim(tc.dim()),
    SpaceQ(D.GetCellFaceQuad(f_id, deg.space)),
    TimeQ(D.GetTimeFaceQuad(f_id, deg.time)) {

  double ww = TC.SpaceCell().LocalFaceArea(fid);
  double a = TC.min();
  double b = TC.max();
  double det_time = b - a;

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  densityST.resize(nQ() * j_dimension());

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
  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          int index = (tq * SpaceQ.size() + sq) * j_dimension() + ti * SS.size() + si;
          densityST[index] = space_values[sq][si] * time_values[tq][ti];
        }
      }
    }
  }
}

SpaceTimeTransportDGTFaceElement::SpaceTimeTransportDGTFaceElement
    (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc, int f_id,
     string dgdg) :
    deg(g.GetDoF().GetDegree(*tc)),
    shape(GetSpaceTimeShape(g, *tc)),
    faceqpoints(D.FaceQPoints()),
    TC(tc), f(g.find_face(tc.Face(f_id))), r(g.find_row(tc())),
    fid(f_id), dim(tc.dim()),
    SpaceQ(D.GetCellQuad(deg.space)),
    TimeQ(GetQuadrature("Qint1")) {

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  densityST.resize(nQ() * j_dimension());

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

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          int index = (tq * SpaceQ.size() + sq) * j_dimension() + ti * SS.size() + si;
          double timeValue = 0.0;
          if (fid == TC.Faces() - 2) timeValue = TS(Point(0.0), ti);
          if (fid == TC.Faces() - 1) timeValue = TS(Point(1.0), ti);
          densityST[index] = space_values[sq][si] * timeValue;
        }
      }
    }
  }
}

SpaceTimeTransportDGTFaceElement::SpaceTimeTransportDGTFaceElement
    (const STDiscretization &D, const VectorMatrixBase &g, const cell &tc, int f_id,
     string dgdg, const cell &tcf) :
    deg(g.GetDoF().GetDegree(*tcf)),
    shape(GetSpaceTimeShape(g, *tc)),
    faceqpoints(D.FaceQPoints()),
    TC(tc), f(g.find_face(tc.Face(f_id))), r(g.find_row(tc())),
    fid(f_id), dim(tc.dim()),
    SpaceQ(D.GetCellQuad(deg.space)),
    TimeQ(GetQuadrature("Qint1")) {

  //Assert(tc.Faces() - 1 == f_id);
  //Assert(tc.Face(tc.Faces() - 1) == tcf.Face(tcf.Faces() - 2));

  qWeight.resize(nQ());
  qPoint.resize(nQ());
  qLocal.resize(nQ());
  qNormal.resize(nQ());
  densityST.resize(nQ() * j_dimension());

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

  for (int tq = 0; tq < TimeQ.size(); ++tq) {
    for (int sq = 0; sq < SpaceQ.size(); ++sq) {
      for (int si = 0; si < SS.size(); ++si) {
        for (int ti = 0; ti < TS.size(); ++ti) {
          int index = (tq * SpaceQ.size() + sq) * j_dimension() + ti * SS.size() + si;
          double timeValue = 0.0;
          if (fid == TC.Faces() - 2) timeValue = TS(Point(0.0), ti);
          if (fid == TC.Faces() - 1) timeValue = TS(Point(1.0), ti);
          densityST[index] = space_values[sq][si] * timeValue;
        }
      }
    }
  }
}
