#ifndef SPACETIME_SPACETIMETRANSFER_HPP
#define SPACETIME_SPACETIMETRANSFER_HPP

#include "Algebra.hpp"

#include "SpaceTimeTools.hpp"
#include "Parallel.hpp"
#include "Functools.hpp"


class SpaceTimeTransfer {
private:

public:
  virtual void Restrict(const Vector &fine, Vector &coarse) const = 0;

  virtual void Prolongate(Vector &fine, const Vector &coarse) const = 0;

  virtual ~SpaceTimeTransfer() {}
};

class SpaceTimeTransferProjectionDGDG : public SpaceTimeTransfer {

  const STDiscretization &disc;
  const int k_max;
public:
  SpaceTimeTransferProjectionDGDG(const STDiscretization &disc) :
      disc(disc), k_max(disc().GetDoF().get_dim()) {};

  void ProlongateOnSameLevel(const Vector &from, Vector &to) const {
    if (from.GetMesh().CellCountGeometry() != to.GetMesh().CellCountGeometry()) {
      THROW("Different cell-counts in same-level transfer!");
    }
    mout << "ProlongateOnSameLevel" << DOUT(from.Level()) << DOUT(to.Level()) << endl;
    to = 0.0;
    std::unordered_map<DegreePair, RMatrix> InvertedMassMatrixCache_t;
    std::unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;
    for (cell c_f = from.cells(); c_f != from.cells_end(); c_f++) {
      row row_f = from.find_row(c_f());
      DegreePair deg_f = from.GetDoF().GetDegree(*c_f);
      auto &space_shape_f = disc.GetCellShape(deg_f.space);
      auto &time_shape_f = disc.GetTimeShape(deg_f.time);
      int N_f = time_shape_f.size() * k_max * space_shape_f.size();

      cell c_t = to.find_cell(c_f());
      row row_t = to.find_row(c_t());
      DegreePair deg_t = to.GetDoF().GetDegree(*c_t);
      auto &space_shape_t = disc.GetCellShape(deg_t.space);
      auto &time_shape_t = disc.GetTimeShape(deg_t.time);
      int N_t = time_shape_t.size() * k_max * space_shape_t.size();

      SpaceTimeCellQuadrature quad_f(c_f, disc.GetCellQuad(deg_f.space),
                                     disc.GetTimeQuad(deg_f.time));

      SpaceTimeCellQuadrature quad_t(c_t, disc.GetCellQuad(deg_t.space),
                                     disc.GetTimeQuad(deg_t.time));


      auto elem = InvertedMassMatrixCache_t.find(deg_t);

      if (elem == InvertedMassMatrixCache_t.end()) {
        RMatrix InvertedMassMatrix(N_t, N_t);
        for (int q = 0; q < quad_t.size(); q++) {
          Point local = quad_t.qPointLocal[q];
          for (int tc = 0; tc < time_shape_t.size(); ++tc) {
            double ts_t = time_shape_t(Point(local.t()), tc);
            for (int k = 0; k < k_max; ++k) {
              for (int sc = 0; sc < space_shape_t.size(); ++sc) {
                double ss_t = space_shape_t(local, sc);
                double v_t = ss_t * ts_t;
                int row = sc + (k + tc * k_max) * space_shape_t.size();
                for (int tc2 = 0; tc2 < time_shape_t.size(); ++tc2) {
                  double ts_t2 = time_shape_t(Point(local.t()), tc2);
                  for (int sc2 = 0; sc2 < space_shape_t.size(); ++sc2) {
                    int col = sc2 + (k + tc2 * k_max) * space_shape_t.size();
                    double ss_t2 = space_shape_t(local, sc2);
                    double v_t2 = ss_t2 * ts_t2;
                    InvertedMassMatrix(row, col) += quad_t.qWeight[q] * v_t * v_t2;
                  }
                }
              }
            }
          }
        }
        InvertedMassMatrix.Invert();
        elem = InvertedMassMatrixCache_t.insert({deg_t, InvertedMassMatrix}).first;
      }
      RMatrix InvertedMassMatrix_t = elem->second;

      RMatrix MassMatrix_ft(N_t, N_f);
      for (int q = 0; q < quad_t.size(); q++) {
        Point local = quad_t.qPointLocal[q];
        for (int tt = 0; tt < time_shape_t.size(); ++tt) {
          double value_time_shape = time_shape_t(Point(local.t()), tt);
          for (int st = 0; st < space_shape_t.size(); ++st) {
            double value_space_shape = space_shape_t(local, st);
            double v_t = value_space_shape * value_time_shape;
            for (int tf = 0; tf < time_shape_f.size(); ++tf) {
              double value_time_shape_f = time_shape_f(Point(local.t()), tf);
              for (int sf = 0; sf < space_shape_f.size(); ++sf) {
                double value_space_shape_f = space_shape_f(local, sf);
                for (int k = 0; k < k_max; ++k) {
                  int t_idx = st + (k + tt * k_max) * space_shape_t.size();
                  int f_idx = sf + (k + tf * k_max) * space_shape_f.size();
                  double v_f = value_space_shape_f * value_time_shape_f;
                  MassMatrix_ft[t_idx][f_idx] += quad_t.qWeight[q] * v_t * v_f;
                }
              }
            }
          }
        }
      }

      RVector value_f(N_f);
      for (int i = 0; i < value_f.size(); i++) {
        value_f[i] = from(row_f)[i];
      }
      //return;
      /*pout << DOUT(value_f.size()) << endl;
      pout << DOUT(MassMatrix_ft.rows()) << endl;
      pout << DOUT(MassMatrix_ft.cols()) << endl;
      pout << DOUT(InvertedMassMatrix_t.rows()) << endl;
      pout << DOUT(deg_f.space) << DOUT(deg_f.time) << endl;
      pout << DOUT(deg_t.space) << DOUT(deg_t.time) << endl;*/

      RVector integral_t(MassMatrix_ft * value_f);
      RVector value_t = InvertedMassMatrix_t * integral_t;

      for (int i = 0; i < value_t.size(); i++) {
        to(row_t)[i] += value_t[i];
      }
    }
    //mout << from << endl;
    //mout << to << endl;
    to.MakeAdditive();
    to.Accumulate();
  }


  void RestrictOnSameLevel(const Vector &F, Vector &C) const {
    std::unordered_map<DegreePair, RMatrix> MassMatrixCache_f;
    std::unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;

    for (cell cc = C.cells(); cc != C.cells_end(); cc++) {
      row rc = C.find_row(cc());
      DegreePair deg_c = C.GetDoF().GetDegree(*cc);
      auto &space_shape_c = disc.GetCellShape(deg_c.space);
      auto &time_shape_c = disc.GetTimeShape(deg_c.time);

      int N_c = time_shape_c.size() * k_max * space_shape_c.size();


      RVector integral_c(N_c);
      row row_f = F.find_row(cc());

      DegreePair deg_f = F.GetDoF().GetDegree(*cc);
      auto elem_q = SpaceTimeQuadCache.find(deg_f);
      if (elem_q == SpaceTimeQuadCache.end()) {
        const Quadrature &s_quad = disc.GetCellQuad(deg_f.space);
        const Quadrature &t_quad = disc.GetTimeQuad(deg_f.time);
        elem_q = SpaceTimeQuadCache.insert(
            {deg_f, SpaceTimeCellQuadrature(cc, s_quad, t_quad)}).first;
      }
      SpaceTimeCellQuadrature &quad_f = elem_q->second;

      auto &space_shape_f = disc.GetCellShape(deg_f.space);
      auto &time_shape_f = disc.GetTimeShape(deg_f.time);

      int N_f = time_shape_f.size() * k_max * space_shape_f.size();

      // ################################
      auto MassMatrixIterator_f = MassMatrixCache_f.find(deg_f);

      if (MassMatrixIterator_f == MassMatrixCache_f.end()) {
        RMatrix InvertedMassMatrix_f(N_f, N_f);

        for (int q = 0; q < quad_f.size(); q++) {
          Point local_f = quad_f.qPointLocal[q];
          for (int tc = 0; tc < time_shape_f.size(); ++tc) {
            double value_time_shape = time_shape_f(Point(local_f.t()), tc);
            for (int k = 0; k < k_max; ++k) {
              for (int sc = 0; sc < space_shape_f.size(); ++sc) {
                double value_space_shape = space_shape_f(local_f, sc);
                double v1 = value_space_shape * value_time_shape;
                int row = sc + (k + tc * k_max) * space_shape_f.size();
                for (int tc2 = 0; tc2 < time_shape_f.size(); ++tc2) {
                  double value_time_shape2 = time_shape_f(Point(local_f.t()), tc2);
                  for (int sc2 = 0; sc2 < space_shape_f.size(); ++sc2) {
                    int col = sc2 + (k + tc2 * k_max) * space_shape_f.size();
                    double value_space_shape2 = space_shape_f(local_f, sc2);
                    double v2 = value_space_shape2 * value_time_shape2;
                    InvertedMassMatrix_f[row][col] += quad_f.qWeight[q] * v1 * v2;
                  }
                }
              }
            }
          }
        }
        InvertedMassMatrix_f.Invert();
        MassMatrixIterator_f = MassMatrixCache_f.insert({deg_f, InvertedMassMatrix_f}).first;
      }
      RMatrix &InvertedMassMatrix_f = MassMatrixIterator_f->second;

      //###############################



      RMatrix MassMatrix_fc(N_f, N_c);

      for (int q = 0; q < quad_f.size(); q++) {
        Point local_f = quad_f.qPointLocal[q];
        Point global = cc.LocalToGlobal(local_f);
        Point local_c = cc.GlobalToLocal(global);
        for (int tc = 0; tc < time_shape_c.size(); ++tc) {
          double value_time_shape = time_shape_c(Point(local_c.t()), tc);
          for (int sc = 0; sc < space_shape_c.size(); ++sc) {
            double value_space_shape = space_shape_c(local_c, sc);
            double v_coarse = value_space_shape * value_time_shape;
            for (int tf = 0; tf < time_shape_f.size(); ++tf) {
              double value_time_shape_f = time_shape_f(Point(local_f.t()), tf);
              for (int sf = 0; sf < space_shape_f.size(); ++sf) {
                double value_space_shape_f = space_shape_f(local_f, sf);
                for (int k = 0; k < k_max; ++k) {
                  int c_idx = sc + (k + tc * k_max) * space_shape_c.size();
                  int f_idx = sf + (k + tf * k_max) * space_shape_f.size();
                  double v_fine = value_space_shape_f * value_time_shape_f;
                  MassMatrix_fc[c_idx][f_idx] += quad_f.qWeight[q] * v_coarse * v_fine;
                }
              }
            }
          }
        }
      }
      RVector integral_f(N_f);
      for (int i = 0; i < integral_f.size(); i++) {
        integral_f[i] = F(row_f)[i];
      }
      RVector value_f(InvertedMassMatrix_f * integral_f);
      integral_c += MassMatrix_fc * value_f;

      for (int i = 0; i < integral_c.size(); i++) {
        C(rc)[i] = integral_c[i];
      }
    }
  }

  void Restrict(const Vector &F, Vector &C) const override {
    TransferInfo info(F, C);
    C.Clear();
    if (!info.in_space && !info.in_time) {
      RestrictOnSameLevel(F, C);
      return;
    }

    //std::unordered_map<DegreePairProduct, RMatrix> InvertedMassMatrixCache;
    std::unordered_map<DegreePair, RMatrix> MassMatrixCache_f;
    std::unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;

    for (cell cc = C.cells(); cc != C.cells_end(); cc++) {
      row rc = C.find_row(cc());
      DegreePair deg_c = C.GetDoF().GetDegree(*cc);
      auto &space_shape_c = disc.GetCellShape(deg_c.space);
      auto &time_shape_c = disc.GetTimeShape(deg_c.time);

      int N_c = time_shape_c.size() * k_max * space_shape_c.size();

      vector<Point> Children(info.GetChildren(cc));
      RVector integral_c(N_c);
      for (int i = 0; i < Children.size(); ++i) {
        row row_f = F.find_row(Children[i]);

        DegreePair deg_f = C.GetDoF().GetDegree(*cc);
        cell child_cell = F.find_cell(Children[i]);
        auto elem_q = SpaceTimeQuadCache.find(deg_f);
        if (elem_q == SpaceTimeQuadCache.end()) {
          const Quadrature &s_quad = disc.GetCellQuad(deg_f.space);
          const Quadrature &t_quad = disc.GetTimeQuad(deg_f.time);
          elem_q = SpaceTimeQuadCache.insert(
              {deg_f, SpaceTimeCellQuadrature(child_cell, s_quad, t_quad)}).first;
        }
        SpaceTimeCellQuadrature &quad_f = elem_q->second;

        auto &space_shape_f = disc.GetCellShape(deg_f.space);
        auto &time_shape_f = disc.GetTimeShape(deg_f.time);

        int N_f = time_shape_f.size() * k_max * space_shape_f.size();

        // ################################
        auto MassMatrixIterator_f = MassMatrixCache_f.find(deg_f);

        if (MassMatrixIterator_f == MassMatrixCache_f.end()) {
          RMatrix InvertedMassMatrix_f(N_f, N_f);

          for (int q = 0; q < quad_f.size(); q++) {
            Point local_f = quad_f.qPointLocal[q];
            for (int tc = 0; tc < time_shape_f.size(); ++tc) {
              double value_time_shape = time_shape_f(Point(local_f.t()), tc);
              for (int k = 0; k < k_max; ++k) {
                for (int sc = 0; sc < space_shape_f.size(); ++sc) {
                  double value_space_shape = space_shape_f(local_f, sc);
                  double v1 = value_space_shape * value_time_shape;
                  int row = sc + (k + tc * k_max) * space_shape_f.size();
                  for (int tc2 = 0; tc2 < time_shape_f.size(); ++tc2) {
                    double value_time_shape2 = time_shape_f(Point(local_f.t()), tc2);
                    for (int sc2 = 0; sc2 < space_shape_f.size(); ++sc2) {
                      int col = sc2 + (k + tc2 * k_max) * space_shape_f.size();
                      double value_space_shape2 = space_shape_f(local_f, sc2);
                      double v2 = value_space_shape2 * value_time_shape2;
                      InvertedMassMatrix_f[row][col] += quad_f.qWeight[q] * v1 * v2;
                    }
                  }
                }
              }
            }
          }
          InvertedMassMatrix_f.Invert();
          MassMatrixIterator_f = MassMatrixCache_f.insert({deg_f, InvertedMassMatrix_f}).first;
        }
        RMatrix &InvertedMassMatrix_f = MassMatrixIterator_f->second;

        //###############################



        RMatrix MassMatrix_fc(N_f, N_c);

        for (int q = 0; q < quad_f.size(); q++) {
          Point local_f = quad_f.qPointLocal[q];
          Point global = child_cell.LocalToGlobal(local_f);
          Point local_c = cc.GlobalToLocal(global);
          for (int tc = 0; tc < time_shape_c.size(); ++tc) {
            double value_time_shape = time_shape_c(Point(local_c.t()), tc);
            for (int sc = 0; sc < space_shape_c.size(); ++sc) {
              double value_space_shape = space_shape_c(local_c, sc);
              double v_coarse = value_space_shape * value_time_shape;
              for (int tf = 0; tf < time_shape_f.size(); ++tf) {
                double value_time_shape_f = time_shape_f(Point(local_f.t()), tf);
                for (int sf = 0; sf < space_shape_f.size(); ++sf) {
                  double value_space_shape_f = space_shape_f(local_f, sf);
                  for (int k = 0; k < k_max; ++k) {
                    int c_idx = sc + (k + tc * k_max) * space_shape_c.size();
                    int f_idx = sf + (k + tf * k_max) * space_shape_f.size();
                    double v_fine = value_space_shape_f * value_time_shape_f;
                    MassMatrix_fc[c_idx][f_idx] += quad_f.qWeight[q] * v_coarse * v_fine;
                  }
                }
              }
            }
          }
        }
        RVector integral_f(N_f);
        for (int i = 0; i < integral_f.size(); i++) {
          integral_f[i] = F(row_f)[i];
        }
        RVector value_f(InvertedMassMatrix_f * integral_f);
        integral_c += MassMatrix_fc * value_f;
      }
      for (int i = 0; i < integral_c.size(); i++) {
        C(rc)[i] = integral_c[i];
      }
    }
  }

  void Prolongate(Vector &F, const Vector &C) const override {
    F = 0.0;
    TransferInfo info(F, C);
    if (!info.in_space && !info.in_time) {
      ProlongateOnSameLevel(C, F);
      return;
    }


    std::unordered_map<DegreePair, RMatrix> InvertedMassMatrixCache_f;
    std::unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;
    for (cell c_c = C.cells(); c_c != C.cells_end(); c_c++) {
      row row_c = C.find_row(c_c());
      DegreePair deg_c = C.GetDoF().GetDegree(*c_c);
      auto &space_shape_c = disc.GetCellShape(deg_c.space);
      auto &time_shape_c = disc.GetTimeShape(deg_c.time);

      int N_c = time_shape_c.size() * k_max * space_shape_c.size();

      vector<Point> Children(info.GetChildren(c_c));
      for (int i = 0; i < Children.size(); ++i) {
        row row_f = F.find_row(Children[i]);
        DegreePair deg_f = F.GetDoF().GetDegree(Children[i]);
        auto &space_shape_f = disc.GetCellShape(deg_f.space);
        auto &time_shape_f = disc.GetTimeShape(deg_f.time);
        int N_f = time_shape_f.size() * k_max * space_shape_f.size();
        cell c_f = F.find_cell(Children[i]);
        auto elem_q = SpaceTimeQuadCache.find(deg_f);
        if (elem_q == SpaceTimeQuadCache.end()) {
          const Quadrature &s_quad_f = disc.GetCellQuad(deg_f.space);
          const Quadrature &t_quad_f = disc.GetTimeQuad(deg_f.time);
          elem_q = SpaceTimeQuadCache.insert(
              {deg_f, SpaceTimeCellQuadrature(c_f, s_quad_f, t_quad_f)}).first;
        }
        SpaceTimeCellQuadrature &quad_f = elem_q->second;

        auto elem = InvertedMassMatrixCache_f.find({deg_f});

        if (elem == InvertedMassMatrixCache_f.end()) {
          RMatrix InvertedMassMatrix(N_f, N_f);
          for (int q = 0; q < quad_f.size(); q++) {
            Point local_f = quad_f.qPointLocal[q];
            for (int tc = 0; tc < time_shape_f.size(); ++tc) {
              double ts_f = time_shape_f(Point(local_f.t()), tc);
              for (int k = 0; k < k_max; ++k) {
                for (int sc = 0; sc < space_shape_f.size(); ++sc) {
                  double ss_f = space_shape_f(local_f, sc);
                  double v_f = ss_f * ts_f;
                  int row = sc + (k + tc * k_max) * space_shape_f.size();
                  for (int tc2 = 0; tc2 < time_shape_f.size(); ++tc2) {
                    double ts_f2 = time_shape_f(Point(local_f.t()), tc2);
                    for (int sc2 = 0; sc2 < space_shape_f.size(); ++sc2) {
                      int col = sc2 + (k + tc2 * k_max) * space_shape_f.size();
                      double ss_f2 = space_shape_f(local_f, sc2);
                      double v_f2 = ss_f2 * ts_f2;
                      InvertedMassMatrix(row, col) += quad_f.qWeight[q] * v_f * v_f2;
                    }
                  }
                }
              }
            }
          }
          InvertedMassMatrix.Invert();
          elem = InvertedMassMatrixCache_f.insert({deg_f, InvertedMassMatrix}).first;
        }
        RMatrix InvertedMassMatrix_f = elem->second;

        RMatrix MassMatrix_cf(N_f, N_c);
        for (int q = 0; q < quad_f.size(); q++) {
          Point local_f = quad_f.qPointLocal[q];
          Point global = c_f.LocalToGlobal(local_f);
          Point local_c = c_c.GlobalToLocal(global);
          for (int tc = 0; tc < time_shape_c.size(); ++tc) {
            double value_time_shape = time_shape_c(Point(local_c.t()), tc);
            for (int sc = 0; sc < space_shape_c.size(); ++sc) {
              double value_space_shape = space_shape_c(local_c, sc);
              double v_coarse = value_space_shape * value_time_shape;
              for (int tf = 0; tf < time_shape_f.size(); ++tf) {
                double value_time_shape_f = time_shape_f(Point(local_f.t()), tf);
                for (int sf = 0; sf < space_shape_f.size(); ++sf) {
                  double value_space_shape_f = space_shape_f(local_f, sf);
                  for (int k = 0; k < k_max; ++k) {
                    int c_idx = sc + (k + tc * k_max) * space_shape_c.size();
                    int f_idx = sf + (k + tf * k_max) * space_shape_f.size();
                    double v_fine = value_space_shape_f * value_time_shape_f;
                    MassMatrix_cf[f_idx][c_idx] += quad_f.qWeight[q] * v_coarse * v_fine;
                  }
                }
              }
            }
          }
        }

        RVector value_c(N_c);
        for (int i = 0; i < value_c.size(); i++) {
          value_c[i] = C(row_c)[i];
        }
        /*mout << DOUT(value_c.size()) << endl;
        mout << DOUT(MassMatrix_cf.rows()) << endl;
        mout << DOUT(InvertedMassMatrix_f.rows()) << endl;
        mout << DOUT(deg_f.space) <<DOUT(deg_f.time) <<  endl;
        mout << DOUT(deg_c.space) <<DOUT(deg_c.time) <<  endl;*/
        RVector integral_f(MassMatrix_cf * value_c);
        RVector value_f = InvertedMassMatrix_f * integral_f;

        for (int i = 0; i < value_f.size(); i++) {
          F(row_f)[i] += value_f[i];
        }
      }
    }


    F.MakeAdditive();
    F.Accumulate();
  }
};

class SpaceTimeTransferInterpolateDGDG : public SpaceTimeTransfer {

  const int maxDegSpace = 6;
  const int maxDegTime = 6;

  const STDiscretization &disc;
  const int k_max;

  vector<const Shape *> space_Shape;
  vector<const Shape *> time_Shape;
  vector<vector<vector<Point>>> z;

  vector<int> DofToDeg;
public:
  SpaceTimeTransferInterpolateDGDG(const STDiscretization &disc) :
      disc(disc), k_max(disc().GetDoF().get_dim()) {

    space_Shape.resize(maxDegSpace + 1);
    for (int sDeg = 0; sDeg < maxDegSpace + 1; sDeg++) {
      const Shape &SpaceShape(disc.GetCellShape(sDeg));
      space_Shape[sDeg] = &SpaceShape;
    }
    time_Shape.resize(maxDegTime + 1);
    for (int tDeg = 0; tDeg < maxDegTime + 1; tDeg++) {
      const Shape &TimeShape(disc.GetTimeShape(tDeg));
      time_Shape[tDeg] = &TimeShape;
    }
    z.resize(maxDegSpace + 1);
    for (int s = 0; s < maxDegSpace + 1; ++s) {
      z[s].resize(k_max);
      for (int k = 0; k < k_max; k++) {
        disc().GetDoF().NodalPointsLocal(s, z[s][k], k);
      }
    }
    DofToDeg.resize(disc().GetDoF().NodalPointsLocal(maxDegSpace) + 1);
    for (int s = 0; s < maxDegSpace + 1; ++s) {
      DofToDeg[disc().GetDoF().NodalPointsLocal(s)] = s;
    }
  };

  void Restrict(const Vector &F, Vector &C) const override {
    C.Clear();
    TransferInfo info(F, C);

    for (cell cc = C.cells(); cc != C.cells_end(); ++cc) {
      row rc = C.find_row(cc());
      DegreePair deg_c = C.GetDoF().GetDegree(*cc);
      vector<Point> Children{info.GetChildren(cc)};
      vector<row> row_f(Children.size());
      vector<DegreePair> deg_f(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_f[i] = F.find_row(Children[i]);
        deg_f[i] = F.GetDoF().GetDegree(Children[i]);
      }
      int idx_c = 0;
      for (int tc = 0; tc <= deg_c.time; ++tc) {
        for (int k = 0; k < k_max; ++k) {
          for (int sc = 0; sc < z[deg_c.space][k].size(); ++sc) {
            Point local_c;
            if (deg_c.time > 0) {
              local_c = z[deg_c.space][k][sc].CopyWithT(tc * 1.0 / deg_c.time);
            } else {
              local_c = z[deg_c.space][k][sc].CopyWithT(0.5);
            }
            vector<double> value_child(Children.size());
            vector<bool> z_in_child(Children.size());
            for (int i = 0; i < Children.size(); ++i) {
              cell cf = F.find_cell(Children[i]);
              value_child[i] = 0.0;
              Point local_f(local_c);

              z_in_child[i] = transformPointLocaltoLocal(local_c, cc, local_f, cf);

              if (z_in_child[i]) {
                int tnf_i = row_f[i].n() / (deg_f[i].time + 1);

                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[deg_f[i].space][kk].size();
                for (int sf = 0; sf < z[deg_f[i].space][k].size(); ++sf) {
                  for (int tf = 0; tf <= deg_f[i].time; ++tf) {
                    auto space_shape = space_Shape[DofToDeg[z[deg_f[i].space][k].size()]];
                    double value_space_shape = (*space_shape)(local_f, sf);
                    auto time_shape = time_Shape[deg_f[i].time];
                    double value_time_shape = (*time_shape)(Point(local_f.t()), tf);
                    value_child[i] += F(row_f[i])[tf * tnf_i + shift + sf]
                                      * value_space_shape * value_time_shape;
                  }
                }
              }
            }
            int cnt = 0;
            double value = 0.0;
            for (int i = 0; i < Children.size(); ++i) {
              value += value_child[i];
              if (z_in_child[i]) cnt++;
            }
            C(rc)[idx_c++] = value / cnt;
          }
        }
      }
    }
  }

  void Prolongate(Vector &F, const Vector &C) const override {
    F = 0.0;
    TransferInfo info(F, C);
    for (cell cf = C.cells(); cf != C.cells_end(); ++cf) {
      row rf = C.find_row(cf());
      DegreePair deg_f = C.GetDoF().GetDegree(*cf);

      vector<Point> Children{info.GetChildren(cf)};
      vector<row> row_c(Children.size());
      vector<int> time_deg_c(Children.size());
      vector<int> space_deg_c(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_c[i] = F.find_row(Children[i]);
        time_deg_c[i] = F.GetDoF().get_time_deg(Children[i]);
        space_deg_c[i] = F.GetDoF().get_space_deg(Children[i]);
      }

      int tnf = rf.n() / (deg_f.time + 1);
      for (int i = 0; i < Children.size(); ++i) {
        cell cc = F.find_cell(Children[i]);
        int tnc_i = row_c[i].n() / (time_deg_c[i] + 1);
        int idx_c = 0;
        for (int tc = 0; tc <= time_deg_c[i]; ++tc) {
          for (int k = 0; k < k_max; ++k) {
            for (int sc = 0; sc < z[space_deg_c[i]][k].size(); ++sc) {
              Point local_c;
              if (time_deg_c[i] > 0)
                local_c = z[space_deg_c[i]][k][sc].CopyWithT(0.5);
              else
                local_c = z[space_deg_c[i]][k][sc].CopyWithT(tc * 1.0 / time_deg_c[i]);

              Point local_f(local_c);
              transformPointLocaltoLocal(local_c, cc, local_f, cf);

              for (int tf = 0; tf <= deg_f.time; ++tf) {
                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[deg_f.space][kk].size();
                for (int sf = 0; sf < z[deg_f.space][k].size(); ++sf) {
                  double value_space_shape = (*space_Shape[DofToDeg[z[deg_f.space][k].size()]])(
                      local_f, sf);
                  double value_time_shape = (*time_Shape[deg_f.time])(Point(local_f.t()), tf);
                  F(row_c[i])[idx_c] +=
                      C(rf)[tf * tnf + shift + sf] * value_space_shape * value_time_shape;
                }
              }
              idx_c++;
            }
          }
        }
      }
    }
    F.MakeAdditive();
    F.Accumulate();
  }
};


class SpaceTimeTransferInterpolation : public SpaceTimeTransfer {

  const int maxDegSpace = 6;
  const int maxDegTime = 6;

  const STDiscretization &disc;
  const int k_max;

  vector<const Shape *> space_Shape;
  vector<const Shape *> time_Shape;
  vector<vector<vector<Point>>> z;

  vector<int> DofToDeg;

  bool dual = false;

public:
  explicit SpaceTimeTransferInterpolation(const STDiscretization &disc) :
      disc(disc), k_max(disc().GetDoF().get_dim()) {
    space_Shape.resize(maxDegSpace + 1);
    for (int sDeg = 0; sDeg < maxDegSpace + 1; sDeg++) {
      const Shape &SpaceShape(disc.GetCellShape(sDeg));
      space_Shape[sDeg] = &SpaceShape;
    }
    time_Shape.resize(maxDegTime + 1);
    for (int tDeg = 0; tDeg < maxDegTime + 1; tDeg++) {
      const Shape &TimeShape(disc.GetTimeShape(tDeg));
      time_Shape[tDeg] = &TimeShape;
    }
    z.resize(maxDegSpace + 1);
    for (int s = 0; s < maxDegSpace + 1; ++s) {
      z[s].resize(k_max);
      for (int k = 0; k < k_max; k++) {
        disc().GetDoF().NodalPointsLocal(s, z[s][k], k);
      }
    }
    DofToDeg.resize(disc().GetDoF().NodalPointsLocal(maxDegSpace) + 1);
    for (int s = 0; s < maxDegSpace + 1; ++s) {
      DofToDeg[disc().GetDoF().NodalPointsLocal(s)] = s;
    }
  }

  void Restrict(const Vector &F, Vector &C) const override {
    C = 0.;
    bool testspace = !dual;
    TransferInfo info(F, C);

    for (cell cc = C.cells(); cc != C.cells_end(); ++cc) {
      row rc = C.find_row(cc());
      DegreePair deg_c = C.GetDoF().GetDegree(*cc);
      vector<Point> Children{info.GetChildren(cc)};
      vector<row> row_f(Children.size());
      vector<DegreePair> deg_f(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_f[i] = F.find_row(Children[i]);
        deg_f[i] = F.GetDoF().GetDegree(Children[i]);
      }


      int tnc = rc.n() / deg_c.time;
      int idx_c = 0;
      for (int tc = 1; tc <= deg_c.time; ++tc) {
        for (int k = 0; k < k_max; ++k) {
          for (int sc = 0; sc < z[deg_c.space][k].size(); ++sc) {
            const Point z_space_c = z[deg_c.space][k][sc];
            if (false)
              for (int sc2 = 0; sc2 < z[deg_c.space].size(); ++sc2) {
                const Point z_space_c = z[deg_c.space][k][sc2];
                if ((*space_Shape[deg_c.space])(z_space_c, sc) < 1e-12)
                  continue;
                if (sc != sc2) Exit("Basisfkt und NP nicht gleich geordnet!");
              }
            Point z_time = tc * 1.0 / deg_c.time;
            vector<double> value_child(Children.size());
            vector<bool> z_in_child(Children.size());
            for (int i = 0; i < Children.size(); ++i) {
              cell cf = F.find_cell(Children[i]);
              int tnf_i = row_f[i].n() / deg_f[i].time;
              value_child[i] = 0.0;
              z_in_child[i] = true;
              Point z_space(z_space_c);
              if (info.in_space) {
                z_in_child[i] = transformPointLocaltoLocal(z_space_c, cc, z_space, cf);
              }
              if (info.in_time) {
                if (cf().t() < cc().t()) {
                  z_time = tc * 1.0 / deg_c.time * 2;
                  if (z_time[0] > 1) z_in_child[i] = false;
                } else if (cf().t() > cc().t()) {
                  z_time = 1 - (1 - tc * 1.0 / deg_c.time) * 2;
                  if (z_time[0] < 0) z_in_child[i] = false;
                }
              }
              if (z_in_child[i]) {
                if (cf.min() == 0.0 && !testspace) {
                } else if (!testspace) {
                  face ff = F.find_face(cf.Face(cf.Faces() - 2));
                  cell cf_prev = F.find_cell_or_overlap_cell(ff.Left());
                  row rf_prev = F.find_row(cf_prev());
                  DegreePair deg_f_prev = F.GetDoF().GetDegree(*cf_prev);
                  int shift = 0;
                  for (int kk = 0; kk < k; kk++)
                    shift += z[deg_f_prev.space][kk].size();

                  for (int sf = 0; sf < z[deg_f_prev.space][k].size(); ++sf) {
                    double value_space_shape
                        = (*space_Shape[DofToDeg[z[deg_f_prev.space][k].size()]])(z_space, sf);
                    //check next line!
                    double value_time_shape = (*time_Shape[deg_f[i].time])(z_time, 0);
                    value_child[i] += F(rf_prev)[rf_prev.n() / deg_f_prev.time + shift + sf]
                                      * value_space_shape * value_time_shape;
                  }
                }

                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[deg_f[i].space][kk].size();
                for (int sf = 0; sf < z[deg_f[i].space][k].size(); ++sf) {
                  for (int tf = 1; tf <= deg_f[i].time; ++tf) {
                    auto space_shape = space_Shape[DofToDeg[z[deg_f[i].space][k].size()]];
                    double value_space_shape = (*space_shape)(z_space, sf);
                    double value_time_shape = (*time_Shape[deg_f[i].time])(z_time, tf);
                    if (testspace) {
                      if (deg_f[i].time == 1) {
                        value_time_shape = 1.0;
                      } else {
                        value_time_shape = (*time_Shape[deg_f[i].time - 1])(z_time, tf - 1);
                      }
                    }
                    value_child[i] += F(row_f[i])[(tf - 1) * tnf_i + shift + sf]
                                      * value_space_shape * value_time_shape;
                  }
                }
              }
            }
            int cnt = 0;
            double value = 0.0;
            for (int i = 0; i < Children.size(); ++i) {
              value += value_child[i];
              if (z_in_child[i]) cnt++;
            }
            C(rc)[idx_c++] = value / cnt;
          }
        }
      }
    }
  }


  void Prolongate(Vector &F, const Vector &C) const override {
    F = 0.0;
    bool testspace = dual;
    TransferInfo info(F, C);

    for (cell cf = C.cells(); cf != C.cells_end(); ++cf) {
      row rf = C.find_row(cf());
      DegreePair deg_f = C.GetDoF().GetDegree(*cf);

      // find prev coarse cell in time
      face ff = C.find_face(cf.Face(cf.Faces() - 2));
      cell c_prev = C.find_cell_or_overlap_cell(ff.Left());
      row rf_prev = C.find_row(c_prev());

      // get corresponding polynomial degrees
      DegreePair deg_f_prev = C.GetDoF().GetDegree(*c_prev);
      vector<Point> Children{info.GetChildren(cf)};
      vector<row> row_c(Children.size());
      vector<DegreePair> deg_c(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_c[i] = F.find_row(Children[i]);
        deg_c[i] = F.GetDoF().GetDegree(Children[i]);
      }

      int tnf = rf.n() / deg_f.time;
      for (int i = 0; i < Children.size(); ++i) {
        cell cc = F.find_cell(Children[i]);
        int tnc_i = row_c[i].n() / deg_c[i].time;
        int idx_c = 0;
        for (int tc = 1; tc <= deg_c[i].time; ++tc) {
          for (int k = 0; k < k_max; ++k) {
            for (int sc = 0; sc < z[deg_c[i].space][k].size(); ++sc) {
              Point z_space(z[deg_c[i].space][k][sc]);
              Point z_time = tc * 1.0 / deg_c[i].time;
              if (info.in_space)
                transformPointLocaltoLocal(z[deg_c[i].space][k][sc], cc, z_space, cf);
              if (info.in_time) {
                if (info.in_space) {
                  if (2 * i >= Children.size())
                    z_time = tc * 0.5 / (deg_c[i].time) + i * 0.5;
                } else {
                  z_time = tc * 0.5 / (deg_c[i].time) + i * 0.5;
                }
              }

              if (testspace == false) {
                if (cf.min() == 0.0) {
                  //for (int sf=0; sf<snf; ++sf){
                  //    double value_space_shape = (*space_Shape[space_deg_f])(z_space,sf);
                  //    double value_time_shape = (*time_Shape[time_deg_f])(z_time,0);
                  //    //Point z_space_glob = (1-z_space[0]-z_space[1]) * cf.Corner(0) + z_space[0] * cf.Corner(1) + z_space[1] * cf.Corner(2);
                  //    for (int k=0; k<k_max; ++k){
                  //        C(row_c[i])[(tc-1)*tnc_i + k*snc_i + sc] += F(rf)[0 + k*snf + sf] * value_space_shape * value_time_shape; //verwende folgenden zeitpunkt
                  //    }
                  //}
                } else {
                  int shift = 0;
                  for (int kk = 0; kk < k; kk++)
                    shift += z[deg_f_prev.space][kk].size();
                  for (int sf = 0; sf < z[deg_f_prev.space][k].size(); ++sf) {
                    //double value_space_shape = (*space_Shape[space_deg_f_prev])(z_space,sf);
                    auto space_shape = space_Shape[DofToDeg[z[deg_f_prev.space][k].size()]];
                    double value_space_shape = (*space_shape)(z_space, sf);
                    double value_time_shape = (*time_Shape[deg_f.time])(z_time, 0);
                    F(row_c[i])[idx_c] += C(rf_prev)[rf_prev.n() / deg_f_prev.time + shift + sf]
                                          * value_space_shape * value_time_shape;
                  }
                }
              }
              for (int tf = 1; tf <= deg_f.time; ++tf) {
                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[deg_f.space][kk].size();
                for (int sf = 0; sf < z[deg_f.space][k].size(); ++sf) {
                  double value_space_shape = (*space_Shape[DofToDeg[z[deg_f.space][k].size()]])(
                      z_space, sf);
                  double value_time_shape = (*time_Shape[deg_f.time])(z_time, tf);
                  if (testspace) {
                    if (deg_f.time == 1) {
                      value_time_shape = 1.0;
                    } else {
                      value_time_shape = (*time_Shape[deg_f.time - 1])(z_time, tf - 1);
                    }
                  }
                  F(row_c[i])[idx_c] += C(rf)[(tf - 1) * tnf + shift + sf]
                                        * value_space_shape * value_time_shape;
                }
              }
              idx_c++;
            }
          }
        }
      }
    }
    F.MakeAdditive();
    F.Accumulate();
  }
};

std::unique_ptr<SpaceTimeTransfer> GetSpaceTimeTransfer(const STDiscretization &disc,
                                                        const std::string &name);


#endif //SPACETIME_SPACETIMETRANSFER_HPP
