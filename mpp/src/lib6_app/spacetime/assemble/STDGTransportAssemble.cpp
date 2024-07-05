#include "STDGTransportAssemble.hpp"
#include "STDGDGTransportElement.hpp"
#include "STDGDGTransportFaceElement.hpp"
#include "DebuggingTools.hpp"
#include "SymRMatrix.hpp"
#include "LinearSolver.hpp"

using namespace std;

void STDGTransportAssemble::MassMatrix(Matrix &M) const {
  M = 0.0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(M, c);
    DGRowEntries M_c(M, *c, *c, false);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        double v_i = elem.Density(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          double v_j = elem.Density(q, j);
          M_c(i, j) += w * v_j * v_i;
        }
      }
    }
  }
}

void STDGTransportAssemble::System(Matrix &B, Vector &RHS) const {
  Date Start;
  mout.PrintInfo("TDGTransportAssemble_DGT", verbose,
                 PrintInfoEntry<int>("Problem Size", B.pSize(), 2),
                 PrintInfoEntry<Date>("start assemble", Start, 2));
  Time t_cell;
  Time t_face;
  const Mesh &M = B.GetMesh();
  const double T = M.GetEndTime();
  std::function<double(Point)> d_T = [T](Point QP) { return (T - QP.t()); };
  std::function<double(Point)> constOne_func = [](Point QP) { return 1.0; };

  B = 0;
  RHS = 0;
  for (cell c = B.cells(); c != B.cells_end(); ++c) {
    Date Start_cell;
    SpaceTimeTransportDGTElement elem(B, c);
    row r_c = B.find_row(c());
    DGRowEntries B_c(B, *c, *c, false);
    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double t = QP.t();
      double w = elem.QWeight(q);
      VectorField flux = prob->Flux(M.Level(), *c, QP);
      for (int i = 0; i < elem.i_dimension(); ++i) {
        double v_i = elem.Density(q, i);
        for (int j = 0; j < elem.j_dimension(); ++j) {
          double v_j = elem.Density(q, j);
          VectorField Dv_j = elem.GradDensity(q, j);
          double Dtv_j = elem.DtDensity(q, j);
          B_c(i, j) += w * (Dtv_j + flux * Dv_j) * v_i;
        }
      }
    }

    t_cell += Date() - Start_cell;
    Date Start_face;

    for (int f = 0; f < c.Faces() - 2; ++f) {
      SpaceTimeTransportDGTFaceElement felem(*disc, B, c, f);
      if (prob->HasOutflow(M.Level(), *c, felem.QPoint(0), felem.QNormal(0))) continue;
      if (B.OnBoundary(*c, f)) {
        for (int q = 0; q < felem.nQ(); ++q) {
          Point QP = felem.QPoint(q);
          VectorField Nq = felem.QNormal(q);
          double fn = prob->FaceNormalFlux(M.Level(), *c, f, Nq, QP);
          if (fn >= 0) continue;
          double w = felem.QWeight(q);
          double g_in = prob->InflowData(M.Level(), *c, QP.t(), QP, Nq);
          for (int i = 0; i < felem.i_dimension(); ++i) {
            double v_i = felem.Density(q, i);
            RHS(c(), i) -= w * g_in * v_i;
            for (int j = 0; j < felem.j_dimension(); ++j) {
              double v_j = felem.Density(q, j);
              B_c(i, j) -= w * v_i * fn * v_j;
            }
          }
        }
      } else {
        cell cf = B.find_neighbour_cell(c, f);
        int f1 = B.find_neighbour_face_id(c.Face(f), cf);
        SpaceTimeTransportDGTFaceElement felem_1(*disc, B, cf, f1);
        DGRowEntries B_cf(B, *c, *cf, false);
        for (int q = 0; q < felem.nQ(); ++q) {
          Point QP = felem.QPoint(q);
          VectorField Nq = felem.QNormal(q);
          double fn = prob->FaceNormalFlux(M.Level(), *c, f, Nq, QP);
          if (fn >= 0) continue;

          int q1 = find_q_id(felem.QPoint(q), felem_1);
          Point QP1 = felem_1.QPoint(q1);

          double w = felem.QWeight(q);
          for (int i = 0; i < felem.i_dimension(); ++i) {
            double v_i = felem.Density(q, i);
            for (int j = 0; j < felem.j_dimension(); ++j) {
              double v_j = felem.Density(q, j);
              B_c(i, j) -= w * fn * v_i * v_j;
            }
            for (int j = 0; j < felem_1.j_dimension(); ++j) {
              double v_j = felem_1.Density(q1, j);
              B_cf(i, j) += w * fn * v_i * v_j;
            }
          }
        }
      }
    }
    SpaceTimeTransportDGTFaceElement felem(*disc, B, c, c.Faces() - 2, "dgdg");
    if (c.min() == 0.0) {
      for (int q = 0; q < felem.nQ(); ++q) {
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
	double u0 = prob->Solution(0.0, *c, small_move * c() + (1 - small_move) * QP);
        for (int i = 0; i < felem.i_dimension(); ++i) {
          double v_i = felem.Density(q, i);	  
          RHS(c(), i) += w * u0 * v_i;
          for (int j = 0; j < felem.j_dimension(); ++j) {
            double v_j = felem.Density(q, j);
            B_c(i, j) += w * v_i * v_j;
          }
        }
      }
    } else {
      cell c_prev = B.find_previous_cell(c);
      if (c_prev != c) {
        //cell c_prev = find_previous_cell(B.GetMesh(),c);
        //if (c_prev != B.cells_end()) {
        SpaceTimeTransportDGTFaceElement felem_1(*disc, B, c_prev, c_prev.Faces() - 1, "dgdg");
        DGRowEntries B_c_prev(B, *c, *c_prev, false);
        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          for (int i = 0; i < felem.i_dimension(); ++i) {
            double v_i = felem.Density(q, i);
            for (int j = 0; j < felem.j_dimension(); ++j) {
              double v_j = felem.Density(q, j);
              B_c(i, j) += w * v_i * v_j;
            }
            for (int j = 0; j < felem_1.j_dimension(); ++j) {
              double v_j = felem_1.Density(q1, j);
              B_c_prev(i, j) -= w * v_i * v_j;
            }
          }
        }
      } else {
        for (int q = 0; q < felem.nQ(); ++q) {
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          for (int i = 0; i < felem.i_dimension(); ++i) {
            double v_i = felem.Density(q, i);
            for (int j = 0; j < felem.j_dimension(); ++j) {
              double v_j = felem.Density(q, j);
              B_c(i, j) += w * v_i * v_j;
            }
          }
        }
      }
    }
    t_face += Date() - Start_face;
  }


  Time t_cell_max = PPM->Max(t_cell.t);
  Time t_cell_min = PPM->Min(t_cell.t);
  Time t_face_max = PPM->Max(t_face.t);
  Time t_face_min = PPM->Min(t_face.t);

  mout.PrintInfo("TDGTransportAssemble_DGT", verbose,
                 PrintInfoEntry<Time>("min cell assemble time", t_cell_min),
                 PrintInfoEntry<Time>("max cell assemble time", t_cell_max),
                 PrintInfoEntry<Time>("min face assemble time", t_face_min),
                 PrintInfoEntry<Time>("max face assemble time", t_face_max),
                 PrintInfoEntry<Time>("assemble time", Date() - Start),
                 PrintInfoEntry<Date>("finish assemble", Date()));
}

pair<double, double> STDGTransportAssemble::DiscNorm(const Matrix &L, const Vector &u) const {
  Vector Lu(L * u);
  Matrix M(u);
  MassMatrix(M);
  Vector Mu(M * u);

  Preconditioner *PC = GetPC("PointBlockGaussSeidel");
  LinearSolver *S = GetLinearSolver("GMRES", PC, "NormSolver");
  (*S)(M);

  Vector mLu(u);
  mLu = (*S) * Lu;

  double discNormW = u * Mu;
  double discNormV = Lu * mLu;

  return pair<double, double>(std::sqrt(discNormW), std::sqrt(discNormW + discNormV));
}

double STDGTransportAssemble::L1Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      double U =
          elem.Density(q, u) - prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
      nrm += w * abs(U);
    }
  }
  return PPM->SumOnCommSplit(nrm, u.CommSplit());
}

double STDGTransportAssemble::L2Error(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      double U =
          elem.Density(q, u) - prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
      nrm += w * U * U;
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::GNError(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point QP = elem.QPoint(q);
      double U =
          elem.Density(q, u) - prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
      double dtU = elem.DtDensity(q, u);
      VectorField DU = elem.GradDensity(q, u);
      double F = dtU + prob->DivFlux(u.Level(), *c, QP, U, DU);
      nrm += w * (U * U + F * F);
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::LInfError(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      Point QP = elem.QPoint(q);
      double U =
          elem.Density(q, u) - prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
      nrm = max(nrm, abs(U));
    }
  }
  return PPM->Max(nrm, u.CommSplit());
};

double STDGTransportAssemble::L1Norm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      double U = elem.Density(q, u);
      nrm += w * abs(U);
    }
  }
  return PPM->SumOnCommSplit(nrm, u.CommSplit());
}

double STDGTransportAssemble::L2Norm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      double U = elem.Density(q, u);
      nrm += w * U * U;
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::L2SpaceNormAtTime(const Vector &u, double time) const {
  double nrm = 0.0;
  if (time == u.GetMesh().GetEndTime()) {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != time) continue;
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 1, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        double U = felem.Density(q, u);
        nrm += w * U * U;
      }
    }
  } else {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.min() != time) continue;
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 2, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        double U = felem.Density(q, u);
        nrm += w * U * U;
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::L2SpaceNormAtTimeError(const Vector &u, double time) const {
  double nrm = 0.0;
  if (time == u.GetMesh().GetEndTime()) {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.max() != time) continue;
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 1, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        double U = felem.Density(q, u) -
                   prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
        nrm += w * U * U;
      }
    }
  } else {
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      if (c.min() != time) continue;
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 2, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        double U = felem.Density(q, u) -
                   prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
        nrm += w * U * U;
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::ResidualErrorEstimateJumpCell(const cell &c,
                                                            const Vector &u) const {
  double eta = 0.0;
  const Mesh &M = u.GetMesh();
  double h = M.MaxMeshWidth();
  SpaceTimeTransportDGTElement elem(u, c);
  for (int q = 0; q < elem.nQ(); q++) {
    Point QP = elem.QPoint(q);
    double w = elem.QWeight(q);
    double dtU = elem.DtDensity(q, u);
    double U = elem.Density(q, u);
    VectorField DU = elem.GradDensity(q, u);
    double F = dtU + prob->DivFlux(u.Level(), *c, QP, U, DU);
    eta += h * w * F * F;
  }
  for (int f = 0; f < c.Faces() - 2; ++f) {
    if (u.OnBoundary(*c, f)) {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, f);
      if (prob->HasOutflow(u.Level(), *c, felem.QPoint(0), felem.QNormal(0))) continue;
      for (int q = 0; q < felem.nQ(); ++q) {
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);
        double fn = prob->FaceNormalFlux(M.Level(), *c, f, Nq, QP);
        if (fn == 0) continue;
        double w = felem.QWeight(q);
        double U = felem.Density(q, u);
        if (prob->HasOutflow(u.Level(), *c, felem.QPoint(0), felem.QNormal(0)))
          eta += 0.5 * w * U * U / abs(fn);
        else {
          double g =
              prob->InflowData(u.Level(), *c, QP.t(), QP, Nq) - U * prob->Flux(u.Level(), *c, QP) * Nq;
          eta += 0.5 * w * g * g / abs(fn);
        }
      }
    } else {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, f);
      cell cf = u.find_neighbour_cell(c, f);
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      SpaceTimeTransportDGTFaceElement felem_1(*disc, u, cf, f1);
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        VectorField Nq = felem.QNormal(q);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
        double U = felem.Density(q, u) - felem_1.Density(q1, u);
        double fn = prob->FaceNormalFlux(M.Level(), *c, f, Nq, QP);
        eta += 0.25 * w * abs(fn) * U * U;
      }
    }
  }
  SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 2, "dgdg");
  if (c.min() == u.t(0)) {
    for (int q = 0; q < felem.nQ(); ++q) {
      double w = felem.QWeight(q);
      Point QP = felem.QPoint(q);
      double U = felem.Density(q, u) -
                 prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
      eta += 0.5 * w * U * U;
    }
  } else {
    cell c_prev = u.find_previous_cell(c);
    if (c_prev != c) {
      SpaceTimeTransportDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        int q1 = find_q_id(felem.QPoint(q), felem_1);
        Point QP = felem.QPoint(q);
        double w = felem.QWeight(q);
        double U = felem.Density(q, u) - felem_1.Density(q1, u);
        eta += 0.25 * w * U * U;
      }
    }
  }
  if (c.max() != u.GetMesh().GetEndTime()) {
    face ff = u.find_face(c.Face(c.Faces() - 1));
    if (ff.Right() != Infty) {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 1, "dgdg");
      cell cf = u.find_neighbour_cell(c, c.Faces() - 1);
      if ((cf != c) && cf != u.cells_end()) {
        SpaceTimeTransportDGTFaceElement felem_1(*disc, u, cf, cf.Faces() - 2, "dgdg");
        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          double w = felem.QWeight(q);
          Point QP = felem.QPoint(q);
          double U = felem.Density(q, u) - felem_1.Density(q1, u);
          eta += 0.25 * w * U * U;
        }
      }
    }
  }
  return sqrt(eta);
}

double STDGTransportAssemble::DGError(const Vector &u) const {
  double nrm = 0;
  double h = u.GetMesh().MaxMeshWidth();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); q++) {
      Point QP = elem.QPoint(q);
      double w = elem.QWeight(q);
      double U = elem.Density(q, u);
      double dtU = elem.DtDensity(q, u);
      VectorField DU = elem.GradDensity(q, u);
      double F = dtU + prob->DivFlux(u.Level(), *c, QP, U, DU);
      nrm += h * w * F * F;
    }
    for (int f = 0; f < c.Faces() - 2; ++f) {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, f);
      if (u.OnBoundary(*c, f)) {
        for (int q = 0; q < felem.nQ(); ++q) {
          VectorField Nq = felem.QNormal(q);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);

          double fn = prob->FaceNormalFlux(u.Level(), *c, f, Nq, QP);
          double U = felem.Density(q, u) -
                     prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);

          nrm += 0.5 * w * abs(fn) * U * U;
        }
      } else {
        cell cf = u.find_neighbour_cell(c, f);
        if (c() < cf()) continue;
        int f1 = u.find_neighbour_face_id(c.Face(f), cf);
        SpaceTimeTransportDGTFaceElement felem_1(*disc, u, cf, f1);
        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          VectorField Nq = felem.QNormal(q);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          double U = felem.Density(q, u) - felem_1.Density(q1, u);
          double fn = prob->FaceNormalFlux(u.Level(), *c, f, Nq, QP);
          nrm += 0.5 * w * abs(fn) * U * U;
        }
      }
    }
    SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 2, "dgdg");
    if (c.min() == u.t(0)) {
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        double U = felem.Density(q, u) -
                   prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
        nrm += 0.5 * w * U * U;
      }
    } else {
      cell c_prev = u.find_previous_cell(c);
      if (c_prev != c) {
        SpaceTimeTransportDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, "dgdg");
        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          double U = felem.Density(q, u) - felem_1.Density(q1, u);
          nrm += 0.5 * w * U * U;
        }
      }
    }
    if (c.max() == u.GetMesh().GetEndTime()) {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 1, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        double U = felem.Density(q, u) -
                   prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
        nrm += 0.5 * w * U * U;
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::DGSemiNorm(const Vector &u) const {
  double nrm = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 2, "dgdg");
    if (c.min() == u.t(0)) {
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        double U = felem.Density(q, u);
        nrm += w * U * U;
      }
    } else {
      cell c_prev = u.find_previous_cell(c);
      if (c_prev != c) {
        SpaceTimeTransportDGTFaceElement felem_1(*disc, u, c_prev, c_prev.Faces() - 1, "dgdg");
        for (int q = 0; q < felem.nQ(); ++q) {
          int q1 = find_q_id(felem.QPoint(q), felem_1);
          Point QP = felem.QPoint(q);
          double w = felem.QWeight(q);
          double U = felem.Density(q, u) - felem_1.Density(q1, u);
          nrm += 0.5 * w * U * U;
        }
      }
    }
    if (c.max() == u.GetMesh().GetEndTime()) {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 1, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Point QP = felem.QPoint(q);
        double U = felem.Density(q, u);
        nrm += w * U * U;
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, u.CommSplit()));
}

double STDGTransportAssemble::DGNormError(const Vector &u) {
  return DGError(u);
}

double STDGTransportAssemble::DGNorm(const Vector &u, const Matrix &L) const {
  double normsq = 0;
  Vector Lu = L * u;
  double buu = Lu * u;

  for (cell c = u.cells(); c != u.cells_end(); c++) {
    double h = max(0.0, dist(c[0], c[c.Corners() - 1]));
    SpaceTimeTransportDGTElement elem(u, c);
    for (int q = 0; q < elem.nQ(); q++) {
      Point QP = elem.QPoint(q);
      double w = elem.QWeight(q);
      double dtU = elem.DtDensity(q, u);
      double U = elem.Density(q, u);
      VectorField DU = elem.GradDensity(q, u);
      double F = dtU + prob->DivFlux(u.Level(), *c, QP, U, DU);
      normsq += h * w * F * F;
    }
  }
  return sqrt(buu + PPM->SumOnCommSplit(normsq, u.CommSplit()));
}

double STDGTransportAssemble::L2ScalarProduct(const Vector &a, const Vector &b) const {
  double nrm = 0;
  for (cell c = a.cells(); c != a.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(a, c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      nrm += w * elem.Density(q, q) * elem.Density(q, b);
    }
  }
  return sqrt(PPM->SumOnCommSplit(nrm, a.CommSplit()));
}

void STDGTransportAssemble::Mass(const Vector &u,
                                 double T,
                                 vector<double> &Mass,
                                 vector<double> &inflow,
                                 vector<double> &outflow,
                                 vector<double> &outflow_bnd) const {
  int N = Mass.size() - 1;
  const Mesh &M = u.GetMesh();
  double dt = T / N;
  double initialMass = 0;
  for (int n = 0; n <= N; ++n) Mass[n] = inflow[n] = outflow[n] = outflow_bnd[n] = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    if (c.min() == 0.0) {
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 2, "dgdg");
      for (int q = 0; q < felem.nQ(); ++q) {
        double w = felem.QWeight(q);
        Mass[0] += w * felem.Density(q, u);
        Point QP = felem.QPoint(q);
        initialMass += w * prob->Solution(QP.t(), *c, small_move * c() + (1 - small_move) * QP);
      }
    }    
    int n = -1;
    for(int i = 0; i < M.GetTimesteps().size(); i++) {
      if (M.GetTimesteps()[i] == c.max()) {
	n = i;
	break;
      }
    }
    if (n < 0) THROW("Timestep not found!")

    SpaceTimeTransportDGTFaceElement felem(*disc, u, c, c.Faces() - 1, "dgdg");
    for (int q = 0; q < felem.nQ(); ++q) {
      double U = felem.Density(q, u);
      double w = felem.QWeight(q);
      Mass[n] += w * U;
    }
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int f = 0; f < c.Faces() - 2; ++f) {
      if (bnd[f] < 0) continue;
      SpaceTimeTransportDGTFaceElement felem(*disc, u, c, f);
      for (int q = 0; q < felem.nQ(); ++q) {
        Point QP = felem.QPoint(q);
        VectorField Nq = felem.QNormal(q);
        double w = felem.QWeight(q);
        double fn = prob->FaceNormalFlux(M.Level(), *c, f, Nq, QP);
        if (fn > 0) {
          outflow[n] += w * fn * felem.Density(q, u);
          outflow_bnd[n] += w;
        } else inflow[n] -= w * fn * felem.Density(q, u);
      }
    }
  }
  initialMass = PPM->SumOnCommSplit(initialMass, u.CommSplit());
  mout << "initial Mass " << initialMass << " used for scaling" << endl;
  Mass[0] = initialMass;
  inflow[0] = outflow[0] = outflow_bnd[0] = 0.0;  
  for (int n = 0; n <= N; ++n) {
    Mass[n] = PPM->SumOnCommSplit(Mass[n], u.CommSplit()) / initialMass;
    inflow[n] = PPM->SumOnCommSplit(inflow[n], u.CommSplit()) / initialMass;
    outflow[n] = PPM->SumOnCommSplit(outflow[n], u.CommSplit()) / initialMass;
    outflow_bnd[n] = PPM->SumOnCommSplit(outflow_bnd[n], u.CommSplit()) / initialMass;
  }
}

void STDGTransportAssemble::PlotSingleSolution(const Vector &u, const string &filename) const {
  auto vtuDisc = std::make_shared<STDiscretizationT_DGDG<>>(u.GetDisc().GetMeshes(), DegreePair{0, 0}, 1);
  Vector U_vtu(0.0, vtuDisc);

  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    SpaceTimeTransportDGTElement elem(u, c);
    U_vtu(c(), 0) = elem.Density(c(), u);
  }

  VtuPlot plot{filename};
  plot.AddData("Density", U_vtu, 0);
  plot.PlotFile();
}


