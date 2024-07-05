#include "LagrangeVectorValuedAssemble.hpp"
#include "VectorFieldElement.hpp"


const char *LagrangeVectorValuedAssemble::Name() const {
  return "LagrangeVectorValuedAssemble";
}

void LagrangeVectorValuedAssemble::Initialize(Vector &u) const {
  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;
    ScalarElement elem(u, *c);
    for (int face = 0; face < c.Faces(); ++face) {
      if (u_c.bc(face) == 1) {
        for (int j = 0; j < u.NumberOfNodalPointsOnFace(*c, face); ++j) {
          int k = u.IdOfNodalPointOnFace(*c, face, j);
          VectorField D(problem.Solution(elem.NodalPoint(k)));
          for (int d = 0; d < c.dim(); ++d) {
            u_c(k, d) = D[d];
            u_c.D(k, d) = true;
          }
        }
      }
    }
  }
  u.DirichletConsistent();
}


double LagrangeVectorValuedAssemble::Energy(const Vector &u) const {
  double energy = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Tensor DU = elem.VectorGradient(q, u);
      Tensor K = problem.Permeability(elem.QPoint(q));
      Tensor K_DU = K * DU;
      energy += w * Frobenius(K_DU, DU);
    }
  }
  return sqrt(PPM->SumOnCommSplit(energy, u.CommSplit()));
}

void LagrangeVectorValuedAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {
  VectorFieldElement elem(u, *c);
  RowBndValues r_c(r, *c);

  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);
    VectorField f = problem.Load(z, c);
    Tensor DU = elem.VectorGradient(q, u);
    Tensor K_DU = K * DU;
    for (int i = 0; i < elem.NodalPoints(); ++i) {
      for (int k = 0; k < c.dim(); ++k) {
        auto U_i = elem.VectorComponentValue(q, i,k);
        auto DU_i = elem.VectorRowGradient(q, i,k);
        r_c(i,k) += w * ( Frobenius(DU_i,K_DU) - (U_i * f));
      }
    }
  }
  if (!r_c.onBnd()) return;
  for (int f = 0; f < c.Faces(); ++f) {
    int bnd = r_c.bc(f);
    if (bnd != 2) continue;
    VectorFieldFaceElement faceElem(u, *c, f);
    RowValues r_f(r, faceElem);
    for (int q = 0; q < faceElem.nQ(); ++q) {
      double w = faceElem.QWeight(q);
      VectorField F = problem.Flux(faceElem.QPoint(q)) * faceElem.QNormal(q);
      for (int i = 0; i < faceElem.NodalPoints(); ++i) {
        for (int k = 0; k < faceElem.Dim(); ++k) {
          auto U_i = faceElem.VectorComponentValue(q, i,k);
          r_f(i, k) -= w * (U_i * F);
        }
      }
    }
  }
}

void LagrangeVectorValuedAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {
  VectorFieldElement elem(u, *c);
  RowEntries A_c(A, elem);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Tensor K = problem.Permeability(elem.QPoint(q));
    for (unsigned i = 0; i < elem.NodalPoints(); ++i) {
      for (int k = 0; k < c.dim(); ++k) {
        auto Dv = K * elem.VectorRowGradient(q, i, k);
        for (unsigned j = 0; j < elem.NodalPoints(); ++j) {
          for (int l = 0; l < c.dim(); ++l) {
            auto Dw = elem.VectorRowGradient(q, j, l);
            A_c(i, j, k, l) += w * (Frobenius(Dv, Dw));
          }
        }
      }
    }
  }
}

double LagrangeVectorValuedAssemble::EnergyError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      Tensor K = problem.Permeability(z);
      Tensor IK = Invert(K);
      Tensor diff = (K * elem.VectorGradient(q, u) -
                     problem.Flux(z));
      err += w * Frobenius(IK * diff , diff); // (IK * diff) * diff
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double LagrangeVectorValuedAssemble::L2(const Vector &u) const {
  double l2 = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      l2 += w * U * U;
    }
  }
  return sqrt(PPM->SumOnCommSplit(l2, u.GetMesh().CommSplit()));
}

double LagrangeVectorValuedAssemble::H1(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      Tensor DU = elem.VectorGradient(q, u);
      err += w * (U * U) + w * Frobenius(DU, DU);
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double LagrangeVectorValuedAssemble::L2Error(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      auto Sol = problem.Solution(elem.QPoint(q));
      auto diff = U - Sol;
      err += w * ((diff) * (diff));
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double LagrangeVectorValuedAssemble::L2CellAvgError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    double w = 0.0;
    VectorField U{};
    auto Sol = zero;
    for (int q = 0; q < elem.nQ(); ++q) {
      w += elem.QWeight(q);
      U += elem.VectorValue(q, u);
      Sol += problem.Solution(elem.QPoint(q));
    }
    err += w * (U - Sol) * (U - Sol);
  }
  return sqrt(PPM->SumOnCommSplit(err, u.GetMesh().CommSplit()));
}

double LagrangeVectorValuedAssemble::MaxError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    VectorFieldElement elem(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField U = elem.VectorValue(q, u);
      auto solution = problem.Solution(elem.QPoint(q));
      for (int k = 0; k < c.dim(); ++k) {
        auto diff = abs(U[k] - solution[k]);
        err = max(err, diff);
      }
    }
  }
  return PPM->Max(err, u.GetMesh().CommSplit());
}

double LagrangeVectorValuedAssemble::FluxError(const Vector &u) const {
  double flux_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] == 2) {
        VectorFieldFaceElement faceElem(u, *c, face);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Tensor G = problem.Flux(faceElem.QPoint(q));
          Tensor K = problem.Permeability(faceElem.QPoint(q));
          Tensor DU = faceElem.VectorGradient(q, u);
          auto F = (K * DU - G) * faceElem.QNormal(q);
          flux_error += w * F * F;
        }
      } else if (bnd[face] == 0) {
        VectorFieldFaceElement faceElem(u, *c, face);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Tensor K = problem.Permeability(faceElem.QPoint(q));
          Tensor DU = faceElem.VectorGradient(q, u);
          auto F = K * DU * faceElem.QNormal(q);
          flux_error += w * F * F;
        }
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(flux_error, u.CommSplit()));
}

double LagrangeVectorValuedAssemble::FaceError(const Vector &u) const {
  double face_error = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    for (int face = 0; face < c.Faces(); ++face) {
      VectorFieldFaceElement faceElem(u, *c, face);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        VectorField U = faceElem.VectorValue(q, u);
        auto Sol = problem.Solution(faceElem.QPoint(q));
        face_error += w * (U - Sol) * (U - Sol);
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(face_error, u.CommSplit()));
}


void LagrangeVectorValuedAssemble::SetExactSolution(Vector &uEx) const {
  for (row r = uEx.rows(); r != uEx.rows_end(); ++r) {
    VectorField solution = problem.Solution(r());
    for (int k = 0; k < uEx.dim(); ++k) {
      uEx(r, k) = solution[k];
    }
  }
}

