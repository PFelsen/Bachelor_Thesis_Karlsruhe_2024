//
// Created by lstengel on 27.07.22.
//

#include "DGVectorValuedAssemble.hpp"

#include "DGElement.hpp"


const char *DGVectorValuedAssemble::Name() const {
  return "DGVectorValuedAssemble";
}

void DGVectorValuedAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);
}

double DGVectorValuedAssemble::Energy(const Vector &u) const {
  double energy = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Tensor DU = elem.VectorGradient(q, u);
      Tensor K = problem.Permeability(elem.QPoint(q));
      Tensor K_DU = K * DU;
      energy += w *Frobenius(K_DU,DU);
    }
  }
  return sqrt(PPM->SumOnCommSplit(energy, u.CommSplit()));
}

void DGVectorValuedAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {
  
  const DGVectorFieldElement &elem = elementPool.Get(u, *c);

  DGSizedRowBndValues r_c(r, *c, elem.NodalPoints());

  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);
    VectorField f = problem.Load(z, c);
    Tensor K_DU = K * elem.VectorGradient(q, u);
    for (int i = 0; i < elem.NodalPoints(); ++i) {
      for (int k = 0; k < elem.Dim(); ++k) {
        auto U_i = elem.VectorComponentValue(q, i, k);
        auto DU_i = elem.VectorRowGradient(q, i, k);
        r_c(i,k) += w * (Frobenius(DU_i,K_DU) - U_i * f);
      }
    }
  }
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar scaledPenalty = (pow(degree,2)*penalty) / u.GetMesh().MaxMeshWidth();

    const DGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, f);
    
    if (r.OnBoundary(*c, f)) {
      int bnd = r_c.bc(f);
      if (bnd == 2) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          VectorField F = problem.Flux(z) * faceElem.QNormal(q);
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            for (int k = 0; k < elem.Dim(); ++k) {
              auto U_i = faceElem.VectorComponentValue(q, i,k);
              r_c(i,k) -= w * (U_i * F);
            }
          }
        }
      } else if (bnd == 1) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          const Point &N = faceElem.QNormal(q);
          Tensor K = problem.Permeability(c()); //?
          VectorField U_diff = faceElem.VectorValue(q, u) - problem.Solution(z);
          Tensor K_DU = K * faceElem.VectorGradient(q, u);
          VectorField F = K_DU * N;
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            for (int k = 0; k < elem.Dim(); ++k) {
              auto phi_i = faceElem.VectorComponentValue(q, i, k);
              auto NDphi_i = (K * faceElem.VectorRowGradient(q, i, k)) * N;
              r_c(i, k) -= w * (F * phi_i +
                                U_diff * NDphi_i * sign - U_diff * phi_i * scaledPenalty);
            }
          }
        }
      }
    } else {
      cell cf = r.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      int f1 = r.find_neighbour_face_id(c.Face(f), cf);

      const DGVectorFieldFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);
      
      DGSizedRowBndValues r_cf(r, *cf, otherFaceElem.NodalPoints());

      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        const Point &N = faceElem.QNormal(q);
        Tensor K = problem.Permeability(c());//?
        auto U = faceElem.VectorValue(q, u);
        auto K_DU = K * faceElem.VectorGradient(q, u);
        const Point &Qf_c = faceElem.QPoint(q);
        int q1 = faceElem.findQPointID(otherFaceElem, Qf_c);
        Tensor K_1 = problem.Permeability(cf());//?
        auto U_1 = otherFaceElem.VectorValue(q1, u);
        auto K_DU_1 = K_1 * otherFaceElem.VectorGradient(q1, u);
        for (int i = 0; i < faceElem.NodalPoints(); ++i) {
          for (int k = 0; k < elem.Dim(); ++k) {
            auto phi_i = faceElem.VectorComponentValue(q, i, k);
            auto NDphi_i = (K * faceElem.VectorRowGradient(q, i, k)) * N;
            r_c(i, k) += w * ((scaledPenalty * U - 0.5 * (K_DU * N)) * phi_i
                              - 0.5 * U * NDphi_i * sign);
            r_c(i, k) += w * ((-scaledPenalty * U_1 - 0.5 * (K_DU_1 * N)) * phi_i
                              + 0.5 * U_1 * NDphi_i * sign);
          }
        }
        for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
          for (int l = 0; l < elem.Dim(); ++l){
          auto phi_j = otherFaceElem.VectorComponentValue(q1, j,l);
          auto NDphi_j = (K_1 * otherFaceElem.VectorRowGradient(q1, j,l)) * N;
          r_cf(j,l) += w * ((0.5 * K_DU_1 * N + scaledPenalty * U_1) * phi_j
                             + 0.5 * U_1 * NDphi_j * sign);
          r_cf(j,l)  += w * ((0.5 * K_DU * N - scaledPenalty * U) * phi_j
                             - 0.5 * U * NDphi_j * sign);
          }
        }
      }
    }
  }
}

void DGVectorValuedAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {
  
  const DGVectorFieldElement &elem = elementPool.Get(u, *c);
  
  DGSizedRowEntries A_c(A, *c, elem.NodalPoints(), *c, elem.NodalPoints());

  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Tensor K = problem.Permeability(elem.QPoint(q));
    for (int i = 0; i < elem.NodalPoints(); ++i) {
      for (int k = 0; k < elem.Dim(); ++k) {
        auto K_DU_i = K * elem.VectorRowGradient(q, i,k);
        for (int j = 0; j < elem.NodalPoints(); ++j) {
          for (int l = 0; l < elem.Dim(); ++l) {
            auto DU_j = elem.VectorRowGradient(q, j,l);
            A_c(i, j, k, l) += w * Frobenius(K_DU_i,DU_j);
          }
        }
      }
    }
  }
  BFParts bnd(u.GetMesh(), *c);
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar s =   (pow(degree,2)*penalty)/  u.GetMesh().MaxMeshWidth();

    const DGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (u.OnBoundary(*c, f)) {
      if (bnd[f] == 1) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          const Point &N = faceElem.QNormal(q);
          Tensor K = problem.Permeability(c());//?
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            for (int k = 0; k < elem.Dim(); ++k) {
              auto phi_i = faceElem.VectorComponentValue(q, i, k);
              VectorField NDphi_i = (K * faceElem.VectorRowGradient(q, i, k)) * N;
              for (int j = 0; j < faceElem.NodalPoints(); ++j) {
                for (int l = 0; l < elem.Dim(); ++l) {
                  auto phi_j = faceElem.VectorComponentValue(q, j, l);
                  VectorField NDphi_j = (K * faceElem.VectorRowGradient(q, j, l)) * N;
                  A_c(i, j, k, l) -= w * (NDphi_i * phi_j * sign +
                                          phi_i * NDphi_j - s * phi_i * phi_j);
                }
              }
            }
          }
        }
      }
    } else {
      cell cf = u.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);

      const DGVectorFieldFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

      DGSizedRowEntries A_cf(A, *c, elem.NodalPoints(), *cf, otherFaceElem.NodalPoints());
      DGSizedRowEntries A_fc(A, *cf, otherFaceElem.NodalPoints(), *c, elem.NodalPoints());
      DGSizedRowEntries A_ff(A, *cf, otherFaceElem.NodalPoints(), *cf, otherFaceElem.NodalPoints());


      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        const Point &N = faceElem.QNormal(q);
        Tensor K = problem.Permeability(c());//?
        // Tensor K = problem.VectorPermeability(*cf());
        const Point &Qf_c = faceElem.QPoint(q);
        int q1 = otherFaceElem.findQPointID(faceElem, Qf_c);
        Tensor K_1 = problem.Permeability(cf());
        for (int i = 0; i < faceElem.NodalPoints(); ++i) {
          for (int k = 0; k < elem.Dim(); ++k) {
            auto phi_i = faceElem.VectorComponentValue(q, i, k);
            VectorField NDphi_i = (K * faceElem.VectorRowGradient(q, i, k)) * N;
            for (int j = 0; j < faceElem.NodalPoints(); ++j) {
              for (int l = 0; l < elem.Dim(); ++l) {
                auto phi_j = faceElem.VectorComponentValue(q, j, l);
                VectorField NDphi_j = (K * faceElem.VectorRowGradient(q, j, l)) * N;
                A_c(i, j, k, l) += w * (-0.5 * NDphi_i * phi_j * sign
                                  - 0.5 * phi_i * NDphi_j
                                  + s * phi_i * phi_j);
              }
            }
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              for (int l = 0; l < elem.Dim(); ++l) {
                auto phi_j = otherFaceElem.VectorComponentValue(q1, j, l);
                VectorField NDphi_j = (K_1 * otherFaceElem.VectorRowGradient(q1, j, l)) * N;
                A_cf(i, j, k, l) += w * (0.5 * NDphi_i * phi_j * sign
                                   - phi_i * 0.5 * NDphi_j
                                   - s * phi_i * phi_j);
                A_fc(j, i, l, k) += w * (0.5 * NDphi_i * phi_j
                                   - phi_i * 0.5 * NDphi_j * sign
                                   - s * phi_i * phi_j);
              }
            }
          }
        }
        for (int i = 0; i < otherFaceElem.NodalPoints(); ++i) {
          for (int k = 0; k < elem.Dim(); ++k) {
            auto phi_i = otherFaceElem.VectorComponentValue(q1, i, k);
            VectorField NDphi_i = (K_1 * otherFaceElem.VectorRowGradient(q1, i, k)) * N;
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              for (int l = 0; l < elem.Dim(); ++l) {
                auto phi_j = otherFaceElem.VectorComponentValue(q1, j, l);
                VectorField NDphi_j = (K_1 * otherFaceElem.VectorRowGradient(q1, j, l)) * N;
                A_ff(i, j, k, l) += w * (0.5 * NDphi_i * phi_j * sign
                                   + phi_i * 0.5 * NDphi_j
                                   + s * phi_i * phi_j);
              }
            }
          }
        }
      }
    }
  }
}

double DGVectorValuedAssemble::EnergyError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      Tensor K = problem.Permeability(z);
      Tensor IK = Invert(K);
      Tensor diff = (K * elem.VectorGradient(q, u) - problem.Flux(z));
      err += w * Frobenius(IK * diff , diff);
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double DGVectorValuedAssemble::L2(const Vector &u) const {
  double l2 = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    
    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      l2 += w * U * U;
    }
  }
  return sqrt(PPM->SumOnCommSplit(l2, u.CommSplit()));
}

double DGVectorValuedAssemble::H1(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    
    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      Tensor DU = elem.VectorGradient(q, u);
      err += w * (U * U) + w * Frobenius(DU,DU);
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double DGVectorValuedAssemble::L2Error(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    
    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      auto Sol = problem.Solution(elem.QPoint(q));
      auto diff = U -Sol;
      err += w * ((diff) * (diff));
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double DGVectorValuedAssemble::L2CellAvgError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    
    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

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
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double DGVectorValuedAssemble::MaxError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    
    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField p = elem.VectorValue(q, u);
      auto solution = problem.Solution(elem.QPoint(q));
      for ( int k = 0; k< elem.Dim(); ++ k){
        auto diff = abs(p[k] - solution[k]);
        err = max(err,diff);
      }
    }
  }
  return PPM->Max(err, u.CommSplit());
}

double DGVectorValuedAssemble::FluxError(const Vector &u) const {
  double flux_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] == 2) {

        const DGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, face);

        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Tensor G = problem.Flux(faceElem.QPoint(q));
          Tensor K = problem.Permeability(faceElem.QPoint(q));
          Tensor DU = faceElem.VectorGradient(q, u);
          auto F = (K * DU - G) * faceElem.QNormal(q);
          flux_error += w * F * F;
        }
      } else if (bnd[face] == 0) {
        DGVectorFieldFaceElement faceElem(u, *c, face);
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

double DGVectorValuedAssemble::FaceError(const Vector &u) const {
  double face_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    for (int face = 0; face < c.Faces(); ++face) {

      const DGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

void DGVectorValuedAssemble::SetExactSolution(Vector &u) const {
  u = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    auto index = u.find_row(c());
    
    const DGVectorFieldElement &elem = elementPool.Get(u, *c);

    for(int i=0;i<elem.NodalPoints(); ++i){
      auto nodalPointValue = problem.Solution(elem.NodalPoint(i));
      for(int k=0;k<elem.Dim();++k){
        u(index, k * elem.NodalPoints() + i) = nodalPointValue[k];
      }
    }
  }
  u.Accumulate();
}
