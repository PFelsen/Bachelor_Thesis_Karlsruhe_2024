#include "EGVectorValuedAssemble.hpp"
#include "MixedEGVectorFieldElement.hpp"
#include "MixedEGElement.hpp"


const char *EGVectorValuedAssemble::Name() const {
  return "EGVectorValuedAssemble";
}

void EGVectorValuedAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);

  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

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


double EGVectorValuedAssemble::Energy(const Vector &u) const {
  double energy = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

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

void EGVectorValuedAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {

  const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

  MixedRowValues r_c(r, *c);

  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Point z = elem.QPoint(q);
    auto K = problem.Permeability(z);
    auto f = problem.Load(z, c);
    auto K_DU = K * elem.VectorGradient(q, u);
    for (int i = 0; i < elem.ValueSize(); ++i) {
      for (int k = 0; k < elem.Dim(); ++k) {
        auto phi_ik = elem.VectorComponentValue(q, i, k);
        auto Dphi_ik = elem.VectorRowGradient(q, i, k);
        r_c(elem.ValueIndex(), i, k) +=
            w * (Frobenius(K_DU, Dphi_ik) - f * phi_ik);
      }
    }
    auto theta = elem.PenaltyVectorValue(q);
    auto Dtheta = elem.PenaltyVectorGradient(q);
    r_c(elem.PenaltyIndex(), 0, 0) +=
        w * (Frobenius(K_DU, Dtheta) - f * theta);
  }
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar s = (pow(degree,2)*penalty)/  u.GetMesh().MaxMeshWidth();

    const MixedEGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (r.OnBoundary(*c, f)) {
      int bnd = elem.Bnd(f);
      if (bnd == 2) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Point z = faceElem.QPoint(q);
          VectorField N = faceElem.QNormal(q);
          VectorField F = problem.Flux(z) * N;
          for (int i = 0; i < faceElem.ValueSize(); ++i) {
            for (int k = 0; k < faceElem.Dim(); ++k) {
              auto phi_ik = faceElem.VectorComponentValue(q, i, k);
              r_c(faceElem.ValueIndex(), i, k) -= w * phi_ik * F;
            }
          }
          auto theta = faceElem.PenaltyVectorValue(q);
          r_c(faceElem.PenaltyIndex(), 0, 0) -= w * theta * F;
        }
      } else if (bnd == 1) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          VectorField N = faceElem.QNormal(q);
          Tensor K = problem.Permeability(c());
          VectorField u_diff = faceElem.VectorValue(q, u) - problem.Solution(z);
          Tensor K_DU = K * faceElem.VectorGradient(q, u);
          VectorField F = K_DU * N;

          VectorField theta = faceElem.PenaltyVectorValue(q);
          VectorField Dtheta = (K * faceElem.PenaltyVectorGradient(q)) * N;
          for (int i = 0; i < faceElem.ValueSize(); ++i) {
            for (int k = 0; k < faceElem.Dim(); ++k) {
              auto phi_ik = faceElem.VectorComponentValue(q, i, k);
              auto NDphi_ik = (K * faceElem.VectorRowGradient(q, i, k)) * N;
              r_c(faceElem.ValueIndex(), i, k) -= w * (
                  phi_ik * F -
                  u_diff * (sign * NDphi_ik +
                            s * (K * N) * N * phi_ik));
            }
          }

          r_c(faceElem.PenaltyIndex(), 0, 0) -= w * (theta * F - u_diff * (sign * Dtheta +
                                                                           s *
                                                                           (K * N) * N * theta));
        }
      }
    }else {
        cell cf = r.find_neighbour_cell(c, f);
        if (cf() < c()) continue;
        int f1 = r.find_neighbour_face_id(c.Face(f), cf);

        const MixedEGVectorFieldFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

        MixedRowValues r_cf(r, *cf);

        Tensor K = problem.Permeability(c()); // kappa_+
        Tensor K_1 = problem.Permeability(cf()); //kappa_-/ kappa_e
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Point z = faceElem.QPoint(q);
          VectorField N = faceElem.QNormal(q);

          int q1 = faceElem.FindQPointID(otherFaceElem, z);
          double K_e = 2 * ((K * N) * N) * ((K_1 * N) * N) / ((K * N) * N + (K_1 * N) * N);

          VectorField K_DU = faceElem.VectorGradient(q, u) * N; // kappa grad v * N
          VectorField U = faceElem.VectorValue(q, u); // v
          VectorField K_DU_1 = (otherFaceElem.VectorGradient(q1, u)) * N; // kappa grad u_nghbr
          VectorField U_1 = otherFaceElem.VectorValue(q1, u);


          for (int i = 0; i < faceElem.ValueSize(); ++i) {
            for (int k = 0; k < faceElem.Dim(); ++k) {
              auto phi_ik = faceElem.VectorComponentValue(q, i, k);
              auto NDphi_ik = (faceElem.VectorRowGradient(q, i, k)) * N;

              r_c(faceElem.ValueIndex(), i, k) +=
                  w * K_e * (s * U * phi_ik - 0.5 * K_DU * phi_ik + 0.5 * sign * U * NDphi_ik);
              r_c(faceElem.ValueIndex(), i, k) +=
                  w * K_e  *
                  (-s * U_1 * phi_ik - 0.5 * K_DU_1 * phi_ik - 0.5 * sign * U_1 * NDphi_ik);
            }
          }

          for (int j = 0; j < otherFaceElem.ValueSize(); ++j) {
            for (int l = 0; l < otherFaceElem.Dim(); ++l) {
              auto phi_jl = otherFaceElem.VectorComponentValue(q1, j, l);
              auto NDphi_jl = (otherFaceElem.VectorRowGradient(q1, j, l) ) * N;
              r_cf(otherFaceElem.ValueIndex(), j,l) +=
                  w * K_e  *
                  (s * U_1 * phi_jl + 0.5 * K_DU_1 * phi_jl - 0.5 * sign * U_1 * NDphi_jl);
              r_cf(otherFaceElem.ValueIndex(), j,l) +=
                  w * K_e  * (-s * U * phi_jl + 0.5 * K_DU * phi_jl + 0.5 * sign * U * NDphi_jl);
            }
          }

          auto theta_i = faceElem.PenaltyVectorValue(q);
          auto Dtheta_i = (faceElem.PenaltyVectorGradient(q)) * N;
          r_c(faceElem.PenaltyIndex(), 0, 0) +=
              w * K_e * (s * U * theta_i  - 0.5 * K_DU * theta_i + 0.5 * sign * U * Dtheta_i);
          r_c(faceElem.PenaltyIndex(), 0, 0) += w * K_e  * (-s * U_1 * theta_i - 0.5 * K_DU_1 *
                                                                                theta_i -
                                                           0.5 * sign * U_1 * Dtheta_i);

          auto theta_j = otherFaceElem.PenaltyVectorValue(q1);
          auto Dtheta_j = (otherFaceElem.PenaltyVectorGradient(q1)) * N;
          r_cf(otherFaceElem.PenaltyIndex(), 0, 0) += w * K_e * (s * U_1 * theta_j  + 0.5 * K_DU_1 *
                                                                                     theta_j -
                                                                 0.5 * sign * U_1 * Dtheta_j);
          r_cf(otherFaceElem.PenaltyIndex(), 0, 0) +=
              w * K_e  * (-s * U * theta_j + 0.5 * K_DU * theta_j + 0.5 * sign * U * Dtheta_j);

        }
      }
    }
}

void EGVectorValuedAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {

  const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

  MixedRowEntries A_c(A, *c);

  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Point z = elem.QPoint(q);
    Tensor K = problem.Permeability(z);

    auto KDtheta = K * elem.PenaltyVectorGradient(q);
    auto Dtheta = elem.PenaltyVectorGradient(q);

    for (int i = 0; i < elem.ValueSize(); ++i) {
      for (int k = 0; k < elem.Dim(); ++k) {
        auto Dphi_ik = K * elem.VectorRowGradient(q, i, k);
        for (int j = 0; j < elem.ValueSize(); ++j) {
          for (int l = 0; l < elem.Dim(); ++l) {
            auto Dphi_jl = elem.VectorRowGradient(q, j, l);
            A_c(elem.ValueIndex(), elem.ValueIndex(), i, j, k, l) +=
                w * Frobenius(Dphi_ik, Dphi_jl);
          }
        }
        A_c(elem.PenaltyIndex(), elem.ValueIndex(), 0, i, 0, k) +=
            w * Frobenius(KDtheta, Dphi_ik);

        A_c(elem.ValueIndex(), elem.PenaltyIndex(), i, 0, k, 0) +=
            w * Frobenius(K * Dphi_ik, Dtheta);
      }
    }

    A_c(elem.PenaltyIndex(), elem.PenaltyIndex(), 0, 0, 0, 0) +=
        w * Frobenius(KDtheta, Dtheta);
  }
  BFParts bnd(u.GetMesh(), *c);
  for (int f = 0; f < c.Faces(); ++f) {
    Scalar s = ((pow(degree,2)*penalty)/  u.GetMesh().MaxMeshWidth());

    const MixedEGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, f);

    if (u.OnBoundary(*c, f)) {
    if (bnd[f] == 1) {
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        VectorField N = faceElem.QNormal(q);
        Tensor K = problem.Permeability(c());

        auto theta = faceElem.PenaltyVectorValue(q);
        auto Dtheta = ( K *  faceElem.PenaltyVectorGradient(q)) * N;

        for (int i = 0; i < faceElem.ValueSize(); ++i) {
          for (int k = 0; k < faceElem.Dim(); ++k) {
            auto phi_ik = faceElem.VectorComponentValue(q, i, k);
            auto NDphi_ik = (K * faceElem.VectorRowGradient(q, i, k)) * N;
            for (int j = 0; j < faceElem.ValueSize(); ++j) {
              for (int l = 0; l < faceElem.Dim(); ++l) {
                auto phi_jl = faceElem.VectorComponentValue(q, j, l);
                auto NDphi_jl = (K * faceElem.VectorRowGradient(q, j, l)) * N;

                A_c(faceElem.ValueIndex(), faceElem.ValueIndex(), i, j, k, l) -=
                    w * (phi_ik * NDphi_jl -
                         sign * phi_jl * NDphi_ik -
                         s * ((K * N) * N) * phi_jl *
                          phi_ik);
              }
            }

            A_c(faceElem.ValueIndex(), faceElem.PenaltyIndex(), i, 0, k, 0) -=
                w * (phi_ik * Dtheta -
                     sign * theta * NDphi_ik -
                     s * ((K * N) * N) * theta * phi_ik);
            A_c(faceElem.PenaltyIndex(), faceElem.ValueIndex(), 0, i, 0, k) -=
                w * (theta * NDphi_ik -
                     sign * phi_ik * Dtheta -
                     s * ((K * N) * N) * phi_ik * theta);
          }
        }

        A_c(faceElem.PenaltyIndex(), faceElem.PenaltyIndex(), 0, 0, 0, 0) -=
            w * (theta * Dtheta -
                 sign * theta * Dtheta -
                 s * ((K * N) * N) * theta * theta);
      }
    }
    } else {
      cell cf = u.find_neighbour_cell(c, f);
      if (cf() < c()) continue;
      int f1 = u.find_neighbour_face_id(c.Face(f), cf);
      MixedEGVectorFieldElement elem_1(u, *cf);

      const MixedEGVectorFieldFaceElement &otherFaceElem = faceElementPool.Get(u, *cf, f1);

      MixedRowEntries A_cf(A, *c, *cf);
      MixedRowEntries A_fc(A, *cf, *c);
      MixedRowEntries A_ff(A, *cf);

      Tensor K = problem.Permeability(c());
      Tensor K_1 = problem.Permeability(cf());
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        Point z = faceElem.QPoint(q);
        VectorField N = faceElem.QNormal(q);

        int q1 = otherFaceElem.FindQPointID(faceElem, z);

        double K_e = 2 * ((K * N) * N) * ((K_1 * N) * N) / ((K * N) * N + (K_1 * N) * N);
        Scalar s = (pow(degree,2)*penalty)/  u.GetMesh().MaxMeshWidth();

        auto theta = faceElem.PenaltyVectorValue(q);
        auto Dtheta = (faceElem.PenaltyVectorGradient(q)) * N;

        auto theta_nghbr = otherFaceElem.PenaltyVectorValue(q1);
        auto Dtheta_nghbr = (otherFaceElem.PenaltyVectorGradient(q1)) * N;

        for (int i = 0; i < faceElem.ValueSize(); ++i) {
          for (int k = 0; k < faceElem.Dim(); ++k) {
            auto phi_ik = faceElem.VectorComponentValue(q, i, k);
            auto NDphi_ik = (faceElem.VectorRowGradient(q, i, k)) * N;
            for (int j = 0; j < faceElem.ValueSize(); ++j) {
              for (int l = 0; l < faceElem.Dim(); ++l) {
                auto phi_jl = faceElem.VectorComponentValue(q, j, l);
                auto NDphi_jl = (faceElem.VectorRowGradient(q, j, l)) * N;

                A_c(faceElem.ValueIndex(), faceElem.ValueIndex(), i, j, k, l) +=
                    w * K_e  * (s * phi_ik * phi_jl - 0.5 * phi_ik * NDphi_jl
                                   + 0.5 * sign * NDphi_ik * phi_jl);
              }
            }

            A_c(faceElem.ValueIndex(), faceElem.PenaltyIndex(), i, 0, k, 0) +=
                w * K_e  * (s * phi_ik * theta - 0.5 * phi_ik * Dtheta
                               + 0.5 * sign * NDphi_ik * theta);
            A_c(faceElem.PenaltyIndex(), faceElem.ValueIndex(), 0, i, 0, k) +=
                w * K_e * (s * theta * phi_ik  - 0.5 * theta * NDphi_ik
                              + 0.5 * sign * Dtheta * phi_ik);

            for (int j = 0; j < otherFaceElem.ValueSize(); ++j) {
              for (int l = 0; l < otherFaceElem.Dim(); ++l) {
                auto phi_jl = otherFaceElem.VectorComponentValue(q1, j, l);
                auto NDphi_jl = (otherFaceElem.VectorRowGradient(q1, j, l)) * N;

                A_cf(faceElem.ValueIndex(), otherFaceElem.ValueIndex(), i, j, k, l) +=
                    w * K_e  * (-s * phi_ik * phi_jl
                  - phi_ik * 0.5 * NDphi_jl
                                   - 0.5 * NDphi_ik * phi_jl * sign);
                A_fc(otherFaceElem.ValueIndex(), faceElem.ValueIndex(), j, i, l, k) +=
                    w * K_e  * (-s * phi_ik * phi_jl + 0.5 * NDphi_ik * phi_jl
                                   + phi_ik * 0.5 * NDphi_jl * sign);
              }
            }

            A_cf(faceElem.ValueIndex(), otherFaceElem.PenaltyIndex(), i, 0, k, 0) +=
                w * K_e  * (-s * phi_ik * theta_nghbr
              - phi_ik * 0.5 * Dtheta_nghbr
                                - 0.5 * NDphi_ik * theta_nghbr * sign);
            A_fc(otherFaceElem.PenaltyIndex(), faceElem.ValueIndex(), 0, i, 0, k) +=
                w * K_e  * (-s * theta_nghbr * phi_ik + 0.5 * NDphi_ik * theta_nghbr
                               + phi_ik * 0.5 * Dtheta_nghbr * sign);
          }
        }


        for (int i = 0; i < otherFaceElem.ValueSize(); ++i) {
          for (int k = 0; k < otherFaceElem.Dim(); ++k) {
            auto phi_ik = otherFaceElem.VectorComponentValue(q1, i, k);
            auto NDphi_ik = (otherFaceElem.VectorRowGradient(q1, i, k)) * N;

            A_cf(faceElem.PenaltyIndex(), otherFaceElem.ValueIndex(), 0, i, 0, k) +=
                w * K_e  *
                (-s * theta * phi_ik - theta * 0.5 * NDphi_ik - 0.5 * Dtheta * phi_ik * sign);
            A_fc(otherFaceElem.ValueIndex(), faceElem.PenaltyIndex(), i, 0, k, 0) +=
                w * K_e  *
                (-s * theta * phi_ik + 0.5 * Dtheta * phi_ik + 0.5 * theta * NDphi_ik * sign);


            for (int j = 0; j < otherFaceElem.ValueSize(); ++j)
              for (int l = 0; l < otherFaceElem.Dim(); ++l) {
                auto phi_jl = otherFaceElem.VectorComponentValue(q1, j, l);
                auto NDphi_jl = (otherFaceElem.VectorRowGradient(q1, j, l)) * N;

                A_ff(otherFaceElem.ValueIndex(), otherFaceElem.ValueIndex(), i, j, k, l) +=
                    w * K_e * (s * phi_ik * phi_jl  + phi_ik * 0.5 * NDphi_jl
                                  - 0.5 * NDphi_ik * phi_jl * sign);
              }

            A_ff(otherFaceElem.ValueIndex(), otherFaceElem.PenaltyIndex(), i, 0, k, 0) +=
                w * K_e  * (s * phi_ik * theta_nghbr + phi_ik * 0.5 * Dtheta_nghbr
                               - 0.5 * NDphi_ik * theta_nghbr * sign);
            A_ff(otherFaceElem.PenaltyIndex(), otherFaceElem.ValueIndex(), 0, i, 0, k) +=
                w * K_e  * (s * theta_nghbr * phi_ik + theta_nghbr * 0.5 * NDphi_ik
                               - 0.5 * Dtheta_nghbr * phi_ik * sign);
          }
        }


        A_c(faceElem.PenaltyIndex(), faceElem.PenaltyIndex(), 0, 0, 0, 0) +=
            w * K_e * (s * theta * theta  - 0.5 * theta * Dtheta
                          + 0.5 * sign * Dtheta * theta);

        A_cf(faceElem.PenaltyIndex(), otherFaceElem.PenaltyIndex(), 0, 0, 0, 0) +=
            w * K_e  *
            (-s * theta * theta_nghbr - theta * 0.5 * Dtheta_nghbr -
             0.5 * Dtheta * theta_nghbr * sign);
        A_fc(otherFaceElem.PenaltyIndex(), faceElem.PenaltyIndex(), 0, 0, 0, 0) +=
            w * K_e  *
            (-s * theta * theta_nghbr + 0.5 * Dtheta * theta_nghbr +
             0.5 * theta * Dtheta_nghbr * sign);

        A_ff(otherFaceElem.PenaltyIndex(), otherFaceElem.PenaltyIndex(), 0, 0, 0, 0) +=
            w * K_e * (s * theta_nghbr * theta_nghbr  + theta_nghbr * 0.5 * Dtheta_nghbr
                          - 0.5 * Dtheta_nghbr * theta_nghbr * sign);
      }
    }
  }
}

double EGVectorValuedAssemble::EnergyError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      Tensor K = problem.Permeability(z);
      Tensor IK = Invert(K);
      Tensor diff = (K * elem.VectorGradient(q, u) -
                     problem.Flux(z));
      err += w * Frobenius(IK * diff, diff); // (IK * diff) * diff
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double EGVectorValuedAssemble::L2(const Vector &u) const {
  double l2 = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      l2 += w * U * U;
    }
  }
  return sqrt(PPM->SumOnCommSplit(l2, u.CommSplit()));
}

double EGVectorValuedAssemble::H1(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      Tensor DU = elem.VectorGradient(q, u);
      err += w * (U * U) + w * Frobenius(DU, DU);
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double EGVectorValuedAssemble::L2Error(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField U = elem.VectorValue(q, u);
      auto Sol = problem.Solution(elem.QPoint(q));
      auto diff = U - Sol;
      err += w * ((diff) * (diff));
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

double EGVectorValuedAssemble::L2CellAvgError(const Vector &u) const {
  double err = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

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

double EGVectorValuedAssemble::MaxError(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const MixedEGVectorFieldElement &elem = elementPool.Get(u, *c);

    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField U = elem.VectorValue(q, u);
      auto solution = problem.Solution(elem.QPoint(q));
      for (int k = 0; k < c.dim(); ++k) {
        auto diff = abs(U[k] - solution[k]);
        err = max(err, diff);
      }
    }
  }
  return PPM->Max(err, u.CommSplit());
}

double EGVectorValuedAssemble::FluxError(const Vector &u) const {
  double flux_error = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), c);
    if (!bnd.onBnd()) continue;
    for (int face = 0; face < c.Faces(); ++face) {
      if (bnd[face] == 2) {

        const MixedEGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, face);

        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          Tensor G = problem.Flux(faceElem.QPoint(q));
          Tensor K = problem.Permeability(faceElem.QPoint(q));
          Tensor DU = faceElem.VectorGradient(q, u);
          auto F = (K * DU - G) * faceElem.QNormal(q);
          flux_error += w * F * F;
        }
      } else if (bnd[face] == 0) {

        const MixedEGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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

double EGVectorValuedAssemble::FaceError(const Vector &u) const {
  double face_error = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    for (int face = 0; face < c.Faces(); ++face) {

      const MixedEGVectorFieldFaceElement &faceElem = faceElementPool.Get(u, *c, face);

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


void EGVectorValuedAssemble::SetExactSolution(Vector &uEx) const {
  for (row r = uEx.rows(); r != uEx.rows_end(); ++r) {
    if (r.size() < uEx.dim()) continue;
    VectorField solution = problem.Solution(r());
    for (int k = 0; k < uEx.dim(); ++k) {
      uEx(r, k) = solution[k];
    }
  }
}
