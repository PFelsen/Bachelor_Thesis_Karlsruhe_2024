#include "DGAcousticAssemble.hpp"
#include "Plotting.hpp"
#include "LagrangeDiscretization.hpp"
#include <functional>

std::unique_ptr<DGAcousticElement>
createAcousticElement(const VectorMatrixBase &u, const Cell &c, int damping) {
  return std::make_unique<DGAcousticElement>(u, c, damping);
}

std::function<std::unique_ptr<DGAcousticElement>(const VectorMatrixBase &, const Cell &)>
DGAcousticElementGenerator(int damping) {
  return std::bind(createAcousticElement, std::placeholders::_1, std::placeholders::_2, damping);
}

std::unique_ptr<DGAcousticFaceElement>
createAcousticFaceElement(const VectorMatrixBase &u, const Cell &c, int face, int damping) {
  return std::make_unique<DGAcousticFaceElement>(u, c, face, damping);
}

std::function<std::unique_ptr<DGAcousticFaceElement>(const VectorMatrixBase &, const Cell &, int)>
DGAcousticFaceElementGenerator(int damping) {
  return std::bind(createAcousticFaceElement, std::placeholders::_1, std::placeholders::_2,
                   std::placeholders::_3, damping);
}

void DGAcousticAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u, DGAcousticElementGenerator(problem.NumDamping()));
  faceElementPool.Initialize(u, DGAcousticFaceElementGenerator(problem.NumDamping()));
}

void DGAcousticAssemble::BoundaryFlux(DGRowEntries &A_c, const DGAcousticFaceElement &felem) const {
  switch (problem.BndID(felem.GetCell()())) {
    case 0:
    case 2:
      NeumannBoundaryFlux(A_c, felem);
      break;
    case 1:
      DirichletBoundaryFlux(A_c, felem);
      break;
    case 3:
      RobinBoundaryFlux(A_c, felem);
      break;
  }
}

void DGAcousticAssemble::RobinBoundaryFlux(DGRowEntries &A_c,
                                           const DGAcousticFaceElement &felem) const {
  for (int q = 0; q < felem.nQ(); ++q) {
    double w = felem.QWeight(q);
    VectorField N = felem.QNormal(q);
    Point z = felem.QPoint(q);
    int nV = felem.V_dimension();
    int nP = felem.shape_size();
    double Z_K = sqrt(problem.Kappa(felem.GetCell(), z)
                      * problem.Rho(felem.GetCell(), z));
    double alpha = upwind ? 1.0 / Z_K : 0.0;
    for (int j = 0; j < nV; ++j) {
      Scalar VN_j = felem.Velocity(q, j) * N;
      for (int i = 0; i < nV; ++i) {
        Scalar VN_i = felem.Velocity(q, i) * N;
        A_c(i, j) -= w * Z_K * VN_i * VN_j;
      }
      for (int i = 0; i < felem.shape_size(); ++i) {
        double P_i = felem.Pressure(q, i);
        for (int n = 0; n < felem.p_dimension(); ++n)
          A_c(nV + n * nP + i, j) -= (fwd ? 1.0 : (-1.0)) * w * VN_j * P_i;
      }
    }
    for (int j = 0; j < nP; ++j) {
      double P_j = felem.Pressure(q, j);
      for (int i = 0; i < nV; ++i) {
        Scalar VN_i = felem.Velocity(q, i) * N;
        for (int n = 0; n < felem.p_dimension(); ++n)
          A_c(i, nV + n * nP + j) += w * VN_i * P_j;
      }
      for (int i = 0; i < nP; ++i) {
        double P_i = felem.Pressure(q, i);
        for (int n = 0; n < felem.p_dimension(); ++n)
          for (int k = 0; k < felem.p_dimension(); ++k)
            A_c(nV + n * nP + i, nV + k * nP + j) -=
              w * alpha * P_i * P_j;
      }
    }
  }
}

void DGAcousticAssemble::NeumannBoundaryFlux(DGRowEntries &A_c,
                                             const DGAcousticFaceElement &felem) const {
  for (int q = 0; q < felem.nQ(); ++q) {
    int nV = felem.V_dimension();
    int nP = felem.shape_size();
    double w = felem.QWeight(q);
    VectorField N = felem.QNormal(q);
    Point z = felem.QPoint(q);
    double Z_K = sqrt(problem.Kappa(felem.GetCell(), z) *
                      problem.Rho(felem.GetCell(), z));
    for (int j = 0; j < nV; ++j) {
      Scalar VN_j = felem.Velocity(q, j) * N;
      for (int i = 0; i < nV; ++i) {
        Scalar VN_i = felem.Velocity(q, i) * N;
        A_c(i, j) -= w * Z_K * VN_i * VN_j;
      }
      for (int i = 0; i < felem.shape_size(); ++i) {
        double P_i = felem.Pressure(q, i);
        for (int n = 0; n < felem.p_dimension(); ++n) {
          A_c(nV + n * nP + i, j) -= (fwd ? 1.0 : (-1.0)) * w * VN_j * P_i;
        }
      }
    }
  }
}

void DGAcousticAssemble::DirichletBoundaryFlux(DGRowEntries &A_c,
                                               const DGAcousticFaceElement &felem) const {
  for (int q = 0; q < felem.nQ(); ++q) {
    double w = felem.QWeight(q);
    VectorField N = felem.QNormal(q);
    Point z = felem.QPoint(q);
    int nV = felem.V_dimension();
    int nP = felem.shape_size();
    double alpha =
        upwind ? 1.0 / sqrt(problem.Kappa(felem.GetCell(), z) *
                            problem.Rho(felem.GetCell(), z)) : 0.0;
    for (int j = 0; j < nP; ++j) {
      double P_j = felem.Pressure(q, j);
      for (int i = 0; i < nV; ++i) {
        Scalar VN_i = felem.Velocity(q, i) * N;
        for (int n = 0; n < felem.p_dimension(); ++n) {
          A_c(i, nV + n * nP + j) -= (fwd ? 1.0 : (-1.0)) * w * VN_i * P_j;
        }
      }
      for (int i = 0; i < nP; ++i) {
        double P_i = felem.Pressure(q, i);
        for (int n = 0; n < felem.p_dimension(); ++n) {
          for (int k = 0; k < felem.p_dimension(); ++k) {
            A_c(nV + n * nP + i, nV + k * nP + j) -=
              w * alpha * P_i * P_j;
          }
        }
      }
    }
  }
}

void DGAcousticAssemble::assembleFlux(Matrix &A_loc, DGRowEntries &A_c,
                                      const cell &c) const {
  for (int f = 0; f < c.Faces(); ++f) {
    
    const DGAcousticFaceElement &faceElem = faceElementPool.Get(A_loc, *c, f, DGAcousticFaceElementGenerator(problem.NumDamping()));

    int nV = faceElem.V_dimension();
    int nP = faceElem.shape_size();
    int pComponents = problem.NumDamping() + 1;
    if (A_loc.OnBoundary(*c, f)) {
      if (!dampingFlux) { BoundaryFlux(A_c, faceElem); }
    } else {
      cell cf = A_loc.find_neighbour_cell(c, f);
      int f1 = A_loc.find_neighbour_face_id(c.Face(f), cf);
      
      const DGAcousticFaceElement &otherFaceElem = faceElementPool.Get(A_loc, *cf, f1, DGAcousticFaceElementGenerator(problem.NumDamping()));

      int nV_1 = otherFaceElem.V_dimension();
      int nP_1 = otherFaceElem.shape_size();
      DGRowEntries A_cf(A_loc, *c, *cf);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        VectorField N = faceElem.QNormal(q);
        Point z = faceElem.QPoint(q);
        double Z_K = sqrt(problem.Kappa(*c, z) * problem.Rho(*c, z));
        double Z_Kf = sqrt(problem.Kappa(*cf, z) * problem.Rho(*cf, z));
        double alpha1 = 1.0 / (Z_K + Z_Kf);
        double alpha2 = Z_Kf * Z_K / (Z_K + Z_Kf);
        double alpha3 = (fwd ? 1.0 : -1.0) * Z_Kf / (Z_K + Z_Kf);
        double alpha4 = (fwd ? 1.0 : -1.0) * Z_K / (Z_K + Z_Kf);
        if (!upwind) {
          alpha1 = alpha2 = 0;
          alpha3 = alpha4 = -0.5;
        }
        for (int i = 0; i < nV; ++i) {
          Scalar VN_i = faceElem.Velocity(q, i) * N;
          for (int j = 0; j < nV; ++j) {
            Scalar VN_j = faceElem.Velocity(q, j) * N;
            A_c(i, j) -= w * alpha2 * VN_i * VN_j;
          }
          for (int j = 0; j < nV_1; ++j) {
            Scalar V1N_j = otherFaceElem.Velocity(q, j) * N;
            A_cf(i, j) += w * alpha2 * VN_i * V1N_j;
          }
          for (int j = 0; j < nP; ++j) {
            double P_j = faceElem.Pressure(q, j);
            for (int n = 0; n < pComponents; ++n) {
              A_c(i, nV + n * nP + j) -= w * alpha4 * P_j * VN_i;
            }
          }
          for (int j = 0; j < nP_1; ++j) {
            double P1_j = otherFaceElem.Pressure(q, j);
            for (int n = 0; n < pComponents; ++n) {
              A_cf(i, nV_1 + n * nP_1 + j) += w * alpha4 * P1_j * VN_i;
            }
          }
        }
        for (int i = 0; i < nP; ++i) {
          double P_i = faceElem.Pressure(q, i);
          for (int j = 0; j < nP; ++j) {
            double P_j = faceElem.Pressure(q, j);
            for (int n = 0; n < pComponents; ++n) {
              for (int k = 0; k < pComponents; ++k) {
                A_c(nV + n * nP + i, nV + k * nP + j) -=
                  w * alpha1 * P_i * P_j;
              }
            }
          }
          for (int j = 0; j < nP_1; ++j) {
            double P1_j = otherFaceElem.Pressure(q, j);
            for (int n = 0; n < pComponents; ++n) {
              for (int k = 0; k < pComponents; ++k) {
                A_cf(nV + n * nP + i, nV_1 + k * nP + j) +=
                  w * alpha1 * P_i * P1_j;
              }
            }
          }
          for (int j = 0; j < nV; ++j) {
            Scalar VN_j = faceElem.Velocity(q, j) * N;
            for (int n = 0; n < pComponents; ++n) {
              A_c(nV + n * nP + i, j) -= w * alpha3 * VN_j * P_i;
            }
          }
          for (int j = 0; j < nV_1; ++j) {
            Scalar V1N_j = otherFaceElem.Velocity(q, j) * N;
            for (int n = 0; n < pComponents; ++n) {
              A_cf(nV + n * nP + i, j) += w * alpha3 * V1N_j * P_i;
            }
          }
        }
      }
    }
  }
}

void DGAcousticAssemble::MassMatrix(Matrix &M) const {
  M = 0;
  for (cell c = M.cells(); c != M.cells_end(); ++c) {

    const DGAcousticElement &elem = elementPool.Get(M, *c, DGAcousticElementGenerator(problem.NumDamping()));

    DGRowEntries M_c(M, *c, *c);
    int nV = elem.V_dimension();
    int nP = elem.shape_size();
    int p_dim = elem.p_dimension();
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      double rho = problem.Rho(*c, z);
      for (int i = 0; i < nV; ++i) {
        VectorField V_i = elem.Velocity(q, i);
        for (int j = 0; j < nV; ++j) {
          VectorField V_j = elem.Velocity(q, j);
          M_c(i, j) += w * rho * (V_i * V_j);
        }
      }
      for (int i = 0; i < nP; ++i) {
        double P_i = elem.Pressure(q, i);
        for (int j = 0; j < nP; ++j) {
          double P_j = elem.Pressure(q, j);
          for (int n = 0; n < p_dim; ++n) {
            M_c(i + nV + n * nP, j + nV + n * nP) +=
                w * 1.0 / problem.Kappa_i(*c, z, n)* P_i * P_j;
          }
        }
      }
    }
  }
}

void DGAcousticAssemble::assembleDifferentialPartFlux(const Matrix &A,
                                                      const cell &c,
                                                      DGRowEntries &A_c) const {

  const DGAcousticElement &elem = elementPool.Get(A, *c, DGAcousticElementGenerator(problem.NumDamping()));

  int nV = elem.V_dimension();
  int nP = elem.shape_size();
  int p_dim = elem.p_dimension();
  if (problem.Damping()) Damping(elem, A_c);
  GradPressure(elem, A_c);
  DivVelocity(elem, A_c);
}

void DGAcousticAssemble::GradPressure(const DGAcousticElement &elem,
                                      DGRowEntries &A_c) const {
  int nV = elem.V_dimension();
  int nP = elem.shape_size();
  int p_dim = elem.p_dimension();
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    for (int i = 0; i < nP; ++i) {
      const VectorField &GradP_i = elem.Gradpressure(q, i);
      for (int j = 0; j < nV; ++j) {
        const VectorField &V_j = elem.Velocity(q, j);
        for (int n = 0; n < p_dim; ++n)
          A_c(j, nV + n * nP + i) += (fwd ? 1.0 : (-1.0)) * w * (GradP_i * V_j);
      }
    }
  }
}

void DGAcousticAssemble::DivVelocity(const DGAcousticElement &elem,
                                     DGRowEntries &A_c) const {
  int nV = elem.V_dimension();
  int nP = elem.shape_size();
  int p_dim = elem.p_dimension();
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    for (int i = 0; i < nV; ++i) {
      double Div_i = elem.Divergence(q, i);
      for (int j = 0; j < nP; ++j) {
        double P_j = elem.Pressure(q, j);
        for (int n = 0; n < p_dim; ++n)
          A_c(nV + n * nP + j, i) += (fwd ? 1.0 : (-1.0)) * w * Div_i * P_j;
      }
    }
  }
}

void DGAcousticAssemble::Damping(const DGAcousticElement &elem, DGRowEntries &A_c) const {
  int nV = elem.V_dimension();
  int nP = elem.shape_size();
  int p_dim = elem.p_dimension();
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    for (int i = 0; i < elem.shape_size(); ++i) {
      for (int j = 0; j < nP; ++j) {
        for (int n = 1; n < p_dim; ++n) {
          A_c(nV + n * nP + j, nV + n * nP + i) -= w * elem.Pressure(q, i) *
            elem.Pressure(q, j) * 1.0 / problem.Kappa_i(elem.GetCell(), z, n)
            * 1.0 / problem.Tau_i(elem.GetCell(), z, n);
        }
      }
    }
  }
}

void DGAcousticAssemble::RHS(double t, Vector &b) const {
  b.Clear();
  if (!problem.HasRHS()) return;
  for (cell c = b.cells(); c != b.cells_end(); ++c) {

    const DGAcousticElement &elem = elementPool.Get(b, *c, DGAcousticElementGenerator(problem.NumDamping()));

    int p_dim = elem.p_dimension();
    int Row = b.find_row(c()).Id();
    int nV = elem.V_dimension();
    int nP = elem.shape_size();
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      VectorField V_rhs(0.0);
      RVector ForceP(0.0, p_dim);
      problem.Force(V_rhs, ForceP, t, *c, z);
      for (int i = 0; i < nV; ++i) {
        VectorField V_i = elem.Velocity(q, i);
        b(Row, i) += w * V_i * V_rhs;
      }
      Point globalPointSourceLocation = Infty;
      if (problem.hasPointSource(globalPointSourceLocation)) {
        if (globalPointSourceLocation == Infty){
          THROW("globalPointSourceLocation not set in Problem!")
        }
        if (c.PointInCell(globalPointSourceLocation)) {
          for (int i = 0; i < nP; ++i) {
            b(Row, nV + i) += elem.Pressure(q, i) *
                              problem.F_PointSource(t, globalPointSourceLocation);
          }
        }
      }
      for (int i = 0; i < nP; ++i) {
        double P_i = elem.Pressure(q, i);
        for (int n = 0; n < p_dim; ++n) {
          b(Row, nV + n * nP + i) += w * P_i * ForceP[n];
        }
      }
    }
    if (problem.HasExactSolution()) {
      for (int f = 0; f < c.Faces(); ++f) {
        if (!b.OnBoundary(*c, f)) continue;

        const DGAcousticFaceElement &faceElem = faceElementPool.Get(b, *c, f, DGAcousticFaceElementGenerator(problem.NumDamping()));

        int nV = faceElem.V_dimension();
        int nP = faceElem.shape_size();
        BFParts bnd(b.GetMesh(), *c);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          VectorField N = faceElem.QNormal(q);
          Point z = faceElem.QPoint(q);
          double Z_L = sqrt(problem.Kappa(*c, z) * problem.Rho(*c, z));
          double alpha = 1 / Z_L;
          if (problem.BndID(z) == 2) {
            Scalar VN_bnd = problem.Neumann(t, *c, z, N);
            for (int i = 0; i < nV; ++i) {
              Scalar VN_i = faceElem.Velocity(q, i) * N;
              b(Row, i) += w * Z_L * VN_i * VN_bnd;
            }
            for (int i = 0; i < nP; ++i) {
              double P_i = faceElem.Pressure(q, i);
              for (int n = 0; n < p_dim; ++n) {
                b(Row, nV + n * nP + i) += w * P_i * VN_bnd;
              }
            }
          } else if (problem.BndID(z) == 1) {
            for (int i = 0; i < nV; ++i) {
              Scalar VN_i = faceElem.Velocity(q, i) * N;
              for (int n = 0; n < p_dim; ++n)
                b(Row, i) += w * VN_i * problem.Dirichlet(t, *c, z, n);
            }
            for (int i = 0; i < nP; ++i) {
              double P_i = faceElem.Pressure(q, i);
              for (int n = 0; n < p_dim; ++n) {
                for (int k = 0; k < p_dim; ++k) {
                  b(c(), nV + n * nP + i) += w * alpha * P_i
                    * problem.Dirichlet(t, *c, z, k);
                }
              }
            }
          }
        }
      }
    }
  }
}

void DGAcousticAssemble::SystemMatrix(Matrix &A) const {
  A = 0;
  for (cell c = A.cells(); c != A.cells_end(); ++c) {
    DGRowEntries A_c(A, *c, *c);
    assembleDifferentialPartFlux(A, c, A_c);
    assembleFlux(A, A_c, c);
  }
}

double DGAcousticAssemble::Energy(const Vector &u) const {
  Scalar energy = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const DGAcousticElement &elem = elementPool.Get(u, *c, DGAcousticElementGenerator(problem.NumDamping()));

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      int p_dim = elem.p_dimension();
      double rho = problem.Rho(*c, z);
      VectorField V = elem.Velocity(q, u);
      energy += w * rho * V * V;
      for (int n = 0; n < p_dim; ++n) {
        double P_i = elem.Pressure_i(q, u, n);
        energy += w * 1.0 / problem.Kappa_i(*c, z, n) * P_i * P_i;
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(energy, u.GetMesh().CommSplit()));
}

MinMaxP DGAcousticAssemble::MinMaxPressure(const Vector &u) const {
  MinMaxP minMaxP = {0.0, 0.0};
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    if (!problem.InRegionOfInterest(c())) continue;

    const DGAcousticElement &elem = elementPool.Get(u, *c, DGAcousticElementGenerator(problem.NumDamping()));
  
    double P = 0;
    double PInt = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      P += w * elem.Pressure(q, u);
    }
    double area = elem.Area();
    P *= (1.0 / area);
    PInt *= (1.0 / area);

    if (minMaxP.second < P) minMaxP.second = P;
    if (minMaxP.first > P) minMaxP.first = P;
  }
  minMaxP.second = PPM->Max(minMaxP.second, u.GetMesh().CommSplit());
  minMaxP.first = PPM->Min(minMaxP.first, u.GetMesh().CommSplit());
  if (abs(minMaxP.second) < 1e-10) minMaxP.second = 0;
  if (abs(minMaxP.first) < 1e-10) minMaxP.first = 0;
  return minMaxP;
}

AcousticNorms DGAcousticAssemble::ComputeNorms(const Vector &u) const {
  AcousticNorms norms;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    if (!problem.InRegionOfInterest(c())) continue;

    const DGAcousticElement &elem = elementPool.Get(u, *c, DGAcousticElementGenerator(problem.NumDamping()));

    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);
      int p_dim = elem.p_dimension();
      double rho = problem.Rho(*c, z);
      VectorField V = elem.Velocity(q, u);
      norms.l2 += w * V * V;
      norms.l1 += w * norm(V);
      norms.energy += w * rho * V * V;
      for (int n = 0; n < p_dim; ++n) {
        double P_i = elem.Pressure_i(q, u, n);
        norms.l2 += w * P_i * P_i;
        norms.l1 += w * abs(P_i);
        norms.energy += w * 1.0 / problem.Kappa_i(*c, z, n) * P_i * P_i;
      }
    }
  }
  norms.l2 = sqrt(PPM->SumOnCommSplit(norms.l2, u.GetMesh().CommSplit()));
  norms.l1 = sqrt(PPM->SumOnCommSplit(norms.l1, u.GetMesh().CommSplit()));
  norms.energy = sqrt(PPM->SumOnCommSplit(norms.energy, u.GetMesh().CommSplit()));
  return norms;
}

AcousticErrors DGAcousticAssemble::ComputeErrors(const Vector &u) const {
  AcousticErrors errors;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const DGAcousticElement &elem = elementPool.Get(u, *c, DGAcousticElementGenerator(problem.NumDamping()));

    int p_dim = elem.p_dimension();
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Point z = elem.QPoint(q);

      VectorField exactV = zero;
      std::vector<double> exactP(p_dim, 0.0);
      for (int i = 0; i < problem.SpaceDim(); ++i)
        exactV[i] = problem.ut(Time(), z, *c, i);
      for (int n = 0; n < p_dim; ++n)
        exactP[n] = problem.ut(Time(), z, *c, problem.SpaceDim() + n);

      VectorField V = elem.Velocity(q, u);
      V -= exactV;
      errors.l2 += w * V * V;
      errors.l1 += w * norm(V);
      errors.inf = max(errors.inf, norm(V));
      for (int n = 0; n < p_dim; ++n) {
        double P_i = elem.Pressure_i(q, u, n) - exactP[n];
        errors.l1 += w * abs(P_i);
        errors.l2 += w * P_i * P_i;
        errors.inf = max(errors.inf, abs(P_i));
      }
    }
  }
  errors.l1 = PPM->SumOnCommSplit(errors.l1, u.GetMesh().CommSplit());
  errors.l2 = sqrt(PPM->SumOnCommSplit(errors.l2, u.GetMesh().CommSplit()));
  errors.inf = PPM->Max(errors.inf, u.GetMesh().CommSplit());
  return errors;
}

void DGAcousticAssemble::InterpolatedExactSolution(Vector &u) const {
  u.Clear();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    Point z_c = c();
    row r = u.find_row(z_c);

    const DGAcousticElement &elem = elementPool.Get(u, *c, DGAcousticElementGenerator(problem.NumDamping()));

    for (int i = 0; i < elem.u_dimension(); ++i) {
      for (int j = 0; j < elem.shape_size(); ++j) {
        Point z = (1 - 1e-10) * elem.NodalPoint(j) + 1e-10 * z_c;
        u(r)[i * elem.shape_size() + j] = problem.ut(Time(), z, *c, i);
      }
    }
  }
}

void DGAcousticAssemble::SetInitialValue(Vector &u) {
  InterpolatedExactSolution(u);
}

void DGAcousticAssemble::FinishTimeStep(const Vector &u) {
  PlotIteration(u);
  PrintIteration(u);
}

void DGAcousticAssemble::PrintIteration(const Vector &u) const {
  if (verbose == 0 || Step() % verbose != 0) return;
  MinMaxP minMaxP = MinMaxPressure(u);
  if (problem.HasExactSolution()) {
    Vector u_int(u);
    InterpolatedExactSolution(u_int);
    AcousticErrors err_h = ComputeErrors(u);
    AcousticErrors err_int = ComputeErrors(u_int);

    mout.PrintIteration(
      verbose,
      PrintIterEntry("n", Step(), 5, 1),
      PrintIterEntry("t", Time(), 13, 1),
      PrintIterEntry("Energy", Energy(u), 13, 1),
      PrintIterEntry("IErr_L1", err_int.l1, 13, 1),
      PrintIterEntry("IErr_L2", err_int.l2, 13, 1),
      PrintIterEntry("IErr_infty", err_int.inf, 13, 1),
      PrintIterEntry("Err_L1", err_h.l1, 13, 1),
      PrintIterEntry("Err_L2", err_h.l2, 13, 1),
      PrintIterEntry("Err_infty", err_h.inf, 13, 1),
      PrintIterEntry("MaxP", minMaxP.second, 13, 1),
      PrintIterEntry("MinP", minMaxP.first, 13, 1)

    );
  } else {

    AcousticNorms norms = ComputeNorms(u);

    mout.PrintIteration(
      verbose,
      PrintIterEntry("n", Step(), 5, 1),
      PrintIterEntry("t", Time(), 13, 1),
      PrintIterEntry("Energy", Energy(u), 13, 1),
      PrintIterEntry("MaxP", minMaxP.second, 13, 1),
      PrintIterEntry("MinP", minMaxP.first, 13, 1),
      PrintIterEntry("L1", norms.l1, 13, 1),
      PrintIterEntry("L2", norms.l2, 13, 1)
    );
  }
}

void DGAcousticAssemble::PlotIteration(const Vector &u) const {
  if (plotting == 0 || Step() % plotting != 0) return;
  if (Step() == 0) PlotParams(u);
  PlotSolution(u);
}

void DGAcousticAssemble::PlotParams(const Vector &u) const {
  if (plotting == 0 || Step() % plotting != 0) return;
  auto matDisc = std::make_shared<const DGDiscretization>(u.GetDisc().GetMeshes(), 0, 3);
  Vector tmp(0.0, matDisc, u.Level());

  Vector RhoKappaVelocity(0.0, matDisc);
  std::vector<double> rho = {infty, -infty};
  std::vector<double> kappa = {infty, -infty};
  std::vector<double> velocity = {infty, -infty};

  for (cell c = RhoKappaVelocity.cells(); c != RhoKappaVelocity.cells_end(); ++c) {
    RhoKappaVelocity(c(), 0) = problem.Rho(*c, c());
    RhoKappaVelocity(c(), 1) = problem.Kappa(*c, c());
    RhoKappaVelocity(c(), 2) = sqrt(problem.Kappa(*c, c()) / problem.Rho(*c, c()));
    rho[0] = min(rho[0], RhoKappaVelocity(c(), 0));
    rho[1] = max(rho[1], RhoKappaVelocity(c(), 0));
    kappa[0] = min(kappa[0], RhoKappaVelocity(c(), 1));
    kappa[1] = max(kappa[1], RhoKappaVelocity(c(), 1));
    velocity[0] = min(velocity[0], RhoKappaVelocity(c(), 2));
    velocity[1] = max(velocity[1], RhoKappaVelocity(c(), 2));
  }

  rho[0] = PPM->Min(rho[0], u.GetMesh().CommSplit());
  rho[1] = PPM->Max(rho[1], u.GetMesh().CommSplit());
  kappa[0] = PPM->Min(kappa[0], u.GetMesh().CommSplit());
  kappa[1] = PPM->Max(kappa[1], u.GetMesh().CommSplit());
  velocity[0] = PPM->Min(velocity[0], u.GetMesh().CommSplit());
  velocity[1] = PPM->Max(velocity[1], u.GetMesh().CommSplit());

  vout(2) << "Rho€[" << rho[0] << ", " << rho[1] << "]" << " ----- "
          << "Kappa€[" << kappa[0] << ", " << kappa[1] << "]" << " ----- "
          << "Velocity€[" << velocity[0] << ", " << velocity[1] << "]" << endl;
  if ((kappa[0] <= 0) || (rho[0] <= 0)) Exit("not admissible material")

  std::string filename = "Params." + problem.Name();
  mpp::plot(to_string(tmp.Level().space)).AddData("Rho", RhoKappaVelocity, 0);
  mpp::plot(to_string(tmp.Level().space)).AddData("Kappa", RhoKappaVelocity, 1);
  mpp::plot(to_string(tmp.Level().space)).AddData("Velocity", RhoKappaVelocity, 2);
  mpp::plot(to_string(tmp.Level().space)).PlotFile(filename);
  Plotting::Instance().Clear();
}

void DGAcousticAssemble::PlotSolution(const Vector &u) const {
  if (plotting == 0 || Step() % plotting != 0) return;
  auto cellDisc = std::make_shared<const DGDiscretization>(GetDisc().GetMeshes(), 0, 1);
  auto velocityDisc = std::make_shared<const DGDiscretization>(GetDisc().GetMeshes(), 0, problem.SpaceDim());


  Vector VComponents(0.0, velocityDisc, u.Level());
  std::vector<Vector> PComponents;
  for (size_t i = 0; i < problem.NumDamping() + 1; ++i)
    PComponents.emplace_back(0.0, cellDisc, u.Level());
  Vector PSumComponent(0.0, cellDisc, u.Level());

  for (cell c = PSumComponent.cells(); c != PSumComponent.cells_end(); ++c) {

    const DGAcousticElement &elem = elementPool.Get(u, *c,DGAcousticElementGenerator(problem.NumDamping()));

    double Area = elem.Area();
    VectorField V = zero;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      for (size_t i = 0; i < PComponents.size(); ++i) {
        PComponents[i](c(), 0) += w * elem.Pressure_i(q, u, i) * (1.0 / Area);
      }
      V += w * elem.Velocity(q, u) * (1.0 / Area);
      PSumComponent(c(), 0) += w * elem.Pressure(q, u) * (1.0 / Area);
    }
    for (int i = 0; i < problem.SpaceDim(); ++i)
      VComponents(c(), i) = V[i];
  }
  std::string spaceLevel = to_string(u.Level().space);
  VtuPlot &plot = mpp::plot(spaceLevel);

  plot.AddData("v_x", VComponents, 0);
  plot.AddData("v_y", VComponents, 1);
  plot.AddData("p_sum", PSumComponent, 0);
  for (int i = 0; i < PComponents.size(); i++) {
    plot.AddData("p"+std::to_string(i), PComponents[i], 0);
  }
  plot.AddData("p_sum", PSumComponent, 0);
  for (int i = 0; i < PComponents.size(); i++) {
    plot.AddData("p"+std::to_string(i), PComponents[i], 0);
  }
  std::string filename = problem.Name() + "." + to_string(Step());
  plot.PlotFile(filename);

  Plotting::Instance().Clear();
}

void DGAcousticAssemble::Interpolate(const Vector &u, Vector &U) const {
  THROW("Think about numL for Elements")
  U = 0;
  int NN = 1 + u.dim();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {

    const DGAcousticElement &elem = elementPool.Get(u, *c);

    row r = u.find_row(c());
    for (int i = 0; i < c.Children(); ++i) {
      cell fc = U.find_cell(c.Child(i));

      // TODO: ElementPool needs to store second element if necessary
      // const DGAcousticElement &elem_fc = elementPool.Get(u, *c);
      DGAcousticElement elem_fc = DGAcousticElement(u, *c);

      row fr = U.find_row(fc());
      for (int iN = 0; iN < NN; ++iN) {
        for (int i = 0; i < elem.shape_size(); ++i)
          U(fr, iN * elem.shape_size() + i) =
              elem.SVValue(c.GlobalToLocal(elem_fc.NodalPoint(i)), u, iN);
      }
    }
  }
  U.Accumulate();
}
