#ifndef EIGENVALUEEXAMPLEASSEMBLE_HPP
#define EIGENVALUEEXAMPLEASSEMBLE_HPP

#include "IAEigenvalueAssemble.hpp"
#include "LagrangeDiscretization.hpp"
#include "VectorFieldElement.hpp"

template<bool VERIFIED>
class EigenvalueExampleAssemble : public IAEigenvalueAssemble<VERIFIED> {
  int lagrange_degree;

  int computeMaxPolynomialDegree() const {
    return 2 * lagrange_degree + 4;
  }
public:
  EigenvalueExampleAssemble(double shift, int lagrange_degree = 3)
      : IAEigenvalueAssemble<VERIFIED>(shift), lagrange_degree(lagrange_degree) {}

  const char *Name() const override {
    return "Eigenvalue methods example assemble";
  }

  void BoundaryConditionsEigenSolver(const cell &c, Eigenfct &u) const override {
    if (!u.OnBoundary(*c)) return;
    for (int face = 0; face < c.Faces(); ++face) {
      switch (u.GetMesh().BoundaryFacePart(c.Face(face))) {
        case 1: {
          rows R(u.GetMatrixGraph(), *c, face);
          for (const row &r: R) {
            for (int k = 0; k < r.NumberOfDofs(); ++k) {
              u(r, k) = 0;
              u.D(r, k) = true;
            }
          }
          break;
        }
        case -1:
          break;
        default: Exit("Boundary condition not implemented;\n");
      }
    }
  }

  void MatrixentriesEigenSolver(const cell &c, Matrix &A, Matrix &B,
                                const double &shift, const double &t) override {
    ScalarElement E(A, *c);
    RowEntries A_c(A, E);
    RowEntries B_c(B, E);
    for (int q = 0; q < E.nQ(); ++q) {
      double tmp = shift - 1.0 + t * (1.0 - 4.0 * E.QPoint(q)[0] * (1 - E.QPoint(q)[0]));
      for (int i = 0; i < E.size(); ++i) {
        for (int j = 0; j < E.size(); ++j) {
          A_c(i, j) += E.QWeight(q) * (E.Derivative(q, i) * E.Derivative(q, j) + tmp * E.Value(q, i) * E.Value(q, j));
          B_c(i, j) += E.QWeight(q) * E.Value(q, i) * E.Value(q, j);
        }
      }
    }
  }

  void MatrixentriesRayleighRitz(const cell &c,
                                 const Eigenfct &U_i, const Eigenfct &U_j,
                                 IAEigenvalueType<VERIFIED> &a, IAEigenvalueType<VERIFIED> &b,
                                 const double &shift, const double &t) override {
    using T = IAEigenvalueType<VERIFIED>;
    ScalarElementT<T> E(U_i, *c);
    for (int q = 0; q < E.nQ(); ++q) {
      T tmp = shift - 1.0 + t * (1.0 - 4.0 * E.QPoint(q)[0] * (1 - E.QPoint(q)[0]));
      a += E.QWeight(q) * (E.Derivative(q, U_i) * E.Derivative(q, U_j) + tmp * E.Value(q, U_i) * E.Value(q, U_j));
      b += E.QWeight(q) * E.Value(q, U_i) * E.Value(q, U_j);
    }
  }

  std::unique_ptr<IDiscretization>
  CreateEVDisc(const Meshes &meshes) const {
    return std::make_unique<LagrangeDiscretization>(meshes, lagrange_degree, computeMaxPolynomialDegree(), 1);
  }

  std::unique_ptr<IAIDiscretization>
  CreateIAEVDisc(const Meshes &meshes) const {
    return std::make_unique<IALagrangeDiscretization>(meshes, lagrange_degree, computeMaxPolynomialDegree(), 1);
  }

  std::unique_ptr<IDiscretization>
  CreateGoerischDisc(const IDiscretization &EVdisc) const override {
    return std::make_unique<LagrangeDiscretization>(
        dynamic_cast<const NonAdaptiveIDiscretization &>(EVdisc), lagrange_degree, 3);
  }

  std::unique_ptr<IAIDiscretization>
  CreateIAGoerischDisc(const IAIDiscretization &iaEVdisc) const override {
    return std::make_unique<IALagrangeDiscretization>(
        dynamic_cast<const IANonAdaptiveIDiscretization &>(iaEVdisc), lagrange_degree, 3);
  }

  void BoundaryConditionsGoerisch(const cell &c, Vector &w) const override {}

  void AssembleGoerischSystem(const cell &c, Matrix &systemMatrix, Goerischfcts &w,
                              double shift, const Eigenpairs &ep, double t) const override {
    VectorFieldElement E(w, *c);
    RowValuesVector rowValues(w, E);
    RowEntries rowEntries(systemMatrix, E);
    for (int q = 0; q < E.nQ(); ++q) {
      double tmp = shift - 2.0 + t * (1.0 - 4.0 * E.QPoint(q)[0] * (1.0 - E.QPoint(q)[0]));

      //assembling rhs
      for (int n = 0; n < w.size(); ++n) {
        double V = E.Value(q, ep(n));
        for (int i = 0; i < E.size(); ++i)
          for (int k = 0; k < 2; ++k) {
            rowValues[n](i, k) -= E.QWeight(q) * V * E.Divergence(q, i, k);
          }
        for (int i = 0; i < E.size(); ++i)
          rowValues[n](i, 2) += E.QWeight(q) * tmp * V * E.Value(q, i);
      }

      //assembling system matrix
      for (int i = 0; i < E.size(); ++i) {
        for (int k = 0; k < 2; ++k) {
          for (int j = 0; j < E.size(); ++j) {
            for (int l = 0; l < 2; ++l) {
              rowEntries(i, j, k, l) += E.QWeight(q) * (E.VectorComponentValue(q, j, l) * E.VectorComponentValue(q, i, k) + E.Divergence(q, j, l) * E.Divergence(q, i, k));
            }
          }
          for (int j = 0; j < E.size(); ++j) {
            rowEntries(i, j, k, 2) -= E.QWeight(q) * tmp * E.Value(q, j) * E.Divergence(q, i, k);
            rowEntries(j, i, 2, k) -= E.QWeight(q) * tmp * E.Value(q, j) * E.Divergence(q, i, k);
          }
        }
      }
      for (int i = 0; i < E.size(); ++i)
        for (int j = 0; j < E.size(); ++j)
          rowEntries(i, j, 2, 2) += E.QWeight(q) * tmp * (1.0 + tmp) * E.Value(q, j) * E.Value(q, i);
    }
  }

  void MatrixentriesGoerisch(const cell &c,
                             const Eigenpair &EP_i, const Eigenpair &EP_j,
                             const Goerischfct &W_i, const Goerischfct &W_j,
                             IAEigenvalueType<VERIFIED> &g,
                             const double &shift, const double &t) override{
    using T = IAEigenvalueType<VERIFIED>;
    VectorFieldElementT<T> E(W_i, *c);
    for (int q = 0; q < E.nQ(); ++q) {
      T tmp = shift - 2.0 + t * (1.0 - 4.0 * E.QPoint(q)[0] * (1.0 - E.QPoint(q)[0]));
      VectorFieldT<T> W1_i = E.VectorValue(q, W_i);
      VectorFieldT<T> W1_j = E.VectorValue(q, W_j);
      T W2_i = E.Value(q, W_i, 2);
      T W2_j = E.Value(q, W_j, 2);
      T W3_i = E.Value(q, EP_i.getEigenfct()) - tmp * W2_i + E.Divergence(q, W_i);
      T W3_j = E.Value(q, EP_j.getEigenfct()) - tmp * W2_j + E.Divergence(q, W_j);

      g += E.QWeight(q) * (W1_i * W1_j + tmp * W2_i * W2_j + W3_i * W3_j);
    }
  }
};

#endif //EIGENVALUEEXAMPLEASSEMBLE_HPP
