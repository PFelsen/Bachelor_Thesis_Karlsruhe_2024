#ifndef IAEIGENVALUEASSEMBLE_H
#define IAEIGENVALUEASSEMBLE_H

#include "Eigenpair.hpp"

#include "Newton.hpp"
#include "Parallel.hpp"


typedef Vector Goerischfct;

typedef Vectors Goerischfcts;

/**
 * This class represents the assemble for the eigenvalue methods dealing with the following
 * eigenvalue problem
 *
 *      \f[ M_t(u,v) = \lambda N_t(u,v) \qquad (\forall v), \f]
 *
 * where t is a parameter which e.g. can be used in a homotopy. The assemble is devided into three
 * parts:
 *
 * 1. Assembling the stiffness matrix A and mass matrix B for the approximative computation of
 * eigenvalues and eigenfunctions using the FEM ansatz functions and an iterative eigensolver.
 *
 * 2. Assembling the matrices for the RayleighRitz method to compute upper eigenvalue bounds.
 *
 * 3. Assembling the matrices needed for the Newton calculation in the LehmannGoerisch method and
 * assembling the additional matrix needed for the computation of lower eigenvalue bounds via the
 * LehmannGoerisch method.
 */

template<bool VERIFIED>
using IAEigenvalueType = typename std::conditional<VERIFIED, IAInterval, double>::type;

template<bool VERIFIED>
class IAEigenvalueAssemble : public ILinearAssemble {
private:
  const Eigenpairs *goerischEigenpairs;
  double goerisch_t;
  double shift;

  template<bool verified>
  friend class IALehmannGoerischMethod;

  void SetGoerischData(const Eigenpairs &_goerischEigenpairs, double _goerisch_t) {
    goerischEigenpairs = &_goerischEigenpairs;
    goerisch_t = _goerisch_t;
  }

  void ClearGoerischData() { goerischEigenpairs = nullptr; }
public:
  explicit IAEigenvalueAssemble(double shift = 0.0) :
      goerischEigenpairs(nullptr), goerisch_t(0.0), shift(shift) {
    verbose = 1;
    Config::Get("EigenvalueAssembleVerbose", verbose);
  }

  virtual ~IAEigenvalueAssemble() = default;

  virtual void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", this->Name(), 0),
                   PrintInfoEntry("Shift", shift, 0));
  }

  double Shift() const { return shift; }

  /// Sets the boundary conditions for the given eigenfunctions
  virtual void BoundaryConditionsEigenSolver(Eigenfcts &U) const {
    for (int i = 0; i < U.size(); ++i)
      BoundaryConditionsEigenSolver(U[i]);
  }

  /// Sets the boundary condition for the given eigenfunction
  virtual void BoundaryConditionsEigenSolver(Eigenfct &u) const {
    u.ClearDirichletFlags();
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      BoundaryConditionsEigenSolver(c, u);
    u.DirichletConsistent();
  }

  /// Sets the boundary condition for the given eigenfunction on a specified cell
  virtual void BoundaryConditionsEigenSolver(const cell &c, Eigenfct &u) const {
    THROW("BoundaryConditionsEigenSolver not implemented!")
  }

  /// Sets the Matrices for the approximative computation of the eigenfunctions
  virtual void MatricesEigenSolver(Matrix &A, Matrix &B, const double &t = 0.0) {
    A = 0;
    B = 0;
    for (cell c = A.cells(); c != A.cells_end(); ++c)
      MatrixentriesEigenSolver(c, A, B, shift, t);
    A.ClearDirichletValues();
    B.ClearDirichletValues();
  }

  /// Sets the Matrixentries for the approximative computation of the eigenfunctions on
  /// a specified cell
  virtual void MatrixentriesEigenSolver(const cell &c, Matrix &A, Matrix &B, const double &shift,
                                        const double &t) {
    THROW("MatrixentriesEigenSolver not implemented!")
  }

  virtual std::unique_ptr<IDiscretization> CreateGoerischDisc(const IDiscretization &EVdisc) const {
    THROW("CreateGoerischDisc not implemented")
  }

  virtual std::unique_ptr<IAIDiscretization>
  CreateIAGoerischDisc(const IAIDiscretization &iaEVdisc) const {
    THROW("CreateIAGoerischDisc not implemented")
  }

  void Initialize(Vectors &u) const override {
    for (int n = 0; n < u.size(); ++n) {
      u[n].ClearDirichletFlags();
      for (cell c = u[n].cells(); c != u[n].cells_end(); ++c)
        BoundaryConditionsGoerisch(c, u[n]);
      u[n].DirichletConsistent();
    }
  }

  void AssembleSystem(const cell &c, Matrix &systemMatrix, Vectors &rhs) const override {
    AssembleGoerischSystem(c, systemMatrix, rhs, shift, *goerischEigenpairs, goerisch_t);
  }

  virtual void BoundaryConditionsGoerisch(const cell &c, Vector &u) const {
    THROW("BoundaryConditionsGoerisch(const cell &, Vector &) not implemented!")
  }

  virtual void AssembleGoerischSystem(const cell &c, Matrix &systemMatrix, Goerischfcts &w,
                                      double shift, const Eigenpairs &ep, double t) const {
    Warning("AssembleGoerischSystem(const cell &, Matrix &, Goerischfcts &, double, const "
            "Eigenpairs &, double)  not implemented!")
  }

  /// Sets the Matrices for the RayleighRitz computation
  virtual void MatricesRayleighRitz(const Eigenfcts &U, SymRMatrixT<IAEigenvalueType<VERIFIED>> &a,
                                    SymRMatrixT<IAEigenvalueType<VERIFIED>> &b,
                                    const double &t = 0.0) {
    a = IAEigenvalueType<VERIFIED>{};
    b = IAEigenvalueType<VERIFIED>{};
    for (cell c = U.cells(); c != U.cells_end(); ++c)
      for (int i = 0; i < U.size(); ++i)
        for (int j = i; j < U.size(); ++j)
          MatrixentriesRayleighRitz(c, U[i], U[j], a(i, j), b(i, j), shift, t);
    a.Accumulate();
    b.Accumulate();
  }

  /// Sets the 1d-Matrices for the RayleighRitz computation (e.g. needed for the
  /// computation of the Rayleigh quotient
  virtual void MatricesRayleighRitz(const Eigenfct &u, IAEigenvalueType<VERIFIED> &a,
                                    IAEigenvalueType<VERIFIED> &b, const double &t = 0.0) {
    a = IAEigenvalueType<VERIFIED>{};
    b = IAEigenvalueType<VERIFIED>{};
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      MatrixentriesRayleighRitz(c, u, u, a, b, shift, t);
    a = PPM->template SumOnCommSplit(a, 0);
    b = PPM->template SumOnCommSplit(b, 0);
  }

  /// Sets the Matrixentries for the RayleighRitz computation on a specified cell
  virtual void MatrixentriesRayleighRitz(const cell &c, const Eigenfct &U_i, const Eigenfct &U_j,
                                         IAEigenvalueType<VERIFIED> &a,
                                         IAEigenvalueType<VERIFIED> &b, const double &shift,
                                         const double &t) {
    THROW("MatrixentriesRayleighRitz not implemented!")
  }

  /// Sets the Matrix for the LehmannGoerisch computation
  virtual void MatrixGoerisch(const Eigenpairs &EP, const Goerischfcts &W,
                              SymRMatrixT<IAEigenvalueType<VERIFIED>> &g, const double &t = 0.0) {
    g = 0.0;
    for (cell c = EP.cells(); c != EP.cells_end(); ++c)
      for (int i = 0; i < EP.size(); ++i)
        for (int j = i; j < EP.size(); ++j)
          MatrixentriesGoerisch(c, EP.getEigenpair(i), EP.getEigenpair(j), W[i], W[j], g(i, j),
                                shift, t);
    g.Accumulate();
  }

  /// Sets the Matrixentires for the LehmannGoerisch computation on a specified cell
  virtual void MatrixentriesGoerisch(const cell &c, const Eigenpair &EP_i, const Eigenpair &EP_j,
                                     const Goerischfct &W_i, const Goerischfct &W_j,
                                     IAEigenvalueType<VERIFIED> &g, const double &shift,
                                     const double &t) {
    THROW("MatrixentriesGoerisch not implemented!")
  }
};

#endif // IAEIGENVALUEASSEMBLE_H
