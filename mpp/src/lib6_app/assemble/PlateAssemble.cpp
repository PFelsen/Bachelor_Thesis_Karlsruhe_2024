#include "PlateAssemble.hpp"
#include "ArgyrisElement.hpp"

void PlateAssemble::Initialize(Vector &u) const {
  u = 0;
  u.ClearDirichletFlags();
  ArgyrisElement::H20BC(problem.corners, u);
}

double PlateAssemble::Energy(const Vector &u) const {
  double energy = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    ArgyrisElement Elem(u, *c);
    for (int q = 0; q < Elem.nQ(); ++q) {
      double w = Elem.QWeight(q);
      Scalar U = Elem.Value(q, u);
      VectorField DU = Elem.Derivative(q, u);
      double DeltaU = Elem.Laplace(q,u);
      energy += w * DeltaU * DeltaU;
    }
  }
  return sqrt(PPM->Sum(energy));
}

void PlateAssemble::Residual(const cell &c, const Vector &u, Vector &defect) const {
    ArgyrisElement E(defect, *c);
    RowValues defect_c(defect, E);
    for (int q = 0; q < E.nQ(); ++q) {
      double f_q = problem.Load(E.QPoint(q));
      double Lu_q = E.Laplace(q, u);
      for (int i = 0; i < E.size(); ++i)
        for (int k = 0; k < E.get_maxk(i); ++k)
          defect_c(i, k) += E.QWeight(q) * (Lu_q * E.Laplace(q, i, k) - f_q * E.Value(q, i, k));
    }
}

void PlateAssemble::Jacobi(const cell &c, const Vector &u, Matrix &jacobi) const {
  ArgyrisElement E(jacobi, *c);
  RowEntries jacobi_c(jacobi, E);
  for (int q = 0; q < E.nQ(); ++q)
    for (int i = 0; i < E.size(); ++i)
      for (int j = 0; j < E.size(); ++j)
	for (int k = 0; k < E.get_maxk(i); ++k)
	  for (int l = 0; l < E.get_maxk(j); ++l)
	    jacobi_c(i, j, k, l) += E.QWeight(q) * E.Laplace(q, i, k) * E.Laplace(q, j, l);
}

void PlateAssemble::Matrices(Matrix &A, Matrix& B) const {
  A = 0;
  B = 0;
  for (cell c = A.cells(); c != A.cells_end(); ++c) {
    ArgyrisElement E(A, *c);
    RowEntries A_c(A, *c);
    RowEntries B_c(B, *c);
    
    for (int q = 0; q < E.nQ(); ++q)
      for (int i = 0; i < E.size(); ++i)
	for (int j = 0; j < E.size(); ++j)
	  for (int k = 0; k < E.get_maxk(i); ++k)
	    for (int l = 0; l < E.get_maxk(j); ++l) {
	      A_c(i, j, k, l) += E.QWeight(q) * E.Laplace(q, i, k) * E.Laplace(q, j, l);
	      B_c(i, j, k, l) += E.QWeight(q) * E.Value(q, i, k) * E.Value(q, j, l);
	    }
  }
  A.ClearDirichletValues();
  B.ClearDirichletValues();
}

PlateAssemble *CreatePlateAssemble(const PlateProblem &problem, const PDESolverConfig &conf) {

  if (conf.modelName.find("Plate") != std::string::npos) { return new PlateAssemble(problem); }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<PlateAssemble> CreatePlateAssembleUnique(const PlateProblem &problem,
                                                         const PDESolverConfig &conf) {
  return std::unique_ptr<PlateAssemble>(CreatePlateAssemble(problem, conf));
}

std::shared_ptr<PlateAssemble> CreatePlateAssembleShared(const PlateProblem &problem,
                                                         const PDESolverConfig &conf) {
  return std::shared_ptr<PlateAssemble>(CreatePlateAssemble(problem, conf));
}
