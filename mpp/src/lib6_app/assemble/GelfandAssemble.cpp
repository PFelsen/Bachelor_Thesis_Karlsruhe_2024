#include "GelfandAssemble.hpp"
#include "ScalarElement.hpp"

void GelfandAssemble::Initialize(Vector &u) const {
  for (row r = u.rows(); r!= u.rows_end(); ++r)
    u(r,0) = problem.InitialValue(r());
  u.ClearDirichletFlags();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    RowBndValues u_c(u, *c);
    if (!u_c.onBnd()) continue;
    ScalarElement elem(u, *c); 
    for (int face = 0; face < c.Faces(); ++face) {
      if (u_c.bc(face) == 1) {
        for (int j = 0; j < u.NumberOfNodalPointsOnFace(*c, face); ++j) {
          int k = u.IdOfNodalPointOnFace(*c, face, j);
          u_c(k) = problem.DirichletValue(elem.NodalPoint(k));
          u_c.D(k) = true;
        }
      }
    }
  }
  u.DirichletConsistent();
}

double GelfandAssemble::Energy(const Vector &u) const {
  double energy = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    ScalarElement elem(u, *c);    
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      VectorField DU = elem.Derivative(q, u);
      energy += w * (DU * DU - problem.F(U));
    }
  }
  return sqrt(PPM->Sum(energy));
}

void GelfandAssemble::Residual(const cell &c, const Vector &u, Vector &r) const {
  ScalarElement elem(u, *c);
  RowBndValues r_c(r, *c);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    const Point &z = elem.QPoint(q);
    Scalar U = elem.Value(q, u);
    VectorField DU = elem.Derivative(q, u);
    for (int i = 0; i < elem.size(); ++i) {
      Scalar U_i = elem.Value(q, i);
      VectorField DU_i = elem.Derivative(q, i);
      r_c(i) += w * (DU * DU_i - problem.F(U) * U_i);
    }
  }
}

void GelfandAssemble::Jacobi(const cell &c, const Vector &u, Matrix &A) const {
  ScalarElement elem(u, *c);
  RowEntries A_c(A, elem);
  for (int q = 0; q < elem.nQ(); ++q) {
    double w = elem.QWeight(q);
    Scalar U = elem.Value(q, u);
    for (unsigned int i = 0; i < elem.size(); ++i) {
      Scalar U_i = elem.Value(q, i);
      VectorField DU_i = elem.Derivative(q, i);
      for (unsigned int j = 0; j < elem.size(); ++j) {
	Scalar U_j = elem.Value(q, j);
        VectorField DU_j = elem.Derivative(q, j);
        A_c(j, i) += w * (DU_i * DU_j - problem.DF(U) * U_i * U_j);
      }
    }
  }
}

GelfandAssemble *CreateGelfandAssemble(const GelfandProblem &problem, const PDESolverConfig &conf) {

  if (conf.modelName.find("Gelfand") != std::string::npos) {
    return new GelfandAssemble(problem, conf.degree);
  }

  Exit(conf.modelName + " not found")
}

std::unique_ptr<GelfandAssemble> CreateGelfandAssembleUnique(const GelfandProblem &problem,
                                                             const PDESolverConfig &conf) {
  return std::unique_ptr<GelfandAssemble>(CreateGelfandAssemble(problem, conf));
}

std::shared_ptr<GelfandAssemble> CreateGelfandAssembleShared(const GelfandProblem &problem,
                                                             const PDESolverConfig &conf) {
  return std::shared_ptr<GelfandAssemble>(CreateGelfandAssemble(problem, conf));
}
