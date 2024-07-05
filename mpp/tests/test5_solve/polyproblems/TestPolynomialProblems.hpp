#ifndef TESTPOLYNOMIALPROBLEMS_HPP
#define TESTPOLYNOMIALPROBLEMS_HPP

#include <string>
#include <vector>
#include "Algebra.hpp"
#include "ArgyrisElement.hpp"
#include "Assemble.hpp"
#include "ScalarElement.hpp"

class TestProblem {
protected:
  int degree;
public:
  TestProblem(int degree) : degree(degree){};

  virtual ~TestProblem(){};

  virtual Scalar Solution(const Point &x) const = 0;

  virtual Scalar Load(const Point &x) const = 0;
};

class TestAssemble : public IAssemble {
protected:
  const TestProblem &problem;
public:
  TestAssemble(const TestProblem &problem) : problem(problem) {}

  virtual double L2Error(const Vector &u) const {
    double err = 0.0;
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      const Cell &C = *c;
      ScalarElement elem(u, C);
      for (int q = 0; q < elem.nQ(); ++q) {
        double w = elem.QWeight(q);
        Scalar U = elem.Value(q, u);
        Scalar UEx = problem.Solution(elem.QPoint(q));
        err += w * (U - UEx) * (U - UEx);
      }
    }
    return sqrt(PPM->SumOnCommSplit(err, 0));
  }
};

/*
 * Laplace Problems
 */
template<int i>
class TestProblemLaplaceMonomialT : public TestProblem {
public:
  TestProblemLaplaceMonomialT(int degree) : TestProblem(degree) {}

  Scalar Solution(const Point &x) const override { return pow(x[i], degree); }

  Scalar Load(const Point &x) const override {
    return -degree * (degree - 1) * pow(x[i], degree - 2);
  }
};

class LaplaceAssemble : public TestAssemble {
public:
  LaplaceAssemble(const TestProblem &problem) : TestAssemble(problem) {}

  const char *Name() const override { return "Laplace Test Assemble"; }

  void Initialize(Vector &u) const override {
    u.ClearDirichletFlags();
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      const Cell &C = *c;
      RowBndValues u_c(u, C);
      if (!u_c.onBnd()) continue;
      ScalarElement elem(u, C);
      for (int face = 0; face < c.Faces(); ++face) {
        if (u_c.bc(face) == 1) {
          for (int j = 0; j < u.NumberOfNodalPointsOnFace(C, face); ++j) {
            int k = u.IdOfNodalPointOnFace(C, face, j);
            u_c(k) = problem.Solution(elem.NodalPoint(k));
            u_c.D(k) = true;
          }
        }
      }
    }
    u.DirichletConsistent();
  }

  void Residual(const cell &c, const Vector &u, Vector &r) const override {
    const Cell &C = *c;
    ScalarElement elem(u, C);
    RowBndValues r_c(r, C);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      const Point &z = elem.QPoint(q);
      Scalar f = problem.Load(z);
      VectorField DU = elem.Derivative(q, u);
      for (int i = 0; i < elem.size(); ++i) {
        Scalar U_i = elem.Value(q, i);
        VectorField DU_i = elem.Derivative(q, i);
        r_c(i) += w * (DU * DU_i - f * U_i);
      }
    }
  }

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override {
    const Cell &C = *c;
    ScalarElement E(u, C);
    RowEntries A_c(A, E);
    for (int q = 0; q < E.nQ(); ++q) {
      double w = E.QWeight(q);
      for (unsigned int i = 0; i < E.size(); ++i) {
        VectorField Dv = E.Derivative(q, i);
        for (unsigned int j = 0; j < E.size(); ++j) {
          VectorField Dw = E.Derivative(q, j);
          A_c(i, j) += w * Dv * Dw;
        }
      }
    }
  }

  double Energy(const Vector &u) const override {
    double energy = 0.0;
    for (cell c = u.cells(); c != u.cells_end(); ++c) {
      const Cell &C = *c;
      ScalarElement elem(u, C);
      for (int q = 0; q < elem.nQ(); ++q) {
        double w = elem.QWeight(q);
        VectorField DU = elem.Derivative(q, u);
        energy += w * DU * DU;
      }
    }
    return sqrt(PPM->SumOnCommSplit(energy, 0));
  }
};

/*
 * BiLaplace Problems
 */
class TestProblemBiLaplace : public TestProblem {
public:
  TestProblemBiLaplace() : TestProblem(0) {}

  Scalar Solution(const Point &x) const override {
    return pow(x[0] * (1 - x[0]) * x[1] * (1 - x[1]), 2);
  }

  Scalar Load(const Point &x) const override {
    return 24 * (pow(x[0] * (1 - x[0]), 2) + pow(x[1] * (1 - x[1]), 2))
           + 8 * (1 - 6 * x[0] * (1 - x[0])) * (1 - 6 * x[1] * (1 - x[1]));
  }
};

class TestProblemBiLaplaceTurned : public TestProblem {
public:
  TestProblemBiLaplaceTurned() : TestProblem(0) {}

  Scalar Solution(const Point &Q) const override {
    double t1 = sqrt(2) * (Q[1] - Q[0]);
    double t2 = sqrt(2) * (Q[1] + Q[0]);
    return pow((t1 - 1.0) * (t1 + 1.0) * (t2 - 1.0) * (t2 + 1.0), 2) / 256;
  }

  Scalar Load(const Point &x) const override {
    return 84 * (pow(x[0], 4) + pow(x[1], 4)) - 36 * (pow(x[0], 2) + pow(x[1], 2))
           - 72 * pow(x[0] * x[1], 2) + 5;
  }
};

class BiLaplaceAssemble : public TestAssemble {
protected:
  const std::list<Point> &corners;
public:
  BiLaplaceAssemble(const TestProblem &problem, const std::list<Point> &corners) :
      TestAssemble(problem), corners(corners) {}

  const char *Name() const override { return "BiLaplace Test Assemble"; }

  void Initialize(Vector &u) const override { ArgyrisElement::H20BC(corners, u); }

  void Residual(const cell &c, const Vector &u, Vector &defect) const override {
    ArgyrisElement E(defect, *c);
    RowValues defect_c(defect, E);
    for (int q = 0; q < E.nQ(); ++q) {
      double f_q = problem.Load(E.QPoint(q));
      double Lu_q = E.Laplace(q, u);
      for (int i = 0; i < E.size(); ++i)
        for (int k = 0; k < E.get_maxk(i); ++k)
          defect_c(i, k) += E.QWeight(q) * (Lu_q * E.Laplace(q, i, k) - f_q * E.Value(q, i, k));
    }
  };

  void Jacobi(const cell &c, const Vector &u, Matrix &jacobi) const override {
    ArgyrisElement E(jacobi, *c);
    RowEntries jacobi_c(jacobi, E);
    for (int q = 0; q < E.nQ(); ++q)
      for (int i = 0; i < E.size(); ++i)
        for (int j = 0; j < E.size(); ++j)
          for (int k = 0; k < E.get_maxk(i); ++k)
            for (int l = 0; l < E.get_maxk(j); ++l)
              jacobi_c(i, j, k, l) += E.QWeight(q) * E.Laplace(q, i, k) * E.Laplace(q, j, l);
  };
};


#endif // TESTPOLYNOMIALPROBLEMS_HPP
