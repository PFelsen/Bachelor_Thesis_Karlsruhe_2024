#ifndef SUPERLU_HPP
#define SUPERLU_HPP

#include "BasicSparseSolver.hpp"
#include "PreconditionerBase.hpp"

#include <set>

class SuperLU : public Preconditioner {
  SparseMatrix *S;
  SparseSolver *Sol;
public:
  SuperLU() : S(0), Sol(0) {}

  void Construct(const Matrix &_A) override;

  void Destruct() override;

  virtual ~SuperLU() { Destruct(); }

  void multiply(Vector &u, const Vector &b) const override;

  string Name() const { return "SuperLU"; }

  friend std::ostream &operator<<(std::ostream &s, const SuperLU &SLU) { return s << "SuperLU"; }
};

class SuperLU_local : public Preconditioner {
  SparseMatrix *S;
  SparseSolver *Sol;
  typedef std::unordered_map<Point, int>::iterator Iterator;
  typedef std::unordered_map<Point, int>::const_iterator ConstIterator;
  std::unordered_map<Point, std::set<short>> *IProc;
  Vector *rhs;
  bool shift_special;
  int verbose = 0;
public:
  SuperLU_local();
private:
  void CommunicateMatrix(Matrix &A);

  void CreateIProc(const Vector &u);

  void CollectResidual(const Vector &b) const;

  void DistributeSolution(Vector &u) const;
public:
  void Construct(const Matrix &_A) override;

  void Destruct() override;

  virtual ~SuperLU_local() { Destruct(); }

  void multiply(Vector &u, const Vector &b) const override;

  string Name() const { return "SuperLU_local"; }

  friend std::ostream &operator<<(std::ostream &s, const SuperLU_local &SLU) {
    return s << "SuperLU_local";
  }
};

#endif // SUPERLU_HPP
