#ifndef GMRES_HPP
#define GMRES_HPP

#include "LinearSolver.hpp"

class GMRES : public LinearSolver {
private:
  int M = 300;
public:
  GMRES(Preconditioner *preCond, const std::string &prefix = "Linear") :
      LinearSolver(preCond, prefix) {
    Config::Get(prefix + "Restart", M);
  }

  GMRES(std::unique_ptr<Preconditioner> preCond, const std::string &prefix = "Linear") :
      LinearSolver(std::move(preCond), prefix) {
    Config::Get(prefix + "Restart", M);
  }

  void solve(const Operator &A, const Operator &B, Vector &u, Vector &r, double d0,
             double epsilon) const;

  void solve(const Operator &A, const Operator &B, Vectors &u, Vectors &r, double d0,
             double eps) const;

  std::string Name() const { return "GMRES"; }
};

#endif // GMRES_HPP
