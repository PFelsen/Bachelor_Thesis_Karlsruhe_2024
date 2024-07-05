#ifndef PROBLEMBASE_HPP
#define PROBLEMBASE_HPP

#include "IProblem.hpp"
#include "Vector.hpp"


class ProblemBase : virtual public IProblem {
protected:
  int verbose = 0;

public:
  ProblemBase() {
    Config::Get("ProblemVerbose", verbose);
  }

  virtual ~ProblemBase() = default;

  virtual double Rho(const Point &p) const { return 1; }

  virtual double Kappa(const Point &p) const { return 1; }

  virtual double Kappa(const Point &p, int var) const {
    Exit("Virtual function not implemented")
  }

  virtual double Tau(const Point &p, int var) const {
    Exit("Virtual function not implemented")
  }

  virtual int Dim() const = 0;

  virtual int nL() const { return 0; };

  bool HasExactSolution() const override { return false; }

  virtual bool HasRHS() const { return false; }

  virtual int BndID(const Point &z) const { return 2; }

  //virtual double ut(double t, const Point &x, COMPONENT comp) const { return 0; }
  /*virtual double ut(const Point &p) const { return 0; }

  virtual double ut(const Point &p, COMPONENT comp) const { return ut(p); }

  virtual double ut(const Point &p, const Cell &c, COMPONENT comp) const { return ut(p,comp); }*/

  virtual std::vector<COMPONENT> GetComponents() const { Exit("no components"); }
};

#endif //PROBLEMBASE_HPP
