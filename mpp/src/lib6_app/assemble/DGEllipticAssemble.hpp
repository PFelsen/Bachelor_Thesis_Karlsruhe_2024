#ifndef _DGLAPLACE_H_
#define _DGLAPLACE_H_

#include "IEllipticAssemble.hpp"
#include "EllipticAssemble.hpp"
#include "DGElement.hpp"
#include "DGDiscretization.hpp"
#include "ElementPool.hpp"


class DGEllipticAssemble : public EllipticAssemble<DGElement, DGFaceElement> {
private:
  std::shared_ptr<const DGDiscretization> disc;

  int degree = disc->Degree();

  double penalty = 20.0;

  int sign = 1;

public:
  DGEllipticAssemble(const IEllipticProblem &problem, int degree) :
      disc(std::make_shared<DGDiscretization>(problem.GetMeshes(), degree)),
      EllipticAssemble(problem) {
    Config::Get("penalty", penalty);
    Config::Get("sign", sign);
  }

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc; }

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose,
                   PrintInfoEntry("Assemble Name", Name()),
                   PrintInfoEntry("penalty", penalty),
                   PrintInfoEntry("sign", sign),
                   PrintInfoEntry("Problem", problem.Name()));
  }

  const char *Name() const override;

  void Initialize(Vector &u) const override;

  void Residual(const cell &c, const Vector &u, Vector &r) const override;

  void Jacobi(const cell &c, const Vector &u, Matrix &A) const override;

  FluxPair InflowOutflow(const Vector &u) const override;

  FluxPair PrescribedInflowOutflow(const Vector &u) const override;

  FluxPair OutflowLeftRight(const Vector &u) const override;

  double GoalFunctional(const Vector &u) const override;

  void SetExactSolution(Vector &uEx) const override;

  void SetFlux(const Vector &u, Vector &flux) override;

  void SetPressure(const Vector &u, Vector &p) const override;
};

#endif
