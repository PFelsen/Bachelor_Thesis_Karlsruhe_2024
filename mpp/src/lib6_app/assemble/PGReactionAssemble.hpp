#ifndef _REACTION_H_
#define _REACTION_H_

#include "LagrangeDiscretization.hpp"
#include "ctools.hpp"
#include "ReactionProblems.hpp"
#include "IReactionAssemble.hpp"
#include "ScalarElement.hpp"
#include "ElementPool.hpp"


class PGReactionAssemble : public IReactionAssemble {
private:
  mutable ElementPool<ScalarElement> elementPool;
  mutable FaceElementPool<ScalarFaceElement> faceElementPool;
  
  std::shared_ptr<const LagrangeDiscretization> disc;

  double delta = 0.0;

public:
  PGReactionAssemble(const IReactionProblem &problem, int degree) :
      disc(std::make_shared<LagrangeDiscretization>(problem.GetMeshes(), degree)),
      IReactionAssemble(problem) {
    Config::Get("delta", delta);
  }

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return disc;
  }

  const char *Name() const override { return "PGReactionAssemble"; }

  void PrintIteration(const Vector &u) const override;;

  void FinishTimeStep(const Vector &u) override;

  void Initialize(Vector &u) const override;

  double Residual(const Vector &u, Vector &r) const override;

  void Jacobi(const Vector &u, Matrix &A) const override;

  double Mass(const Vector &u) const override;

  RatePair InflowOutflow(const Vector &u) const override;

  void SetInitialValue(Vector &u) override;

  void PlotIteration(const Vector &u) const override;

  void PrintInfo() const override;

private:
  double Delta(const cell &c) const;

  double diam(const cell &c) const;
};

#endif
