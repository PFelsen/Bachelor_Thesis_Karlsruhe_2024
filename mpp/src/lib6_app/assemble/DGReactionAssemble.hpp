#ifndef _DG_REACTION_H_
#define _DG_REACTION_H_

#include "IReactionAssemble.hpp"
#include "DGElement.hpp"
#include "ctools.hpp"
#include "DGDiscretization.hpp"
#include "ElementPool.hpp"


class DGReactionAssemble : public IReactionAssemble {
private:
  mutable ElementPool<DGElement> elementPool;
  mutable FaceElementPool<DGFaceElement> faceElementPool;

  std::shared_ptr<const DGDiscretization> disc;

  double flux_alpha = 1.0;

  double penalty = 1.0;

  int sign = 1;

public:
  DGReactionAssemble(const IReactionProblem &problem, int degree) :
    disc(std::make_shared<DGDiscretization>(problem.GetMeshes(), degree)),
    IReactionAssemble(problem) {
    Config::Get("flux_alpha", flux_alpha);
    Config::Get("penalty", penalty);
    Config::Get("sign", sign);
  }

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return disc;
  }

  const char *Name() const override { return "DGReactionAssemble"; }

  void PrintInfo() const override {
    int verbose = 1;
    mout.PrintInfo("Assemble", verbose,
                   PrintInfoEntry("Name", Name()),
                   PrintInfoEntry("Problem", problem.Name()),
                   PrintInfoEntry("Convection", problem.GetConvection()),
                   PrintInfoEntry("Diffusion", problem.GetDiffusion()),
                   PrintInfoEntry("Reaction", problem.GetReaction()),
                   PrintInfoEntry("sign", sign),
                   PrintInfoEntry("penalty", penalty));
  }

  void Initialize(Vector &u) const override;

  double Energy(const Vector &u) const override { return 0.0; }

  double Residual(const Vector &u, Vector &r) const override;

  void Jacobi(const Vector &u, Matrix &A) const override;

  double Mass(const Vector &u) const override;

  std::pair<double, double> InflowOutflow(const Vector &u) const override;

  void SetInitialValue(Vector &u) override;

  void PlotIteration(const Vector &u) const override;

  void FinishTimeStep(const Vector &u);
};

#endif
