#ifndef TUTORIAL_DGTRANSPORT_HPP
#define TUTORIAL_DGTRANSPORT_HPP

#include "ITransportAssemble.hpp"
#include "DGElement.hpp"
#include "DGDiscretization.hpp"
#include "ElementPool.hpp"


class DGTransportAssemble : public ITransportAssemble {
private:
  mutable ElementPool<DGElement> elementPool;
  mutable FaceElementPool<DGFaceElement> faceElementPool;

  std::shared_ptr<const DGDiscretization> disc;

  double flux_alpha = 1.0;

  double diffusion = 0.0;

  double IFR = 0;

  double OFR = 0;

  double IFR_old = 0;

  double OFR_old = 0;

  double IFR_sum = 0;

  double OFR_sum = 0;

  const double small_move = 1e-12;

  int rkorder = -1;
  
public:
  DGTransportAssemble(const ITransportProblem &problem, int degree) :
      disc(std::make_shared<DGDiscretization>(problem.GetMeshes(), degree)),
      ITransportAssemble(problem) {
    Config::Get("flux_alpha", flux_alpha);
    Config::Get("rkorder", rkorder);
  }

  const IDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override { return disc;  }

  const char *Name() const override { return "DGTransportAssemble"; }

  void PrintIteration(const Vector &u) const override;

  void Initialize(Vector &u) const override;

  void FinishTimeStep(const Vector &u) override;

  void MassMatrix(Matrix &massMatrix) const override;

  void SystemMatrix(Matrix &systemMatrix) const override;

  void RHS(double t, Vector &rhs) const override;

  double Energy(const Vector &u) const override;

  double Mass(const Vector &u) const override;

  double Error(const Vector &u) const override;

  RatePair InflowOutflow() const override;

  RatePair InflowOutflow(const Vector &u) const override;

  double MaxFlux(const Vector &u) const;

  double Flow(const Vector &u) const;
  
  void SetInitialValue(Vector &u) override;

  void PlotIteration(const Vector &u) const override;

  void SetExactSolution(Vector &u_ex) const;
};

#endif //TUTORIAL_DGTRANSPORT_HPP
