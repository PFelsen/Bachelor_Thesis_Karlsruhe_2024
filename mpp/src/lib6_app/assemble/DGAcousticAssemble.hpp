#ifndef DGACOUSTICASSEMBLE_HPP
#define DGACOUSTICASSEMBLE_HPP

#include <utility>

#include "AcousticProblems.hpp"
#include "DGAcousticElement.hpp"
#include "DGDiscretization.hpp"
#include "ElementPool.hpp"
#include "IAcousticAssemble.hpp"
#include "ctools.hpp"

class DGAcousticAssemble : public IAcousticAssemble {
protected:
  mutable FaceElementPool<DGAcousticFaceElement> faceElementPool;

  mutable ElementPool<DGAcousticElement> elementPool;

  std::shared_ptr<const DGDiscretization> disc;

  bool fwd = true;

  bool upwind = true;

  bool dampingFlux = false;

  void Damping(const DGAcousticElement &elem,
               DGRowEntries &A_c) const;

  void GradPressure(const DGAcousticElement &elem, DGRowEntries &A_c) const;

  void DivVelocity(const DGAcousticElement &elem, DGRowEntries &A_c) const;

  void assembleDifferentialPartFlux(const Matrix &A,
                                    const cell &c,
                                    DGRowEntries &A_c) const;

  void DirichletBoundaryFlux(DGRowEntries &A_c,
                             const DGAcousticFaceElement &felem) const;

  void NeumannBoundaryFlux(DGRowEntries &A_c,
                           const DGAcousticFaceElement &felem) const;

  void RobinBoundaryFlux(DGRowEntries &A_c,
                         const DGAcousticFaceElement &felem) const;

  void BoundaryFlux(DGRowEntries &A_c, const DGAcousticFaceElement &felem) const;

  void assembleFlux(Matrix &, DGRowEntries &, const cell &) const;

public:
  DGAcousticAssemble(const AcousticProblem &problem, int degree) :
      disc(std::make_shared<const DGDiscretization>(problem.GetMeshes(), degree,
                                                    problem.GetMeshes().dim() +
                                                    problem.NumDamping() + 1)),
      IAcousticAssemble(problem) {
    Config::Get("PDESolverPlotting", plotting);
    Config::Get("upwind", upwind);
    Config::Get("dampingFlux", dampingFlux);
  }

  const IDiscretization &GetDisc() const override { return *disc; };

  std::shared_ptr<const IDiscretization> GetSharedDisc() const override {
    return disc;
  }

  const char *Name() const override { return "DGAcousticAssemble"; }

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose,
                   PrintInfoEntry("Name", Name(), 1),
//                   PrintInfoEntry("Problem Name", problem.Name(), 1),
                   PrintInfoEntry("Start Time", FirstTStep(), 1),
                   PrintInfoEntry("End Time", LastTStep(), 1),
                   PrintInfoEntry("Step size", StepSize(), 1));
  }


  void Initialize(Vector &u) const override;

  void MassMatrix(Matrix &M) const override;

  void SystemMatrix(Matrix &A) const override;

  void RHS(double t, Vector &b) const override;

  double Energy(const Vector &u) const override;

  MinMaxP MinMaxPressure(const Vector &u) const override;

  AcousticNorms ComputeNorms(const Vector &u) const override;

  AcousticErrors ComputeErrors(const Vector &u) const override;

  void InterpolatedExactSolution(Vector &u) const;

  void SetInitialValue(Vector &u) override;

  void FinishTimeStep(const Vector &u) override;

  void PrintIteration(const Vector &u) const override;

  void PlotIteration(const Vector &u) const override;

  void PlotParams(const Vector &u) const;

  void PlotSolution(const Vector &u) const;

  void Interpolate(const Vector &u, Vector &U) const;
};

#endif //DGACOUSTICASSEMBLE_HPP
