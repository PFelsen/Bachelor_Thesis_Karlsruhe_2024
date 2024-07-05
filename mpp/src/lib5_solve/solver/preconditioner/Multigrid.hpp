#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

#include "LinearSolver.hpp"
#include "Transfers.hpp"

class Multigrid : public Preconditioner {
  class StepData {
  public:
    std::unique_ptr<ITransfer> transfer = nullptr;
    std::unique_ptr<Preconditioner> smoother = nullptr;
    const Matrix *A = nullptr; // gets deleted in ~Multigrid if not on finelevel
    std::shared_ptr<Vector> v_coarse = nullptr;
    int presmooth = 1;
    int postsmooth = 1;
    int cycle = 1;
  public:
    StepData() = default;

    StepData(const Matrix *a, const IAssemble &assemble);

    StepData(const Vector &v, const IAssemble &assemble);
  private:
    const Matrix *createMatrix(const Vector &v, const IAssemble &assemble);
  };

  LevelPair fineLevel{};
  LevelPair baseLevel{};

  std::unordered_map<LevelPair, StepData> stepData{};

  std::unique_ptr<LinearSolver> baseSolver = nullptr;
  std::unique_ptr<Vector> baseRHS;
  std::unique_ptr<Matrix> baseA;
  double theta = 1;
public:
  Multigrid() : Preconditioner() {
    Config::Get("MultigridVerbose", verbose);
    Config::Get("SmootherDamp", theta);
  }

  void Construct(const Matrix &) override;

  void Construct(const Matrix &a, const IAssemble &assemble) override;

  void Destruct() override;

  void Cycle(Vector &u, Vector &r) const;

  void multiply(Vector &u, const Vector &b) const override;

  void multiply_transpose(Vector &u, Vector &v) const;

  void Cycle(Vectors &u, Vectors &r) const;

  void multiply(Vectors &u, const Vectors &b) const override;

  void multiply_transpose(Vectors &u, Vectors &v) const;

  ~Multigrid() override { Destruct(); }

  string Name() const { return "Multigrid"; }
};

#endif // MULTIGRID_HPP
