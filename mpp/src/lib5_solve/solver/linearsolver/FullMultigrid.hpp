#ifndef TUTORIAL_FULLMULTIGRID_HPP
#define TUTORIAL_FULLMULTIGRID_HPP

#include "LinearSolver.hpp"
#include "Transfers.hpp"

class FullMultigrid : public LinearSolver {
protected:
  class StepData {
  public:
    const IAssemble &assemble;
    std::unique_ptr<ITransfer> transfer = nullptr;
    std::unique_ptr<Preconditioner> smoother = nullptr;
    std::unique_ptr<Matrix> A = nullptr;
    std::unique_ptr<Vector> solution = nullptr;
    std::unique_ptr<Vector> rhs = nullptr;
    std::unique_ptr<Vector> v_coarse = nullptr;
    std::unique_ptr<LinearSolver> solver = nullptr;

    int presmooth = 1;
    int postsmooth = 1;
    int cycle = 1;
    double theta = 1;

    StepData(FullMultigrid *fm, const IAssemble &assemble,
             std::shared_ptr<const IDiscretization> disc, LevelPair level);
  };

  LevelPair fineLevel{};
  LevelPair baseLevel{};

  std::unordered_map<LevelPair, StepData> stepData{};

  double theta = 1;
  int verbose = -1;

  class MultigridPreconditioner : public Preconditioner {
  public:
    FullMultigrid *fm;

    MultigridPreconditioner(FullMultigrid *fm) : fm(fm) {
      Config::Get("MultigridVerbose", verbose);
    };

    void Destruct() {
      // Nothing to do.
    }

    void Construct(const Matrix &M) override {
      // Do nothing, already initialized by FullMultigrid.
    }

    void Cycle(Vector &u, Vector &r) const;

    void multiply(Vector &u, const Vector &b) const {
      Vector r = b;
      u = 0;
      Cycle(u, r);
    }

    virtual string Name() const { return "MultigridPreconditioner"; }
  };
public:
  FullMultigrid() : LinearSolver(nullptr) { Config::Get("MultigridVerbose", verbose); }

  void multiply(Vector &u, const Vector &b) const override;

  void operator()(const ILinearAssemble &assemble,
                  Vector &u) override{THROW("Not implemented for FullMultigrid")}

  LinearSolver &
  operator()(const Matrix &A, const IAssemble &assemble, int reassemblePC) override {

    if (reassemblePC) {
      mout << "Initialize FM" << endl;
      fineLevel = A.Level();
      baseLevel = A.GetDisc().GetMeshes().CoarseLevel();
      if (!stepData.empty()) { stepData.clear(); }
      for (LevelPair level = baseLevel; level != fineLevel.NextInSpace();
           level = level.NextInSpace()) {
        mout << "Initialize StepData on level " << level.space << endl;
        stepData.emplace(level, StepData(this, assemble, A.GetSharedDisc(), level));
      }
    }
    return *this;
  }

  std::string Name() const override { return "FullMultigrid"; };

  void Solve(Vector &u, const Operator &A, const Operator &B, Vector &r) const override {
    THROW("not implemented for FullMultigrid");
  }

  void Solve(Vectors &u, const Operator &A, const Operator &B, Vectors &r) const override {
    THROW("not implemented for FullMultigrid");
  }

  void Solve(Vector &u, const Operator &A, const Operator &B, Vector &&r) const override {
    THROW("not implemented for FullMultigrid");
  }

  void SolveWithInitialValue(Vector &u, const Vector &b) const override {
    THROW("not implemented for FullMultigrid");
  }

  void SolveWithInitialValue(Vectors &u, const Vectors &b) const override {
    THROW("not implemented for FullMultigrid");
  }

  const Preconditioner &GetPreconditioner() override { THROW("not implemented for FullMultigrid"); }

  const int &GetLinearVerbose() override { THROW("not implemented for FullMultigrid"); }

  const int &GetLinearMinIter() override { THROW("not implemented for FullMultigrid"); }

  const int &GetLinearMaxIter() override { THROW("not implemented for FullMultigrid"); }

  const double &GetLinearDefaultEps() override { THROW("not implemented for FullMultigrid"); }

  const double &GetLinearDefaultReduction() override { THROW("not implemented for FullMultigrid"); }

  const std::string &GetLinearPrefix() override { THROW("not implemented for FullMultigrid"); }

  void operator()(const ILinearAssemble &assemble, Vectors &u) override {
    THROW("not implemented for FullMultigrid");
  }

  LinearSolver &operator()(Preconditioner *pc) override {
    THROW("not implemented for FullMultigrid");
  }

  LinearSolver &operator()(const Matrix &A, int reassemblePC = 1) override {
    THROW("not implemented for FullMultigrid");
  }

  LinearSolver &operator()(const Matrix &AA, const Matrix &WW, int reassemblePC = 1) override {
    THROW("not implemented for FullMultigrid");
  }

  void multiply(Vectors &u, const Vectors &b) const override {
    THROW("not implemented for FullMultigrid");
  }

  void multiply_plus(Vector &u, const Vector &b) const override {
    THROW("not implemented for FullMultigrid");
  }

  void multiply_plus(Vectors &u, const Vectors &b) const override {
    THROW("not implemented for FullMultigrid");
  }

  int Iter() const override { THROW("not implemented for FullMultigrid"); }

  void PrintInfo() const override { THROW("not implemented for FullMultigrid"); }

  void SetVerbose(int v) override { verbose = v; }

  void SetSteps(int s) override { THROW("not implemented for FullMultigrid"); }

  void SetEpsilon(double s) override { THROW("not implemented for FullMultigrid"); }

  void SetReduction(double s) override { THROW("not implemented for FullMultigrid"); }

  void SetStartingValue(const Vector &v) override { THROW("not implemented for FullMultigrid"); }

  Iteration GetIteration() const override { return stepData.at(fineLevel).solver->GetIteration(); }
};


#endif
