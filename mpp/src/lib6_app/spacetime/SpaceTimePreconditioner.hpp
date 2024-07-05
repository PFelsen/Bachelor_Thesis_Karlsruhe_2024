#ifndef _SPACETIMEPRECONDITIONER_H_
#define _SPACETIMEPRECONDITIONER_H_

#include "m++.hpp"
#include "SpaceTimeTools.hpp"
#include "STAssemble.hpp"
#include "SpaceTimeTransfer.hpp"
#include "LinearSolver.hpp"
#include "STPathManager.hpp"


struct MultiGridData {
  std::unique_ptr<InterpolateWave_dG> CT;
  std::unique_ptr<SpaceTimeTransfer> transfer;
  std::unique_ptr<Vector> u;
  std::unique_ptr<Preconditioner> SM;
  const Matrix * A;
  int smoothing_steps;
  int pre;
  int pos;
  int cycle;
  double theta;
};


class STMultiGridPC : public Preconditioner {
  std::shared_ptr<const STDiscretization> disc;
  const STAssemble &assemble;
  std::shared_ptr<Preconditioner> BP;
  std::shared_ptr<LinearSolver> BS;
  int space_cycle;
  int time_cycle;

  double theta_space = 0.9;
  double theta_time = 0.9;
  int pre_space = 5;
  int pre_time = 5;
  int pos_space = 5;
  int pos_time = 5;

  string name_time = "V";
  string name_space = "V";
  string name_transfer = "Interpolation";
  string name_smoother = "PointBlockGaussSeidel";

  bool dual;
  bool useSparseMatrix = false;

  bool skipbasesolver = false;

private:
  mutable std::vector<MultiGridData> data{};
  std::unique_ptr<PathStrategy> pathStrategy;

  int cycleNameToNumber(const std::string &cycle_name);

public:

  STMultiGridPC(std::unique_ptr<PathStrategy> pathStrategy,
                const STAssemble &assemble,
                const std::string &prefix = "",
                std::shared_ptr<Preconditioner> precond = nullptr,
                std::shared_ptr<LinearSolver> base_solver = nullptr);

  string Name() const override;

  void Construct(const Matrix &A) override;

  void Cycle(int level, Vector &u, Vector &r) const;

  void multiply(Vector &u, const Vector &b) const override;

  void transpose() const override;

  void Destruct() override;

  ~STMultiGridPC() override;
};

#endif
