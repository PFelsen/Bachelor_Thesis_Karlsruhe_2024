#ifndef _STDGMATRIX_FREE_VAA_H_
#define _STDGMATRIX_FREE_VAA_H_

#include <utility>

#include "AcousticProblems.hpp"
#include "CImg.hpp"
#include "STAssemble.hpp"
#include "SpaceTimeDiscretization_DGDG.hpp"
#include "SpaceTimeTools.hpp"

#include "m++.hpp"

class STDGMatrixFreeViscoAcousticAssemble : public STAssemble {
protected:
  std::shared_ptr<AcousticProblem> prob;

  std::shared_ptr<STDiscretization> disc;

  double tau_zero_inv = 0;

  double penalty = 1e2;

  STDiscretization &getDisc() override { return *disc; }

public:
  STDGMatrixFreeViscoAcousticAssemble(const Meshes &meshes, DegreePair degree,
                                      std::shared_ptr<AcousticProblem> problem)
      : prob(std::move(problem)), STAssemble() {
    disc = std::make_shared<STDiscretization_DGDG>(meshes, degree, prob->Dim(),
                                                   false);
    Config::Get("tau_zero_inv", tau_zero_inv);
    Config::Get("penalty", penalty);
  }

  void PrintInfo() const override {
    bool WeightedAssemble = false;
    Config::Get("WeightedAssemble", WeightedAssemble);
    mout.PrintInfo("Assemble", verbose, PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Problem Name", GetProblem().Name(), 1),
                   PrintInfoEntry("Discretization Name", GetDisc().DiscName(),
                                  1),
                   PrintInfoEntry("PG(weighted Testfunc)", WeightedAssemble,
                                  1));
  }

  STDGMatrixFreeViscoAcousticAssemble(const Meshes &meshes, DegreePair degree,
                                      const string &problemName)
      : STDGMatrixFreeViscoAcousticAssemble(meshes, degree,
                                            CreateAcousticProblemShared(
                                                problemName)) {}

  ~STDGMatrixFreeViscoAcousticAssemble() override = default;

  const IProblem &GetAbstractProblem() const override { return *prob; }

  void SetProblem(std::shared_ptr<IProblem> problem) override {
    this->prob = std::dynamic_pointer_cast<AcousticProblem>(problem);
  }

  const STDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const STDiscretization> GetSharedDisc() const override { return disc; }

  const AcousticProblem &GetProblem() const { return *prob; }

  static bool isDGinTime() { return true; }

  const char *Name() const { return "STDGMatrixFreeViscoAcousticAssemble"; }

  void get_projected_exact_solution(Vector &Exact_Solution) const override;

  void get_exact_solution(Vector &Exact_Solution) const override {
    Exact_Solution.Clear();
    const IDoF &dofs = Exact_Solution.GetDoF();
    for (cell c = Exact_Solution.cells(); c != Exact_Solution.cells_end();
         ++c) {
      vector<Point> z = Exact_Solution.GetDoF().GetNodalPoints(*c);
      row r = Exact_Solution.find_row(c());
      DegreePair deg = dofs.GetDegree(*c);
      int cnt = 0;
      for (int ti = 0; ti <= deg.time; ti++) {
        for (int pi = 0; pi < prob->Dim(); pi++) {
          COMPONENT comp = prob->GetComponents()[pi];
          for (int s = 0; s < dofs.NodalPointsLocal(deg.space, pi); s++) {
            Exact_Solution(r, cnt) = prob->ut(z[cnt].t(), z[cnt], *c, comp);
            cnt++;
          }
        }
      }
    }
    Exact_Solution.Accumulate();
  }

  void printConformReconstructionL2Error(Vector &u) const;

  void PX_P0(const vector<Scalar> &F, vector<Scalar> &C, const int &f_deg,
             const int &c_deg, const VectorMatrixBase &u) const {
    for (int i = 0; i < C.size(); ++i)
      C[i] = 0.0;

    vector<Point> np;
    u.GetDoF().NodalPointsLocal(c_deg, np);

    const Shape &PX_space = (*disc).GetCellShape(f_deg);

    for (int k = 0; k < prob->Dim(); ++k)
      for (int i = 0; i < np.size(); ++i)
        for (int j = 0; j < u.GetDoF().NodalPointsLocal(f_deg); ++j)
          C[k * np.size() + i] +=
              F[k * u.GetDoF().NodalPointsLocal(f_deg) + j] *
              PX_space(np[i], j);
  }

  void MassMatrix(Matrix &M) const override;

  void System(Matrix &, Vector &) const override;

  void SystemAddDoubleD(Matrix &) const override;

  virtual double MhalfNorm(const Vector &U) const override;

  virtual double MhalfInvLNorm(const Vector &U) const override;

  std::pair<double, double> DiscNorm(const Matrix &,
                                     const Vector &) const override;

  double L1Error(const Vector &) const override;

  double L2Error(const Vector &) const override;

  double DGError(const Vector &) const override;

  double LInfError(const Vector &u) const override;

  double GNError(const Vector &) const override;

  double EnergyNorm(const Vector &u) const override;

  double VNorm(const Matrix &L, const Vector &) const override;

  double L1Norm(const Vector &) const override;

  double L2Norm(const Vector &) const override;

  double L2ScalarProduct(const Vector &A, const Vector &B) const override;

  double L2SpaceNormAtTime(const Vector &u, double time) const override;

  double L2SpaceNormAtTimeError(const Vector &u, double time) const override;

  double DGSemiNorm(const Vector &u) const override;

  double DGNormError(const Vector &u) override;

  double DGNorm(const Vector &u, const Matrix &M) const override;

  void DualRHS_linear(Vector &, const Vector &, Point, Point) const override;

  void DualRHS_quadratic(Vector &, const Vector &, Point, Point) const override;

  double MminushalfLuMinusF(const cell &c, const Vector &U) const override;

  double MhalfLu_exactMinusF(LevelPair levels) const override;

  double DualErrorEstimateJumpCell(const cell &, const Vector &,
                                   const Vector &) const override;

  double ResidualErrorEstimateJumpCell(const cell &c,
                                       const Vector &u) const override;

  double L2ReliableResidualErrorEstimateJumpCell(
      const cell &c, const Vector &u, const Vector &conf) const override;

  double DGErrorEstimateJumpCell(const cell &c, const Vector &u,
                                 const Vector &conf, double eta_res_c,
                                 double eta_conf_c) const override;

  bool HasConformingLagrangeInterpolation() { return true; }

  void ConformingLagrangeInterpolation(const Vector &u, Vector &conf,
                                       int degree = 1) const override;

  void ConformingInterpolation(const Vector &u, Vector &conf) const override;

  void ConformingProjectionWithPenalty(const Vector &u, Vector &u_conf,
                                       Vector &r_conf) const override;

  void RedistributeCells(LevelPair levels) override {
    redistribute_cells(disc, levels);
  }

  double Goal_Functional(const Vector &, Point, Point) const override;

  using STAssemble::Energy;

  double Energy(const Vector &, Point, Point) const override;

  void PlotSingleSolution(const Vector &u,
                          const string &filename) const override;
};

#endif
