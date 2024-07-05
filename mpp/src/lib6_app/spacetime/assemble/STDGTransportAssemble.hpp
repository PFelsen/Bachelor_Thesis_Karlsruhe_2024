#ifndef _DGTRANSPORT_H_
#define _DGTRANSPORT_H_

#include <utility>

#include "m++.hpp"
#include "SpaceTimeTools.hpp"
#include "STAssemble.hpp"
#include "CImg.hpp"
#include "SpaceTimeDiscretization_DGDG.hpp"
#include "TransportProblems.hpp"

class STDGTransportAssemble : public STAssemble {
  const double small_move = 1e-12;
protected:
  std::shared_ptr<ITransportProblem> prob;

  std::shared_ptr<STDiscretization> disc;

  STDiscretization &getDisc() override {
    return *disc;
  }

public:

  STDGTransportAssemble(DegreePair degree, std::shared_ptr<ITransportProblem> problem) :
      prob(std::move(problem)), STAssemble() {
    disc = std::make_shared<STDiscretization_DGDG>(prob->GetMeshes(), degree, 1, false);
  }

  void PrintInfo() const override {
    mout.PrintInfo("Assemble", verbose,
                   PrintInfoEntry("Name", Name(), 1),
                   PrintInfoEntry("Problem Name", GetProblem().Name(), 1),
                   PrintInfoEntry("Discretization Name", GetDisc().DiscName(), 1));
  }

  ~STDGTransportAssemble() override = default;

  const STDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const STDiscretization> GetSharedDisc() const override { return disc; }

  const IProblem &GetAbstractProblem() const override {
    return *prob;
  }

  const ITransportProblem &GetProblem() const { return *prob; }

  static bool isDGinTime() { return true; }

  const char *Name() const { return "TDGTransportAssemble_DGT"; }

  void printConformReconstructionL2Error(Vector &u) const {Exit("not implemented"); }

  void MassMatrix(Matrix &M) const override;

  void System(Matrix &, Vector &) const override;

  virtual double MhalfNorm(const Vector &U) const {
    return 0;
    Exit("not implemented");
  }

  virtual double MhalfInvLNorm(const Vector &U) const {
    return 0;
    Exit("not implemented");
  }

  std::pair<double, double> DiscNorm(const Matrix &, const Vector &) const override;

  double L1Error(const Vector &) const override;

  double L2Error(const Vector &) const override;

  double DGError(const Vector &) const override;

  double LInfError(const Vector &u) const override;

  double GNError(const Vector &) const override;

  double EnergyNorm(const Vector &u) const override {
    return 0;
    Exit("not implemented");
  }


  double VNorm(const Matrix &L, const Vector &) const override {
    Exit("not implemented");
  }

  double L1Norm(const Vector &) const override;

  double L2Norm(const Vector &) const override;

  double L2ScalarProduct(const Vector &A, const Vector &B) const override;

  double L2SpaceNormAtTime(const Vector &u, double time) const override;

  double L2SpaceNormAtTimeError(const Vector &u, double time) const override;

  double DGSemiNorm(const Vector &u) const override;

  double DGNormError(const Vector &u) override;

  double DGNorm(const Vector &u, const Matrix &M) const override;

  double MminushalfLuMinusF(const cell &c, const Vector &U) const override {
    return 0;
    Exit("not implemented");
  }

  double MhalfLu_exactMinusF(LevelPair levels) const override {
    return 0;
    Exit("not implemented");
  }

  double DualErrorEstimateJumpCell(const cell &,
                                   const Vector &,
                                   const Vector &) const override {
    return 0;
    Exit("not implemented");
  }

  double ResidualErrorEstimateJumpCell(const cell &c, const Vector &u) const override;

  double L2ReliableResidualErrorEstimateJumpCell(const cell &c,
                                                 const Vector &u,
                                                 const Vector &conf) const override {
    return 0;
    Exit("not implemented");
  }

  double DGErrorEstimateJumpCell(const cell &c,
                                 const Vector &u,
                                 const Vector &conf,
                                 double eta_res_c,
                                 double eta_conf_c) const override {
    return 0;
    Exit("not implemented");
  }


  bool HasConformingLagrangeInterpolation() { return false; }

  void ConformingLagrangeInterpolation(const Vector &u,
                                       Vector &conf,
                                       int degree = 1) const override {
    return;
    Exit("not implemented");
  }

  void ConformingInterpolation(const Vector &u, Vector &conf) const override {
    return;
    Exit("not implemented");
  }

  void ConformingProjectionWithPenalty(const Vector &u,
                                       Vector &u_conf,
                                       Vector &r_conf) const override {
    return;
    Exit("not implemented");
  }

  void RedistributeCells(LevelPair levels) override {
    redistribute_cells(disc, levels);
  }

  double Goal_Functional(const Vector &, Point, Point) const {
    double e = 0.0;
    return e;
  }

  double Energy(const Vector &, Point, Point) const {
    return 0;
    Exit("Not implemented");
  }

  void Energy(const Vector &u, vector<double> &Energy) const {
    for (int n = 0; n < Energy.size(); ++n)
      Energy[n] = n;
  }

  void Mass(const Vector &u,
            double T,
            vector<double> &Mass,
            vector<double> &inflow,
            vector<double> &outflow,
            vector<double> &outflow_bnd) const;

  void get_exact_solution(Vector &Exact_Solution) const override {
    Exact_Solution.Clear();
    const IDoF &dofs = Exact_Solution.GetDoF();
    for (cell c = Exact_Solution.cells(); c != Exact_Solution.cells_end(); ++c) {
      vector<Point> z = Exact_Solution.GetDoF().GetNodalPoints(*c);
      row r = Exact_Solution.find_row(c());
      DegreePair deg = dofs.GetDegree(*c);
      int cnt = 0;
      for (int ti = 0; ti <= deg.time; ti++) {
        for (int s = 0; s < dofs.NodalPointsLocal(deg.space); s++) {
          Exact_Solution(r, cnt) = prob->Solution(z[cnt].t(), *c, z[cnt]);
          cnt++;
        }
      }
    }
    Exact_Solution.Accumulate();
  }

  void get_projected_exact_solution(Vector &Exact_Solution) const {
    get_exact_solution(Exact_Solution);
  }


  void PlotSingleSolution(const Vector &u, const string &filename) const override;

};

#endif
