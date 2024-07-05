#ifndef _DGELASTICITY_H_
#define _DGELASTICITY_H_

#include "STAssemble.hpp"
#include "ViscoElasticProblems.hpp"
#include "SpaceTimeElasticityElement.hpp"
#include "SpaceTimeDiscretization_PGDG.hpp"

class TDGElasticityAssemble : public STAssemble {
protected:
  std::unique_ptr<TProblemElasticity> prob;

  std::shared_ptr<STDiscretization> disc;

  double flux_alpha = 1;

  vector<vector<Point>> v;

  STDiscretization &getDisc() override {
    return *disc;
  }

public:
  TDGElasticityAssemble(const Meshes &meshes, DegreePair degree, const string &problemName) :
    prob(CreateElasticProblemUnique(problemName)),
    STAssemble() {
    Warning("Please remove this class."
            " All functionality should be in STPGViscoElasticAssemble")
    disc = std::make_shared<STDiscretization_PGDG>(meshes, degree, prob->Dim());
    Config::Get("flux_alpha", flux_alpha);
    for (int i = 0; i < 7; ++i) (*disc)().GetDoF().NodalPointsLocal(i, v[i]);
  }

  const IProblem& GetAbstractProblem() const override {
    return *prob;
  }

  const STDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const STDiscretization> GetSharedDisc() const override { return disc; }

  //const TProblemElasticity &GetProblem() const { return *prob; }

  static bool isDGinTime() {
    return false;
  }

  const char *Name() const { return "TDGElasticityAssemble"; }

  void PlotSingleSolution(const Vector &u, const string &filename) const override {
    THROW("Implement PlotSingleSolution for TDGElasticityAssemble")
  }

  int get_space_deg(int m) const {
    if (m == v[0].size()) return 0;
    else if (m == v[1].size()) return 1;
    else if (m == v[2].size()) return 2;
    else if (m == v[3].size()) return 3;
    else if (m == v[4].size()) return 4;
    else if (m == v[5].size()) return 5;
    else if (m == v[6].size()) return 6;
    Exit("not implemented");
  }

  int get_m(int i) const { return v[i].size() * (*disc)().GetDoF().get_dim(); }

  void get_exact_solution(Vector &Exact_Solution) const {
    Exact_Solution = 0;
    for (cell c = Exact_Solution.cells(); c != Exact_Solution.cells_end(); ++c) {
      vector<Point> z = Exact_Solution.GetDoF().GetNodalPoints(*c);
      row r = Exact_Solution.find_row(c());
      int time_deg = Exact_Solution.GetDoF().get_time_deg(*c);
      int NN = r.n() / time_deg;
      for (int i = 0; i < z.size(); ++i) {
        int p = int(i % NN) / v[Exact_Solution.GetDoF().get_space_deg(*c)].size();
        COMPONENT comp = prob->GetComponents()[p];
        Exact_Solution(r, i) = prob->ut(z[i], comp);
      }
    }
    Exact_Solution.Accumulate();
  }

  void PX_P0(const vector<Scalar> &F, vector<Scalar> &C,
             const int &f_deg, const int &c_deg, const VectorMatrixBase &g) const {
    for (int i = 0; i < C.size(); ++i)
      C[i] = 0.0;

    vector<Point> np;
    g.GetDoF().NodalPointsLocal(c_deg, np);
    const Shape &PX_space = disc->GetCellShape(f_deg);
    const int size_localnp = g.GetDoF().NodalPointsLocal(f_deg);
    for (int k = 0; k < prob->Dim(); ++k)
      for (int i = 0; i < np.size(); ++i)
        for (int j = 0; j < size_localnp; ++j)
          C[k * np.size() + i] += F[k * size_localnp + j] * PX_space(np[i], j);
  }

  void System(Matrix &, Vector &) const;

  void SystemAddDoubleD(Matrix &) const { return; }

  std::pair<double, double> DiscNorm(const Matrix &, const Vector &) const;

  double L2Error(Vector &) const;

  double GNError(Vector &) const;

  double L2Norm(Vector &) const;

  void DualRHS_linear(Vector &, const Vector &, Point, Point) const;

  void DualRHS_quadratic(Vector &, const Vector &, Point, Point) const;

  double DualErrorEstimateCell(const cell &, const Vector &, const Vector &) const;

  double DualErrorEstimate(const Vector &u, const Vector &u_star, Vector &Eta) const {
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      Eta(c(), 0) = DualErrorEstimateCell(c, u, u_star);
    double sum = norm(Eta);
    double max = Eta.Max();
    mout << " || Eta ||_infty = " << max << " || Eta ||_2 = " << sum << endl;
    return max;
  }

  double DualErrorEstimateJumpCell(const cell &, const Vector &, const Vector &) const;

  double DualErrorEstimateCell_quadratic(const cell &, const Vector &,
                                         const Vector &, const Vector &) const;

  double DualErrorEstimate_quadratic(const Vector &u, const Vector &e,
                                     const Vector &e_star, Vector &Eta) const {
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      if (c().t() > 0.0)
        Eta(c(), 0) = DualErrorEstimateCell_quadratic(c, u, e, e_star);
    double sum = norm(Eta);
    double max = Eta.Max();
    mout << " || Eta ||_infty = " << max << " || Eta ||_2 = " << sum << endl;
    return max;
  }

  double ErrorEstimateCell(const cell &, const Vector &) const;

  double ErrorEstimate(const Vector &u, Vector &Eta) const {
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      Eta(c(), 0) = ErrorEstimateCell(c, u);
    double sum = norm(Eta);
    double max = Eta.Max();
    mout << " || Eta ||_infty = " << max << " || Eta ||_2 = " << sum << endl;
    return max;
  }

  double Goal_Functional(const Vector &, Point, Point) const;

  double Energy(const Vector &, Point, Point) const;

  Scalar VFlux_c(const VectorField &V, const VectorField &TN) const {
    Scalar TNV = TN * V;
    return TNV;
  }

  Scalar VFlux_cf(const VectorField &V, const VectorField &TN, bool bnd) const {
    Scalar TNV = TN * V;
    if (!bnd) {
      return TNV;
    } else {
      return TNV;
    }
  }

  Scalar SFlux_c(const Tensor &S, const VectorField &TN, const VectorField &N) const {
    VectorField NS = S * N;
    Scalar TNNS = NS * TN;
    return TNNS;
  }

  Scalar SFlux_cf(const Tensor &S, const VectorField &TN,
                  const VectorField &N, bool bnd) const {
    VectorField NS = S * N;
    Scalar TNNS = NS * TN;
    if (!bnd)
      return TNNS;
    else {
      return -TNNS;
      return TNNS;
    }
  }

  void RedistributeCells(LevelPair levels) override {
    redistribute_cells(disc, levels);
  }

};

#endif
