#ifndef _DGVISCOELASTIC_H_
#define _DGVISCOELASTIC_H_

#include "STAssemble.hpp"
#include "ViscoElasticProblems.hpp"
#include "SpaceTimeViscoElasticElement.hpp"
#include "SpaceTimeDiscretization_PGDG.hpp"

class STPGViscoElasticAssemble : public STAssemble {
protected:
  std::unique_ptr<TProblemViscoElastic> prob;

  std::shared_ptr<STDiscretization> disc;

  STDiscretization &getDisc() override {
    return *disc;
  }

public:
  STPGViscoElasticAssemble(const Meshes &meshes, DegreePair degree, const string &problemName) :
    prob(CreateViscoElasticProblemUnique(problemName)), STAssemble() {
    disc = std::make_shared<STDiscretization_PGDG>(meshes, degree, prob->Dim());
  }

  const STDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const STDiscretization> GetSharedDisc() const override { return disc; }

  const IProblem& GetAbstractProblem() const override {
    return *prob;
  }

  const TProblemViscoElastic &GetProblem() const { return *prob; }

  static bool isDGinTime() {
    return false;
  }

  const char *Name() const { return "STPGViscoElasticAssemble"; }

  void PlotSingleSolution(const Vector &u, const string &filename) const override {
    THROW("Implement PlotSingleSolution for STPGViscoElasticAssemble");
  }

  void get_exact_solution(Vector &Exact_Solution) const override {
    Exact_Solution = 0;
    for (cell c = Exact_Solution.cells(); c != Exact_Solution.cells_end(); ++c) {
      vector<Point> z = Exact_Solution.GetDoF().GetNodalPoints(*c);
      row r = Exact_Solution.find_row(c());
      DegreePair deg = Exact_Solution.GetDoF().GetDegree(*c);
      int cnt = 0;
      for (int pi = 0; pi < prob->Dim(); pi++) {
        COMPONENT comp = prob->GetComponents()[pi];
        for (int s = 0; s < Exact_Solution.GetDoF().NodalPointsLocal(deg.space, pi); s++) {
          Exact_Solution(r, cnt) = prob->ut(z[cnt], comp);
          cnt++;
        }
      }
    }
    Exact_Solution.Accumulate();
  }

  void PX_P0(const vector<Scalar> &F,
             vector<Scalar> &C,
             const int &f_deg,
             const int &c_deg, const VectorMatrixBase &g) const {

    for (int i = 0; i < C.size(); ++i) C[i] = 0.0;

    vector<Point> np;
    g.GetDoF().NodalPointsLocal(c_deg, np);

    const Shape &PX_space = disc->GetCellShape(f_deg);

    int size_localnp = g.GetDoF().NodalPointsLocal(f_deg);

    for (int k = 0; k < prob->Dim(); ++k)
      for (int i = 0; i < np.size(); ++i)
        for (int j = 0; j < size_localnp; ++j)
          C[k * np.size() + i] += F[k * size_localnp + j] * PX_space(np[i], j);
  }

  void System(Matrix &, Vector &) const override;

  void SystemAddDoubleD(Matrix &) const override;

  std::pair<double, double> DiscNorm(const Matrix &, const Vector &) const override;

  double L2Error(const Vector &) const override;

  double GNError(const Vector &) const override;

  double L2Norm(const Vector &) const override;

  //void DualRHS (Vector&, Point, Point) const;
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
    mout << " || Eta ||_infty = " << max << endl
         << " || Eta ||_2     = " << sum << endl;
    return max;
  }

  double ErrorEstimateCell(const cell &, const Vector &) const;

  double ErrorEstimate(const Vector &u, Vector &Eta) const {
    Eta = 0;
    for (cell c = u.cells(); c != u.cells_end(); ++c)
      Eta(c(), 0) = ErrorEstimateCell(c, u);
    double sum = norm(Eta);
    double max = Eta.Max();
    mout << " || Eta ||_infty = " << max << endl
         << " || Eta ||_2     = " << sum << endl;
    return max;
  }

  double Goal_Functional(const Vector &, Point, Point) const;

  double Energy(const Vector &, Point, Point) const;

  Scalar VFlux_c(const VectorField &V, const VectorField &TN) const {
    Scalar TNV = TN * V;
    return TNV;
  }

  Scalar VFlux_cf(const VectorField &V,
                  const VectorField &TN,
                  bool bnd,
                  int bnd_id = 1) const {
    Scalar TNV = TN * V;
    if (!bnd) return TNV;
    else {
      if (bnd_id == 2) return TNV;
      return -TNV;
    }
  }

  Scalar SFlux_c(const Tensor &S, const VectorField &TN, const VectorField &N) const {
    VectorField NS = S * N;
    Scalar TNNS = NS * TN;
    return TNNS;
  }

  Scalar SFlux_cf(const Tensor &S,
                  const VectorField &TN,
                  const VectorField &N,
                  bool bnd,
                  int bnd_id = 1) const {
    VectorField NS = S * N;
    Scalar TNNS = NS * TN;
    if (!bnd) return TNNS;
    else {
      if (bnd_id == 2) return -TNNS; //Neumann bnd condition
      return TNNS;                   // Dirichlet bnd condition
    }
  }

  void measure(const Vector &u) const;

  void RedistributeCells(LevelPair levels) override {
    redistribute_cells(disc, levels);
  }
};

#endif
