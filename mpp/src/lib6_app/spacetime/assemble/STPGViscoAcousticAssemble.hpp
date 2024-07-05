#ifndef _DGVISCOACUSTIC_H_
#define _DGVISCOACUSTIC_H_

#include "AcousticProblems.hpp"
#include "CImg.hpp"
#include "Functools.hpp"
#include "STAssemble.hpp"
#include "SpaceTimeViscoAcousticFaceElement.hpp"
#include "ViscoAcousticElements.hpp"
#include "SpaceTimeDiscretization_PGDG.hpp"
#include "m++.hpp"

class STPGViscoAcousticAssemble : public STAssemble {
protected:
  std::shared_ptr<AcousticProblem> prob;

  std::shared_ptr<STDiscretization> disc;

  mutable Time t_cell;
  mutable Time t_face;

  STDiscretization &getDisc() override {
    return *disc;
  }

public:
  STPGViscoAcousticAssemble(const Meshes &meshes, DegreePair degree,
                            std::shared_ptr<AcousticProblem> problem) :
      prob(std::move(problem)), STAssemble() {
    disc = std::make_shared<STDiscretization_PGDG>(meshes, degree, prob->Dim());
  }

  STPGViscoAcousticAssemble(const Meshes &meshes, DegreePair degree, const string &problemName) :
      STPGViscoAcousticAssemble(meshes, degree, CreateAcousticProblemShared(problemName)) {}

  void SetProblem(std::shared_ptr<IProblem> problem) override {
    this->prob = std::dynamic_pointer_cast<AcousticProblem>(problem);
  }

  using STAssemble::Energy;

  const STDiscretization &GetDisc() const override { return *disc; }

  std::shared_ptr<const STDiscretization> GetSharedDisc() const override { return disc; }

  const IProblem &GetAbstractProblem() const override {
    return *prob;
  }

  const AcousticProblem &GetProblem() const { return *prob; }

  void PlotSingleSolution(const Vector &u, const std::string &filename) const override {
    THROW("Implement PlotSingleSolution for STPGViscoAcousticAssemble");
  }

  static bool isDGinTime() {
    return false;
  }

  const char *Name() const override { return "PGDGViscoAcousticAssemble"; }

  void get_exact_solution(Vector &Exact_Solution) const override {
    Exact_Solution.Clear();
    for (cell c = Exact_Solution.cells(); c != Exact_Solution.cells_end(); ++c) {
      vector<Point> z = Exact_Solution.GetDoF().GetNodalPoints(*c);
      row r = Exact_Solution.find_row(c());
      DegreePair deg = Exact_Solution.GetDoF().GetDegree(*c);
      int NN = r.n() / deg.time;
      for (int i = 0; i < z.size(); ++i) {
        int p = int(i % NN) / Exact_Solution.GetDoF().NodalPointsLocal(deg.space, i);
        COMPONENT comp = prob->GetComponents()[p];
        Exact_Solution(r, i) = prob->ut(z[i].t(), z[i], *c, comp);
      }
    }
    Exact_Solution.Accumulate();
  }

  vector<vector<Scalar>> calc_TMatrix(int mc, int prev_m, int deg, int prev_deg) const;

  void PX_P0(const vector<Scalar> &F, vector<Scalar> &C,
             const int &f_deg, const int &c_deg) const {
    for (int i = 0; i < C.size(); ++i)
      C[i] = 0.0;

    vector<Point> np;
    (*disc)().GetDoF().NodalPointsLocal(c_deg, np);

    const Shape &PX_space = disc->GetCellShape(f_deg);
    const int NodalPointsSize = (*disc)().GetDoF().NodalPointsLocal(f_deg);

    for (int k = 0; k < prob->Dim(); ++k) {
      for (int i = 0; i < np.size(); ++i) {
        for (int j = 0; j < NodalPointsSize; ++j) {
          C[k * np.size() + i] += F[k * NodalPointsSize + j] * PX_space(np[i], j);
        }
      }
    }
  }

  void PX_P0(const Scalar *F, Scalar *C, size_t size_C,
             const int &f_deg, const int &c_deg) const {
    std::fill(C, C + size_C, 0.0);

    vector<Point> np;
    (*disc)().GetDoF().NodalPointsLocal(c_deg, np);

    const Shape &PX_space = disc->GetCellShape(f_deg);

    const int NodalPointsSize = (*disc)().GetDoF().NodalPointsLocal(f_deg);

    for (int k = 0; k < prob->Dim(); ++k) {
      for (int i = 0; i < np.size(); ++i) {
        for (int j = 0; j < NodalPointsSize; ++j) {
          C[k * np.size() + i] += F[k * NodalPointsSize + j] * PX_space(np[i], j);
        }
      }
    }
  }

  void System(Matrix &, Vector &) const override;

  void SystemAddDoubleD(Matrix &) const override;

  std::pair<double, double> DiscNorm(const Matrix &, const Vector &) const override;

  double LInfError(const Vector &) const override;

  double L1Error(const Vector &) const override;

  double L2Error(const Vector &) const override;

  double GNError(const Vector &) const override;

  double L2Norm(const Vector &) const override;

  void DualRHS_linear(Vector &, const Vector &, Point, Point) const override;

  void DualRHS_quadratic(Vector &, const Vector &, Point, Point) const override;

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

  double DualErrorEstimateJumpCell(const cell &, const Vector &, const Vector &) const override;

  double DualErrorEstimateCell_quadratic(const cell &,
                                         const Vector &,
                                         const Vector &,
                                         const Vector &) const;

  double DualErrorEstimate_quadratic(const Vector &u,
                                     const Vector &e,
                                     const Vector &e_star,
                                     Vector &Eta) const {
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

  double Goal_Functional(const Vector &, Point, Point) const override;

  double Energy(const Vector &, Point, Point) const override;

  Scalar VFlux_c(const VectorField &V, const VectorField &TN) const {
    Scalar TNV = TN * V;
    return TNV;
  }

  Scalar VFlux_cf(const VectorField &V, const VectorField &TN,
                  bool bnd, int bnd_id = 2) const {
    Scalar TNV = TN * V;
    if (!bnd) return TNV;
    else {
      if (bnd_id == 1) return TNV;  // Dirichlet bnd
      if (bnd_id == 3) return -TNV; // Robin bnd
      return -TNV;                  // Neumann bnd
    }
  }

  Scalar PFlux_cf(const double &P, bool bnd, int bnd_id = 2) const {
    if (!bnd) return P;
    else {
      if (bnd_id == 1) return -P;  // Dirichlet bnd
      if (bnd_id == 3) return -P;  // Robin bnd
      return P;                    // Neumann bnd
    }
  }

  void fillNodalValuesAtStart(int deg,
                              cell &c,
                              vector<Scalar> &c_nodal_value,
                              const VectorMatrixBase &) const;

  void SystemCell(Vector &RHS, Matrix &M, cell &c) const;

  void RedistributeCells(LevelPair levels) override {
    redistribute_cells(disc, levels);
  }
};

#endif
