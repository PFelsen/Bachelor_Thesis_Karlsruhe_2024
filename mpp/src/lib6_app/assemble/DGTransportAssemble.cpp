#include "DGTransportAssemble.hpp"
#include "ctools.hpp"
#include "DGElement.hpp"
#include <iomanip>

void DGTransportAssemble::Initialize(Vector &u) const {
  elementPool.Initialize(u);
  faceElementPool.Initialize(u);
}

void DGTransportAssemble::MassMatrix(Matrix &massMatrix) const {
  massMatrix = 0;
  for (cell c = massMatrix.cells(); c != massMatrix.cells_end(); ++c) {
    const DGElement &elem = elementPool.Get(massMatrix, *c);
    DGRowEntries M_c(massMatrix, *c, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      for (int i = 0; i < elem.NodalPoints(); ++i) {
        Scalar Phi_i = elem.Value(q, i);
        for (int j = 0; j < elem.NodalPoints(); ++j) {
          Scalar Phi_j = elem.Value(q, j);
          M_c(i, j) += w * Phi_i * Phi_j;
        }
      }
    }
  }
}

void DGTransportAssemble::SystemMatrix(Matrix &systemMatrix) const {
  systemMatrix = 0;
  for (cell c = systemMatrix.cells(); c != systemMatrix.cells_end(); ++c) {
    const DGElement &elem = elementPool.Get(systemMatrix, *c);
    DGRowEntries A_c(systemMatrix, *c, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      VectorField B = problem.Flux(systemMatrix.Level(), *c, elem.QPoint(q));
      for (int i = 0; i < elem.NodalPoints(); ++i) {
        Scalar Phi_i = elem.Value(q, i);
        for (int j = 0; j < elem.NodalPoints(); ++j) {
          VectorField gradPhi_j = elem.Derivative(q, j);
          A_c(i, j) -= w * (gradPhi_j * B * Phi_i);
        }
      }
    }
    for (int f = 0; f < c.Faces(); ++f) {
      const DGFaceElement &faceElem = faceElementPool.Get(systemMatrix, *c, f);
      if (systemMatrix.OnBoundary(*c, f)) {
        for (int q = 0; q < faceElem.nQ(); ++q) {
          const Point &Qf_c = faceElem.QPoint(q);
          VectorField Nq = faceElem.QNormal(q);
          Scalar BN = problem.FaceNormalFlux(systemMatrix.Level(), *c, f, Nq, Qf_c);
          if (BN >= 0) continue;
          double w = faceElem.QWeight(q);
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            for (int j = 0; j < faceElem.NodalPoints(); ++j) {
              Scalar phi_j = faceElem.Value(q, j);
              A_c(i, j) += w * BN * phi_j * phi_i;
            }
          }
        }
      } else {
        cell cf = systemMatrix.find_neighbour_cell(c, f);
        int f1 = systemMatrix.find_neighbour_face_id(c.Face(f), cf);
        const DGFaceElement &otherFaceElem = faceElementPool.Get(systemMatrix, *cf, f1);
        DGRowEntries A_cf(systemMatrix, *c, *cf);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          const Point &Qf_c = faceElem.QPoint(q);
          VectorField Nq = faceElem.QNormal(q);
          Scalar BN = problem.FaceNormalFlux(systemMatrix.Level(), *c, f, Nq, Qf_c);
          if (BN >= 0) continue;
          int q1 = otherFaceElem.findQPointID(faceElem, Qf_c);
          Point QP1 = otherFaceElem.QPoint(q1);
          double w = faceElem.QWeight(q);
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            for (int j = 0; j < faceElem.NodalPoints(); ++j) {
              Scalar phi_j = faceElem.Value(q, j);
              A_c(i, j) += w * BN * phi_j * phi_i;
            }
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              Scalar phi_j = otherFaceElem.Value(q1, j);
              A_cf(i, j) -= w * BN * phi_j * phi_i;
            }
          }
        }
        if (diffusion == 0.0) {
          continue;
        }
        if (cf() < c()) continue;
        DGRowEntries A_fc(systemMatrix, *cf, *c);
        DGRowEntries A_ff(systemMatrix, *cf, *cf);
        for (int q = 0; q < faceElem.nQ(); ++q) {
          double w = faceElem.QWeight(q);
          const Point &z = faceElem.QPoint(q);
          const Point &N = faceElem.QNormal(q);
          const Point &Qf_c = faceElem.QPoint(q);
          int q1 = otherFaceElem.findQPointID(faceElem, Qf_c);
          double s = 1;
          for (int i = 0; i < faceElem.NodalPoints(); ++i) {
            Scalar phi_i = faceElem.Value(q, i);
            Scalar NDphi_i = diffusion * (faceElem.Derivative(q, i) * N);
            for (int j = 0; j < faceElem.NodalPoints(); ++j) {
              Scalar phi_j = faceElem.Value(q, j);
              Scalar NDphi_j = diffusion * (faceElem.Derivative(q, j) * N);
              A_c(i, j) -= w * (-0.5 * NDphi_i * phi_j
                                - 0.5 * phi_i * NDphi_j
                                + s * phi_i * phi_j);
            }
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              Scalar phi_j = otherFaceElem.Value(q1, j);
              Scalar NDphi_j = diffusion * (otherFaceElem.Derivative(q1, j) * N);
              THROW("Take a look at the scalingfactor s in DGTransportAssemble::SystemMatrix.\nMaybe it should be multiplied by diffusion?")
              A_cf(i, j) -= w * (0.5 * NDphi_i * phi_j
                                 - phi_i * 0.5 * NDphi_j
                                 - s * phi_i * phi_j);
              A_fc(j, i) -= w * (0.5 * NDphi_i * phi_j
                                 - phi_i * 0.5 * NDphi_j
                                 - s * phi_i * phi_j);
            }
          }
          for (int i = 0; i < otherFaceElem.NodalPoints(); ++i) {
            Scalar phi_i = otherFaceElem.Value(q1, i);
            Scalar NDphi_i = diffusion
                             * (otherFaceElem.Derivative(q1, i) * N);
            for (int j = 0; j < otherFaceElem.NodalPoints(); ++j) {
              Scalar phi_j = otherFaceElem.Value(q1, j);
              Scalar NDphi_j = diffusion
                               * (otherFaceElem.Derivative(q1, j) * N);
              A_ff(i, j) -= w * (0.5 * NDphi_i * phi_j
                                 + phi_i * 0.5 * NDphi_j
                                 + s * phi_i * phi_j);
            }
          }
        }
      }
    }
  }
}

void DGTransportAssemble::RHS(double t, Vector &rhs) const {
  rhs.Clear();
  if (!problem.RHS()) return;
  const Mesh &M = rhs.GetMesh();    
  for (cell c = rhs.cells(); c != rhs.cells_end(); ++c) {
    row r = rhs.find_row(c());
    if (problem.InternalInflow()) {
      const DGElement &cellElement = elementPool.Get(rhs, *c);
      for (int q = 0; q < cellElement.nQ(); ++q) {
        const Point &z = cellElement.QPoint(q);
        double w = cellElement.QWeight(q);
        Scalar Inflow = problem.Inflow(t, *c, z);
        if (Inflow > 0) {
          for (int i = 0; i < cellElement.NodalPoints(); ++i) {
            Scalar Phi_i = cellElement.Value(q, i);
            rhs(r, i) += w * Inflow * Phi_i;
          }
        }
      }
    }
    for (int f = 0; f < c.Faces(); ++f) {
      if (!rhs.OnBoundary(*c, f)) continue;
      const DGFaceElement &faceElem = faceElementPool.Get(rhs, *c, f);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        const Point &z = faceElem.QPoint(q);
        double w = faceElem.QWeight(q);
        VectorField N = faceElem.QNormal(q);
        Scalar BN = problem.FaceNormalFlux(rhs.Level(), *c, f, N, z);
        if (BN > 0) continue;
	double g_in = problem.InflowData(M.Level(), *c, t, z, N);
        for (int i = 0; i < faceElem.NodalPoints(); ++i) {
          Scalar Phi_i = faceElem.Value(q, i);
          rhs(r, i) -= w * g_in * Phi_i;
        }
      }
    }
  }
}

double DGTransportAssemble::Energy(const Vector &u) const {
  double energy = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    const DGElement &elem = elementPool.Get(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      energy += w * (U * U);
    }
  }
  return 0.5 * PPM->SumOnCommSplit(energy, u.CommSplit());
}

double DGTransportAssemble::MaxFlux(const Vector &u) const {
  double flux = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    const DGElement &elem = elementPool.Get(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      VectorField B = problem.Flux(u.Level(), *c, elem.QPoint(q));
      flux = std::max(flux,maxNorm(B));
    }
    for (int f = 0; f < c.Faces(); ++f) {
      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);
      for (int q = 0; q < faceElem.nQ(); ++q) {
	VectorField B = problem.Flux(u.Level(), *c, faceElem.QPoint(q));
	flux = std::max(flux,maxNorm(B));
      }
    }
  }
  return PPM->Max(flux, u.CommSplit());
}

double DGTransportAssemble::Flow(const Vector &u) const {
  double flow = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    for (int f = 0; f < c.Faces(); ++f) {
      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);
      for (int q = 0; q < faceElem.nQ(); ++q) {
	VectorField Nq = faceElem.QNormal(q);
	Scalar BN = problem.FaceNormalFlux(u.Level(), *c, f, Nq, faceElem.QPoint(q));
	Scalar U = faceElem.Value(q, u);
	flow += faceElem.QWeight(q) * abs(BN) * U * U;
      }
    }
  }
  return sqrt(PPM->SumOnCommSplit(flow, u.CommSplit()));
}

double DGTransportAssemble::Error(const Vector &u) const {
  double err = 0.0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    const DGElement &elem = elementPool.Get(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar Sol = problem.Solution(Time(), *c,  small_move * c() + (1 - small_move) * elem.QPoint(q));
      Scalar E = elem.Value(q, u) - Sol;
      if (abs(E) > 10 * small_move)
        err += w * E * E;
    }
  }
  return sqrt(PPM->SumOnCommSplit(err, u.CommSplit()));
}

RatePair DGTransportAssemble::InflowOutflow() const {
  return {IFR_sum, OFR_sum};
}

RatePair DGTransportAssemble::InflowOutflow(const Vector &u) const {
  double inflow = 0.0;
  double outflow = 0.0;  
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    BFParts bnd(u.GetMesh(), *c);
    if (!bnd.onBnd()) continue;
    for (int f = 0; f < c.Faces(); ++f) {
      if (bnd[f] < 0) continue;
      const DGFaceElement &faceElem = faceElementPool.Get(u, *c, f);
      for (int q = 0; q < faceElem.nQ(); ++q) {
        double w = faceElem.QWeight(q);
        VectorField Nq = faceElem.QNormal(q);
        Scalar BN = problem.FaceNormalFlux(u.Level(), *c, f, Nq, faceElem.QPoint(q));
        if (BN > 0) outflow += w * BN * faceElem.Value(q, u);
        else        inflow -= w * BN * faceElem.Value(q, u);
      }
    }
  }
  inflow = PPM->SumOnCommSplit(inflow, u.CommSplit());
  outflow = PPM->SumOnCommSplit(outflow, u.CommSplit());
  return {inflow, outflow};
}

double DGTransportAssemble::Mass(const Vector &u) const {
  Scalar e = 0;
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    const DGElement &elem = elementPool.Get(u, *c);
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      Scalar U = elem.Value(q, u);
      e += w * U;
    }
  }
  return PPM->SumOnCommSplit(e, u.CommSplit());
}

void DGTransportAssemble::SetInitialValue(Vector &u) {
  u.Clear();
  for (cell c = u.cells(); c != u.cells_end(); ++c) {
    row r = u.find_row(c());

    const DGElement &elem = elementPool.Get(u, *c);

    for (int j = 0; j < elem.NodalPoints(); ++j)
      u(r, j) = problem.Solution(0, *c, small_move * c() + (1 - small_move) * elem.NodalPoint(j));
  }
}

void DGTransportAssemble::FinishTimeStep(const Vector &u) {
  std::pair<double, double> rate = InflowOutflow(u);
  IFR = rate.first;
  OFR = rate.second;
  if (rkorder == -1) {
    IFR_sum += StepSize() * IFR;
    OFR_sum += StepSize() * OFR;
  } else {
    if (Step() == 0) {
      IFR_old = rate.first;
      OFR_old = rate.second;
    } else {
      IFR_sum += 0.5 * StepSize() * (IFR_old + IFR);
      OFR_sum += 0.5 * StepSize() * (OFR_old + OFR);
    }
  }
  IFR_old = IFR;
  OFR_old = OFR;
  PrintIteration(u);
  PlotIteration(u);
}

void DGTransportAssemble::PrintIteration(const Vector &u) const {
  if (verbose == 0 || Step() % verbose != 0) return;
  std::vector<PrintIterEntry<double>> entries{
      {"n",    (double)Step(),  6,  1},
      {"t",    Time(),  11, 1},
      {"Energy", Energy(u), 11, 1},
      {"Mass", Mass(u), 11, 1},
      {"Flow", Flow(u), 11, 1},
      {"IFR",  IFR,     11, 1},
      {"IF_sum",  IFR_sum, 11, 1},
      {"OFR",  OFR,     11, 1},
      {"OF_sum",  OFR_sum, 11, 1}
  };
  if (problem.HasExactSolution())
    entries.emplace_back("Error", Error(u), 11, 1);
  mout.PrintIteration(verbose, entries);
}

void DGTransportAssemble::PlotIteration(const Vector &u) const {
  if (plotting == 0 || Step() % plotting != 0) return;
  auto cellDisc = std::make_shared<const DGDiscretization>(u.GetDisc().GetMeshes(), 0, 1);
  Vector tmp(0.0, cellDisc, u.Level());
  for (cell c = tmp.cells(); c != tmp.cells_end(); ++c) {
    DGElement elem(u, *c);
    Scalar U = 0.0;
    double a = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      U += w * elem.Value(q, u);
      a += w;
    }
    U *= (1 / a);
    row r = u.find_row(c());
    tmp(r)[0] = U;
  }
  std::string filename = problem.Name() + "." + to_string(Step());
  VtuPlot plot(filename);
  plot.AddData("Density", tmp);
  plot.PlotFile();

  if (!problem.HasExactSolution()) return;
  for (cell c = tmp.cells(); c != tmp.cells_end(); ++c) {
    DGElement elem(u, *c);
    Scalar U = 0.0;
    double a = 0;
    for (int q = 0; q < elem.nQ(); ++q) {
      double w = elem.QWeight(q);
      U += w * problem.Solution(Time(), *c, elem.QPoint(q));
      a += w;
    }
    U *= (1 / a);
    row r = u.find_row(c());
    tmp(r)[0] = U;
  }
  filename = problem.Name() + ".exact." + to_string(Step());
  VtuPlot plot2(filename);
  plot2.AddData("Density", tmp);
  plot2.PlotFile();
}

void DGTransportAssemble::SetExactSolution(Vector &u_ex) const {
  u_ex = 0;
  for (cell c = u_ex.cells(); c != u_ex.cells_end(); ++c) {
    row r = u_ex.find_row(c());
    const DGElement &elem = elementPool.Get(u_ex, *c);
    for (int j = 0; j < elem.NodalPoints(); ++j)
      u_ex(r, j) = problem.Solution(Time(), *c, elem.NodalPoint(j));
  }
  u_ex.Accumulate();
}
