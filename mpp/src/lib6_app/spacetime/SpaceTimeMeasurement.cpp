#include "SpaceTimeMeasurement.hpp"

void measurePkt(const Vector &u, const AcousticProblem &problem,
                const STDiscretization &disc, int run) {

  bool dgInTime = disc.isDgInTime();
  int N = -1;
  Config::Get("receiver_count", N, -1);

  if (N == -1) return;

  if (N <= 0) {
    mout
        << "******************************************************************* *"
        << endl;
    mout
        << "* ERROR: invalid number of recievers - WILL NOT STORE MEASUREMENTS! *"
        << endl;
    mout
        << "*********************************************************************"
        << endl;
    return;
  }

  Date Start;
  mout << "Start Measure at " << Start << endl;

  Point r_first = Origin;
  Config::Get("first_receiver", r_first);
  Point r_last = Origin;
  Config::Get("last_receiver", r_last);
  Point time_config = Origin;
  Config::Get("receive_times", time_config);
  std::list<Point> receiver_positions;
  for (int i = 0; i < N; ++i)
    receiver_positions.push_back(r_first + i * (r_last - r_first) / (N - 1));
  FWIObservationSpecification
      spec(time_config[0], time_config[1], time_config[2], receiver_positions);

  //const char* TMP_SPEC_FILE = "test_spec.txt";
  //spec.writeToFile(TMP_SPEC_FILE);
  //FWIObservationSpecification tmp_spec(TMP_SPEC_FILE);
  //FWISeismogram s(tmp_spec, STATE_COMPS);
  //auto& spec = s.getSpecification();

  int probDim = problem.Dim();
  FWISeismogram s(spec, probDim);
  FWISeismogram s_conformReconstructionRadau(spec, probDim);
  vector<vector<double>> ci(6);
  ci[0].resize(2);
  ci[0][0] = 0.0;
  ci[0][1] = 1.0;
  ci[1].resize(3);
  ci[1][0] = 0.0;
  ci[1][1] = 1.0 / 3.0;
  ci[1][2] = 1.0;
  ci[2].resize(4);
  ci[2][0] = 0.0;
  ci[2][1] = 0.4 - sqrt(0.06);
  ci[2][2] = 0.4 + sqrt(0.06);
  ci[2][3] = 1.0;
  ci[3].resize(5);
  ci[3][0] = 0.0;
  ci[3][1] = 0.088588;
  ci[3][2] = 0.409467;
  ci[3][3] = 0.7876595;
  ci[3][4] = 1.0;
  ci[4].resize(6);
  ci[4][0] = 0.0;
  ci[4][1] = 0.057104196114518;
  ci[4][2] = 0.276843013638124;
  ci[4][3] = 0.583590432368917;
  ci[4][4] = 0.860240135656219;
  ci[4][5] = 1.0;
  ci[5].resize(7);
  ci[5][0] = 0.0;
  ci[5][1] = 0.039809857051469;
  ci[5][2] = 0.198013417873608;
  ci[5][3] = 0.437974810247386;
  ci[5][4] = 0.695464273353636;
  ci[5][5] = 0.901464914201174;
  ci[5][6] = 1.0;

  for (auto it = spec.begin(); it != spec.end(); ++it) {
    Point rec(*it);
    for (int i = 0; i < spec.nSteps(); ++i) {
      double t = spec.step(i);
      rec = rec.CopyWithT(t);

      RVector X(0.0, probDim);
      RVector X_conformReconstructionRadau(0.0, probDim);
      int multiple_times = 0;

      for (cell c = u.cells(); c != u.cells_end(); ++c) {
        Point localPT;
        if (transformPointGlobalToLocal(rec, localPT, c)) {
          if (localPT.t() < 1e-14) continue;
          multiple_times++;

          row r = u.find_row(c());
          cell c_prev = c;
          if (c.min() != 0) {
            face f = u.find_face(c.Face(c.Faces() - 2));
            c_prev = u.find_cell(f.Left());
            if (c_prev == u.cells_end())
              c_prev = u.find_overlap_cell(f.Left());
            if (c_prev == u.overlap_end())
              c_prev = c;
          }
          DegreePair deg_c = u.GetDoF().GetDegree(*c);

          if (dgInTime) {
            for (int s = 0; s <= deg_c.time; ++s) {
              int shift = s * r.n() / (deg_c.time + 1);
              int NX = 0;
              for (int p = 0; p < probDim; ++p) {
                const Shape &shape(disc.GetCellShape(deg_c.space));
                const Shape &timeShape(disc.GetTimeShape(deg_c.time));
                for (int n = 0; n < u.GetDoF().NodalPointsLocal(deg_c.space, p); n++) {
                  X[p] += shape(localPT, n) * timeShape(Point(localPT.t()), s)
                          * u(r)[shift + NX + n];
                }
                NX += u.GetDoF().NodalPointsLocal(deg_c.space, p);
              }
            }

            // conform reconstruction Radau
            const Shape &spaceShape(disc.GetCellShape(deg_c.space));
            const Shape &timeShape(disc.GetTimeShape(deg_c.time));
            DegreePair deg_c_prev = u.GetDoF().GetDegree(*c_prev);
            const Shape &spaceShape_prev(disc.GetCellShape(deg_c_prev.space));

            const int dofPerTdeg = r.n() / (deg_c.time + 1);
            for (int pi = 0; pi < probDim; pi++) {
              for (int si = 0; si < spaceShape_prev.size(); si++) {
                if (c.min() != 0.0 && c != c_prev) {
                  row r_prev = u.find_row(c_prev());
                  int shift = deg_c_prev.time * r_prev.n() / (deg_c_prev.time + 1);
                  double tmp = u(r_prev)[shift + pi * spaceShape_prev.size() + si]
                               * spaceShape_prev(localPT, si);
                  for (int tii = 1; tii <= deg_c.time + 1; tii++)
                    tmp *= (localPT.t() - ci[deg_c.time][tii]) / (0.0 - ci[deg_c.time][tii]);
                  X_conformReconstructionRadau[pi] += tmp;
                }
              }
              for (int si = 0; si < spaceShape.size(); si++) {
                for (int ti = 1; ti < deg_c.time + 1; ti++) {
                  double tmp = 0.0;
                  for (int tii = 1; tii <= deg_c.time + 1; tii++) {
                    tmp += u(r)[(tii - 1) * dofPerTdeg + pi * spaceShape.size() + si]
                           * spaceShape(localPT, si) *
                           timeShape(Point(ci[deg_c.time][ti]), tii - 1);
                  }
                  for (int tii = 0; tii <= deg_c.time + 1; tii++) {
                    if (tii == ti) continue;
                    tmp *= (localPT.t() - ci[deg_c.time][tii]) /
                           (ci[deg_c.time][ti] - ci[deg_c.time][tii]);
                  }
                  X_conformReconstructionRadau[pi] += tmp;
                }
                for (int ti = deg_c.time + 1; ti <= deg_c.time + 1; ti++) {
                  double tmp = u(r)[(ti - 1) * dofPerTdeg
                                    + pi * spaceShape.size() + si]
                               * spaceShape(localPT, si)
                               * timeShape(Point(ci[deg_c.time][ti]), ti - 1);
                  for (int tii = 0; tii <= deg_c.time + 1; tii++) {
                    if (tii == ti) continue;
                    tmp *= (localPT.t() - ci[deg_c.time][tii]) /
                           (ci[deg_c.time][ti] - ci[deg_c.time][tii]);
                  }
                  X_conformReconstructionRadau[pi] += tmp;
                }
              }
            }
          } else {
            if (c.min() == 0) {
              vector<Point> c_nodal = u.GetDoF().GetNodalPoints(*c);

              int j = 0;
              int NX = 0;
              for (int p = 0; p < probDim; ++p) {
                COMPONENT comp = problem.GetComponents()[p];
                const Shape &shape(disc.GetCellShape(deg_c.space));
                const Shape &timeShape(disc.GetTimeShape(deg_c.time));
                for (int n = 0; n < u.GetDoF().NodalPointsLocal(deg_c.space, p); n++) {
                  Point &x = c_nodal[j++];
                  X[p] += shape(localPT, n)
                          * timeShape(Point(localPT.t()), 0)
                          * problem.ut(x.t(), x, *c, comp);
                }
                NX += u.GetDoF().NodalPointsLocal(deg_c.space, p);
              }
            } else if (c != c_prev) {
              row r_prev = u.find_row(c_prev());
              DegreePair deg_c_prev = u.GetDoF().GetDegree(*c_prev);
              int shift = (deg_c_prev.time - 1) * r_prev.n() / deg_c_prev.time;
              int NX = 0;
              for (int p = 0; p < probDim; ++p) {
                const Shape &shape(disc.GetCellShape(deg_c_prev.space));
                const Shape &timeShape(disc.GetTimeShape(deg_c_prev.time));
                for (int n = 0; n < u.GetDoF().NodalPointsLocal(deg_c_prev.space, p); n++) {
                  X[p] += shape(localPT, n)
                          * timeShape(Point(localPT.t()), 0)
                          * u(r_prev)[shift + NX + n];
                }
                NX += u.GetDoF().NodalPointsLocal(deg_c_prev.space, p);
              }
            }

            for (int s = 1; s <= deg_c.time; ++s) {
              int shift = (s - 1) * r.n() / deg_c.time;
              int NX = 0;
              for (int p = 0; p < probDim; ++p) {
                const Shape &shape(disc.GetCellShape(deg_c.space));
                const Shape &timeShape(disc.GetTimeShape(deg_c.time));
                for (int n = 0; n < u.GetDoF().NodalPointsLocal(deg_c.space, p); n++) {
                  X[p] += shape(localPT, n)
                          * timeShape(Point(localPT.t()), s)
                          * u(r)[shift + NX + n];
                }
                NX += u.GetDoF().NodalPointsLocal(deg_c.space, p);
              }
            }
          }
        }
      }
      multiple_times = PPM->SumOnCommSplit(multiple_times, u.CommSplit());
      if (multiple_times > 0) {
        X *= 1.0 / double(multiple_times);
        X_conformReconstructionRadau *= 1.0 / double(multiple_times);
      }

      rec = rec.CopyWithT(0);
      s.addMeasurement(rec, t, X);
      if (dgInTime)
        s_conformReconstructionRadau.addMeasurement(rec, t, X_conformReconstructionRadau);
    }
  }

  mout << "     collect measure" << endl;
  s.collect();

  int spaceDeg(0);
  Config::Get("degree", spaceDeg, -1);
  int timeDeg(1);
  Config::Get("time_degree", timeDeg, -1);
  int level(0);
  Config::Get("level", level, -1);
  int plevel(0);
  Config::Get("plevel", plevel, -1);
  int problemLevel(0);
  Config::Get("ProblemLevel", problemLevel, -1);
  int numL(0);
  Config::Get("numL", numL, -1);
  double pml(0.0);
  Config::Get("pml", pml, -1);
  const int InvDtCell = 1.0 / (u.cells().max() - u.cells().min() - 1e-8);
  bool truncateSTMesh = false;
  Config::Get("truncateSTMesh", truncateSTMesh, -1);

  double tau_zero_inv = 0;
  Config::Get<double>("tau_zero_inv", tau_zero_inv);
  double tau_p = 0.1;
  Config::Get<double>("tau_p", tau_p);
  double kappa_factor = 1;
  Config::Get<double>("kappa_factor", kappa_factor);


  string sname = Config::GetDataPath() + "s_ST_l_" + to_string(level)
                 + "_probl_" + to_string(problemLevel)
                 + "_sDeg_" + to_string(spaceDeg)
                 + "_tDeg_" + to_string(timeDeg)
                 + "_nL_" + to_string(numL)
                 + "_pml_" + to_string(pml)
                 + "_dt_" + to_string(InvDtCell)
                 + "_adapt_" + to_string(run)
                 + "_tauz_" + to_string(tau_zero_inv)
                 + "_taup_" + to_string(tau_p)
                 + "_kapaf_" + to_string(kappa_factor);

  if (dgInTime) sname += "_dgInTime";
  if (truncateSTMesh) sname += "_trunc";

  mout << "     write measure to " << sname << endl;
  s.writeToFile(sname);
  mout << "TIME:        measuring time Pkt: " << Date() - Start << endl;

  if (dgInTime) {
    sname += "_confReconstrRadau";
    s_conformReconstructionRadau.collect();
    s_conformReconstructionRadau.writeToFile(sname);
  }
}
