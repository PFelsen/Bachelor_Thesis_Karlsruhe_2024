#include <memory>
#include <utility>

#include "MeshesCreator.hpp"
#include "SpaceTimeTools.hpp"
#include "RVector.hpp"
#include "TimeDate.hpp"
#include "SpaceTimeDiscretization.hpp"
#include "M_IOFiles.hpp"
#include "Functools.hpp"
#include "MeshesCreator.hpp"


using namespace std;

// Todo: Is this useful for kernel ? is it already implemented there ?
bool transformPointGlobalToLocal(const Point &globalPoint, Point &localPoint, const cell &CT) {
  //const double eps = 1.0e-12;

  localPoint = CT.GlobalToLocal(globalPoint);
  return CT.PointInCell(globalPoint);
  /*switch (SpaceCellType(CT.Type())) {
    case TRIANGLE:
      return -eps < localPoint[0] && -eps < localPoint[1] && localPoint[0] + localPoint[1] < 1 + eps
             && -eps < localPoint.t() && localPoint.t() < 1 + eps;
    case QUADRILATERAL:
      return -eps < localPoint[0] && -eps < localPoint[1]
             && localPoint[0] < 1 + eps && localPoint[1] < 1 + eps
             && -eps < localPoint.t() && localPoint.t() < 1 + eps;
    default: Exit("not Implemented!");
  }*/

}

// Todo: Is this useful for kernel ? is it already implemented there ?
bool transformPointLocaltoLocal(const Point PF,
                                const cell &CF,
                                Point &PT,
                                const cell &CT) {
  if (CF.Type() != CT.Type()) Exit("types are different");
  return transformPointGlobalToLocal(CF.LocalToGlobal(PF), PT, CT);
}

void redistribute_cells(std::shared_ptr<STDiscretization> &disc,
                        LevelPair levels) {

  int verbose = -1;
  Config::Get("DistributionVerbose", verbose);

  const Meshes &STM = disc->GetMeshes();
  const Mesh &pMesh = STM[STM.pLevel()];
  const Mesh &mesh = STM[levels];

  mout << " communicate cell degrees for load balancing" << endl;
  Date Start_Comm;
  DegreePair default_degree{-1, -1};
  //int default_space_deg(-1), default_time_deg(-1);
  Config::Get("degree", default_degree.space);
  Config::Get("time_degree", default_degree.time);
  int lb_default_weight = 0;
  Config::Get("lb_default_weight", lb_default_weight, -1);
  int probDim = (*disc)(levels).GetDoF().get_dim();

  std::unordered_map<Point, DegreePair> degs;

  PPM->Barrier(levels.commSplit);

  Date Start_Comm_degs_to_master;
  ExchangeBuffer E_degs_send_to_master(levels.commSplit);
  vout(1) << "   send degrees on fine mesh to master: ";

  for (cell c = (*disc)(levels).cells(); c != (*disc)(levels).cells_end(); ++c) {
    degs[c()] = (*disc)(levels).GetDoF().GetDegree(*c);
  }

  if (PPM->Proc(levels.commSplit) != 0) {
    for (cell c = (*disc)(levels).cells(); c != (*disc)(levels).cells_end(); ++c) {
      DegreePair deg = degs[c()];
      E_degs_send_to_master.Send(0) << c();
      E_degs_send_to_master.Send(0) << deg.space;
      E_degs_send_to_master.Send(0) << deg.time;
    }
  }

  E_degs_send_to_master.Communicate();

  if (PPM->Proc(levels.commSplit) == 0) {
    for (short q = 1; q < PPM->Size(levels.commSplit); ++q) {
      while (E_degs_send_to_master.Receive(q).size() < E_degs_send_to_master.Receive(q).Size()) {
        Point p(0., 0., 0.);
        DegreePair deg;
        E_degs_send_to_master.Receive(q) >> p;
        E_degs_send_to_master.Receive(q) >> deg.space;
        E_degs_send_to_master.Receive(q) >> deg.time;
        degs[p] = deg;
      }
    }
  }

  vout(1) << " comm_time " << Date() - Start_Comm_degs_to_master << endl;
  mout << "      Master Process stores " << degs.size() << " spacetime degrees " << endl;

  E_degs_send_to_master.ClearBuffers();

  PPM->Barrier(levels.commSplit);

  std::unordered_map<Point, short> balance_weights;

  Date Start_Comm_dofs_to_master;
  ExchangeBuffer E_cells_plevel(levels.commSplit);
  vout(1) << "   send DoFs on coarse mesh to master: ";

  for (cell c = (*disc)(STM.PLevel()).cells(); c != (*disc)(STM.PLevel()).cells_end(); ++c) {
    balance_weights[c()] = lb_default_weight;
    for (cell cl = (*disc)(levels).cells(); cl != (*disc)(levels).cells_end(); ++cl) {
      Point P;
      if (transformPointGlobalToLocal(cl(), P, c)) {
        if (balance_weights[c()] < (*disc)(levels).GetDoF().NumberOfNodalPoints(*cl))
          balance_weights[c()] = (*disc)(levels).GetDoF().NumberOfNodalPoints(*cl);
      }
    }
    if (balance_weights[c()] == -1) {
      Exit("NodalPoints not found");
    }
  }

  if (PPM->Proc(levels.commSplit) != 0)
    for (cell c = (*disc)(STM.PLevel()).cells(); c != (*disc)(STM.PLevel()).cells_end(); ++c) {
      E_cells_plevel.Send(0) << c();
      E_cells_plevel.Send(0) << balance_weights[c()];
    }

  E_cells_plevel.Communicate();

  if (PPM->Proc(levels.commSplit) == 0)
    for (short q = 1; q < PPM->Size(levels.commSplit); ++q) {
      while (E_cells_plevel.Receive(q).size() < E_cells_plevel.Receive(q).Size()) {
        Point p(0., 0., 0.);
        short w = -1;
        E_cells_plevel.Receive(q) >> p;
        E_cells_plevel.Receive(q) >> w;
        balance_weights[p] = w;
      }
    }

  E_cells_plevel.ClearBuffers();
  vout(1) << " comm_time " << Date() - Start_Comm_dofs_to_master << endl;

  std::shared_ptr<Meshes> newSTM =
      MeshesCreator()
          .WithBalanceWeights(std::move(balance_weights))
          .WithTimeRefinement(STM.Settings().timeRefinement)
          .CreateShared();

  // Todo: this line is nasty, has to go. find a way to manipulate a ref
  disc = (*disc).create(*newSTM, default_degree,probDim);

  PPM->Barrier(levels.commSplit);

  std::unordered_map<Point, short> proc_id;

  Date Start_Comm_proc_id_to_master;
  ExchangeBuffer E_proc_id(levels.commSplit);
  vout(1) << "   send proc id on coarse mesh to master: ";

  if (PPM->Proc(levels.commSplit) != 0)
    for (cell c = (*disc)(levels).cells(); c != (*disc)(levels).cells_end(); ++c) {
      E_proc_id.Send(0) << c();
      E_proc_id.Send(0) << short(PPM->Proc(levels.commSplit));
    }

  E_proc_id.Communicate();

  if (PPM->Proc(levels.commSplit) == 0)
    for (short q = 1; q < PPM->Size(levels.commSplit); ++q) {
      while (E_proc_id.Receive(q).size() < E_proc_id.Receive(q).Size()) {
        Point p(0., 0., 0.);
        short w = -1;
        E_proc_id.Receive(q) >> p;
        E_proc_id.Receive(q) >> w;
        proc_id[p] = w;
      }
    }

  E_proc_id.ClearBuffers();
  vout(1) << " comm_time " << Date() - Start_Comm_proc_id_to_master << endl;

  PPM->Barrier(levels.commSplit);

  Date Start_Comm_degs_from_master;
  ExchangeBuffer E_degs_get_from_master(levels.commSplit);
  vout(1) << "   get degrees on fine mesh from master: ";

  if (PPM->Proc(levels.commSplit) == 0)
    for (auto &[point, deg] : degs) {
      short cell_proc_id = proc_id[point];
      E_degs_get_from_master.Send(cell_proc_id) << point;
      E_degs_get_from_master.Send(cell_proc_id) << deg.space;
      E_degs_get_from_master.Send(cell_proc_id) << deg.time;
    }

  E_degs_get_from_master.Communicate();

  if (PPM->Proc(levels.commSplit) != 0)
    while (E_degs_get_from_master.Receive(0).size() < E_degs_get_from_master.Receive(0).Size()) {
      Point p(0., 0., 0.);
      DegreePair deg{-1, -1};
      E_degs_get_from_master.Receive(0) >> p;
      E_degs_get_from_master.Receive(0) >> deg.space;
      E_degs_get_from_master.Receive(0) >> deg.time;
      degs[p] = deg;
    }

  E_degs_get_from_master.ClearBuffers();
  vout(1) << " comm_time " << Date() - Start_Comm_degs_from_master << endl;

  PPM->Barrier(levels.commSplit);

  mout << " communicate cell degrees for load balancing done" << endl;
  vout(1) << " comm_time " << Date() - Start_Comm << endl;

  for (cell c = (*disc)(levels).cells(); c != (*disc)(levels).cells_end(); ++c) {
    auto iter = degs.find(c());
    if (iter != degs.end()) {
      DegreePair deg = iter->second;
      (*disc)(levels).GetDoF().SetDegree(c(), deg);
    }
  }

  disc->communicate(levels);
  (*disc).adaptQuadrature(levels, -1, true);
}

void increase_poly_degs_uniform(STDiscretization &discretization) {
  for (cell c = discretization().cells(); c != discretization().cells_end(); ++c) {
    discretization().GetDoF().increase_space_deg(*c);
    discretization().GetDoF().increase_time_deg(*c);
  }
}

std::vector<double> get_cumulative_cells(const Vector &vec) {
  const size_t SIZE = 10000 + 1;
  std::vector<double> histogram(SIZE, 0.0);
  double max_value = abs(vec(vec.cells()(), 0));
  for (cell c = vec.cells(); c != vec.cells_end(); ++c) {
    max_value = max(max_value, abs(vec(c(), 0)));
  }
  max_value = PPM->Max(max_value, vec.CommSplit());
  for (cell c = vec.cells(); c != vec.cells_end(); ++c) {
    double percent = abs(vec(c(), 0)) / max_value;
    int percent_int = int((histogram.size() - 1) * percent);
    histogram[percent_int]++;
  }
  ExchangeBuffer buffer(vec.CommSplit());
  if (!PPM->Master(vec.CommSplit())) {
    buffer.Send(0) << histogram;
  }
  buffer.Communicate();
  if (PPM->Master(vec.CommSplit())) {
    for (int q = 1; q < PPM->Size(vec.CommSplit()); q++) {
      std::vector<double> temp;
      buffer.Receive(q) >> temp;
      for (int i = 0; i < temp.size(); i++) {
        histogram[i] += temp[i];
      }
    }
  }
  buffer.ClearBuffers();
  buffer = ExchangeBuffer(vec.CommSplit());

  std::vector<double> cumulated_histrogram(histogram.size(), 0.0);
  if (PPM->Master(vec.CommSplit())) {
    cumulated_histrogram[histogram.size() - 1] = histogram[histogram.size() - 1];
    for (int i = int(histogram.size()) - 2; i >= 0; i--) {
      cumulated_histrogram[i] += cumulated_histrogram[i + 1] + histogram[i];
    }
    for (int q = 1; q < PPM->Size(vec.CommSplit()); q++) {
      buffer.Send(q) << cumulated_histrogram;
    }
  }
  buffer.Communicate();
  if (!PPM->Master(vec.CommSplit())) {
    buffer.Receive(0) >> cumulated_histrogram;
  }
  buffer.ClearBuffers();
  return cumulated_histrogram;
}

void print_vec_distribution(const Vector &vec) {
  int cell_count = vec.GetMesh().CellCountGeometry();
  std::vector<double> cumulated_histrogram = get_cumulative_cells(vec);
  for (int i = int(cumulated_histrogram.size()) - 1; i >= 0; i--) {
    std::string s = to_string(i) + " / " + to_string(cumulated_histrogram.size() - 1);
    std::string percent_cells = to_string(cumulated_histrogram[i] / cell_count);
    percent_cells = percent_cells.substr(0, 5);
    mout << "theta: " << s << "% ~~ "
         << percent_cells << "% of cells refined" << endl;
  }

}

double calculateEtaCrit(const Vector &Eta, std::string name, double theta){
  int cell_count = Eta.GetMesh().CellCountGeometry();
  double eta_crit = std::numeric_limits<double>::max();

  double eta_sum = 0.0;
  double eta_max = 0.0;
  int n_dof = 0;
  for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
    ++n_dof;
    eta_sum += abs(Eta(c(), 0));
    eta_max = max(eta_max, abs(Eta(c(), 0)));
  }
  n_dof = PPM->SumOnCommSplit(n_dof, Eta.CommSplit());
  eta_sum = PPM->SumOnCommSplit(eta_sum, Eta.CommSplit());
  eta_max = PPM->Max(eta_max, Eta.CommSplit());
  mout << "  | Eta |_1 =      " << eta_sum << endl;
  //	 << " eta_max = " << eta_max << endl;

  if (name == "abs_value") {
    eta_crit = eta_max * theta;
  } else if (name == "my_percentage") {
    if (theta == 1.0) {
      eta_crit = 0;
    }else {
      bool print_eta_dist = false;
      Config::Get("print_eta_dist", print_eta_dist);
      std::vector<double> cumulative_cells = get_cumulative_cells(Eta);
      for (int i = int(cumulative_cells.size()) - 1; i >= 0; i--) {
        if (print_eta_dist)
          mout << i << " " << cumulative_cells[i] / cell_count << " >? " << theta << endl;
        if (cumulative_cells[i] / cell_count > theta) {
          eta_crit = (eta_max * i) / (cumulative_cells.size() - 1);
          mout << "choosen theta = " << i / (cumulative_cells.size() - 1.0) << endl;
          break;
        }

        if (i == 1 || i == 0) {
          eta_crit = eta_max / (cumulative_cells.size() - 1);
          mout << "choosen theta = " << 1.0 / (cumulative_cells.size() - 1) << "(saved)" << endl;
        }
      }
    }

  } else if (name == "percentage") {
    THROW("implement percentage-adaptive strategy.")
  } else if (name == "equidistribution") {
    eta_crit = eta_sum / n_dof;
  }
  int verbose = 0;
  Config::Get("AdaptivityVerbose", verbose, true);

  vout(1) << DOUT(eta_crit) << DOUT(eta_max) << DOUT(eta_sum) << DOUT(theta)<< endl;
  return eta_crit;
}

std::unordered_map<Point, DegreePair>
    createPolynomialDistribution(STDiscretization &discretization,
                        const Vector &Eta,
                        const double theta,
                        const double theta_min,
                        const string name) {
  std::unordered_map<Point, DegreePair> polyDistribution;
  int cell_count = Eta.GetMesh().CellCountGeometry();

  bool refine_neighbours = false;
  Config::Get("refine_neightbours", refine_neighbours);

  int max_space_deg = 6;
  int max_time_deg = 6;

  double eta_crit = calculateEtaCrit(Eta, name, theta);

  int n_ref = 0;
  int n_deref = 0;
  for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
    DegreePair deg = Eta.GetDoF().GetDegree(c());
    polyDistribution[c()] = deg;
    if (abs(Eta(c(), 0)) < eta_crit * theta_min) {
      ++n_deref;
      deg = deg.decreaseBoth();
      polyDistribution[c()] = deg;
      continue;
    }
    if (abs(Eta(c(), 0)) > eta_crit) {
      if (deg.space < max_space_deg) {
        ++n_ref;
        polyDistribution[c()] = {(short)(deg.space + 1),
                                 (short)max(deg.space + 1, 1)};
      }
    } else if (refine_neighbours) {
      if (deg.space < max_space_deg || deg.time < max_time_deg) {
        for (int ff = 0; ff < c.Faces(); ++ff) {
          face f = Eta.find_face(c.Face(ff));
          cell cf = Eta.find_cell(f.Left());
          if (Eta.find_cell(f.Left()) == Eta.cells_end())
            cf = Eta.find_overlap_cell(f.Left());
          if (cf() == c()) {
            if (Eta.find_cell(f.Right()) == Eta.cells_end()) {
              // c has no neighbour at f on the same proc
              auto overlap_right = Eta.find_overlap_cell(f.Right());
              if (overlap_right != Eta.overlap_end())
                cf = overlap_right; // c has neighbour at f on a different proc
              else
                continue; // f is a boundary face
            } else
              cf = Eta.find_cell(f.Right()); // c has a neighbour at f on the same proc
          }

          if (abs(Eta(cf(), 0)) > eta_crit) {
            if (deg.space < max_space_deg) {
              ++n_ref;
              polyDistribution[c()] = {(short)(deg.space + 1),
                                       (short)max(deg.space + 1, 1)};
            }
            break;
          }
        }
      }
    }
  }
  n_ref = PPM->SumOnCommSplit(n_ref, Eta.CommSplit());
  n_deref = PPM->SumOnCommSplit(n_deref, Eta.CommSplit());
  int n_cell = Eta.GetMesh().CellCount();
  n_cell = PPM->SumOnCommSplit(n_cell, Eta.CommSplit());
  mout << " derefine " << n_deref
       << " cells and refine degree on " << n_ref << " cells of "
       << n_cell << " space-time cells " << endl;
  return polyDistribution;
}

void increase_poly_degs(STDiscretization &discretization,
                        const Vector &Eta,
                        const double theta,
                        const double theta_min,
                        const string &name) {
  bool print_eta_dist = false;
  Config::Get("print_eta_dist", print_eta_dist);

  int cell_count = Eta.GetMesh().CellCountGeometry();
  //print_vec_distribution(Eta);

  bool refine_neighbours = false;
  Config::Get("refine_neightbours", refine_neighbours);

  int max_space_deg = 6;
  int max_time_deg = 6;
  double eta_crit = std::numeric_limits<double>::max();

  double eta_sum = 0.0;
  double eta_max = 0.0;
  int n_dof = 0;
  for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
    ++n_dof;
    eta_sum += abs(Eta(c(), 0));
    eta_max = max(eta_max, abs(Eta(c(), 0)));
  }
  n_dof = PPM->SumOnCommSplit(n_dof, Eta.CommSplit());
  eta_sum = PPM->SumOnCommSplit(eta_sum, Eta.CommSplit());
  eta_max = PPM->Max(eta_max, Eta.CommSplit());
  mout << "  | Eta |_1 =      " << eta_sum << endl;
  //	 << " eta_max = " << eta_max << endl;

  if (name == "abs_value") {
    eta_crit = eta_max * theta;
  } else if (name == "my_percentage") {
    if (theta == 1.0) {
      eta_crit = 0;
    }else {
      std::vector<double> cumulative_cells = get_cumulative_cells(Eta);
      for (int i = int(cumulative_cells.size()) - 1; i >= 0; i--) {
        if (print_eta_dist)
          mout << i << " " << cumulative_cells[i] / cell_count << " >? " << theta << endl;
        if (cumulative_cells[i] / cell_count > theta) {
          eta_crit = (eta_max * i) / (cumulative_cells.size() - 1);
          mout << "choosen theta = " << i / (cumulative_cells.size() - 1.0) << endl;
          break;
        }

        if (i == 1 || i == 0) {
          eta_crit = eta_max / (cumulative_cells.size() - 1);
          mout << "choosen theta = " << 1.0 / (cumulative_cells.size() - 1) << "(saved)" << endl;
        }
      }
    }

  } else if (name == "percentage") {
    increase_poly_degs_by_percentage(discretization, Eta, theta);
    return;
  } else if (name == "equidistribution") {
    eta_crit = eta_sum / n_dof;
  } else if (name == "none") {
    return;
  }
  int n_ref = 0;
  int n_deref = 0;
  for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
    if (abs(Eta(c(), 0)) < eta_crit * theta_min) {
      ++n_deref;
      Eta.GetDoF().decrease_space_deg(*c);
      Eta.GetDoF().decrease_time_deg(*c);
      continue;
    }
    DegreePair deg = Eta.GetDoF().GetDegree(*c);
    if (abs(Eta(c(), 0)) > eta_crit) {
      if (deg.space < max_space_deg) {
        ++n_ref;
        Eta.GetDoF().increase_space_deg(*c);
        Eta.GetDoF().set_time_deg(*c, max(deg.space + 1, 1));
      }
    } else if (refine_neighbours) {
      if (deg.space < max_space_deg || deg.time < max_time_deg) {
        for (int ff = 0; ff < c.Faces(); ++ff) {
          face f = Eta.find_face(c.Face(ff));
          cell cf = Eta.find_cell(f.Left());
          if (Eta.find_cell(f.Left()) == Eta.cells_end())
            cf = Eta.find_overlap_cell(f.Left());
          if (cf() == c()) {
            if (Eta.find_cell(f.Right()) == Eta.cells_end()) {
              // c has no neighbour at f on the same proc
              auto overlap_right = Eta.find_overlap_cell(f.Right());
              if (overlap_right != Eta.overlap_end())
                cf = overlap_right; // c has neighbour at f on a different proc
              else
                continue; // f is a boundary face
            } else
              cf = Eta.find_cell(f.Right()); // c has a neighbour at f on the same proc
          }

          if (abs(Eta(cf(), 0)) > eta_crit) {
            if (deg.space < max_space_deg) {
              ++n_ref;
              Eta.GetDoF().increase_space_deg(*c);
              int tdeg = max(deg.space + 1, 1);
              Eta.GetDoF().set_time_deg(*c, tdeg);
            }
            break;
          }
        }
      }
    }
  }
  n_ref = PPM->SumOnCommSplit(n_ref, Eta.CommSplit());
  n_deref = PPM->SumOnCommSplit(n_deref, Eta.CommSplit());
  int n_cell = Eta.GetMesh().CellCount();
  n_cell = PPM->SumOnCommSplit(n_cell, Eta.CommSplit());
  mout << " derefine " << n_deref
       << " cells and refine degree on " << n_ref << " cells of "
       << n_cell << " space-time cells " << endl;
}

void increase_poly_degs_by_percentage(STDiscretization &disc_vector,
                                      const Vector &Eta,
                                      const double per_crit) {
  int max_space_deg = 6;
  int max_time_deg = 6;

  Vector increased(0.0, Eta);
  increased.SetAccumulateFlag(false);
  Vector increased_add(0.0, Eta);

  struct DataContainer {
    int q;
    Point p;
    double eta;

    DataContainer() {
      q = -1;
      p = Point(0, 0);
      eta = 0;
    }
  };

  struct {
    bool operator()(DataContainer const &a, DataContainer const &b) {
      return abs(a.eta) > abs(b.eta);
    }
  } DataBigger;

  vector<DataContainer> eta_sort(PPM->SumOnCommSplit((Eta.GetMesh()).CellCount(), Eta.CommSplit()));
  ExchangeBuffer Eta_send_to_master(Eta.CommSplit());

  unsigned int iter = 0;
  if (PPM->Proc(Eta.CommSplit()) != 0) {
    for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
      DegreePair deg = Eta.GetDoF().GetDegree(*c);
      if (deg.time < max_time_deg && deg.space < max_space_deg) {
        Eta_send_to_master.Send(0) << c();
        Eta_send_to_master.Send(0) << abs(Eta(c(), 0));
      }
    }
  } else {
    for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
      DegreePair deg = Eta.GetDoF().GetDegree(*c);
      if (deg.time < max_time_deg && deg.space < max_space_deg) {
        DataContainer d;
        d.eta = Eta(c(), 0);
        d.q = 0;
        d.p = c();
        eta_sort[iter++] = d;
      }
    }
  }
  Eta_send_to_master.Communicate();

  if (PPM->Proc(Eta.CommSplit()) == 0)
    for (short q = 1; q < PPM->Size(Eta.CommSplit()); ++q) {
      while (Eta_send_to_master.Receive(q).size() < Eta_send_to_master.Receive(q).Size()) {
        Point p(0., 0., 0.);
        DataContainer d;
        d.q = q;
        Eta_send_to_master.Receive(q) >> p;
        Eta_send_to_master.Receive(q) >> d.eta;
        d.p = p;
        eta_sort[iter++] = d;
      }
    }

  sort(eta_sort.begin(), eta_sort.end(), DataBigger);

  ExchangeBuffer Eta_send_to_procs(Eta.CommSplit());

  int cnt = 0;
  for (unsigned int i = 0; i < per_crit * eta_sort.size(); ++i) {
    if (eta_sort[i].q >= 0)
      Eta_send_to_procs.Send(eta_sort[i].q) << eta_sort[i].p;
  }

  Eta_send_to_procs.Communicate();

  while (Eta_send_to_procs.Receive(0).size() < Eta_send_to_procs.Receive(0).Size()) {
    Point p(0., 0., 0.);
    Eta_send_to_procs.Receive(0) >> p;
    cell c = Eta.find_cell(p);
    Eta.GetDoF().increase_time_deg(*c);
    Eta.GetDoF().increase_space_deg(*c);
    increased(p, 0) = 1;
  }

  increased.Accumulate();

  for (cell c = Eta.cells(); c != Eta.cells_end(); ++c) {
    DegreePair deg = Eta.GetDoF().GetDegree(*c);
    if (increased(c(), 0) == 0.0) {
      if (deg.space < max_space_deg || deg.time < max_time_deg)
        for (int ff = 0; ff < c.Faces(); ++ff) {
          face f = Eta.find_face(c.Face(ff));
          cell cf = Eta.find_cell_or_overlap_cell(f.Left());
          if (cf() == c()) {
            if (Eta.find_cell(f.Right()) == Eta.cells_end()) {
              // c has no neighbour at f on the same proc
              if (Eta.find_overlap_cell(f.Right()) != Eta.overlap_end()) {
                // c has neighbour at f on a different proc
                cf = Eta.find_overlap_cell(f.Right());
              } else {
                continue; // f is a boundary face
              }
            } else {
              cf = Eta.find_cell(f.Right()); // c has a neighbour at f on the same proc
            }
          }

          if (increased(cf(), 0) == 1.0) {
            if (deg.time < max_time_deg)
              Eta.GetDoF().increase_time_deg(*c);
            if (deg.space < max_space_deg)
              Eta.GetDoF().increase_space_deg(*c);
            increased_add(c(), 0) = 1;
            break;
          }
        }
    }
  }
  increased.MakeAdditive();
  increased_add.MakeAdditive();
  unsigned int refined_cells = increased * increased + increased_add * increased_add;
  mout << refined_cells << " of " << eta_sort.size() << " ("
       << refined_cells / (double) eta_sort.size() * 100.0
       << "\%) space-time cells refined " << endl;
}

void write_pattern(Vector &U, STDiscretization &disc, double x_pos) {
  int NN = U.GetDoF().get_dim();
  int comp = 0;
  if (NN == 3) comp = 2;

  ExchangeBuffer E_send_to_master(U.CommSplit());

  for (cell c = U.cells(); c != U.cells_end(); ++c) {
    int f = 0;
    for (; f < c.Faces(); ++f)
      if (c.Face(f)[0] == x_pos) break;
    if (f == c.Faces()) continue;

    vector<Point> z = U.GetDoF().GetNodalPoints(*c);
    row r = U.find_row(c());
    int time_deg = U.GetDoF().GetDegree(*c).time;
    int size = r.n() / (NN * time_deg);

    for (int td = 0; td < time_deg; ++td)
      for (int iN = 0; iN < NN; ++iN)
        for (int i = 0; i < size; ++i)
          if (z[td * size * NN + iN * size + i][0] == x_pos && iN == comp) {
            E_send_to_master.Send(0)
                << float(U(r, td * size * NN + iN * size + i));
            E_send_to_master.Send(0)
                << float(z[td * size * NN + iN * size + i][1]);
            E_send_to_master.Send(0)
                << float(z[td * size * NN + iN * size + i].t());
          }
  }

  E_send_to_master.Communicate();

  if (PPM->Master(U.CommSplit())) {
    string DataPathName("data");
    string filename = DataPathName + string("/csv/pattern.csv");
    M_ofstream out(filename.c_str());
    out << "value,y,t" << endl;

    for (short q = 0; q < PPM->Size(U.CommSplit()); ++q) {
      while (E_send_to_master.Receive(q).size() < E_send_to_master.Receive(q).Size()) {
        float value, y_coord, time;
        E_send_to_master.Receive(q) >> value;
        E_send_to_master.Receive(q) >> y_coord;
        E_send_to_master.Receive(q) >> time;
        out << value << "," << y_coord << "," << time << endl;
      }
    }
  }

  E_send_to_master.ClearBuffers();
}

void NewCellTransfer_2D::Interpolate(const Vector &F,
                                     Vector &C,
                                     string s,
                                     bool testspace = false) const {
  C = 0.;
  if (s == "Res") {
    cell cc_test = C.cells();

    vector<Point> Children_in_time_test(2);
    vector<Point> Children_in_space_test(cc_test.Children());
    Children_in_time_test[0] = cc_test().CopyWithT(0.5 * (cc_test().t() + cc_test.min()));
    Children_in_time_test[1] = cc_test().CopyWithT(0.5 * (cc_test().t() + cc_test.max()));
    for (int i = 0; i < cc_test.Children(); ++i) {
      Children_in_space_test[i] = cc_test.Child(i).
          CopyWithT(cc_test().t());
    }
    bool in_time = true;
    bool in_space = true;
    for (int i = 0; i < Children_in_time_test.size(); ++i) {
      cell cc_test = F.find_cell(Children_in_time_test[i]);
      if (cc_test == F.cells_end()) in_time = false;
    }
    for (int i = 0; i < Children_in_space_test.size(); ++i) {
      cell cc_test = F.find_cell(Children_in_space_test[i]);
      if (cc_test == F.cells_end()) in_space = false;
    }
    if (in_time == in_space) {
      vector<Point> Children(2 * cc_test.Children());
      for (int i = 0; i < cc_test.Children(); ++i) {
        double t_min = 0.5 * (cc_test().t() + cc_test.min());
        double t_max = 0.5 * (cc_test().t() + cc_test.max());
        Children[i] = cc_test.Child(i).CopyWithT(t_min);
        Children[i + cc_test.Children()] = cc_test.Child(i).CopyWithT(t_max);
      }
      bool in_st = true;
      for (int i = 0; i < Children.size(); ++i) {
        cell cc_test = F.find_cell(Children[i]);
        if (cc_test == F.cells_end()) in_st = false;
      }
      if (in_st) {
        in_space = true;
        in_time = true;
      } else {
        if (cc_test == F.cells()) mout << "Interpolate on same level" << endl;
        else Exit("Error in SpaceTimeTools.C\n");
      }
    }

    for (cell cc = C.cells(); cc != C.cells_end(); ++cc) {
      row rc = C.find_row(cc());
      DegreePair deg_c = C.GetDoF().GetDegree(*cc);

      vector<Point> Children;
      if (in_time && !in_space) {
        Children.resize(2);
        Children[0] = cc().CopyWithT(0.5 * (cc().t() + cc.min()));
        Children[1] = cc().CopyWithT(0.5 * (cc().t() + cc.max()));
      }
      if (!in_time && in_space) {
        Children.resize(cc.Children());
        for (int i = 0; i < cc.Children(); ++i) {
          Children[i] = cc.Child(i).CopyWithT(cc().t());
        }
      }
      if (in_time && in_space) {
        Children.resize(cc.Children() * 2);
        for (int i = 0; i < cc.Children(); ++i) {
          Children[i] = cc.Child(i).
              CopyWithT(0.5 * (cc().t() + cc.min()));
          Children[i + cc.Children()] = cc.Child(i).
              CopyWithT(0.5 * (cc().t() + cc.max()));
        }
      }
      if (!in_time && !in_space) {
        Children.resize(1);
        Children[0] = cc();
      }

      vector<row> row_f(Children.size());
      vector<DegreePair> deg_f(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_f[i] = F.find_row(Children[i]);
        deg_f[i] = F.GetDoF().GetDegree(Children[i]);
      }

      int tnc = rc.n() / deg_c.time;
      int idx_c = 0;
      for (int tc = 1; tc <= deg_c.time; ++tc) {
        for (int k = 0; k < k_max; ++k) {
          for (int sc = 0; sc < z[deg_c.space][k].size(); ++sc) {
            const Point z_space_c = z[deg_c.space][k][sc];
            Point z_time = tc * 1.0 / deg_c.time;
            vector<double> value_child(Children.size());
            vector<bool> z_in_child(Children.size());
            for (int i = 0; i < Children.size(); ++i) {
              cell cf = F.find_cell(Children[i]);
              int tnf_i = row_f[i].n() / deg_f[i].time;
              value_child[i] = 0.0;
              z_in_child[i] = true;
              Point z_space(z_space_c);
              if (in_space) {
                z_in_child[i] = transformPointLocaltoLocal(z_space_c, cc, z_space, cf);
              }
              if (in_time) {
                if (cf().t() < cc().t()) {
                  z_time = tc * 1.0 / deg_c.time * 2;
                  if (z_time[0] > 1) z_in_child[i] = false;
                } else if (cf().t() > cc().t()) {
                  z_time = 1 - (1 - tc * 1.0 / deg_c.time) * 2;
                  if (z_time[0] < 0) z_in_child[i] = false;
                }
              }
              if (z_in_child[i]) {
                if (cf.min() == 0.0 && !testspace) {
                } else if (!testspace) {
                  face ff = F.find_face(cf.Face(cf.Faces() - 2));
                  cell cf_prev = F.find_cell_or_overlap_cell(ff.Left());
                  row rf_prev = F.find_row(cf_prev());
                  DegreePair deg_f_prev = F.GetDoF().GetDegree(*cf_prev);
                  int shift = 0;
                  for (int kk = 0; kk < k; kk++)
                    shift += z[deg_f_prev.space][kk].size();

                  for (int sf = 0; sf < z[deg_f_prev.space][k].size(); ++sf) {
                    int dofs = DofToDeg[z[deg_f_prev.space][k].size()];
                    double value_space_shape = (*space_Shape[dofs])(z_space, sf);
                    double value_time_shape = (*time_Shape[deg_f[i].time])(z_time, 0);
                    int index = rf_prev.n() / deg_f_prev.time + shift + sf;
                    value_child[i] += F(rf_prev)[index] * value_space_shape * value_time_shape;
                  }
                }

                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[deg_f[i].space][kk].size();
                for (int sf = 0; sf < z[deg_f[i].space][k].size(); ++sf) {
                  for (int tf = 1; tf <= deg_f[i].time; ++tf) {
                    int dofs = DofToDeg[z[deg_f[i].space][k].size()];
                    double value_space_shape = (*space_Shape[dofs])(z_space, sf);
                    double value_time_shape = (*time_Shape[deg_f[i].time])(z_time, tf);
                    if (testspace) {
                      if (deg_f[i].time == 1) {
                        value_time_shape = 1.0;
                      } else {
                        value_time_shape = (*time_Shape[deg_f[i].time - 1])(z_time, tf - 1);
                      }
                    }
                    int index = (tf - 1) * tnf_i + shift + sf;
                    value_child[i] += F(row_f[i])[index] * value_space_shape * value_time_shape;
                  }
                }
              }
            }
            int cnt = 0;
            double value = 0.0;
            for (int i = 0; i < Children.size(); ++i) {
              value += value_child[i];
              if (z_in_child[i]) cnt++;
            }
            C(rc)[idx_c++] = value / cnt;
          }
        }
      }
    }
  } else if (s == "Pro") {
    cell cf_test = F.cells();
    vector<Point> Children_in_time_test(2);
    vector<Point> Children_in_space_test(cf_test.Children());
    Children_in_time_test[0] = cf_test().CopyWithT(0.5 * (cf_test().t() + cf_test.min()));
    Children_in_time_test[1] = cf_test().CopyWithT(0.5 * (cf_test().t() + cf_test.max()));
    for (int i = 0; i < cf_test.Children(); ++i) {
      Children_in_space_test[i] = cf_test.Child(i).CopyWithT(cf_test().t());
    }
    bool in_time = true;
    bool in_space = true;
    for (int i = 0; i < Children_in_time_test.size(); ++i) {
      cell cc_test = C.find_cell(Children_in_time_test[i]);
      if (cc_test == C.cells_end()) in_time = false;
    }
    for (int i = 0; i < Children_in_space_test.size(); ++i) {
      cell cc_test = C.find_cell(Children_in_space_test[i]);
      if (cc_test == C.cells_end()) in_space = false;
    }
    if (in_time == in_space) {
      vector<Point> Children(2 * cf_test.Children());
      for (int i = 0; i < cf_test.Children(); ++i) {
        double t_min = 0.5 * (cf_test().t() + cf_test.min());
        double t_max = 0.5 * (cf_test().t() + cf_test.max());
        Children[i] = cf_test.Child(i).CopyWithT(t_min);
        Children[i + cf_test.Children()] = cf_test.Child(i).CopyWithT(t_max);
      }
      bool in_st = true;
      for (int i = 0; i < Children.size(); ++i) {
        cell cc_test = C.find_cell(Children[i]);
        if (cc_test == C.cells_end()) in_st = false;
      }
      if (in_st) {
        in_space = true;
        in_time = true;
      } else {
        if (cf_test == C.cells()) mout << "Interpolate on same level" << endl;
        else Exit("Error in SpaceTimeTools.C\n");
      }
    }

    for (cell cf = F.cells(); cf != F.cells_end(); ++cf) {
      row rf = F.find_row(cf());
      DegreePair deg_f = C.GetDoF().GetDegree(*cf);
      // find prev coarse cell in time
      face ff = F.find_face(cf.Face(cf.Faces() - 2));
      cell c_prev = F.find_cell_or_overlap_cell(ff.Left());
      row rf_prev = F.find_row(c_prev());

      // get corresponding polynomial degrees
      DegreePair deg_f_prev = C.GetDoF().GetDegree(*c_prev);

      vector<Point> Children;
      if (in_time == true && in_space == false) {
        Children.resize(2);
        Children[0] = cf().CopyWithT(0.5 * (cf().t() + cf.min()));
        Children[1] = cf().CopyWithT(0.5 * (cf().t() + cf.max()));
      }
      if (in_time == false && in_space == true) {
        Children.resize(cf.Children());
        for (int i = 0; i < cf.Children(); ++i) {
          Children[i] = cf.Child(i).CopyWithT(cf().t());
        }
      }
      if (in_time == true && in_space == true) {
        Children.resize(cf.Children() * 2);
        for (int i = 0; i < cf.Children(); ++i) {
          Children[i] = cf.Child(i).
              CopyWithT(0.5 * (cf().t() + cf.min()));
          Children[i + cf.Children()] = cf.Child(i).
              CopyWithT(0.5 * (cf().t() + cf.max()));
        }
      }
      if (in_time == false && in_space == false) {
        Children.resize(1);
        Children[0] = cf();
      }

      vector<row> row_c(Children.size());
      vector<DegreePair> deg_c(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_c[i] = C.find_row(Children[i]);
        deg_c[i] = F.GetDoF().GetDegree(Children[i]);
      }
      int tnf = rf.n() / deg_f.time;
      for (int i = 0; i < Children.size(); ++i) {
        cell cc = C.find_cell(Children[i]);
        int tnc_i = row_c[i].n() / deg_c[i].time;
        int idx_c = 0;
        for (int tc = 1; tc <= deg_c[i].time; ++tc) { //tc=0 wird in der entspr. Zelle interpoliert
          for (int k = 0; k < k_max; ++k) {
            for (int sc = 0; sc < z[deg_c[i].space][k].size(); ++sc) {
              Point z_space(z[deg_c[i].space][k][sc]);
              Point z_time = tc * 1.0 / deg_c[i].time;
              if (in_space)
                transformPointLocaltoLocal(z[deg_c[i].space][k][sc], cc, z_space, cf);
              if (in_time) {
                if (in_space) {
                  if (2 * i >= Children.size())
                    z_time = tc * 0.5 / (deg_c[i].time) + i * 0.5;
                } else {
                  z_time = tc * 0.5 / (deg_c[i].time) + i * 0.5;
                }
              }

              if (testspace == false) {
                if (cf.min() == 0.0) {
                } else {
                  int shift = 0;
                  for (int kk = 0; kk < k; kk++)
                    shift += z[deg_f_prev.space][kk].size();
                  for (int sf = 0; sf < z[deg_f_prev.space][k].size(); ++sf) {
                    const Shape *SS = space_Shape[DofToDeg[z[deg_f_prev.space][k].size()]];
                    double value_space_shape = (*SS)(z_space, sf);
                    double value_time_shape = (*time_Shape[deg_f.time])(z_time, 0);
                    double v = F(rf_prev)[rf_prev.n() / deg_f_prev.time + shift + sf];
                    C(row_c[i])[idx_c] += v * value_space_shape * value_time_shape;
                  }
                }
              }
              for (int tf = 1; tf <= deg_f.time; ++tf) {
                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[deg_f.space][kk].size();
                for (int sf = 0; sf < z[deg_f.space][k].size(); ++sf) {
                  int deg = DofToDeg[z[deg_f.space][k].size()];
                  double value_space_shape = (*space_Shape[deg])(z_space, sf);
                  double value_time_shape = (*time_Shape[deg_f.time])(z_time, tf);
                  if (testspace) {
                    if (deg_f.time == 1) {
                      value_time_shape = 1.0;
                    } else {
                      value_time_shape = (*time_Shape[deg_f.time - 1])(z_time, tf - 1);
                    }
                  }
                  double v = F(rf)[(tf - 1) * tnf + shift + sf];
                  C(row_c[i])[idx_c] += v * value_space_shape * value_time_shape;
                }
              }
              idx_c++;
            }
          }
        }
      }
    }
  }

  if (s == "Pro") {
    C.MakeAdditive();
    C.Accumulate();
  }
};

void NewCellTransfer_2D::L2Projection(const Vector &F,
                                      Vector &C,
                                      string s,
                                      bool testspace = false) const {
  C = 0.;
  if (s == "Res") {
    cell cc_test = C.cells();
    vector<Point> Children_in_time_test(2);
    vector<Point> Children_in_space_test(cc_test.Children());
    Children_in_time_test[0] = cc_test().CopyWithT(0.5 * (cc_test().t() + cc_test.min()));
    Children_in_time_test[1] = cc_test().CopyWithT(0.5 * (cc_test().t() + cc_test.max()));
    for (int i = 0; i < cc_test.Children(); ++i) {
      Children_in_space_test[i] = cc_test.Child(i).CopyWithT(cc_test().t());
    }
    bool in_time = true;
    bool in_space = true;
    for (int i = 0; i < Children_in_time_test.size(); ++i) {
      cell cf_test = F.find_cell(Children_in_time_test[i]);
      if (cf_test == F.cells_end()) in_time = false;
    }
    for (int i = 0; i < Children_in_space_test.size(); ++i) {
      cell cf_test = F.find_cell(Children_in_space_test[i]);
      if (cf_test == F.cells_end()) in_space = false;
    }
    if (in_time == in_space) {
      Exit("Error in SpaceTimePreconditioner.C\n");
    }

    for (cell cc = C.cells(); cc != C.cells_end(); ++cc) {
      row rc = C.find_row(cc());
      DegreePair deg_c = C.GetDoF().GetDegree(*cc);

      vector<Point> Children;
      if (in_time) {
        Children.resize(2);
        Children[0] = cc().CopyWithT(0.5 * (cc().t() + cc.min()));
        Children[1] = cc().CopyWithT(0.5 * (cc().t() + cc.max()));
      }
      if (in_space) {
        Children.resize(cc.Children());
        for (int i = 0; i < cc.Children(); ++i) {
          Children[i] = cc.Child(i).CopyWithT(cc().t());
        }
      }

      vector<row> row_f(Children.size());
      vector<int> space_deg_f(Children.size());
      vector<int> time_deg_f(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_f[i] = F.find_row(Children[i]);
        space_deg_f[i] = F.GetDoF().get_space_deg(Children[i]);
        time_deg_f[i] = F.GetDoF().get_time_deg(Children[i]);
      }

      for (int ci = 0; ci < Children.size(); ++ci) {

        cell cf = F.find_cell(Children[ci]);
        row rf = F.find_row(cf());
        int space_deg_f = F.GetDoF().get_space_deg(cf());
        int time_deg_f = F.GetDoF().get_time_deg(cf());

        int tnf = rf.n() / time_deg_f;
        int probDim = tnf / z[space_deg_f].size();
        vector<vector<std::unique_ptr<RVector>>> Space(time_deg_f);

        for (int tf = 0; tf < time_deg_f; ++tf) {
          Space[tf].resize(probDim);
          for (int p = 0; p < probDim; ++p) {
            RVector tmp(0, z[space_deg_f].size());
            for (int sf = 0; sf < z[space_deg_f].size(); ++sf) {
              tmp[sf] = F(rf)[tf * tnf + p * z[space_deg_f].size() + sf];
            }

            RVector tmp2(IntMatrixSpace[deg_c.space][space_deg_f][ci].multiplyWith(tmp));
            Space[tf][p] = std::make_unique<RVector>(l2SpaceMatrix[deg_c.space].multiplyWith(tmp2));
          }
        }

        if (testspace && in_space) {
          for (int p = 0; p < probDim; ++p) {
            for (int sc = 0; sc < z[deg_c.space].size(); ++sc) {

              RVector tmp(0.0, time_deg_f);
              for (int tf = 0; tf < time_deg_f; ++tf) {
                tmp[tf] = (*Space[tf][p])[sc];
              }

              RVector tmp2(IntMatrixTime[deg_c.time - 1][time_deg_f - 1][0].multiplyWith(tmp));
              RVector tmp3(l2TimeMatrix[deg_c.time - 1].multiplyWith(tmp2));

              for (int i = 0; i < tmp3.size(); ++i) {
                C(rc)[i * rc.n() / deg_c.time + p * z[deg_c.space].size() + sc] += tmp3[i];
              }
            }
          }
        }

        if (testspace && in_time) {
          for (int p = 0; p < probDim; ++p) {
            for (int sc = 0; sc < z[deg_c.space].size(); ++sc) {

              RVector tmp(0.0, time_deg_f);
              for (int tf = 0; tf < time_deg_f; ++tf) {
                tmp[tf] = (*Space[tf][p])[sc];
              }

              RVector tmp2(IntMatrixTime[deg_c.time - 1][time_deg_f - 1][ci].multiplyWith(tmp));
              RVector tmp3(l2TimeMatrix[deg_c.time - 1].multiplyWith(tmp2));

              for (int i = 0; i < tmp3.size(); ++i) {
                C(rc)[i * rc.n() / deg_c.time + p * z[deg_c.space].size() + sc] += tmp3[i];
              }
            }
          }
        }
      }

      if (testspace && in_time) continue;
      if (testspace && in_space) continue;

      int tnc = rc.n() / deg_c.time;
      int idx_c = 0;

      for (int tc = 1; tc <= deg_c.time; ++tc) {
        for (int k = 0; k < k_max; ++k) {
          for (int sc = 0; sc < z[deg_c.space].size(); ++sc) {
            const Point z_space_c = z[deg_c.space][k][sc];
            Point z_time = tc * 1.0 / deg_c.time;
            vector<double> value_child(Children.size());
            vector<bool> z_in_child(Children.size());

            for (int i = 0; i < Children.size(); ++i) {
              cell cf = F.find_cell(Children[i]);
              int tnf_i = row_f[i].n() / time_deg_f[i];

              value_child[i] = 0.0;
              z_in_child[i] = true;
              Point z_space;
              if (in_space) {
                z_in_child[i] = transformPointLocaltoLocal(z_space_c, cc, z_space, cf);
              }
              if (cf().t() < cc().t()) {
                z_time = tc * 1.0 / deg_c.time * 2;
                if (z_time[0] > 1) z_in_child[i] = false;
              } else if (cf().t() > cc().t()) {
                z_time = 1 - (1 - tc * 1.0 / deg_c.time) * 2;
                if (z_time[0] < 0) z_in_child[i] = false;
              }
              if (z_in_child[i]) {
                if (cf.min() == 0.0 && testspace == false) {
                } else if (testspace == false) {
                  face ff = F.find_face(cf.Face(cf.Faces() - 2));
                  cell cf_prev = F.find_cell_or_overlap_cell(ff.Left());
                  row rf_prev = F.find_row(cf_prev());
                  DegreePair deg_f_prev = F.GetDoF().GetDegree(*cf_prev);
                  int shift = 0;
                  if (k < k_max - cc.dim()) {
                    shift = k * z[deg_f_prev.space].size();
                  } else {
                    shift = (k_max - cc.dim()) * z[deg_f_prev.space].size()
                            + (k + cc.dim() - k_max) * z[deg_f_prev.space].size();
                  }

                  int time_deg_f_prev = F.GetDoF().get_time_deg(cf_prev());
                  for (int sf = 0; sf < z[deg_f_prev.space].size(); ++sf) {
                    double value_space_shape = (*space_Shape[deg_f_prev.space])(z_space, sf);
                    double value_time_shape = (*time_Shape[time_deg_f[i]])(z_time, 0);
                    value_child[i] += F(rf_prev)[rf_prev.n() / time_deg_f_prev + shift + sf]
                                      * value_space_shape * value_time_shape;
                  }
                }

                int shift = 0;
                DegreePair deg_cf = F.GetDoF().GetDegree(*cf);
                space_deg_f[i] = deg_cf.space;
                if (k < k_max - cc.dim()) {
                  shift = k * z[space_deg_f[i]].size();
                } else {
                  shift = (k_max - cc.dim()) * z[deg_cf.space].size()
                          + (k + cc.dim() - k_max) * z[space_deg_f[i]].size();
                }

                for (int sf = 0; sf < z[space_deg_f[i]].size(); ++sf) {
                  for (int tf = 1; tf <= time_deg_f[i]; ++tf) {
                    double value_space_shape = (*space_Shape[space_deg_f[i]])(z_space, sf);
                    double value_time_shape = (*time_Shape[time_deg_f[i]])(z_time, tf);
                    if (testspace) {
                      if (time_deg_f[i] == 1) {
                        value_time_shape = 1.0;
                      } else {
                        value_time_shape = (*time_Shape[time_deg_f[i] - 1])(z_time, tf - 1);
                      }
                    }
                    value_child[i] += F(row_f[i])[(tf - 1) * tnf_i + shift + sf]
                                      * value_space_shape * value_time_shape;
                  }
                }
              }
            }
            int cnt = 0;
            double value = 0.0;
            for (int i = 0; i < Children.size(); ++i) {
              value += value_child[i];
              if (z_in_child[i]) cnt++;
            }

            C(rc)[idx_c++] = value / cnt;
          }
        }
      }
    }
  } else if (s == "Pro") Exit("not here");
};

void NewCellTransfer_2D::InterpolateDgInTime(const Vector &F, Vector &C, string s) const {
  C = 0.;
  if (s == "Res") {
    cell cc_test = C.cells();

    bool in_time = true;
    Point Child_in_time_test = cc_test().CopyWithT(0.5 * (cc_test().t() + cc_test.min()));
    cell cf_test = F.find_cell(Child_in_time_test);
    if (cf_test == F.cells_end()) in_time = false;
    Point Child_in_space_test = cc_test.Child(0).CopyWithT(cc_test().t());
    bool in_space = true;
    cf_test = F.find_cell(Child_in_space_test);
    if (cf_test == F.cells_end()) in_space = false;
    if (in_time == in_space) {
      Point Child_in_spaceTime_test = cc_test.Child(0).
          CopyWithT(0.5 * (cc_test().t() + cc_test.min()));
      bool in_st = true;
      cf_test = C.find_cell(Child_in_spaceTime_test);
      if (cf_test == C.cells_end()) in_st = false;
      if (in_st) {
        in_space = true;
        in_time = true;
      } else {
        if (cc_test == F.cells()) mout << "Interpolate on same level" << endl;
        else Exit("Error in SpaceTimeTools.C\n");
      }
    }

    for (cell cc = C.cells(); cc != C.cells_end(); ++cc) {
      row rc = C.find_row(cc());
      int space_deg_c = C.GetDoF().get_space_deg(*cc);
      int time_deg_c = C.GetDoF().get_time_deg(cc());

      vector<Point> Children;
      if (in_time == true && in_space == false) {
        Children.resize(2);
        Children[0] = cc().CopyWithT(0.5 * (cc().t() + cc.min()));
        Children[1] = cc().CopyWithT(0.5 * (cc().t() + cc.max()));
      } else if (in_time == false && in_space == true) {
        Children.resize(cc.Children());
        for (int i = 0; i < cc.Children(); ++i) {
          Children[i] = cc.Child(i).CopyWithT(cc().t());
        }
      } else if (in_time == true && in_space == true) {
        Children.resize(cc.Children() * 2);
        for (int i = 0; i < cc.Children(); ++i) {
          Children[i] = cc.Child(i).CopyWithT(0.5 * (cc().t() + cc.min()));
          Children[i + cc.Children()] =
              cc.Child(i).CopyWithT(0.5 * (cc().t() + cc.max()));
        }
      } else if (in_time == false && in_space == false) {
        Children.resize(1);
        Children[0] = cc();
      }

      vector<row> row_f(Children.size());
      vector<int> space_deg_f(Children.size());
      vector<int> time_deg_f(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_f[i] = F.find_row(Children[i]);
        space_deg_f[i] = F.GetDoF().get_space_deg(Children[i]);
        time_deg_f[i] = F.GetDoF().get_time_deg(Children[i]);
      }
      int idx_c = 0;
      for (int tc = 0; tc <= time_deg_c; ++tc) {
        for (int k = 0; k < k_max; ++k) {
          for (int sc = 0; sc < z[space_deg_c][k].size(); ++sc) {
            Point local_c;
            if (time_deg_c > 0) {
              local_c = z[space_deg_c][k][sc].
                  CopyWithT(tc * 1.0 / time_deg_c);
            } else {
              local_c = z[space_deg_c][k][sc].
                  CopyWithT(0.5);
            }
            vector<double> value_child(Children.size());
            vector<bool> z_in_child(Children.size());
            for (int i = 0; i < Children.size(); ++i) {
              cell cf = F.find_cell(Children[i]);
              value_child[i] = 0.0;
              Point local_f(local_c);

              z_in_child[i] = transformPointLocaltoLocal(local_c, cc, local_f, cf);

              if (z_in_child[i]) {
                int tnf_i = row_f[i].n() / (time_deg_f[i] + 1);

                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[space_deg_f[i]][kk].size();
                for (int sf = 0; sf < z[space_deg_f[i]][k].size(); ++sf) {
                  for (int tf = 0; tf <= time_deg_f[i]; ++tf) {
                    auto space_shape = space_Shape[DofToDeg[z[space_deg_f[i]][k].size()]];
                    double value_space_shape = (*space_shape)(local_f, sf);
                    auto time_shape = time_Shape[time_deg_f[i]];
                    double value_time_shape = (*time_shape)(Point(local_f.t()), tf);
                    value_child[i] += F(row_f[i])[tf * tnf_i + shift + sf]
                                      * value_space_shape * value_time_shape;
                  }
                }
              }
            }
            int cnt = 0;
            double value = 0.0;
            for (int i = 0; i < Children.size(); ++i) {
              value += value_child[i];
              if (z_in_child[i]) cnt++;
            }
            C(rc)[idx_c++] = value / cnt;
          }
        }
      }
    }
  } else if (s == "Pro") {
    cell cf_test = F.cells();

    vector<Point> Children_in_space_test(cf_test.Children());
    Point Child_in_time_test = cf_test().
        CopyWithT(0.5 * (cf_test().t() + cf_test.min()));
    bool in_time = true;
    cell cc_test = C.find_cell(Child_in_time_test);
    if (cc_test == C.cells_end()) in_time = false;
    Point Child_in_space_test = cf_test.Child(0).CopyWithT(cf_test().t());
    bool in_space = true;
    cc_test = C.find_cell(Child_in_space_test);
    if (cc_test == C.cells_end()) in_space = false;
    if (in_time == in_space) {
      Point Child_in_spaceTime_test = cf_test.Child(0).
          CopyWithT(0.5 * (cf_test().t() + cf_test.min()));
      bool in_st = true;
      cc_test = C.find_cell(Child_in_spaceTime_test);
      if (cc_test == C.cells_end()) in_st = false;
      if (in_st) {
        in_space = true;
        in_time = true;
      } else {
        if (cf_test == C.cells()) mout << "Interpolate on same level" << endl;
        else Exit("Error in SpaceTimePreconditioner.C\n");
      }
    }

    for (cell cf = F.cells(); cf != F.cells_end(); ++cf) {
      row rf = F.find_row(cf());
      int time_deg_f = C.GetDoF().get_time_deg(*cf);
      int space_deg_f = C.GetDoF().get_space_deg(*cf);

      vector<Point> Children;
      if (in_time == true && in_space == false) {
        Children.resize(2);
        Children[0] = cf().CopyWithT(0.5 * (cf().t() + cf.min()));
        Children[1] = cf().CopyWithT(0.5 * (cf().t() + cf.max()));
      } else if (in_time == false && in_space == true) {
        Children.resize(cf.Children());
        for (int i = 0; i < cf.Children(); ++i) {
          Children[i] = cf.Child(i).CopyWithT(cf().t());
        }
      } else if (in_time == true && in_space == true) {
        Children.resize(cf.Children() * 2);
        for (int i = 0; i < cf.Children(); ++i) {
          Children[i] = cf.Child(i).
              CopyWithT(0.5 * (cf().t() + cf.min()));
          Children[i + cf.Children()] = cf.Child(i).
              CopyWithT(0.5 * (cf().t() + cf.max()));
        }
      } else if (in_time == false && in_space == false) {
        Children.resize(1);
        Children[0] = cf();
      }

      vector<row> row_c(Children.size());
      vector<int> time_deg_c(Children.size());
      vector<int> space_deg_c(Children.size());
      for (int i = 0; i < Children.size(); ++i) {
        row_c[i] = C.find_row(Children[i]);
        time_deg_c[i] = F.GetDoF().get_time_deg(Children[i]);
        space_deg_c[i] = F.GetDoF().get_space_deg(Children[i]);
      }

      int tnf = rf.n() / (time_deg_f + 1);
      for (int i = 0; i < Children.size(); ++i) {
        cell cc = C.find_cell(Children[i]);
        int tnc_i = row_c[i].n() / (time_deg_c[i] + 1);
        int idx_c = 0;
        for (int tc = 0; tc <= time_deg_c[i]; ++tc) {
          for (int k = 0; k < k_max; ++k) {
            for (int sc = 0; sc < z[space_deg_c[i]][k].size(); ++sc) {
              Point local_c;
              if (time_deg_c[i] > 0) {
                local_c = z[space_deg_c[i]][k][sc].
                    CopyWithT(tc * 1.0 / time_deg_c[i]);
              } else {
                local_c = z[space_deg_c[i]][k][sc].CopyWithT(0.5);
              }

              Point local_f(local_c);
              transformPointLocaltoLocal(local_c, cc, local_f, cf);

              for (int tf = 0; tf <= time_deg_f; ++tf) {
                int shift = 0;
                for (int kk = 0; kk < k; kk++)
                  shift += z[space_deg_f][kk].size();
                for (int sf = 0; sf < z[space_deg_f][k].size(); ++sf) {
                  int deg = DofToDeg[z[space_deg_f][k].size()];
                  double value_space_shape = (*space_Shape[deg])(local_f, sf);
                  double
                      value_time_shape = (*time_Shape[time_deg_f])(Point(local_f.t()), tf);
                  double v = F(rf)[tf * tnf + shift + sf];
                  C(row_c[i])[idx_c] += v * value_space_shape * value_time_shape;
                }
              }
              idx_c++;
            }
          }
        }
      }
    }
  }

  if (s == "Pro") {
    C.MakeAdditive();
    C.Accumulate();
  }
};

void NewCellTransfer_2D::L2ProjectionDgInTime(const Vector &F, Vector &C) const {

  TransferInfo info(F, C);

  unordered_map<DegreePair, RMatrix> InvertedMassMatrixCache;
  unordered_map<DegreePair, SpaceTimeCellQuadrature> SpaceTimeQuadCache;

  for (cell cc = C.cells(); cc != C.cells_end(); cc++) {
    row rc = C.find_row(cc());
    DegreePair deg_c = C.GetDoF().GetDegree(*cc);
    auto &space_shape = disc.GetCellShape(deg_c.space);
    auto &time_shape = disc.GetTimeShape(deg_c.time);

    int N = time_shape.size() * k_max * space_shape.size();
    RMatrix InvertedMassMatrix(N, N);
    auto elem = InvertedMassMatrixCache.find(deg_c);

    if (elem == InvertedMassMatrixCache.end()) {
      Quadrature s_quad = disc.GetCellQuad(deg_c.space);
      Quadrature t_quad = disc.GetTimeQuad(deg_c.time);
      auto quad_c = SpaceTimeCellQuadrature(cc, s_quad, t_quad);

      for (int q = 0; q < quad_c.size(); q++) {
        Point local_c = quad_c.qPointLocal[q];
        for (int tc = 0; tc < time_shape.size(); ++tc) {
          double value_time_shape = time_shape(Point(local_c.t()), tc);
          for (int k = 0; k < k_max; ++k) {
            for (int sc = 0; sc < space_shape.size(); ++sc) {
              double value_space_shape = space_shape(local_c, sc);
              double v1 = value_space_shape * value_time_shape;
              int row = sc + (k + tc * k_max) * space_shape.size();
              for (int tc2 = 0; tc2 < time_shape.size(); ++tc2) {
                double value_time_shape2 = time_shape(Point(local_c.t()), tc2);
                for (int sc2 = 0; sc2 < space_shape.size(); ++sc2) {
                  int col = sc2 + (k + tc2 * k_max) * space_shape.size();
                  double value_space_shape2 = space_shape(local_c, sc2);
                  double v2 = value_space_shape2 * value_time_shape2;
                  InvertedMassMatrix[row][col] += quad_c.qWeight[q] * v1 * v2;
                }
              }
            }
          }
        }
      }
      InvertedMassMatrix.Invert();
      elem = InvertedMassMatrixCache.insert({deg_c, InvertedMassMatrix}).first;
    }

    vector<Point> Children(info.GetChildren(cc));
    vector<row> row_f(Children.size());
    vector<DegreePair> deg_f(Children.size());
    vector<cell> child_cells(Children.size());
    vector<SpaceTimeCellQuadrature *> quad_f;
    for (int i = 0; i < Children.size(); ++i) {
      row_f[i] = F.find_row(Children[i]);
      deg_f[i] = F.GetDoF().GetDegree(Children[i]);
      child_cells[i] = F.find_cell(Children[i]);
      auto elem_q = SpaceTimeQuadCache.find(deg_f[i]);
      if (elem_q == SpaceTimeQuadCache.end()) {
        const Quadrature &s_quad = disc.GetCellQuad(deg_f[i].space);
        const Quadrature &t_quad = disc.GetTimeQuad(deg_f[i].time);
        elem_q = SpaceTimeQuadCache.insert(
            {deg_f[i], SpaceTimeCellQuadrature(child_cells[i], s_quad, t_quad)}).first;
      }
      quad_f.push_back(&(elem_q->second));
    }
    RVector fine_RHS(N);
    for (int i = 0; i < Children.size(); ++i) {
      auto &space_shape_f = disc.GetCellShape(deg_f[i].space);
      auto &time_shape_f = disc.GetTimeShape(deg_f[i].time);
      for (int q = 0; q < quad_f[i]->size(); q++) {
        Point local_f = quad_f[i]->qPointLocal[q];
        Point local_c;
        transformPointLocaltoLocal(local_f, child_cells[i], local_c, cc);
        for (int tc = 0; tc < time_shape.size(); ++tc) {
          double value_time_shape = time_shape(Point(local_c.t()), tc);
          for (int sc = 0; sc < space_shape.size(); ++sc) {
            double value_space_shape = space_shape(local_c, sc);
            double v_coarse = value_space_shape * value_time_shape;
            for (int tf = 0; tf < time_shape_f.size(); ++tf) {
              double value_time_shape_f = time_shape_f(Point(local_f.t()), tf);
              for (int sf = 0; sf < space_shape_f.size(); ++sf) {
                double value_space_shape_f = space_shape_f(local_f, sf);
                for (int k = 0; k < k_max; ++k) {
                  int row = sc + (k + tc * k_max) * space_shape.size();
                  int f_idx = sf + (k + tf * k_max) * space_shape_f.size();
                  double v_fine = value_space_shape_f * value_time_shape_f;
                  fine_RHS[row] += quad_f[i]->qWeight[q] * v_coarse * v_fine * F(row_f[i])[f_idx];
                }
              }
            }
          }
        }
      }
    }
    InvertedMassMatrix = elem->second;
    RVector coarse_solution(InvertedMassMatrix.multiplyWith(fine_RHS));

    for (int i = 0; i < coarse_solution.size(); i++) {
      C(rc)[i] = coarse_solution[i];
    }
  }
}

void NewCellTransfer_2D::L2ProlongationDgInTime(const Vector &F, Vector &C) const {
  TransferInfo info(C, F);
  //TODO
}

void fillSpaceTimeVector(Vector &V, int pdim,
                         const std::function<double(const Point &,
                                                    const cell &,
                                                    int)> &func) {
  for (cell tc = V.cells(); tc != V.cells_end(); tc++) {
    vector<Point> nodalp = V.GetDoF().GetNodalPoints(*tc);
    DegreePair deg = V.GetDoF().GetDegree(*tc);

    int cnt = 0;
    for (int t = 0; t <= deg.time; t++) {
      for (int pi = 0; pi < pdim; ++pi) {
        for (int s = 0; s < V.GetDoF().NodalPointsLocal(deg.space); ++s) {
          V(tc())[cnt] = func(nodalp[cnt], tc, pi);
          cnt++;
        }
      }
    }
  }
}

