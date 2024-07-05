#ifndef _SPACETIMETOOLS_H_
#define _SPACETIMETOOLS_H_

#include "m++.hpp"
#include "SpaceTimeDiscretization.hpp"
#include "Meshes.hpp"
#include "Parallel.hpp"
#include "Algebra.hpp"
#include "RMatrix.hpp"
#include "Shapes.hpp"
#include "RMatrix.hpp"
#include "LagrangeShapesInterval.hpp"
#include <memory>
#include <iosfwd>
#include <iomanip>
#include <fstream>
#include <functional>

// Todo: Is this useful for kernel ? is it already implemented there ?
bool transformPointGlobalToLocal(const Point &P_glob, Point &PT, const cell &CT);

// Todo: Is this useful for kernel ? is it already implemented there ?
Point transformPointLocaltoGlobal(const Point PF, const cell &CF);

// Todo: Is this useful for kernel ? is it already implemented there ?
bool transformPointLocaltoLocal(const Point PF,
                                const cell &CF,
                                Point &PT,
                                const cell &CT);

// Redistribute cells for a better load balance
// old STM, disc_vector and dual_disc_vector will be deleted and replaced by a better balanced version
void redistribute_cells(std::shared_ptr<STDiscretization> &disc,
                        LevelPair levels);

void saveAdapt(std::shared_ptr<Meshes> &STM, STDiscretization &disc, LevelPair level, int run);

void restoreAdapt(std::shared_ptr<Meshes> &STM, STDiscretization &disc, LevelPair level, int &run);

void increase_poly_degs_uniform(STDiscretization &discretization);

void increase_poly_degs(STDiscretization &,
                        const Vector &,
                        const double,
                        const double,
                        const string&);

std::unordered_map<Point, DegreePair>
createPolynomialDistribution(STDiscretization &discretization,
                             const Vector &Eta,
                             const double theta,
                             const double theta_min,
                             const string name);

void increase_poly_degs_by_percentage(STDiscretization &,
                                      const Vector &,
                                      const double);

// Performs a space-time cut at x = x_pos. Data on this slice is written to /data/csv/pattern.csv
void write_pattern(Vector &U, STDiscretization *disc, double x_pos = 6.);

// ----------------------------------------------------
//    InterpolateWave TransferTWave
// ----------------------------------------------------

class InterpolateWave_dG : public Operator {
public:
  InterpolateWave_dG() {}

  virtual void Interpolate(const Vector &f, Vector &c) const {}

  virtual void multiply(Vector &u, const Vector &v) const {}

  virtual void multiply_transpose(Vector &u, const Vector &v) const {}

  virtual int get_deg(int m) const {Exit("Not Implemented"); }

  bool dual;
};

inline constAB<Operator, Vector>
operator*(const InterpolateWave_dG &T, const Vector &v) {
  return constAB<Operator, Vector>(T, v);
}

inline constAB<Vector, Operator>
operator*(const Vector &v, const InterpolateWave_dG &T) {
  return constAB<Vector, Operator>(v, T);
}

struct TransferInfo {
  const bool in_space;
  const bool in_time;

  TransferInfo(const Vector &F, const Vector &C) :
        in_space(F.Level().space > C.Level().space),
        in_time(F.Level().time > C.Level().time){
    if (!in_time && !in_space) {
      mout << "Interpolate on same level" << endl;
    }
  }


  vector<Point> GetChildren(cell cc) {
    vector<Point> Children;
    if (in_time && !in_space) {
      Children.resize(2);
      Children[0] = cc().CopyWithT(0.5 * (cc().t() + cc.min()));
      Children[1] = cc().CopyWithT(0.5 * (cc().t() + cc.max()));
    } else if (!in_time && in_space) {
      Children.resize(cc.SpaceCell().Children());
      for (int i = 0; i < cc.SpaceCell().Children(); ++i) {
        Children[i] = cc.SpaceCell().Child(i).CopyWithT(cc().t());
      }
    } else if (in_time && in_space) {
      Children.resize(cc.SpaceCell().Children() * 2);
      for (int i = 0; i < cc.SpaceCell().Children(); ++i) {
        Children[i] = cc.SpaceCell().Child(i).CopyWithT(0.5 * (cc().t() + cc.min()));
        Children[i + cc.SpaceCell().Children()] = cc.SpaceCell().Child(i).CopyWithT(0.5 * (cc().t() + cc.max()));
      }
    } else if (!in_time && !in_space) {
      Children.resize(1);
      Children[0] = cc();
    }
    return Children;
  }

};


struct SpaceTimeCellQuadrature {
  vector<double> qWeight;
  vector<Point> qPoint;

  vector<Point> qPointLocal;

  SpaceTimeCellQuadrature(const cell tc, const Quadrature &SpaceQ, const Quadrature &TimeQ) :
      qWeight(SpaceQ.size() * TimeQ.size()),
      qPoint(SpaceQ.size() * TimeQ.size()),
      qPointLocal(SpaceQ.size() * TimeQ.size()) {
    for (int tq = 0; tq < TimeQ.size(); ++tq) {
      for (int sq = 0; sq < SpaceQ.size(); ++sq) {
        Transformation T = tc.GetTransformation();
        double det = T.Det();
        qWeight[tq * SpaceQ.size() + sq] = det * TimeQ.Weight(tq) * SpaceQ.Weight(sq);
        qPoint[tq * SpaceQ.size() + sq] = tc.LocalToGlobal(SpaceQ.QPoint(sq).CopyWithT(TimeQ.QPoint(tq)[0]));
        qPointLocal[tq * SpaceQ.size() + sq] = SpaceQ.QPoint(sq).CopyWithT(TimeQ.QPoint(tq)[0]);
      }
    }
  }

  int size() const {
    return qWeight.size();
  }
};


class NewCellTransfer_2D : public InterpolateWave_dG {
private:
  const STDiscretization &disc;
  const int k_max; //problemDim
  vector<const Shape *> space_Shape;
  vector<const Shape *> time_Shape;
  vector<vector<vector<Point>>> z;

  bool useL2Proj;
  vector<vector<vector<RMatrix>>> IntMatrixSpace;
  vector<vector<vector<RMatrix>>> IntMatrixTime;
  vector<RMatrix> l2SpaceMatrix;
  vector<RMatrix> l2TimeMatrix;

  const int maxDegSpace;
  const int maxDegTime;
  vector<int> DofToDeg;
  const bool dgInTime;
public:
  NewCellTransfer_2D(const STDiscretization &disc)
      :
      disc(disc),
      k_max(disc().GetDoF().get_dim()),
      useL2Proj(false),
      maxDegSpace(6), maxDegTime(6),
      dgInTime(disc.isDgInTime()) {
    Config::Get("useL2Projection", useL2Proj);


    space_Shape.resize(maxDegSpace + 1);
    for (int sDeg = 0; sDeg < maxDegSpace + 1; sDeg++) {
      const Shape &SpaceShape(disc.GetCellShape(sDeg));
      space_Shape[sDeg] = &SpaceShape;
    }
    time_Shape.resize(maxDegTime + 1);
    for (int tDeg = 0; tDeg < maxDegTime + 1; tDeg++) {
      const Shape &TimeShape(disc.GetTimeShape(tDeg));
      time_Shape[tDeg] = &TimeShape;
    }

    dual = false;

    z.resize(maxDegSpace + 1);
    for (int s = 0; s < maxDegSpace + 1; ++s) {
      z[s].resize(k_max);
      for (int k = 0; k < k_max; k++) {
        disc().GetDoF().NodalPointsLocal(s, z[s][k], k);
      }
    }

    DofToDeg.resize(disc().GetDoF().NodalPointsLocal(maxDegSpace) + 1);
    for (int s = 0; s < maxDegSpace + 1; ++s) {
      DofToDeg[disc().GetDoF().NodalPointsLocal(s)] = s;
    }
  }

  ~NewCellTransfer_2D() {}

  void Interpolate(const Vector &, Vector &, string, bool) const;

  void L2Projection(const Vector &, Vector &, string, bool) const;

  void InterpolateDgInTime(const Vector &, Vector &, string) const;

  void L2ProjectionDgInTime(const Vector &, Vector &) const;

  void L2ProlongationDgInTime(const Vector &, Vector &) const;

  void multiply(Vector &u, const Vector &v) const {
    bool dgInTime = disc.isDgInTime();
    if (dgInTime) {
      if (false && useL2Proj) {
        L2ProlongationDgInTime(v, u);
      } else {
        InterpolateDgInTime(v, u, "Pro");
      }
      return;
    }
    if (!dual) Interpolate(v, u, "Pro", false); else Interpolate(v, u, "Pro", true);
  }

  void multiply_transpose(Vector &u, const Vector &v) const {
    bool dgInTime = disc.isDgInTime();
    if (dgInTime) {
      if (useL2Proj) {
        L2ProjectionDgInTime(v, u);
      } else {
        InterpolateDgInTime(v, u, "Res");
      }
      return;
    }
    if (!dual) {
      if (useL2Proj) {
        L2Projection(v, u, "Res", true);
        return;
      } else Interpolate(v, u, "Res", true);
    } else Interpolate(v, u, "Res", false);
  }
};

typedef RVector FWIState;

class FWIReceiverTimeseries : public std::vector<std::pair<double, FWIState *> > {
  int state_comps;
public:
  explicit FWIReceiverTimeseries(int state_comps) : state_comps(state_comps) {}

  explicit FWIReceiverTimeseries(const FWIReceiverTimeseries &ts) : state_comps(ts.state_comps) {
    for (const std::pair<double, FWIState *> &e : ts) {
      push_back(std::make_pair(e.first, new FWIState(*e.second)));
    }
  }

  ~FWIReceiverTimeseries() {
    for (std::pair<double, FWIState *> &i : *this)
      delete i.second;
  }

  friend std::ofstream &operator<<(std::ofstream &o, const FWIReceiverTimeseries &ts) {
    o << ts.size() << endl;
    for (auto m = ts.begin(); m != ts.end(); ++m) {
      o << "  " << m->first << '\t';
      const FWIState &state = *m->second;
      for (int i = 0; i < ts.state_comps; ++i)
        o << state[i] << " ";
      o << endl;
    }
    return o;
  }

  friend std::istream &operator>>(std::istream &in, FWIReceiverTimeseries &ts) {
    ts.clear();
    int ts_size = 0;
    in >> ts_size;
    for (int i = 0; i < ts_size; ++i) {
      double t;
      FWIState *state = new FWIState(0.0, ts.state_comps);
      in >> t;
      for (int j = 0; j < ts.state_comps; ++j)
        in >> (*state)[j];
      ts.push_back(std::make_pair(t, state));
    }
    return in;
  }

  void setStateComponents(int comps) { state_comps = comps; }
};

class FWIObservationSpecification : public std::list<Point> {
  double t_start;
  double t_end;
  double t_delta;
  int n_steps;

  void set_n_steps() { n_steps = round((t_end - t_start) / t_delta) + 1; }

public:
  FWIObservationSpecification()
      : t_start(0), t_end(1), t_delta(0.001) { set_n_steps(); }

  FWIObservationSpecification(const FWIObservationSpecification &spec)
      : list<Point>(spec), t_start(spec.t_start), t_end(spec.t_end),
        t_delta(spec.t_delta) { set_n_steps(); }

  FWIObservationSpecification(double start_time, double end_time, double delta_t,
                              const list <Point> &receiver_positions)
      : list<Point>(receiver_positions), t_start(start_time), t_end(end_time),
        t_delta(delta_t) { set_n_steps(); }

  FWIObservationSpecification(const char *fn) { readFromFile(fn); }

  double startTime() const { return t_start; }

  double endTime() const { return t_end; }

  double timeDelta() const { return t_delta; }

  int nSteps() const { return n_steps; }

  double step(int i) const { return t_start + i * t_delta; }

  // -- IO stuff --
  void printInfo(const char *offset = "") const {
    mout << offset << "FWIObservationSpecification:" << endl;
    mout << offset << "  times: [" << t_start << ", " << t_end << "],  delta_t: "
         << t_delta << endl;
    mout << offset << "  receiver count: " << size() << endl;
    mout << offset << "  receiver positions:" << endl;
    for (auto r = begin(); r != end(); ++r)
      mout << offset << "    " << *r << endl;
  }

  void readFromFile(const string &filename) {
    std::ifstream file;
    file.open(filename.c_str());
    file >> *this;
    file.close();
  }

  void writeToFile(const string &filename) {
    if (PPM->Proc(0) != 0) return;
    std::ofstream file;
    file.open(filename.c_str());
    file << *this;
    file.close();
  }

  friend std::ofstream &operator<<(std::ofstream &o, const FWIObservationSpecification &spec) {
    o << spec.t_start << endl;
    o << spec.t_end << endl;
    o << spec.t_delta << endl;
    o << spec.size() << endl;
    for (auto i = spec.begin(); i != spec.end(); ++i) {
      for (int d = 0; d < 3; ++d)
        o << (*i)[d] << "  ";
      o << endl;
    }
    return o;
  }

  friend std::istream &operator>>(std::istream &in, FWIObservationSpecification &spec) {
    spec.clear();
    in >> spec.t_start;
    in >> spec.t_end;
    in >> spec.t_delta;
    int receiver_count = 0;
    in >> receiver_count;
    for (int i = 0; i < receiver_count; ++i) {
      Point p = Origin;
      for (int d = 0; d < 3; ++d)
        in >> p[d];
      spec.push_back(p);
    }
    spec.set_n_steps();
    return in;
  }
};

inline bool comp_timeseries_entry(const std::pair<double, FWIState *> &a,
                                  const std::pair<double, FWIState *> &b) {
  return (a.first < b.first);
}

class FWISeismogram {
protected:
  FWIObservationSpecification spec;
  std::unordered_map<Point, FWIReceiverTimeseries *> timeseries;
  int state_comps;
public:
  FWISeismogram()
      : state_comps(0) {}

  FWISeismogram(const FWIObservationSpecification &spec, int state_comps)
      : spec(spec), state_comps(state_comps) {}

  FWISeismogram(const char *filename, int state_comps)
      : state_comps(state_comps) { readFromFile(filename); }

  FWISeismogram(const FWISeismogram &s)
      : spec(s.spec), state_comps(s.state_comps) {
    for (auto ts = s.ts_begin(); ts != s.ts_end(); ++ts)
      timeseries[ts->first] = new FWIReceiverTimeseries(*ts->second);
  }

  ~FWISeismogram() { ts_clear(); }

  double alpha() const { return spec.timeDelta() * spec.timeDelta(); }

  FWISeismogram &operator=(const FWISeismogram &s) {
    spec = s.spec;
    state_comps = s.state_comps;

    ts_clear();
    for (auto ts = s.ts_begin(); ts != s.ts_end(); ++ts) {
      timeseries[ts->first] = new FWIReceiverTimeseries(*ts->second);
    }
    return *this;
  }

  void setSpecification(const FWIObservationSpecification &spec_new) {
    ts_clear();
    spec = spec_new;
    for (auto r = spec.begin(); r != spec.end(); ++r)
      timeseries[*r] = new FWIReceiverTimeseries(state_comps);
  }

  int
  getStateComponents() const { return state_comps; } // TODO: sollte auch in Datei gespeichert werden
  void setStateComponents(int state_comps_new) {
    state_comps = state_comps_new;
    for (auto ts = ts_begin(); ts != ts_end(); ++ts)
      ts->second->setStateComponents(state_comps_new);
  }

  void addMeasurement(const Point &receiver, double t, const FWIState &measurement) {
    addMeasurement(receiver, t, new FWIState(measurement));
//        if (timeseries.find(receiver) == timeseries.end())
//            timeseries.emplace(std::make_pair(receiver, new FWIReceiverTimeseries(*assemble)));
//        timeseries.at(receiver).push_back(make_pair(t, new FWIState(measurement)));
  }

  void addMeasurement(const Point &receiver, double t, FWIState *measurement) {
    if (timeseries.find(receiver) == timeseries.end())
//            timeseries.emplace(std::make_pair(receiver, new FWIReceiverTimeseries(*assemble)));
      timeseries[receiver] = new FWIReceiverTimeseries(state_comps);
    timeseries[receiver]->push_back(std::make_pair(t, measurement));
  }

//    const FWIReceiverTimeseries& getTimeSeries(const Point& receiver) const { return *timeseries.at(receiver); }
  const FWIObservationSpecification &getSpecification() const { return spec; }

  void collect(bool average = false) {
    ExchangeBuffer E(0);
    if (!PPM->Master(0)) {
      for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts) {
        E.Send(0) << (Point) ts->first;
        E.Send(0) << (int) ts->second->size();
        for (auto i = ts->second->begin(); i != ts->second->end(); ++i) {
          E.Send(0) << (double) i->first;

          const FWIState &state = *i->second;
          for (int k = 0; k < state.size(); ++k)
            E.Send(0) << state[k];
        }
      }
    }
    E.Communicate();

    if (PPM->Master(0)) {
      for (short q = 0; q < PPM->Size(0); ++q) {
        while (E.Receive(q).size() < E.ReceiveSize(q)) {
          Point receiver;
          E.Receive(q) >> receiver;
          int ts_size;
          E.Receive(q) >> ts_size;

          if (timeseries.find(receiver) == timeseries.end())
            timeseries[receiver] = new FWIReceiverTimeseries(state_comps);

          FWIReceiverTimeseries *current_ts = timeseries[receiver];
          for (int i = 0; i < ts_size; ++i) {
            double t;
            E.Receive(q) >> t;

            FWIState *state = new FWIState(0.0, state_comps);
            for (int k = 0; k < state->size(); ++k)
              E.Receive(q) >> (*state)[k];

            current_ts->push_back(std::make_pair(t, state));
          }
        }
      }
      for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts) {
        FWIReceiverTimeseries *current_ts = ts->second;

        vector<std::pair<double, FWIState *> >
            tmp(current_ts->begin(), current_ts->end());
        std::sort(tmp.begin(), tmp.end(), comp_timeseries_entry);
        current_ts->clear();

        for (auto e = tmp.begin(); e != tmp.end(); ++e) {
          auto e1 = e;
          FWIState *sum = e->second;
          ++e1;
          int counter = 1;
          while (e1 != tmp.end() && e1->first == e->first) {
            *sum += *e1->second;
            ++counter;
            ++e1;
          }
          if (average)
            *sum /= counter;

          current_ts->push_back(std::make_pair(e->first, sum));
          e += (counter - 1);
        }
      }
    }
  }

  void operator=(double val) {
    for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts) {
      FWIReceiverTimeseries *current_ts = ts->second;
      for (auto m = current_ts->begin(); m != current_ts->end(); ++m)
        *m->second = val;
    }
  }

  void operator*=(double val) {
    for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts) {
      FWIReceiverTimeseries *current_ts = ts->second;
      for (auto m = current_ts->begin(); m != current_ts->end(); ++m)
        *m->second *= val;
    }
  }

  void setComponent(int comp, double val) {
    for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts) {
      FWIReceiverTimeseries *current_ts = ts->second;
      for (auto m = current_ts->begin(); m != current_ts->end(); ++m)
        (*m->second)[comp] = val;
    }
  }

  std::unordered_map<Point, FWIReceiverTimeseries *>::iterator ts_begin() {
    return timeseries.begin();
  }

  std::unordered_map<Point, FWIReceiverTimeseries *>::iterator ts_end() {
    return timeseries.end();
  }

  void ts_clear() {
    for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts)
      delete ts->second;
    timeseries.clear();
  }

  std::unordered_map<Point,
      FWIReceiverTimeseries *>::const_iterator ts_begin() const { return timeseries.begin(); }

  std::unordered_map<Point,
      FWIReceiverTimeseries *>::const_iterator ts_end() const { return timeseries.end(); }

  // -- IO stuff --
  void printSpecification() const { spec.printInfo(); }

  void printObservation() const {
    mout << "FWIObservationRecording:" << endl;
    spec.printInfo("  ");
    mout << "  measurements:" << endl;
    for (auto ts = timeseries.begin(); ts != timeseries.end(); ++ts) {
      FWIReceiverTimeseries *current_ts = ts->second;
      mout << "    receiver: " << ts->first << endl;
      for (auto m = current_ts->begin(); m != current_ts->end(); ++m)
        mout << "      t: " << std::setw(6) << m->first << "  -  " << *m->second
             << endl;
    }
  }

  void readFromFile(const string &filename) {
    std::ifstream file;
    file.open(filename.c_str());
    file >> *this;
    file.close();
  }

  void writeToFile(const string &filename) const {
    if (PPM->Proc(0) != 0) return;
    std::ofstream file;
    file.open(filename.c_str());
    if (file.fail()) Exit("Could not write seismogram")
    file << *this;
    file.close();
  }

  friend std::ofstream &operator<<(std::ofstream &o, const FWISeismogram &rec) {
    o << rec.spec;
    o << rec.state_comps << endl;
    o << rec.timeseries.size() << endl;
    for (auto ts = rec.timeseries.begin(); ts != rec.timeseries.end(); ++ts) {
      for (int d = 0; d < 3; ++d)
        o << ts->first[d] << " ";
      o << endl;
//            o << ts->first << endl;
      o << *ts->second;
    }
    return o;
  }

  friend std::istream &operator>>(std::istream &in, FWISeismogram &rec) {
    rec.spec.clear();
    in >> rec.spec;
    in >> rec.state_comps;
    int ts_count = 0;
    in >> ts_count;
    for (int i = 0; i < ts_count; ++i) {
      Point p = Origin;
      for (int d = 0; d < 3; ++d)
        in >> p[d];
      auto *current_ts = new FWIReceiverTimeseries(rec.state_comps);
      in >> *current_ts;
      rec.timeseries[p] = current_ts;
    }
    return in;
  }
};

// L1-Integral = 1
class FWIMeasurementKernel_L1 {
  Point c;
  double r;
  double scaleFactor;
  const double to_zero_factor = 7;
public:
  FWIMeasurementKernel_L1(int dim, const Point &center, double radius)
      : c(center), r(radius), scaleFactor(r) {
    int STdim = dim + 1;
    scaleFactor = 1.0 / pow(2 * M_PI * (r * r) / (2 * to_zero_factor), (STdim) / 2.0);
  }

  virtual double oneDbasis(double t) const {
    if (t > 1) return 0.0;
    else return exp(-to_zero_factor * t * t);
  }

  double operator()(const Point &z) const {
    Point d = c - z;
    double delta = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + d.t() * d.t());
    if (delta >= r) return 0;
    delta /= r;
    return oneDbasis(delta) * scaleFactor;
  }
};

inline double squared(double x) { return x * x; }

static bool cellAndSphereOverlap(const cell &C, const Point &P, double radius) {
  Point minCor = Infty.CopyWithT(C.min());
  Point maxCor = (-Infty).CopyWithT(C.max());
  for (int i = 0; i < C.Corners(); ++i) {
    Point corner = C.Corner(i);
    for (int d = 0; d < C.dim(); ++d) {
      minCor[d] = min(minCor[d], corner[d]);
      maxCor[d] = max(maxCor[d], corner[d]);
    }
  }

  double dist_squared = 0.0;
  for (int d = 0; d < C.dim(); ++d) {
    if (P[d] < minCor[d]) dist_squared += squared(P[d] - minCor[d]);
    else if (P[d] > maxCor[d]) dist_squared += squared(P[d] - maxCor[d]);
  }
  if (P.t() < minCor.t()) dist_squared += squared(P.t() - minCor.t());
  else if (P.t() > maxCor.t()) dist_squared += squared(P.t() - maxCor.t());

  return (dist_squared < squared(radius));
}

void fillSpaceTimeVector(Vector &V, int pdim,
                         const std::function<double(const Point &, const cell &, int)> &func);


#endif
