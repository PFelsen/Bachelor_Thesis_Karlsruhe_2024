#include "PrintUtil.hpp"

void PrintScalar(std::ostream &s, const Scalar &a) {
#ifdef NDOUBLE
  s << " " << a;
#else
  char buf[128];
  buf[0] = 0;
  if (a == 0)
    s << "         0";
  else if (isZero(a))
    s << "         0";
  else if ((abs(a) > 1000) || (abs(a) < 0.001)) {
    sprintf(buf, "%9.3e", a);
    if (a > 0) s << " ";
    s << buf;
  } else {
    sprintf(buf, "%10.5f", a);
    s << buf;
  }
#endif
}

std::string format_double(const double v) {
  std::stringstream ss;
  if (v >= 0) {
    ss << " ";
  }
  int precision = 6;
  Config::Get("precision", precision, true);
  ss.precision(precision);
  bool scientific = false;
  Config::Get("scientific", scientific, true);
  if (scientific){
    ss << std::scientific;
  }else{
    ss << std::fixed;
  }
  ss << v;
  return ss.str();
}

std::string to_string(const vector<int> &V) {
  int precision = 6;
  Config::Get("precision", precision, true);
  std::stringstream ss;
  ss << "[ ";
  for (int i = 0; i < V.size(); i++) {
    std::string s = std::to_string(V[i]);
    ss << s;
    if (i < V.size() - 1) {
      ss << "," << std::string(MAX((8 + precision) - s.size(),0), ' ');
    } else {
      ss << " ]";
    }
  }
  return ss.str();
}

std::string to_string(const RVector &V) {
  int precision = 6;
  Config::Get("precision", precision, true);
  std::stringstream ss;
  ss << "[ ";
  for (int i = 0; i < V.size(); i++) {
    std::string s = format_double(V[i]);
    ss << s;
    if (i < V.size() - 1) {
      ss << ",";
      if ((8 + precision) > s.size()) {
        ss << std::string((8 + precision) - s.size(), ' ');
      }
    }
  }
  ss << " ]";
  return ss.str();
}

std::ostream &operator<<(std::ostream &s, const vector<int> &v) {
  char buf[128];
  for (int i = 0; i < v.size(); ++i) {
    sprintf(buf, "%10d", v[i]);
    s << buf;
    if ((i + 1) % 25 == 0)
      s << endl;
  }
  return s << endl;
}

void PrintRVectorHorizontal(std::ostream &s, const RVector &v) {
  for (int i = 0; i < v.size(); ++i) {
    PrintScalar(s, v[i]);
    if ((i + 1) % 25 == 0)
      s << endl;
  }
  s << endl;
}

std::ostream &operator<<(std::ostream &s, const std::pair<RVector, double> &va) {
  for (int i = 0; i < va.first.size(); ++i) {
    PrintScalar(s, va.first[i]);
    if ((i + 1) % 25 == 0)
      s << endl;
  }
  s << "          ";
  PrintScalar(s, va.second);
  return s << endl;
}

void PrintValues2(const string &name, const RVector &V) {
  for (int i = 0; i < V.size(); ++i) {
    mout.printEntryRaw(name, V[i]);
  }
  if (V.size() < 2) {
    return;
  }

  RVector Difference(0.0, V.size() - 1);
  for (int i = 0; i < Difference.size(); ++i)
    Difference[i] = V[i + 1] - V[i];
  if (V.size() <= 2) {
    return;
  }
  RVector Rate(0.0, V.size() - 2);
  for (int i = 0; i < Rate.size(); ++i) {
    if (abs(Difference[i]) > 1e-9) {
      Rate[i] = Difference[i] / Difference[i + 1];
      mout.printEntryRaw(name + "_Rate", Rate[i]);
    }
  }
  if (Rate[Rate.size() - 1] > 1) {
    double limit = V[V.size() - 1] + Difference[Difference.size() - 1] / (Rate[Rate.size()] - 1);
    mout.printEntryRaw(name + "_LimitRate", limit);
  }
}

void PrintValues(const string &name, const RVector &V, int run, int verbose) {
  vout(1).printEntryRaw(name, to_string(V));
  int N = V.size();
  if (run + 1 < N) N = run + 1;
  if (N < 2) return;
  RVector Difference(0.0, N - 1);
  RVectorT Order(0.0, N - 1);
  for (int i = 0; i < N - 1; ++i) {
    Difference[i] = V[i + 1] - V[i];
    Order[i] = abs(V[i + 1]) > 1e-9 ? log(V[i] / V[i + 1]) / log(2) : -1;
    if (Order[i] == -1 || abs(Order[i]) < 1e-6){
      return;
    }
  }
  vout(2).printEntryRaw(name + " Order", to_string(Order));
  if (N < 3) return;
  RVector Rate(0.0, N - 2);
  for (int i = 0; i < N - 2; ++i) {
    if (abs(Difference[i]) > 1e-9) {
      Rate[i] = Difference[i] / Difference[i + 1];
    }
  }
  vout(2).printEntryRaw(name + " Rates", to_string(Rate));

  if (Rate[N - 3] > 1) {
    double limit = V[N - 1] + Difference[N - 2] / (Rate[N - 3] - 1);
    vout(5).printEntryRaw(name + " Limit", format_double(limit));
  }
}

void PrintMemoryInfo(const IMatrixGraph &mg_fine) {
  if (mg_fine.Size() < 0) {
    THROW("Problem too big. ABORT.")
  }
  int memMB = mg_fine.Size() / (128 * 1024);
  int memMax = PPM->Max(memMB, 0);
  int memTotal = PPM->SumOnCommSplit(memMB, 0);
  int memMaxGB = memMax / 1024;
  int memTotalGB = memTotal / 1024;
  int dof_sum = 0;
  auto &STM = mg_fine.GetMesh();
  for (row r = mg_fine.rows(); r != mg_fine.rows_end(); ++r) {
    if (mg_fine.find_cell(r()) != mg_fine.cells_end()) {
      dof_sum += r.n();
    }
  }
  double psize = mg_fine.pSize();
  double balancingFactor = PPM->Max(dof_sum, 0) / double(PPM->Min(dof_sum, 0));

  mout.PrintInfo("Memory", 1,
                 PrintInfoEntry("Problem size", psize),
                 PrintInfoEntry("Mem/proc in MB", memMax),
                 PrintInfoEntry("Mem/proc in GB", memMaxGB),
                 PrintInfoEntry("Total Memory MB", memTotal),
                 PrintInfoEntry("Total Memory GB", memTotalGB),
                 PrintInfoEntry("LoadBalancingFactor", balancingFactor));

}
