#include "ProcSet.hpp"
#include "Parallel.hpp"

#include <unordered_set>
#include <ext/numeric>

/*
 * ProcSet
 */

void ProcSet::Add(intproc q) {
  long unsigned int m = size();
  if (m == 0) {
    resize(1);
    (*this)[0] = q;
  }
  for (long unsigned int i = 0; i < m; ++i)
    if ((*this)[i] == q) return;
  resize(m + 1);
  if (q < (*this)[0]) {
    (*this)[m] = (*this)[0];
    (*this)[0] = q;
  } else (*this)[m] = q;
}

void ProcSet::Append(intproc q) {
  long unsigned int m = size();
  for (long unsigned int i = 0; i < m; ++i)
    if ((*this)[i] == q) return;
  resize(m + 1);
  (*this)[m] = q;
}

void ProcSet::Append(const ProcSet &PS) {
  for (int i = 0; i < PS.size(); ++i)
    Append(PS[i]);
}

void ProcSet::Add(const ProcSet &PS) {
  for (long unsigned int i = 0; i < PS.size(); ++i)
    Add(PS[i]);
}

bool ProcSet::Erase() {
  if ((size() == 1) && ((*this)[0] == PPM->Proc(commSplit))) return true;
  for (long unsigned int i = 0; i < size(); ++i)
    if ((*this)[i] == PPM->Proc(commSplit)) return false;
  return false;
}

void ProcSet::erase(intproc q) {
  long unsigned int m = size();
  long unsigned i = 0;
  for (; i < m; ++i)
    if ((*this)[i] == q) break;
  if (i == m) return;
  for (++i; i < m; ++i)
    (*this)[i - 1] = (*this)[i];
  resize(m - 1);
}

intproc ProcSet::master() const { return (*this)[0]; }

void ProcSet::SetMaster(intproc q) {
  if ((*this)[0] == q) return;
  for (long unsigned int i = 1; i < size(); ++i)
    if ((*this)[i] == q) {
      (*this)[i] = (*this)[0];
      (*this)[0] = q;
      return;
    }
  THROW("not valid proc id")
}

bool const ProcSet::equalProcSet(const ProcSet &P) const {
  if (P.size() != (*this).size()) return false;
  for (long unsigned int i = 0; i < (*this).size(); ++i) {
    bool equal = false;
    for (int j = 0; j < P.size(); ++j)
      if ((*this)[i] == P[j]) {
        equal = true;
        break;
      }
    if (!equal) return false;
  }
  return true;
}

bool ProcSet::existselementof(const ProcSet &Q) const {
  for (long unsigned int i = 0; i < (*this).size(); ++i) {
    for (long unsigned int j = 0; j < Q.size(); ++j) {
      if ((*this)[i] == Q[j]) { return true; }
    }
  }
  return false;
}

bool ProcSet::existselementof(const intproc &q) const {
  for (int i = 0; i < (*this).size(); ++i)
    if ((*this)[i] == q) return true;
  return false;
}

bool ProcSet::subset(const ProcSet &Q) const {
  for (long unsigned int i = 0; i < (*this).size(); ++i) {
    bool goon = false;
    for (int j = 0; j < Q.size(); ++j)
      if ((*this)[i] == Q[j]) {
        goon = true;
        break;
      }
    if (goon) continue;
    return false;
  }
  return true;
}

void ProcSet::Sort() { std::sort((*this).begin(), (*this).end()); }

void ProcSet::InsertFront(const ProcSet &PS) {
  if (commSplit != PS.CommSplit()) { THROW("Comsplit do not coincide in ProcSet::InsertFront"); }
  ProcSet copy = *this;
  *this = PS;
  this->Append(copy);
}

std::string ProcSet::ToString() const {
  std::stringstream ss;
  ss << "[ ";
  for (auto q : *this) {
    ss << q << " ";
  }
  ss << "]";
  return ss.str();
}

/*
 * ProcSets
 */

void ProcSets::Copy(const procset &p, const Point &z) { (*this)[z] = *p; }

void ProcSets::Copy(const procset &p) { (*this)[p()] = *p; }

void ProcSets::Add(const Point &z, intproc q) {
  auto p = find(z);
  if (p == end()) (*this)[z] = ProcSet(q, commSplit);
  else p->second.Add(q);
}

void ProcSets::Add(const Point &z, const procset &q) {
  std::unordered_map<Point, ProcSet>::iterator p = find(z);
  if (p == end()) (*this)[z] = *q;
  else
    for (int i = 0; i < q.size(); ++i)
      p->second.Add(q[i]);
}

void ProcSets::AddInfty() {
  auto p = find(Infty);
  if (p == end()) (*this)[Infty] = ProcSet(0);
  for (intproc q = 1; q < PPM->Size(commSplit); ++q)
    p->second.Add(q);
}

void ProcSets::Append(const Point &z, intproc q) {
  auto p = find(z);
  if (p == end()) (*this)[z] = ProcSet(intproc(PPM->Proc(commSplit)), q, commSplit);
  // could be removed if "int" instead of "short".
  else p->second.Append(q);
}

void ProcSets::Append(const Point &z, const procset &q) { Insert(z, *q); }

void ProcSets::Insert(const Point &z, const ProcSet &PS) {
  auto p = find(z);
  if (p == end()) (*this)[z] = PS;
  else p->second.Append(PS);
}

void ProcSets::InsertFront(const Point &z, const ProcSet &PS) {
  auto p = find(z);
  if (p == end()) {
    (*this)[z] = PS;
  } else {
    p->second.InsertFront(PS);
  }
}

void ProcSets::Replace(const Point &z, const ProcSet &PS) { (*this)[z] = PS; }

bool ProcSets::on(std::unordered_map<Point, ProcSet>::iterator p, intproc q) {
  for (int i = 0; i < p->second.size(); ++i)
    if (p->second[i] == q) return true;
  return false;
}

void ProcSets::RemoveIfOn(const Point &z, intproc q) {
  std::unordered_map<Point, ProcSet>::iterator p = find(z);
  if (p == end()) return;
  if (!on(p, q)) return;
  p->second.erase(q);
}

void ProcSets::RemoveSingle() {
  for (std::unordered_map<Point, ProcSet>::iterator p = begin(); p != end();) {
    std::unordered_map<Point, ProcSet>::iterator q = p++;
    if (q->second.size() == 1) erase(q);
  }
}

void ProcSets::Clean() {
  for (auto p = begin(); p != end();) {
    auto q = p++;
    if (!on(q, intproc(PPM->Proc(commSplit)))) erase(q);
    // could be removed if int instead of short
  }
}

bool ProcSets::master(const Point &z) const {
  /* Answers Question:
   * Is proc on which this function is called
   * the master proc of this point ?
   */
  procset p = procset(this->Find(z));
  if (p == procset(this->End())) return true;
  return (p.master() == PPM->Proc(commSplit));
}

void ProcSets::CheckConsistency() const {
  pout << "Running consistency check for ProcSets" << endl;

  bool error = false;

  ExchangeBuffer buffer(commSplit);
  for (procset ps = procsets(); ps != procsets_end(); ++ps) {
    buffer.Send(0) << ps;
  }
  buffer.Communicate();

  if (PPM->Master(commSplit)) {
    std::map<Point, std::map<short, ProcSet>> collectedData;
    for (short q = 0; q < PPM->Size(commSplit); ++q) {
      while (buffer.Receive(q).size() < buffer.Receive(q).Size()) {
        std::pair<Point, ProcSet> ps;
        buffer.Receive(q) >> ps;
        if (collectedData.find(ps.first) == collectedData.end()) {
          collectedData[ps.first] = std::map<short, ProcSet>();
        }
        collectedData[ps.first][q] = ps.second;
      }
    }

    for (auto const &[point, data_p] : collectedData) {
      const ProcSet &firstPs = data_p.begin()->second;
      short first_proc = data_p.begin()->first;
      for (auto const &[q, ps] : data_p) {
        if (ps[0] != firstPs[0]) {
          mout << "master for ProcSet at " << point << " for proc " << first_proc << " is "
               << firstPs[0] << endl;
          mout << "master for ProcSet at " << point << " for proc " << q << " is " << ps[0] << endl
               << endl;
          error = true;
        }
        if (std::find(ps.begin(), ps.end(), q) == ps.end()) {
          mout << "proc not contained in own ProcSet at " << point << " for proc " << q << " " << ps
               << endl
               << endl;
          error = true;
        }
      }

      std::unordered_set<short> first_set(firstPs.begin(), firstPs.end());
      for (auto const &[q, ps] : data_p) {
        std::unordered_set<short> set(ps.begin(), ps.end());
        if (first_set != set) {
          mout << "ProcSet at " << point << " for proc " << first_proc << " is " << firstPs << endl;
          mout << "ProcSet at " << point << " for proc " << q << " is " << ps << endl << endl;
          error = true;
        }
      }
    }
  }

  if (PPM->Or(error, commSplit)) { THROW("ProcSets not consistent") }
}

Buffer &operator<<(Buffer &b, const ProcSet &P) {
  b << P.commSplit << short(P.size());
  for (int i = 0; i < P.size(); ++i)
    b << P[i];
  return b;
}

Buffer &operator>>(Buffer &b, ProcSet &P) {
  short n, q;
  b >> P.commSplit >> n;
  // TODO: The received ProcSet should not be appended, but replace the old one!
  for (int i = 0; i < n; ++i) {
    b >> q;
    P.Append(q);
  }
  return b;
}
