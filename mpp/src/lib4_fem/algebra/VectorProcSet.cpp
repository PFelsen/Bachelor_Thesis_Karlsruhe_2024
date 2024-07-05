#include "VectorProcSet.hpp"
#include "Vector.hpp"

vector_ps::vector_ps() : invIND(0) {
  local_size = 0;
  global_size = 0;
  position = -1;
}

vector_ps::vector_ps(int n) : invIND(n, -1) { local_size = n; }

void VectorProcSet::Add(const ProcSet &Q, int local, int global) {
  ProcSet P(Q);
  P.Sort();
  if (findequalProcSet(P) == -1) {
    int m = vecPS.size();

    vecPS.resize(m + 1);
    vecPS[m] = new vector_ps();
    vecPS[m]->getps().Add(P);
    if (local >= 0) vecPS[m]->set_local_size(local);
    if (global >= 0) vecPS[m]->set_global_size(global);

    vecPS[m]->set_position(m);
    ++num_sorted_elements;
  }
}

void VectorProcSet::Add2(const ProcSet &Q, int &m, int local, int global) {
  ProcSet P(Q);
  P.Sort();
  std::vector<vector_ps *>::iterator it;
  it = vecPS.begin();
  it += m;

  while (it != vecPS.end() && (*it)->getps().size() < P.size()) {
    ++it;
    ++m;
  }
  while (it != vecPS.end() && (*it)->getps().size() == P.size()) {
    int oldm = m;
    unsigned int i = 0;
    while (i < P.size() && (*it)->getps()[i] == P[i])
      ++i;
    if (i == P.size()) { return; }
    if (it != vecPS.end() && (*it)->getps().size() == P.size() && (*it)->getps()[i] < P[i]) {
      ++it;
      ++m;
    }
    if (oldm == m) break;
  }
  vector_ps *v = new vector_ps();
  v->getps().Add(P);
  v->set_local_size(local);
  v->set_global_size(global);
  v->set_position(-1);
  ++num_sorted_elements;
  vecPS.insert(it, v);
}

int VectorProcSet::findequalProcSet(const ProcSet &P) {
  for (int i = 0; i < (int)vecPS.size(); ++i)
    if (vecPS[i]->getps().equalProcSet(P)) return i;
  return -1;
}

void VectorProcSet::swap(int i, int j) {
  vector_ps *vp = vecPS[j];
  vecPS[j] = vecPS[i];
  vecPS[i] = vp;
}

void VectorProcSet::set_position(int index, int pos) {
  if (vecPS[index]->get_position() == -1) ++num_sorted_elements;
  vecPS[index]->set_position(pos);
}

int VectorProcSet::min_PS(const ProcSet &P1) {
  int m = P1[0];
  for (unsigned int i = 1; i < P1.size(); ++i)
    if (P1[i] < m) m = P1[i];
  return m;
}

VectorProcSet::VectorProcSet(const Vector &U) {
  own_proc = PPM->Proc(U.CommSplit());
  total_proc = PPM->Size(U.CommSplit());

  int size = U.GetMatrixGraph().GetRows().SumOfDoFs();

  VectorProcSet::Add(ProcSet(own_proc, 0));

  vecPS.resize(1);
  vecPS[0]->create(size);
  vecPS[0]->set_local_size(0);
  vecPS[0]->set_global_size(0);

  std::vector<std::vector<row>> _R(0);
  for (row r = U.rows(); r != U.rows_end(); ++r) {
    procset _ps = U.find_procset(r());
    if (_ps != U.procsets_end()) {
      int i = findequalProcSet(_ps->second);
      if (i != -1) {
        _R[i - 1].push_back(r);
      } else {
        _R.push_back(std::vector<row>(1, r));
        VectorProcSet::Add(_ps->second);
      }
    } else {
      int idx = U.Index(r.Id());
      for (int k = 0; k < r.n(); ++k)
        vecPS[0]->Add(idx + k);
    }
  }

  vecPS[0]->set_real_local_size();
  vecPS[0]->set_global_size(vecPS[0]->get_local_size());

  for (int b = 0; b < _R.size(); ++b) {
    std::sort(_R[b].begin(), _R[b].end(), [](const row &r0, const row &r1) { return r0() < r1(); });

    int vecb = b + 1;
    vecPS[vecb]->create(size);
    vecPS[vecb]->set_local_size(0);
    int bsort = b;
    for (int n = 0; n < _R[bsort].size(); ++n) {
      int idx = U.Index(_R[b][n].Id());
      for (int k = 0; k < _R[b][n].NumberOfDofs(); ++k)
        vecPS[vecb]->Add(idx + k);
    }
    vecPS[vecb]->set_real_local_size();
    vecPS[vecb]->set_global_size(vecPS[vecb]->get_local_size());
  }
  sort_global();
}

void VectorProcSet::sort_global() {
  for (unsigned int i = 0; i < vecPS.size(); ++i)
    vecPS[i]->getps().Sort();
  for (unsigned int i = 0; i < vecPS.size(); ++i)
    for (unsigned int j = i + 1; j < vecPS.size(); ++j)
      if (vecPS[j]->getps().size() < vecPS[i]->getps().size()) { swap(i, j); }
  for (unsigned int i = 0; i < vecPS.size(); ++i)
    for (unsigned int j = i + 1; j < vecPS.size(); ++j)
      if (vecPS[i]->getps().size() == vecPS[j]->getps().size()) {
        ProcSet Phelp1;
        Phelp1.Add(vecPS[i]->getps());
        ProcSet Phelp2;
        Phelp2.Add(vecPS[j]->getps());
        while (1) {
          int minP1 = min_PS(Phelp1);
          int minP2 = min_PS(Phelp2);
          if (minP1 == minP2) {
            Phelp1.erase(minP1);
            Phelp2.erase(minP2);
          } else if (minP1 > minP2) {
            swap(i, j);
            break;
          } else break;
        }
      }
  num_sorted_elements = vecPS.size();
}

void VectorProcSet::combine_procs() {
  // Bestimme jeweilige Groesse fuer die Kommunikation
  int own_sz = 0;
  int sendsize = 1;
  for (unsigned int i = 0; i < vecPS.size(); ++i)
    own_sz += vecPS[i]->getps().size() + 2;
  int *sending = new int[own_sz];
  int *all_sz = new int[total_proc];
  PPM->Allgather(&own_sz, sendsize, all_sz, vecPS[0]->getps().CommSplit());

  size_t total_rcv = 0;
  int *displs = new int[total_proc];
  displs[0] = 0;
  for (int i = 0; i < total_proc; ++i) {
    total_rcv += all_sz[i];
    if (i < total_proc - 1) displs[i + 1] = displs[i] + all_sz[i];
  }

  int *receiving = new int[total_rcv];
  int *pos = sending;
  for (unsigned int i = 0; i < vecPS.size(); ++i) {
    *(pos++) = int(vecPS[i]->getps().size());
    for (unsigned int j = 0; j < vecPS[i]->getps().size(); ++j)
      *(pos++) = int(vecPS[i]->getps()[j]);
    *(pos++) = int(GetLocalSize(i));
  }
  PPM->Allgatherv(sending, own_sz, receiving, all_sz, displs, vecPS[0]->getps().CommSplit());

  int i = 0;
  int *rcv = receiving;
  int tmps, tmpq, sz;
  while (i < (int)total_rcv) {
    ProcSet PS;
    PS.resize(0);
    tmps = *(rcv++);
    for (int s = 0; s < tmps; ++s) {
      tmpq = *(rcv++);
      PS.Add(tmpq);
    }
    sz = *(rcv++);
    i += tmps + 2;
    int m = 0;
    Add2(PS, m, -1, sz);
  }

  delete[] sending;
  delete[] receiving;
  delete[] all_sz;
  delete[] displs;
}
