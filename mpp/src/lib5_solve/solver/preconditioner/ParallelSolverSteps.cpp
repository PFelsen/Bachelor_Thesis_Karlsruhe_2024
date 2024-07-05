#include "ParallelSolverSteps.hpp"

void ParallelSolverClusterStep::add_to_Solve(int k) {
  int s = Sol.size();
  Sol.resize(s + 1);
  Sol[s] = k;
}

void ParallelSolverClusterStep::add_to_Schur(int k) {
  int s = Schur.size();
  Schur.resize(s + 1);
  Schur[s] = k;
}

void ParallelSolverClusterStep::set(ProcSet PS1,
                                    ProcSet PS2) { // combined_PS = PS1 \cup PS2
  combined_PS.Add(PS1);
  combined_PS.Add(PS2);
  combined_PS.Sort();
}

void ParallelSolverClusterStep::add(ProcSet PS) {
  combined_PS.Add(PS);
  combined_PS.Sort();
}

void ParallelSolverClusterStep::set_Solve_and_Schur_elements() {
  for (int i = 0; i < vps.size(); ++i)
    if (vps.get_position(i) == -1) {
      if (vps.GetPS(i).subset(combined_PS)) {
        vps.add_to_position(i);
        add_to_Solve(i);
      } else if (vps.GetPS(i).existselementof(combined_PS)) add_to_Schur(i);
    }
}

int ParallelSolverClusterStep::get_Sol_element(int i) const {
  if (i >= (int)Sol.size()) return -1;
  return Sol[i];
}

int ParallelSolverClusterStep::get_Sol_shift(int j) const {
  int size = 0;
  if (j >= (int)Sol.size()) { return -1; }
  for (int i = 0; i < j; ++i)
    size += vps.GetGlobalSize(Sol[i]);
  return size;
}

int ParallelSolverClusterStep::get_Schur_element(int i) const {
  if (i >= (int)Schur.size()) return -1;
  return Schur[i];
}

int ParallelSolverClusterStep::get_Schur_shift(int j) const {
  int size = 0;
  if (j >= (int)Schur.size()) { return -1; }
  for (int i = 0; i < j; ++i)
    size += vps.GetGlobalSize(Schur[i]);
  return size;
}

void ParallelSolverClusterStep::Destruct() {
  Sol.clear();
  Schur.clear();
  if (CommunicationModule) delete CommunicationModule;
  CommunicationModule = NULL;
}

int ParallelSolverClusterStep::matrix_size_Sol() const {
  int size = 0;
  for (unsigned int i = 0; i < Sol.size(); ++i)
    size += vps.GetGlobalSize(Sol[i]);
  return size;
}

int ParallelSolverClusterStep::matrix_size_Schur() const {
  int size = 0;
  for (unsigned int i = 0; i < Schur.size(); ++i)
    size += vps.GetGlobalSize(Schur[i]);
  return size;
}

void ParallelSolverClusterStep::Sort_Schur(std::vector<int> SORT) {
  std::vector<int> INVSORT(Schur);
  for (unsigned int i = 0; i < INVSORT.size(); ++i)
    INVSORT[i] = vps.get_position(Schur[i]);
  sort(INVSORT.begin(), INVSORT.end());
  for (unsigned int i = 0; i < INVSORT.size(); ++i)
    Schur[i] = SORT[INVSORT[i]];
}

void ParallelSolverOneStep::add(ParallelSolverClusterStep *cl1, ParallelSolverClusterStep *cl2) {
  int T = cluster.size();
  cluster.resize(T + 1);
  cluster[T] = new ParallelSolverClusterStep(vps);
  cluster[T]->set(cl1->cPS(), cl2->cPS());

  if (cluster[T]->cPS().existselementof(vps.get_own_proc())) _t = T;
}

void ParallelSolverOneStep::add(ParallelSolverClusterStep *cl1) {
  int T = cluster.size();
  cluster.resize(T + 1);
  cluster[T] = new ParallelSolverClusterStep(vps);
  cluster[T]->add(cl1->cPS());
  if (cluster[T]->cPS().existselementof(vps.get_own_proc())) _t = T;
}

void ParallelSolverOneStep::add(int p) {
  int T = cluster.size();
  cluster.resize(T + 1);
  cluster[T] = new ParallelSolverClusterStep(vps);
  // TODO: change commSplit
  cluster[T]->add(ProcSet(p, 0));
  if (p == vps.get_own_proc()) _t = T;
}

void ParallelSolverOneStep::set_position() {
  for (unsigned int t = 0; t < cluster.size(); ++t)
    cluster[t]->set_Solve_and_Schur_elements();
}

void ParallelSolverOneStep::Destruct() {
  for (unsigned int t = 0; t < cluster.size(); ++t) {
    delete cluster[t];
  }
  cluster.clear();
}

void ParallelSolverOneStep::CreateCommunicationModule(const Communicator &PC) {
  if (_t == -1) { exit(0); }
  cluster[_t]->CreateCommunicationModule(PC, _t);
}

void ParallelSolverOneStep::CreateCommunicationModuleGlobal() {
  if (_t == -1) { exit(0); }
  cluster[_t]->CreateCommunicationModuleGlobal();
}

int ParallelSolverOneStep::total_Solsize() {
  int sz = 0;
  for (unsigned int t = 0; t < cluster.size(); ++t)
    sz += cluster[t]->matrix_size_Sol();
  return sz;
}

ParallelSolverAllSteps::ParallelSolverAllSteps(const Matrix &A) : s(0) {
  vps = new VectorProcSet(A.GetVector());
  vps->combine_procs();

  int P = PPM->Size(A.CommSplit());
  vps->clear_position();
  s.resize(1);
  s[0] = new ParallelSolverOneStep(*vps);
  for (int p = 0; p < P; ++p)
    s[0]->add(p);
  // First step always starts with one proc.

  s[0]->set_position();
  int S = 1;

  while (s[S - 1]->size() > 1) {
    s.resize(S + 1);
    s[S] = new ParallelSolverOneStep(*vps);
    int sz = s[S - 1]->size();
    for (int i = 0; i < sz - 2; i = i + 2)
      s[S]->add(s[S - 1]->get_cluster(i),
                s[S - 1]->get_cluster(
                    i + 1)); // new cluster is the merging of cluster 0,1; 2,3; 4,5 ...
    if (sz % 2 == 1) {
      s[S]->add(s[S - 1]->get_cluster(sz - 1));
    } else s[S]->add(s[S - 1]->get_cluster(sz - 1), s[S - 1]->get_cluster(sz - 2));
    s[S]->set_position();
    S++;
  }

  std::vector<int> SORT;
  for (int i = s.size() - 1; i >= 0; --i)
    for (int t = s[i]->size() - 1; t >= 0; --t)
      for (int k = s[i]->get_cluster(t)->size_Sol() - 1; k >= 0; --k)
        SORT.insert(SORT.begin(), s[i]->get_cluster(t)->get_Sol_element(k));

  for (unsigned int i = 0; i < s.size(); ++i)
    for (int t = 0; t < s[i]->size(); ++t)
      s[i]->get_cluster(t)->Sort_Schur(SORT);

  s[S - 1]->CreateCommunicationModuleGlobal();
  for (int i = S - 2; i >= 0; --i)
    s[i]->CreateCommunicationModule(s[i + 1]->get_own_cluster().getCommunicationModule());
}

void ParallelSolverAllSteps::Add_element(vector<short> ps, int q) {
  int m = int(ps.size());
  if (m == 0) {
    ps.resize(1);
    ps[0] = q;
  }
  for (int i = 0; i < m; ++i)
    if (ps[i] == q) return;
  ps.resize(m + 1);
  if (q < ps[0]) {
    ps[m] = ps[0];
    ps[0] = q;
  } else ps[m] = q;
}

int ParallelSolverAllSteps::parallel_size() {
  int sz = 0;
  for (unsigned int t = 0; t < s.size(); ++t)
    sz += s[t]->total_Solsize();
  return sz;
}

void ParallelSolverAllSteps::Destruct() {
  for (unsigned int i = 0; i < s.size(); ++i) {
    delete s[i];
    s[i] = NULL;
  }
  s.clear();
  if (vps) delete vps;
  vps = NULL;
}