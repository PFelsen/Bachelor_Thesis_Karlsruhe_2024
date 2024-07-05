#ifndef _PARALLELSOLVERSTEPS_H_
#define _PARALLELSOLVERSTEPS_H_

#include "VectorProcSet.hpp"

#include "Communicator.hpp"

#include "Algebra.hpp"

#include <algorithm>

class ParallelSolverClusterStep {
  VectorProcSet &vps; // reference to the vecprocset
  ProcSet combined_PS;
  Communicator *CommunicationModule;
  std::vector<int> Sol;   // elements to Solve
  std::vector<int> Schur; // elements for Schur complement

  void add_to_Solve(int k);
  // if changing the numbering of Sol -> change also vps->position!

  void add_to_Schur(int k);
public:
  ParallelSolverClusterStep(VectorProcSet &v) : vps(v), CommunicationModule(0), Sol(0), Schur(0) {}

  void set(ProcSet PS1, ProcSet PS2);

  void add(ProcSet PS);

  ProcSet cPS() { return combined_PS; }

  void set_Solve_and_Schur_elements();

  int get_Sol_element(int i) const;

  int get_Sol_size(int j) const { return vps.GetGlobalSize(Sol[j]); }

  int get_Sol_shift(int j) const;

  int get_Schur_element(int i) const;

  int get_Schur_size(int j) const { return vps.GetGlobalSize(Schur[j]); }

  int get_Schur_shift(int j) const;

  void Destruct();

  ~ParallelSolverClusterStep() { Destruct(); }

  int size_Sol() const { return Sol.size(); }

  int size_Schur() const { return Schur.size(); }

  int matrix_size_Sol() const;

  int matrix_size_Schur() const;

  void CreateCommunicationModule(const Communicator &PC, int col) {
    CommunicationModule = new Communicator(PC, col);
  }

  void CreateCommunicationModuleGlobal() { CommunicationModule = PPM->Copy(0); }

  Communicator &getCommunicationModule() const { return *CommunicationModule; }

  VectorProcSet &get_vps() const { return vps; }

  void Sort_Schur(std::vector<int> SORT);
};

class ParallelSolverOneStep {
  VectorProcSet &vps; // reference to the vecprocset
  int _t; // defining on which cluster element this processor works on (p \in combined_PS)
  std::vector<ParallelSolverClusterStep *> cluster;
public:
  ParallelSolverOneStep(VectorProcSet &v) : vps(v), _t(-1), cluster(0) {}

  void add(ParallelSolverClusterStep *cl1, ParallelSolverClusterStep *cl2);

  void add(ParallelSolverClusterStep *cl1);

  void add(int p);

  int cluster_number() const { return _t; }

  void set_position();

  ParallelSolverClusterStep *get_cluster(int t) { return cluster[t]; }

  ParallelSolverClusterStep *get_cluster(int t) const { return cluster[t]; }

  ParallelSolverClusterStep &get_own_cluster() const { return *cluster[_t]; }

  void Destruct();

  ~ParallelSolverOneStep() { Destruct(); }

  int size() { return cluster.size(); }

  void CreateCommunicationModule(const Communicator &PC);

  void CreateCommunicationModuleGlobal();

  VectorProcSet &get_vps() const { return vps; }

  int total_Solsize();
};

class ParallelSolverAllSteps {
  VectorProcSet *vps;
  std::vector<ParallelSolverOneStep *> s;
public:
  ParallelSolverAllSteps(const Matrix &A);

  ParallelSolverClusterStep &get_cluster(int step, int t) { return *s[step]->get_cluster(t); }

  ParallelSolverClusterStep &get_own_cluster(int step) {
    return *s[step]->get_cluster(s[step]->cluster_number());
  }

  ParallelSolverOneStep &get_step(int step) { return *s[step]; }

  void Add_element(vector<short> ps, int q);

  int parallel_size();

  int parallel_size(int t) { return s[t]->total_Solsize(); }

  void Destruct();

  int size() { return s.size(); }

  ~ParallelSolverAllSteps() { Destruct(); }
};

#endif // _PARALLELSOLVERSTEPS_H_
