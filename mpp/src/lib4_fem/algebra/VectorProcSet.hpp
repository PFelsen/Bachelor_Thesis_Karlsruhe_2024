#ifndef VECTORPROCSET_HPP
#define VECTORPROCSET_HPP

#include "ProcSet.hpp"

class Vector;

class vector_ps {
  ProcSet ps;
  int local_size;
  int global_size;
  std::vector<int> invIND;
  int position;
public:
  vector_ps();

  vector_ps(int n);

  ~vector_ps() { Destruct(); }

  void create(int n) { invIND.resize(n); }

  void set_local_size(int n) { local_size = n; }

  void set_global_size(int n) { global_size = n; }

  int size() const { return invIND.size(); }

  int operator[](int n) const { return invIND[n]; }

  ProcSet &getps() { return ps; }

  int const get_local_size() const { return local_size; }

  int const get_global_size() const { return global_size; }

  int const get_position() const { return position; }

  void set_position(int pos) { position = pos; }

  void Add(int k) {
    invIND[local_size] = k;
    local_size++;
  }

  void set_real_local_size() { invIND.resize(local_size); }

  void Destruct() { invIND.clear(); }
};

class VectorProcSet {
  int own_proc;
  int total_proc;
  int num_sorted_elements;

  std::vector<vector_ps *> vecPS;

  void Add(const ProcSet &Q, int local = -1, int global = -1);

  void Add2(const ProcSet &Q, int &m, int local = -1, int global = -1);

  int findequalProcSet(const ProcSet &P);

  void swap(int i, int j);

  void set_position(int index, int pos);

  int min_PS(const ProcSet &P1);
public:
  VectorProcSet(){};

  VectorProcSet(const Vector &A);

  ProcSet &GetPS(int i) const { return vecPS[i]->getps(); }

  int const GetGlobalSize(int i) const { return vecPS[i]->get_global_size(); }

  int const GetLocalSize(int i) const { return vecPS[i]->get_local_size(); }

  int const Get_Position(int i) const { return vecPS[i]->get_position(); }

  void SetGlobalSize(int i, int s) { vecPS[i]->set_global_size(s); }

  void SetLocalSize(int i, int s) { vecPS[i]->set_local_size(s); }

  void add_to_position(int index) { set_position(index, num_sorted_elements); }

  void clear_position() {
    for (unsigned int i = 0; i < vecPS.size(); ++i)
      vecPS[i]->set_position(-1);
    num_sorted_elements = 0;
  }

  int get_position(int index) { return vecPS[index]->get_position(); }

  int number_of_sorted_elements() { return num_sorted_elements; }

  int size() { return vecPS.size(); }

  void sort_global();

  void combine_procs();

  vector_ps *getvps(int i) { return vecPS[i]; }

  const vector_ps *getvps(int i) const { return vecPS[i]; }

  void Destruct() {
    for (unsigned int i = 0; i < vecPS.size(); ++i)
      delete vecPS[i];
    vecPS.clear();
  }

  ~VectorProcSet() { Destruct(); }

  int get_own_proc() { return own_proc; }

  int get_total_proc() { return total_proc; }

  std::vector<vector_ps *> getallvps() { return vecPS; }
};

#endif // VECTORPROCSET_HPP
