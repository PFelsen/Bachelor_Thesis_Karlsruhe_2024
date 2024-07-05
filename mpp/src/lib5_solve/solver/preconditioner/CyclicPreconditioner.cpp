#include "CyclicPreconditioner.hpp"


#include <memory>

#include "Preconditioner.hpp"

CyclicPreconditioner::CyclicPreconditioner(const std::string &prefix) {
  vector<string> pc_names;
  Config::Get(prefix + "CyclicPCNames", pc_names);
  Config::Get(prefix + "CyclicPCIndices", cycle_indices);
  for (const string &name : pc_names) {
    preconditioners.push_back(std::shared_ptr<Preconditioner>(GetPC(name)));
  }
}

CyclicPreconditioner::CyclicPreconditioner(
    const vector<std::shared_ptr<Preconditioner>> &preconditioners,
    const vector<int> &cycle_indices) :
    preconditioners(preconditioners), cycle_indices(cycle_indices) {}

void CyclicPreconditioner::Construct(const Matrix &M) {
  for (auto &PC : preconditioners) {
    PC->Construct(M);
  }
};

void CyclicPreconditioner::multiply(Vector &u, const Vector &b) const {
  auto PC = preconditioners[cycle_indices[current_cycle_index]];

  Date start;
  u *= 0.975;
  PC->multiply(u, b);

  Time t = Date() - start;


  if (current_cycle_index < cycle_indices.size()) {
    vout(1) << "[CPC] using " << PC->Name() << " t=" << t << endl;
  }

  current_cycle_index++;
  current_cycle_index %= cycle_indices.size();
}

void CyclicPreconditioner::Destruct() {
  for (auto &PC : preconditioners) {
    PC->Destruct();
  }
}

void CyclicPreconditioner::reset() { current_cycle_index = 0; }

string CyclicPreconditioner::Name() const {

  std::stringstream ss;
  ss << "CyclicPreconditioner[" << preconditioners[0]->Name();
  for (int i = 1; i < preconditioners.size(); i++) {
    ss << ", " << preconditioners[i]->Name();
  }
  ss << "](" << cycle_indices[0];
  for (int i = 1; i < cycle_indices.size(); i++) {
    ss << ", " << cycle_indices[i];
  }
  ss << ")";

  return ss.str();
}

void CyclicPreconditioner::transpose() const {
  for (auto &PC : preconditioners) {
    PC->transpose();
  }
}
