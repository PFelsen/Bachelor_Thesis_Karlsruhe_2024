#ifndef CYCLICPRECONDITIONER_HPP
#define CYCLICPRECONDITIONER_HPP

#include <memory>
#include "PreconditionerBase.hpp"

class CyclicPreconditioner : public Preconditioner {
private:
  vector<std::shared_ptr<Preconditioner>> preconditioners;
  vector<int> cycle_indices;
  mutable int current_cycle_index = 0;
public:
  void multiply(Vector &u, const Vector &b) const override;

  explicit CyclicPreconditioner(const std::string &prefix = "");

  CyclicPreconditioner(const vector<std::shared_ptr<Preconditioner>> &preconditioners,
                       const vector<int> &cycle_indices);

  void Construct(const Matrix &M) override;

  void Destruct() override;

  void reset() override;

  string Name() const override;

  void transpose() const override;
};

#endif // CYCLICPRECONDITIONER_HPP
