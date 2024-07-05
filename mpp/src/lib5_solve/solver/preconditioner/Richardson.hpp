#ifndef RICHARDSON_HPP
#define RICHARDSON_HPP

#include "PreconditionerBase.hpp"

class Richardson : public Preconditioner {
  double damp;
public:
  Richardson() {
    damp = 1.0;
    Config::Get("PreconditionerDamp", damp);
  }

  void Construct(const Matrix &A) override;

  void multiply(Vector &u, const Vector &b) const override;

  virtual string Name() const { return "Richardson[Damp=" + std::to_string(damp) + "]"; }

  void Destruct() override;

  ~Richardson() override { Destruct(); }
};

#endif
