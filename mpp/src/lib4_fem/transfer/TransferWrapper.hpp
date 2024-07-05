#ifndef _TRANSFER_H_
#define _TRANSFER_H_

#include "Algebra.hpp"

class Transfer {
protected:
  int verbose;
  vector<bool> dirichlet;
  vector<vector<int>> I;
public:
  [[deprecated("Use matching TransferT instead")]]
  Transfer() : verbose(0) {
    Config::Get("TransferVerbose", verbose);
  }

  virtual ~Transfer() {}

  virtual void Destruct() { I = vector<vector<int>>(); }

  virtual void Construct(const IMatrixGraph &, const IMatrixGraph &) = 0;

  virtual void loop(Vector &, const Vector &, int, int) const {}

  virtual void multiply(Vector &, const Vector &) const = 0;

  virtual void loop_transpose(Vector &, const Vector &, int, int) const {}

  virtual void multiply_transpose(Vector &, const Vector &) const = 0;

  virtual void Project(const Vector &f, Vector &c) const = 0;

  virtual Transfer *GetTransferPointer() const { THROW("Not implemented") }

  const vector<vector<int>> &Get_I() const { return I; }
};

class TransferWrapper {
  Transfer *transfer;

  Transfer *GetTransfer(const string &, Point kkk = Origin);

  TransferWrapper(Transfer *transfer, const Vector &coarse, const Vector &fine);
public:
  TransferWrapper(const std::string &name, const Vector &coarse, const Vector &fine);

  ~TransferWrapper() { delete transfer; }

  void Prolongate(const Vector &coarse, Vector &fine) const { transfer->multiply(fine, coarse); }

  void ProlongateTransposed(Vector &coarse, const Vector &fine) const {
    transfer->multiply_transpose(coarse, fine);
  }

  void Restrict(Vector &coarse, const Vector &fine) const {
    THROW("Not implemented in TransferWrapper")
  }

  void Project(Vector &coarse, const Vector &fine) const { transfer->Project(fine, coarse); }

  RMatrix AsMatrix() const { THROW("Not implemented in TransferWrapper") }
};

#endif
