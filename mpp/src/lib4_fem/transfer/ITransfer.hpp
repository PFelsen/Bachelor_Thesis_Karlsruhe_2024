#ifndef ITRANSFER_HPP
#define ITRANSFER_HPP

#include <memory>
#include "Config.hpp"
#include "Operator.hpp"
#include "RMatrix.hpp"
#include "Vector.hpp"

/*
 * ITransfer uses an idiom called "TypeErasure" to ensure usablity within code without relying on
 * Polymorphism. Good examples can be found here: http://www.cplusplus.com/articles/oz18T05o/
 * https://davekilian.com/cpp-type-erasure.html
 */

class ITransfer : public Operator {
  struct TransferConcept {
    virtual ~TransferConcept() = default;

    virtual void Prolongate(const Vector &coarse, Vector &fine) const = 0;

    virtual void ProlongateTransposed(Vector &coarse, const Vector &fine) const = 0;

    virtual void Restrict(Vector &coarse, const Vector &fine) const = 0;

    virtual void Project(Vector &coarse, const Vector &fine) const = 0;

    virtual RMatrix AsMatrix() const = 0;
  };

  template<typename T>
  struct TransferModel : TransferConcept {
    explicit TransferModel(T *transfer) : transfer(transfer) {}

    void Prolongate(const Vector &coarse, Vector &fine) const override {
      transfer->Prolongate(coarse, fine);
    }

    void ProlongateTransposed(Vector &coarse, const Vector &fine) const override {
      transfer->ProlongateTransposed(coarse, fine);
    }

    void Restrict(Vector &coarse, const Vector &fine) const override {
      transfer->Restrict(coarse, fine);
    }

    void Project(Vector &coarse, const Vector &fine) const override {
      transfer->Project(coarse, fine);
    }

    RMatrix AsMatrix() const override { return transfer->AsMatrix(); }

    ~TransferModel() override = default;
  private:
    std::unique_ptr<T> transfer;
  };

  int verbose{-1};
  std::unique_ptr<TransferConcept> transfer;
public:
  template<typename T>
  explicit ITransfer(T *_transfer) : transfer(new TransferModel<T>(_transfer)) {
    Config::Get("TransferVerbose", verbose);
  };

  /**
   * Linear interpolation from coarse vector to fine vector.
   * @param coarse Coarse Vector with same discretization as in constructor.
   * @param fine Fine Vector with same discretization as in constructor.
   */
  void Prolongate(const Vector &coarse, Vector &fine) const { transfer->Prolongate(coarse, fine); }

  /**
   * Restriction from fine vector to coarse vector.
   * In algebraic form, this is the transpose of Prolongate(coarse, fine).
   * @param coarse Coarse Vector with same discretization as in constructor.
   * @param fine Fine Vector with same discretization as in constructor.
   */
  void ProlongateTransposed(Vector &coarse, const Vector &fine) const {
    transfer->ProlongateTransposed(coarse, fine);
  }

  /**
   * Restriction from fine vector to coarse vector.
   * Similar to ProlongateTransposed(coarse, fine) but with additional weight factor.
   * This is NOT supposed to be the transpose of Prolongate(coarse, fine).
   * @param coarse Coarse Vector with same discretization as in constructor.
   * @param fine Fine Vector with same discretization as in constructor.
   */
  void Restrict(Vector &coarse, const Vector &fine) const { transfer->Restrict(coarse, fine); }

  /**
   * Projection from fine vector to coarse vector.
   * @param coarse Coarse Vector with same discretization as in constructor.
   * @param fine Fine Vector with same discretization as in constructor.
   */
  void Project(Vector &coarse, const Vector &fine) const { transfer->Project(coarse, fine); }

  /**
   * Constructs the elongation matrix A, such that A * coarse = fine.
   * @return Elongation matrix with fine.nR() rows and coarse.nR() columns.
   */
  RMatrix AsMatrix() const { return transfer->AsMatrix(); }

  /**
   * Used methods in Operator which have to be implemented
   */
  void multiply(Vector &fine, const Vector &coarse) const override { Prolongate(coarse, fine); }

  void multiply_transpose(Vector &coarse, const Vector &fine) const override {
    ProlongateTransposed(coarse, fine);
  }
};

static constAB<Operator, Vector> operator*(const ITransfer &T, const Vector &v) {
  return constAB<Operator, Vector>(T, v);
}

static constAB<Vector, Operator> operator*(const Vector &v, const ITransfer &T) {
  return constAB<Vector, Operator>(v, T);
}

static constAB<Operator, Vectors> operator*(const ITransfer &T, const Vectors &v) {
  return constAB<Operator, Vectors>(T, v);
}

static constAB<Vectors, Operator> operator*(const Vectors &v, const ITransfer &T) {
  return constAB<Vectors, Operator>(v, T);
}

#endif // ITRANSFER_HPP
