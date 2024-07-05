#ifndef ALGEBRAFWD_HPP
#define ALGEBRAFWD_HPP

/**
 * This class is used to avoid unnecessary copies of large data types cf. e.g. "Vector * Operator"
 */
template<class A, class B>
class constAB {
  const A &a;
  const B &b;
public:
  constAB(const A &_a, const B &_b) : a(_a), b(_b) {}

  const A &first() const { return a; }

  const B &second() const { return b; }
};

class DirichletFlags;

class Vector;

class Vectors;

class RowValues;

class RowBndValues;

class MixedRowValues;

class MixedRowBndValues;

class RowBndFaceValues;

class Operator;

class Matrix;

class SparseMatrix;

class DGSizedRowEntries;

class EGRowEntries;

class RowEntries;

class MixedRowEntries;

class PartRowBndValues;

class PartRowEntries;

#endif // ALGEBRAFWD_HPP
