#ifndef MATRIXTRANSFER_HPP
#define MATRIXTRANSFER_HPP

#include "ITransfer.hpp"

class MatrixTransfer {
private:
  class TransferMatrixGraph : public Rows {
    const VectorMatrixBase &cg;
    const VectorMatrixBase &fg;
    bool simple;
    int d;
    int m;
    int *index;
    int *diag;
    int *column;
    int *matentry;

    void addCell(const Cell &cc, const Cell &fc);

    void numbering();
  public:
    const VectorMatrixBase &CoarseMatrixGraph() const { return cg; }

    const VectorMatrixBase &FineMatrixGraph() const { return fg; }

    bool Simple() const { return simple; }

    int Size() const { return m; }

    int nR() const { return size(); }

    int Index(int i) const { return index[i]; }

    int Diag(int i) const { return diag[i]; }

    int Column(int i) const { return column[i]; }

    int Entry(int i) const { return matentry[i]; }

    TransferMatrixGraph(const VectorMatrixBase &, const VectorMatrixBase &);

    template<typename S>
    friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const TransferMatrixGraph &TG) {
      if (TG.Simple())
        s << endl << "Simple Matrix Transfer" << endl << "Rows: " << endl << (Rows)(TG) << "diag ";
      for (int i = 0; i <= TG.nR(); ++i)
        s << TG.Diag(i) << " ";
      s << endl << "col ";
      for (int i = 0; i < TG.Diag(TG.nR()); ++i)
        s << TG.Column(i) << " ";
      s << endl << "entry ";
      for (int i = 0; i < TG.Diag(TG.nR()); ++i)
        s << TG.Entry(i) << " ";
      return s << endl;
    };

    row rows() const { return row(Rows::begin()); }

    row rows_end() const { return row(Rows::end()); }

    row find_row(const Point &z) const { return row(Rows::find(z)); }
  };

  class TransferMatrix : public Operator {
    BasicVector data;
    const TransferMatrixGraph &TG;
  public:
    double Max() const { return data.MaxAccumulated(); }

    double Min() const { return data.MinAccumulated(); }

    double norm() const { return data.normAccumulated(); }

    bool Simple() const { return TG.Simple(); }

    const VectorMatrixBase &CoarseMatrixGraph() const { return TG.CoarseMatrixGraph(); }

    const VectorMatrixBase &FineMatrixGraph() const { return TG.FineMatrixGraph(); }

    TransferMatrix(const TransferMatrixGraph &_TG) : data(_TG.Size()), TG(_TG) {}

    TransferMatrix &operator=(Scalar b) {
      data = b;
      return *this;
    }

    const Scalar *operator()() const { return data(); }

    Scalar *operator()() { return data(); }

    Scalar *operator()(const Point &x, const Point &y);

    const Scalar *operator()(const Point &x, const Point &y) const;

    void multiply_plus(Vector &fv, const Vector &cv) const override;

    void multiply_transpose_plus(Vector &cv, const Vector &fv) const override;

    template<typename S>
    friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const TransferMatrix &TM) {
      return s << TM.TG;
    }

    friend constAB<Operator, Vector> operator*(const TransferMatrix &A, const Vector &v) {
      return {A, v};
    };

    friend constAB<Vector, Operator> operator*(const Vector &v, const TransferMatrix &A) {
      return {v, A};
    };
  };

  TransferMatrixGraph tMatrixGraph;
  TransferMatrix tMatrix;

  void constructLagrangeTransferMatrix();

  void constructDGTransferMatrix();

  void constructEGTransferMatrix();
public:
  MatrixTransfer(const Vector &coarse, const Vector &fine);

  void Prolongate(const Vector &coarse, Vector &fine) const;

  void ProlongateTransposed(Vector &coarse, const Vector &fine) const;

  void Restrict(Vector &coarse, const Vector &fine) const;

  void Project(Vector &coarse, const Vector &fine) const;

  RMatrix AsMatrix() const;
};


#endif // MATRIXTRANSFER_HPP
