#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "AlgebraFwd.hpp"
#include "RVector.hpp"
#include "VectorMatrixBase.hpp"


void SetQuasiperiodic(const Point &);

const Point &GetQuasiperiodic();

using BCEquation = std::unordered_map<int, double>;

/*
 * Contains equations for BC where not only single DoFs are effected but entire equations have to be
 * considered for the correct BC. Note that the (system of) equations need(s) to be solved for
 * specific DoF, i.e., the i-th equation needs to be of the form
 *    \varphi_j(p) + \sum_{l\neq j} \lambda_l \varphi_j(p) = 0
 * where p denotes the nodal point under consideration. In particular, the coefficient of the i-th
 * shape function is supposed to be 1.0. Note that only conditions within DoFs of one single node
 * can be used (cf. equation above holds for fixed p, i.e., for fixed node)
 */
class BCEquations {
  // key (int) corresponds to row id; value contains (k, \lambda_k) as presented above
  std::unordered_map<int, BCEquation> equations{};
public:
  // j denotes the index of node with respect to entire vector size (will be mapped to numbering of
  // row in vector) k denotes the number of DoF at node i
  void AddCondition(int i, int l, double lambda);

  bool Exists(int i) const { return equations.contains(i); }

  bool Exists(int i, int l) const { return equations.at(i).contains(l); }

  double Lambda(int i, int l) const { return equations.at(i).at(l); }

  const BCEquation &Get(int i) const { return equations.at(i); }
};

class DirichletFlags {
  bool *f;
  int n;
public:
  void ClearDirichletFlags();

  DirichletFlags(int N);

  ~DirichletFlags();

  const bool *D() const { return f; }

  bool D(int i) const { return f[i]; }

  bool *D() { return f; }
};

class Vector : public VectorMatrixBase {
  bool accumulateFlag = false;
  BasicVector data;
  std::shared_ptr<DirichletFlags> dirichletFlags;
  std::shared_ptr<BCEquations> bcEquations;
public:
  double Max() const { return data.MaxAccumulated(CommSplit()); }

  double Min() const { return data.MinAccumulated(CommSplit()); }

  // not parallel for consistent Vectors
  double norm() const { return data.normAccumulated(CommSplit()); }

  // abs(sqrt((*this)*(*this))); }

  // Default value provides Vector on fine level
  Vector(std::shared_ptr<const IDiscretization> disc, int spaceLevel = -1, int timeLevel = -1);

  Vector(std::shared_ptr<const IDiscretization> disc, LevelPair levels);

  Vector(Scalar b, std::shared_ptr<const IDiscretization> disc, int spaceLevel = -1,
         int timeLevel = -1);

  Vector(Scalar b, std::shared_ptr<const IDiscretization> disc, LevelPair levels);

  Vector(const Vector &u);

  Vector(Vector &&u);

  Vector(Scalar b, const Vector &u);

  Vector(const constAB<Operator, Vector> &Ov);

  int size() const { return data.size(); }

  void PrintInfo() const { graph.PrintVectorMemoryInfo(); }

#ifdef BUILD_IA

  void SetIADisc(std::shared_ptr<const IAIDiscretization> disc) { iadisc = disc; }

#endif

  Vector &operator=(const Vector &u);

  Vector &operator=(Scalar b);

  Vector &operator*=(Scalar b);

  Vector &operator*=(int b);

  Vector &operator/=(Scalar b);

  Vector &operator+=(const Vector &u);

  Vector &operator-=(const Vector &u);

  Vector &operator*=(const Vector &u);

  Vector &operator/=(const Vector &u);

  void ConsistentAddition(const Vector &u);

  void Clear();

  Vector &operator=(const constAB<Scalar, Vector> &au);

  Vector &operator+=(const constAB<Scalar, Vector> &au);

  Vector &operator-=(const constAB<Scalar, Vector> &au);

  Vector &operator=(const constAB<int, Vector> &au);

  Vector &operator=(const constAB<Vector, Operator> &vO);

  Vector &operator+=(const constAB<Vector, Operator> &vO);

  Vector &operator-=(const constAB<Vector, Operator> &vO);

  Vector &operator=(const constAB<Operator, Vector> &Ov);

  Vector &operator+=(const constAB<Operator, Vector> &Ov);

  Vector &operator-=(const constAB<Operator, Vector> &Ov);

  const Scalar &operator[](int i) const { return data[i]; }

  Scalar &operator[](int i) { return data[i]; }

  const Scalar *operator()() const { return data(); }

  Scalar *operator()() { return data(); }

  const BasicVector &GetData() const { return data; }

  BasicVector &GetData() { return data; }

  const Scalar *operator()(int i) const { return (*this)() + Index(i); }

  Scalar *operator()(int i) { return (*this)() + Index(i); }

  const Scalar *operator()(const row &r) const { return (*this)(r.Id()); }

  Scalar *operator()(const row &r) { return (*this)(r.Id()); }

  const Scalar *operator()(const Point &z) const { return (*this)(find_row(z)); }

  Scalar *operator()(const Point &z) { return (*this)(find_row(z)); }

  Scalar operator()(int i, int k) const { return (*this)(i)[k]; }

  Scalar &operator()(int i, int k) { return (*this)(i)[k]; }

  Scalar operator()(const row &r, int k) const { return (*this)(r)[k]; }

  Scalar &operator()(const row &r, int k) { return (*this)(r)[k]; }

  Scalar operator()(const Point &z, int k) const { return (*this)(z)[k]; }

  Scalar &operator()(const Point &z, int k) { return (*this)(z)[k]; }

  const DirichletFlags &GetDirichletFlags() const { return *dirichletFlags; }

  std::shared_ptr<DirichletFlags> GetDirichletFlagsPtr() const { return dirichletFlags; }

  void ClearDirichletFlags() { dirichletFlags->ClearDirichletFlags(); }

  const bool *D() const { return dirichletFlags->D(); }

  bool *D() { return dirichletFlags->D(); }

  const bool *D(int i) const { return dirichletFlags->D() + Index(i); }

  bool *D(int i) { return dirichletFlags->D() + Index(i); }

  const bool *D(const row &r) const { return D(r.Id()); }

  bool *D(const row &r) { return D(r.Id()); }

  const bool *D(const Point &z) const { return D(find_row(z)); }

  bool *D(const Point &z) { return D(find_row(z)); }

  bool D(int i, int k) const { return D(i)[k]; }

  bool &D(int i, int k) { return D(i)[k]; }

  bool D(const row &r, int k) const { return D(r)[k]; }

  bool &D(const row &r, int k) { return D(r)[k]; }

  bool D(const Point &z, int k) const { return D(z)[k]; }

  bool &D(const Point &z, int k) { return D(z)[k]; }

  const BCEquations &GetBCEquations() const { return *bcEquations; };

  bool BC() const { return bcEquations.operator bool(); }

  bool BC(int i, int k) const { return bcEquations->Exists(Index(i) + k); }

  bool BC(int i, int k, int l) { return bcEquations->Exists(Index(i) + k, l); }

  double BCLambda(int i, int k, int l) const { return bcEquations->Lambda(Index(i) + k, l); }

  const BCEquation &BCGet(int i, int k) const { return bcEquations->Get(Index(i) + k); }

  void BCadd(int i, int k, int l, double lambda);

  bool BC(const row &r) const;

  bool BC(const row &r, int k) const { return BC(r.Id(), k); }

  bool BC(const row &r, int k, int l) { return BC(r.Id(), k, l); }

  double BCLambda(const row &r, int k, int l) const { return BCLambda(r.Id(), k, l); }

  const BCEquation &BCGet(const row &r, int k) const { return BCGet(r.Id(), k); }

  void BCadd(const row &r, int k, int l, double lambda) { BCadd(r.Id(), k, l, lambda); }

  bool BC(const Point &z) const { return BC(find_row(z)); }

  bool BC(const Point &z, int k) const { return BC(find_row(z), k); }

  bool BC(const Point &z, int k, int l) { return BC(find_row(z), k, l); }

  double BCLambda(const Point &z, int k, int l) const { return BCLambda(find_row(z), k, l); }

  const BCEquation &BCGet(const Point &z, int k) const { return BCGet(find_row(z), k); }

  void BCadd(const Point &z, int k, int l, double lambda) { BCadd(find_row(z), k, l, lambda); }

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);

  void print() const;

  friend Scalar operator*(const Vector &u, const Vector &v);

  friend Scalar ConsistentScalarProduct(const Vector &u, const Vector &v);

  Scalar ConsistentNorm() const { return sqrt(ConsistentScalarProduct(*this, *this)); }

  void ClearDirichletValues();

  void DirichletConsistent();

  void Consistent2Additive();

  void Accumulate();

  void Collect();

  void MakeAdditive();

  void Average();

  bool GetAccumulateFlag() const { return accumulateFlag; }

  void SetAccumulateFlag(bool flag) { accumulateFlag = flag; }

  bool IsAccumulated();

  bool IsAccumulatedParallel();

  bool IsAccumulatedIdentify();

  bool IsAdditive();

  void CommunicateVector(ExchangeBuffer &exBuffer);

  template<typename S>
  friend LogTextStream<S> &operator<<(LogTextStream<S> &s, const Vector &u);
private:
  void AccumulateIdentify();

  void AccumulateParallel();

  void CollectIdentify();

  void CollectParallel();
};

inline Saver &operator<<(Saver &saver, const Vector &u) { return u.save(saver); }

inline Loader &operator>>(Loader &loader, Vector &u) { return u.load(loader); }

constAB<Vector, Operator> operator*(const Vector &, const Operator &);

constAB<Operator, Vector> operator*(const Operator &, const Vector &);

constAB<Scalar, Vector> operator*(const Scalar &, const Vector &);

constAB<int, Vector> operator*(const int &, const Vector &);

Vector operator+(const Vector &u, const Vector &v);

Vector operator-(const Vector &u, const Vector &v);

Vector HadamardProductVector(const Vector &u, const Vector &v);

Vector ComponentSqrt(const Vector &u);

Vector ComponentDivide(const Vector &u, const Vector &v);

double norm(const Vector &u);

double norm(const double &u);

double norm(double &u);

class Vectors : public VectorMatrixBase {
  bool sharesFlags;
protected:
  vector<std::unique_ptr<Vector>> V;
public:
  Vectors(int n, std::shared_ptr<const IDiscretization> disc, int spaceLevel = -1,
          int timeLevel = -1, bool shared = true);

  Vectors(int n, Scalar b, std::shared_ptr<const IDiscretization> disc, int spaceLevel = -1,
          int timeLevel = -1, bool shared = true);

  Vectors(int n, std::shared_ptr<const IDiscretization> disc, LevelPair levels, bool shared = true);

  Vectors(int n, Scalar b, std::shared_ptr<const IDiscretization> disc, LevelPair levels,
          bool shared = true);

  Vectors(int n, const Vector &u, bool shared = true);

  Vectors(int n, const Vectors &u, bool shared = true);

  Vectors(const Vectors &u);

  Vectors(Vectors &&u);

  /// selects part of u (indices included)
  Vectors(const Vectors &u, int startIdx, int endIdx = -1);

  Vectors(const constAB<Operator, Vectors> &Ov);

  int size() const { return V.size(); }

  void resize(int N);

  /// Data remains while resizing the object
  void resizeData(int N);

  /// Sets the data, the vectors object will not be resized (cf assigning operator)
  void set(const Vectors &u);

  const Vector &operator[](int i) const { return *(V[i]); }

  Vector &operator[](int i) { return *(V[i]); }

  void erase(int i);

  void reverse() { std::reverse(V.begin(), V.end()); }

  Vector &back() { return *V.back(); }

  const Vector &back() const { return *V.back(); }

#ifdef BUILD_IA

  void SetIADisc(std::shared_ptr<const IAIDiscretization> disc) {
    iadisc = disc;
    for (int i = 0; i < V.size(); ++i)
      V[i]->SetIADisc(disc);
  }

#endif

  Vectors &operator=(const Vector &u);

  Vectors &operator=(const Vectors &u);

  Vectors &operator=(Scalar b);

  Vectors &operator*=(Scalar b);

  Vectors &operator*=(int b);

  Vectors &operator/=(Scalar b);

  Vectors &operator+=(const Vector &u);

  Vectors &operator+=(const Vectors &u);

  Vectors &operator-=(const Vector &u);

  Vectors &operator-=(const Vectors &u);

  Vectors &operator=(const constAB<Scalar, Vectors> &au);

  Vectors &operator+=(const constAB<Scalar, Vectors> &au);

  Vectors &operator-=(const constAB<Scalar, Vectors> &au);

  Vectors &operator=(const constAB<int, Vectors> &au);

  Vectors &operator=(const constAB<Vectors, Operator> &vO);

  Vectors &operator+=(const constAB<Vectors, Operator> &vO);

  Vectors &operator-=(const constAB<Vectors, Operator> &vO);

  Vectors &operator=(const constAB<Operator, Vectors> &Ov);

  Vectors &operator+=(const constAB<Operator, Vectors> &Ov);

  Vectors &operator-=(const constAB<Operator, Vectors> &Ov);

  const Scalar *operator()(int n, int i) const { return (*V[n])(i); }

  Scalar *operator()(int n, int i) { return (*V[n])(i); }

  const Scalar *operator()(int n, const row &r) const { return (*V[n])(r.Id()); }

  Scalar *operator()(int n, const row &r) { return (*V[n])(r.Id()); }

  Scalar operator()(int n, int i, int k) const { return (*V[n])(i, k); }

  Scalar &operator()(int n, int i, int k) { return (*V[n])(i, k); }

  Scalar operator()(int n, const row &r, int k) const { return (*V[n])(r, k); }

  Scalar &operator()(int n, const row &r, int k) { return (*V[n])(r, k); }

  vector<double> norm() const;

  double normScalar() const;

  void ConsistentAddition(const Vectors &u);

  friend double ConsistentScalarProduct(const Vectors &u, const Vectors &v);

  double ConsistentScalarNorm() const { return std::sqrt(ConsistentScalarProduct(*this, *this)); }

  bool GetSharesFlags() const;

  void SetSharesFlags(bool sharesFlags);

  Saver &save(Saver &saver) const;

  Loader &load(Loader &loader);

  void ClearDirichletValues();

  void DirichletConsistent();

  void Accumulate();

  void Collect();

  void MakeAdditive();

  void Average();

  void SetAccumulateFlags(bool flag);

  void Clear();
};

inline Saver &operator<<(Saver &saver, const Vectors &u) { return u.save(saver); }

inline Loader &operator>>(Loader &loader, Vectors &u) { return u.load(loader); }

constAB<Vectors, Operator> operator*(const Vectors &, const Operator &);

constAB<Operator, Vectors> operator*(const Operator &, const Vectors &);

constAB<Scalar, Vectors> operator*(const Scalar &, const Vectors &);

Scalar operator*(const Vectors &u, const Vectors &v);

constAB<int, Vectors> operator*(const int &, const Vectors &);

Vectors operator+(const Vectors &u, const Vectors &v);

Vectors operator-(const Vectors &u, const Vectors &v);

Vector SumVectorOnCommSplit(const Vector &u, int commSplit, int commSplitOfReturn);

Vector SumVectorAcrossCommSplit(const Vector &u, int commSplit, int commSplitOfReturn);

/*
 *  Allows lazy-memory allocation of Vectors.
 */

class LazyVectors {

  vector<std::unique_ptr<Vector>> v;
  std::unique_ptr<Vector> u;
public:
  LazyVectors(int size, const Vector &u) : v(size), u(new Vector(u)){};

  Vector &operator[](int i);
};

template<typename S>
LogTextStream<S> &operator<<(LogTextStream<S> &s, const Vector &u) {
  std::vector<row> R(u.nR());

  for (row r = u.rows(); r != u.rows_end(); ++r)
    R[r.Id()] = r;

  for (int k = 0; k < R.size(); ++k) {
    identifyset identifySet = u.find_identifyset(R[k]());
    procset procSet = u.find_procset(R[k]());

    s << R[k]() << " ";
    if (identifySet == u.identifysets_end()) s << " ";
    else s << "I";

    for (int i = 0; i < R[k].n(); ++i)
      s << " " << u.D(R[k], i);

    s << " :";
    for (int i = 0; i < R[k].n(); ++i) {
      Scalar a = u(R[k], i);
      if (abs(a) < 1e-10) s << " 0";
      else s << " " << a;
    }

    if (u.BC()) {
      s << " :";
      if (u.BC(R[k])) {
        for (int i = 0; i < R[k].n(); ++i) {
          if (u.BC(k, i)) {
            s << " n_" << i;
            const BCEquation &eq_i = u.BCGet(R[k], i);
            for (auto part : eq_i) {
              if (part.second < 0.0) s << part.second << "*n_" << part.first;
              else s << "+" << part.second << "*n_" << part.first;
            }
          }
        }
      } else {
        s << " empty";
      }
    }

    if (procSet != u.procsets_end()) {
      s << " p "
        << "[ ";
      for (int i = 0; i < procSet.size(); i++) {
        s << procSet[i] << " ";
      }
      s << "]";
    }
    s << "\n";

    if (identifySet == u.identifysets_end()) continue;

    for (int i = 0; i < identifySet.size(); ++i) {
      int l = u.Idx(identifySet[i]);
      if (l == -1) continue;
      if ((i == 0) && (!identifySet.master())) s << R[l]() << " M";
      else s << R[l]() << " i";
      for (int i = 0; i < R[l].n(); ++i)
        s << " " << u.D(R[l], i);
      s << " :";
      for (int i = 0; i < R[l].n(); ++i) {
        Scalar a = u(R[l], i);
        if (abs(a) < 1e-10) s << " 0";
        else s << " " << a;
      }
      s << "\n";
    }
  }

  return s;
}

#endif // VECTOR_HPP
