#include "Vector.hpp"
#include "Operator.hpp"
#include "Row.hpp"
#include "VectorAccess.hpp"
Point K = Origin;

void SetQuasiperiodic(const Point &k) { K = k; }

const Point &GetQuasiperiodic() { return K; }

void BCEquations::AddCondition(int i, int k, double lambda) {
  auto entry = equations.find(i);
  if (entry != equations.end()) {
    entry->second[k] = lambda;
  } else {
    std::unordered_map<int, double> cd;
    cd[k] = lambda;
    equations[i] = cd;
  }
}

void DirichletFlags::ClearDirichletFlags() {
  for (int i = 0; i < n; ++i)
    f[i] = false;
}

DirichletFlags::DirichletFlags(int N) : n(N) {
  f = new bool[N];
  ClearDirichletFlags();
}

DirichletFlags::~DirichletFlags() { delete[] f; }

constAB<Vector, Operator> operator*(const Vector &v, const Operator &A) {
  return constAB<Vector, Operator>(v, A);
}

constAB<Operator, Vector> operator*(const Operator &A, const Vector &v) {
  return constAB<Operator, Vector>(A, v);
}

constAB<Scalar, Vector> operator*(const Scalar &a, const Vector &v) {
  return constAB<Scalar, Vector>(a, v);
}

constAB<int, Vector> operator*(const int &a, const Vector &v) { return constAB<int, Vector>(a, v); }

Vector::Vector(std::shared_ptr<const IDiscretization> disc, int spaceLevel, int timeLevel) :
    Vector(std::move(disc), {spaceLevel, timeLevel}) {}

Vector::Vector(std::shared_ptr<const IDiscretization> disc, LevelPair levels) :
    VectorMatrixBase(std::move(disc), levels), data(VectorMatrixBase::size()),
    dirichletFlags(std::make_shared<DirichletFlags>(VectorMatrixBase::size())),
    bcEquations(nullptr) {}

Vector::Vector(Scalar b, std::shared_ptr<const IDiscretization> disc, int spaceLevel,
               int timeLevel) : Vector(b, std::move(disc), {spaceLevel, timeLevel}) {}

Vector::Vector(Scalar b, std::shared_ptr<const IDiscretization> disc, LevelPair levels) :
    VectorMatrixBase(std::move(disc), levels), data(b, VectorMatrixBase::size()),
    dirichletFlags(std::make_shared<DirichletFlags>(VectorMatrixBase::size())),
    bcEquations(nullptr) {
  (*this).accumulateFlag = true;
}

Vector::Vector(const Vector &u) :
    VectorMatrixBase(u), data(u.data), dirichletFlags(u.dirichletFlags),
    bcEquations(u.bcEquations) {
  (*this).accumulateFlag = u.GetAccumulateFlag();
}

Vector::Vector(Vector &&u) :
    VectorMatrixBase(u), data(std::move(u.data)), dirichletFlags(u.dirichletFlags),
    bcEquations(u.bcEquations) {
  (*this).accumulateFlag = u.GetAccumulateFlag();
}

Vector::Vector(Scalar b, const Vector &u) :
    VectorMatrixBase(u), data(b, u.size()), dirichletFlags(u.dirichletFlags),
    bcEquations(u.bcEquations) {
  (*this).accumulateFlag = true;
}

Vector::Vector(const constAB<Operator, Vector> &Ov) :
    VectorMatrixBase(Ov.second()), data(Ov.second().data),
    dirichletFlags(Ov.second().dirichletFlags), bcEquations(Ov.second().bcEquations) {
  Ov.first().multiply(*this, Ov.second());
}

void Vector::Clear() {
  (*this).accumulateFlag = false;
  data = 0.0;
}

Vector &Vector::operator=(const Vector &u) {
#ifdef BUILD_IA
  iadisc = u.iadisc;
#endif
  (*this).accumulateFlag = u.GetAccumulateFlag();
  if (u.size() == size()) {
    data = u.data;
    return *this;
  } else {
    THROW("Sizes of Vector do not fit");
  }
  return *this;
}

Vector &Vector::operator=(Scalar b) {
  (*this).accumulateFlag = true;
  data = b;
  return *this;
}

Vector &Vector::operator*=(Scalar b) {
  data *= b;
  return *this;
}

Vector &Vector::operator*=(int b) {
  data *= b;
  return *this;
}

Vector &Vector::operator/=(Scalar b) {
  data /= b;
  return *this;
}

Vector &Vector::operator+=(const Vector &u) {
  //  data += u.data;
  //  return *this;
  if (u.GetAccumulateFlag() != (*this).GetAccumulateFlag()) {
    Warning("Adding Vectors with different accumulate Flags: this("
            + std::to_string((*this).GetAccumulateFlag()) + ") and u("
            + std::to_string(u.GetAccumulateFlag()) + ")");
  }
  if (u.size() == size()) {
    data += u.data;
    return *this;
  } else {
    THROW("Sizes of Vector do not fit");
  }
  return *this;
}

void Vector::ConsistentAddition(const Vector &u) {
  for (row r = u.rows(); r != u.rows_end(); ++r) {
    procset p = u.find_procset(r());
    if (p == u.procsets_end()) {
      for (int i = 0; i < r.n(); ++i) {
        (*this)(r, i) += u(r, i);
      }
    } else {
      if (p.master() == PPM->Proc(CommSplit())) {
        for (int i = 0; i < r.n(); ++i) {
          (*this)(r, i) += u(r, i);
        }
      }
    }
  }
}

Scalar ConsistentScalarProduct(const Vector &u, const Vector &v) {
  double s = 0.0;
  for (row r = u.rows(); r != u.rows_end(); ++r) {
    procset p = u.find_procset(r());
    if (p == u.procsets_end()) {
      for (int i = 0; i < r.n(); ++i) {
        s += u(r, i) * v(r, i);
      }
    } else {
      if (p.master() == PPM->Proc(u.CommSplit())) {
        for (int i = 0; i < r.n(); ++i) {
          s += u(r, i) * v(r, i);
        }
      }
    }
  }
  return PPM->SumOnCommSplit(s, u.CommSplit());
}

Vector &Vector::operator-=(const Vector &u) {
  if (u.GetAccumulateFlag() != (*this).GetAccumulateFlag()) {
    Warning("Subtracting accumulated and not accumulated Vector");
  }
  data -= u.data;
  return *this;
}

Vector &Vector::operator*=(const Vector &u) {
  if (u.GetAccumulateFlag() != (*this).GetAccumulateFlag()) {
    Warning("Subtracting accumulated and not accumulated Vector");
  }
  data *= u.data;
  return *this;
}

Vector &Vector::operator/=(const Vector &u) {
  if (u.GetAccumulateFlag() != (*this).GetAccumulateFlag()) {
    Warning("Subtracting accumulated and not accumulated Vector");
  }
  data /= u.data;
  return *this;
}

Vector &Vector::operator=(const constAB<Scalar, Vector> &au) {
  (*this).accumulateFlag = au.second().GetAccumulateFlag();
  data.Multiply(au.first(), au.second().data);
  return *this;
}

Vector &Vector::operator+=(const constAB<Scalar, Vector> &au) {
  if (au.second().GetAccumulateFlag() != (*this).GetAccumulateFlag()) {
    Warning("Adding accumulated and not accumulated Vector");
  }
  data.MultiplyPlus(au.first(), au.second().data);
  return *this;
}

Vector &Vector::operator-=(const constAB<Scalar, Vector> &au) {
  if (au.second().GetAccumulateFlag() != (*this).GetAccumulateFlag()) {
    Warning("Subtracting accumulated and not accumulated Vector");
  }
  data.MultiplyMinus(au.first(), au.second().data);
  return *this;
}

Vector &Vector::operator=(const constAB<int, Vector> &au) {
  (*this).accumulateFlag = au.second().GetAccumulateFlag();
  data.Multiply(au.first(), au.second().data);
  return *this;
}

Vector &Vector::operator=(const constAB<Vector, Operator> &vO) {
  vO.second().multiply_transpose(*this, vO.first());
  return *this;
}

Vector &Vector::operator+=(const constAB<Vector, Operator> &vO) {
  vO.second().multiply_transpose_plus(*this, vO.first());
  return *this;
}

Vector &Vector::operator-=(const constAB<Vector, Operator> &vO) {
  vO.second().multiply_transpose_minus(*this, vO.first());
  return *this;
}

Vector &Vector::operator=(const constAB<Operator, Vector> &Ov) {
  Ov.first().multiply(*this, Ov.second());
  return *this;
}

Vector &Vector::operator+=(const constAB<Operator, Vector> &Ov) {
  Ov.first().multiply_plus(*this, Ov.second());
  return *this;
}

Vector &Vector::operator-=(const constAB<Operator, Vector> &Ov) {
  Ov.first().multiply_minus(*this, Ov.second());
  return *this;
}

void Vector::ClearDirichletValues() {
  for (int i = 0; i < size(); ++i) {
    if (dirichletFlags->D(i)) { data[i] = Scalar{}; }
  }

  // linear combinations of DoFs at bnd
  if (!BC()) return;
  for (int i = 0; i < size(); ++i) {
    if (bcEquations->Exists(i)) { data[i] = Scalar{}; }
  }
}

void Vector::Consistent2Additive() {
  for (procset p = (*this).procsets(); p != (*this).procsets_end(); ++p) {
    int i = (*this).Id(p());
    identifyset is = (*this).find_identifyset(p());
    if (is != (*this).identifysets_end()) continue;
    if (p.master() != PPM->Proc((*this).CommSplit()))
      for (int k = 0; k < (*this).Dof(i); ++k)
        (*this)(i, k) = 0.0;
  }
  Vector v((*this));
  (*this) = 0.0;
  for (cell c = (*this).cells(); c != (*this).cells_end(); ++c) {
    ::rows R((*this).GetMatrixGraph(), *c);
    RowValues u_c((*this), R);
    RowValues v_c(v, R);
    for (int i = 0; i < R.size(); ++i) {
      for (int k = 0; k < R[i].n(); ++k)
        u_c(i, k) = v_c(i, k);
    }
  }
}

void Vector::DirichletConsistent() {
  ExchangeBuffer &exBuffer = Buffers().DirichletParallelBuffer();
  for (procset p = procsets(); p != procsets_end(); ++p) {
    int i = Id(p());
    for (int j = 0; j < p.size(); ++j) {
      int q = p[j];
      if (q == PPM->Proc(CommSplit())) continue;
      exBuffer.Send(q) << p();
      for (int k = 0; k < Dof(i); ++k)
        if (D(i, k)) exBuffer.Send(q) << (*this)(i, k);
        else exBuffer.Send(q) << Scalar(infty);
    }
  }
  exBuffer.Communicate();
  for (short q = 0; q < PPM->Size(CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point z;
      //      dpout(100) << "s " << exBuffer.Receive(q).size()
      //                 << " S " << exBuffer.Receive(q).Size()
      //                 << " q " << q
      //                 << endl;
      exBuffer.Receive(q) >> z;
      int i = Idx(z);
      if (i == -1) pout << OUT(z) << " not found \n";
      assert(i != -1);
      for (int k = 0; k < Dof(i); ++k) {
        Scalar a;
        exBuffer.Receive(q) >> a;
        if (a != infty) {
          D(i, k) = true;
          (*this)(i, k) = a;
        }
      }
    }
  exBuffer.ClearBuffers();
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    int j = Id(is());
    for (int i = 0; i < is.size(); ++i) {
      int l = Idx(is[i]);
      if (l == -1) continue;
      for (int k = 0; k < Dof(l); ++k)
        if (D(l, k)) {
          D(j, k) = true;
          (*this)(j, k) = (*this)(l, k);
        }
    }
  }
}

void Vector::Accumulate() {
  if (accumulateFlag == true) {
    Warning("Accumulating already accumulated Vector!");
  } else {
    accumulateFlag = true;
  }
  if (identify()) AccumulateIdentify();
  if (parallel()) AccumulateParallel();
  if (!identify()) return;
  const Point &K = GetQuasiperiodic();
  if (K == Origin) return;
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    const Point &z = is();
    int j = Id(z);
    double zK = 0;
    for (int i = 0; i < z.SpaceDim(); ++i)
      if (z[i] == 1) zK += K[i];
    Scalar s = exp((zK)*iUnit);
    for (int k = 0; k < Dof(j); ++k)
      (*this)(j, k) *= s;
  }
}

void Vector::Collect() {
  accumulateFlag = false;
  if (identify()) CollectIdentify();
  if (parallel()) CollectParallel();
}

void Vector::MakeAdditive() {
  accumulateFlag = false;
  for (row r = rows(); r != rows_end(); ++r) {
    procset p = find_procset(r());
    if (p == procsets_end()) continue;
    if (p.master() == PPM->Proc(CommSplit())) continue;
    for (int k = 0; k < r.n(); ++k)
      (*this)(r, k) = 0;
  }
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    row r = find_row(is());
    if (is.master()) continue;
    for (int k = 0; k < r.n(); ++k)
      (*this)(r, k) = 0;
  }
}

void Vector::Average() {
  ExchangeBuffer &exBuffer = Buffers().AverageParallelBuffer();
  for (procset p = procsets(); p != procsets_end(); ++p) {
    int i = Id(p());
    if (p.master() != PPM->Proc(CommSplit())) continue;
    double s = 1.0 / double(p.size());
    for (int k = 0; k < Dof(i); ++k)
      (*this)(i, k) *= s;
    for (int j = 0; j < p.size(); ++j) {
      int q = p[j];
      if (q == PPM->Proc(CommSplit())) continue;
      exBuffer.Send(q) << p();
      for (int k = 0; k < Dof(i); ++k)
        exBuffer.Send(q) << (*this)(i, k);
    }
  }
  CommunicateVector(exBuffer);
}

bool Vector::IsAccumulated() {
  bool accumulated = true;
  if (identify()) accumulated = accumulated && IsAccumulatedIdentify();
  if (parallel()) accumulated = accumulated && IsAccumulatedParallel();
  return accumulated;
}

bool Vector::IsAccumulatedParallel() {
  ExchangeBuffer &exBuffer = Buffers().AccumulateParallelBuffer();
  for (procset p = procsets(); p != procsets_end(); ++p) {
    int i = Id(p());
    for (int j = 0; j < p.size(); ++j) {
      int q = p[j];
      if (q == PPM->Proc(CommSplit())) continue;
      exBuffer.Send(q) << p();
      for (int k = 0; k < Dof(i); ++k)
        exBuffer.Send(q) << (*this)(i, k);
    }
  }
  exBuffer.Communicate();
  bool accumulated = true;
  for (short q = 0; q < PPM->Size(CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point z;
      exBuffer.Receive(q) >> z;
      int i = Idx(z);
      for (int k = 0; k < Dof(i); ++k) {
        Scalar a;
        exBuffer.Receive(q) >> a;
        if (i != -1) accumulated &= ((*this)(i, k) == a);
      }
    }
  exBuffer.ClearBuffers();

  return PPM->And(accumulated, CommSplit());
}

bool Vector::IsAccumulatedIdentify() {
  ExchangeBuffer &exBuffer = Buffers().AccumulateIdentifyBuffer();
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    const Point &z = is();
    int j = Id(z);
    for (int i = 0; i < is.size(); ++i) {
      Point y = is[i];
      exBuffer.Send(PPM->Proc(CommSplit())) << y;
      for (int k = 0; k < Dof(j); ++k)
        exBuffer.Send(PPM->Proc(CommSplit())) << (*this)(j, k);
    }
  }
  exBuffer.Communicate();
  bool accumulated = true;
  for (short q = 0; q < PPM->Size(CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point z;
      exBuffer.Receive(q) >> z;
      int i = Idx(z);
      for (int k = 0; k < Dof(i); ++k) {
        Scalar a;
        exBuffer.Receive(q) >> a;
        if (i != -1) accumulated &= ((*this)(i, k) == a);
      }
    }
  exBuffer.ClearBuffers();

  return PPM->And(accumulated, CommSplit());
}

bool Vector::IsAdditive() {
  bool additive = true;

  for (row r = rows(); r != rows_end(); ++r) {
    procset p = find_procset(r());
    if (p == procsets_end()) continue;
    if (p.master() == PPM->Proc(CommSplit())) continue;
    for (int k = 0; k < r.n(); ++k)
      additive &= (*this)(r, k) == 0.0;
  }
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    row r = find_row(is());
    if (is.master()) continue;
    for (int k = 0; k < r.n(); ++k)
      additive &= (*this)(r, k) == 0.0;
  }
  return PPM->And(additive, CommSplit());
}

void Vector::CommunicateVector(ExchangeBuffer &exBuffer) {
  exBuffer.Communicate();
  for (short q = 0; q < PPM->Size(CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point z;
      //      std::cout << PPM->Proc() << ": s " << exBuffer.Receive(q).size()
      //                << " S " << exBuffer.Receive(q).Size()
      //                << " q " << q
      //                << endl;
      exBuffer.Receive(q) >> z;
      int i = Idx(z);
      if constexpr (DebugLevel > 0) {
        if (exBuffer.Receive(q).size() + Dof(i) >= exBuffer.ReceiveSize(q)) {
          //          std::cout << PPM->Proc()
          //                    << ": s " << exBuffer.Receive(q).size()
          //                    << " S " << exBuffer.Receive(q).Size()
          //                    << " q " << q
          //                    << " i " << i
          //                    << " z " << z
          //                    << " "
          //                    << endl;
          Exit(to_string(PPM->Proc(CommSplit())) + ": Out of bounds.");
        }
      }
      for (int k = 0; k < Dof(i); ++k) {
        Scalar a;
        exBuffer.Receive(q) >> a;
        //        std::cout << PPM->Proc() << ": " << z << " from " << q
        //                  << " recv " << a
        //                  << " with i = " << i
        //                  << " and k = " << k << endl;
        if (i != -1) (*this)(i, k) += a;
      }
    }
  exBuffer.ClearBuffers();
}

void Vector::AccumulateIdentify() {
  ExchangeBuffer &exBuffer = Buffers().AccumulateIdentifyBuffer();
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    const Point &z = is();
    int j = Id(z);
    for (int i = 0; i < is.size(); ++i) {
      Point y = is[i];
      exBuffer.Send(PPM->Proc(CommSplit())) << y;
      for (int k = 0; k < Dof(j); ++k)
        exBuffer.Send(PPM->Proc(CommSplit())) << (*this)(j, k);
    }
  }
  CommunicateVector(exBuffer);
}

void Vector::AccumulateParallel() {
  ExchangeBuffer &exBuffer = Buffers().AccumulateParallelBuffer();
  for (procset p = procsets(); p != procsets_end(); ++p) {
    int i = Id(p());
    for (int j = 0; j < p.size(); ++j) {
      int q = p[j];
      if (q == PPM->Proc(CommSplit())) continue;
      exBuffer.Send(q) << p();
      for (int k = 0; k < Dof(i); ++k)
        exBuffer.Send(q) << (*this)(i, k);
    }
  }
  CommunicateVector(exBuffer);
}

void Vector::CollectIdentify() {
  ExchangeBuffer &exBuffer = Buffers().CollectIdentifyBuffer();
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    if (is.master()) continue;
    const Point &z = is();
    int j = Id(z);
    exBuffer.Send(PPM->Proc(CommSplit())) << is[0];
    for (int k = 0; k < Dof(j); ++k) {
      exBuffer.Send(PPM->Proc(CommSplit())) << (*this)(j, k);
      (*this)(j, k) = 0;
    }
  }
  CommunicateVector(exBuffer);
}

void Vector::CollectParallel() {
  ExchangeBuffer &exBuffer = Buffers().CollectParallelBuffer();
  for (procset p = procsets(); p != procsets_end(); ++p) {
    int i = Id(p());
    int q = p.master();
    if (q == PPM->Proc(CommSplit())) continue;
    exBuffer.Send(q) << p();
    for (int k = 0; k < Dof(i); ++k) {
      exBuffer.Send(q) << (*this)(i, k);
      (*this)(i, k) = 0;
    }
  }
  CommunicateVector(exBuffer);
}

void Vector::print() const {
  vector<row> R(Size());
  for (row r = rows(); r != rows_end(); ++r)
    R[r.Id()] = r;
  for (int k = 0; k < R.size(); ++k) {
    mout << R[k]() << " :";
    for (int i = 0; i < R[k].n(); ++i)
      mout << " " << D(R[k], i);
    mout << " :";
    for (int i = 0; i < R[k].n(); ++i)
      mout << " " << (*this)(R[k], i);
    mout << endl;
  }
  return;

  for (row r = rows(); r != rows_end(); ++r) {
    mout << r() << " :";
    for (int i = 0; i < r.n(); ++i)
      mout << " " << D(r, i);
    mout << " :";
    for (int i = 0; i < r.n(); ++i)
      mout << " " << (*this)(r, i);
    mout << endl;
  }
}

void Vector::BCadd(int i, int k, int l, double lambda) {
  if (!BC()) { bcEquations = std::make_shared<BCEquations>(); }
  bcEquations->AddCondition(Index(i) + k, l, lambda);
}

bool Vector::BC(const row &r) const {
  bool a = BC(r, 0);
  for (int k = 1; k < r.n(); ++k)
    a |= BC(r, k);
  return a;
}

Saver &Vector::save(Saver &saver) const {
  saver << size();
  for (int i = 0; i < size(); ++i)
    saver << (*this)()[i];
  return saver;
}

Loader &Vector::load(Loader &loader) {
  int s;
  loader >> s;
  if (s != size()) THROW("Size does not fit")
  for (int i = 0; i < s; ++i)
    loader >> (*this)()[i];
  return loader;
}

constAB<Vectors, Operator> operator*(const Vectors &v, const Operator &A) {
  return constAB<Vectors, Operator>(v, A);
}

constAB<Operator, Vectors> operator*(const Operator &A, const Vectors &v) {
  return constAB<Operator, Vectors>(A, v);
}

constAB<Scalar, Vectors> operator*(const Scalar &a, const Vectors &v) {
  return constAB<Scalar, Vectors>(a, v);
}

constAB<int, Vectors> operator*(const int &a, const Vectors &v) {
  return constAB<int, Vectors>(a, v);
}

Vectors::Vectors(int n, std::shared_ptr<const IDiscretization> disc, int spaceLevel, int timeLevel,
                 bool shared) : Vectors(n, disc, {spaceLevel, timeLevel}, shared) {}

Vectors::Vectors(int n, Scalar b, std::shared_ptr<const IDiscretization> disc, int spaceLevel,
                 int timeLevel, bool shared) :
    Vectors(n, b, disc, {spaceLevel, timeLevel}, shared) {}

Vectors::Vectors(int n, std::shared_ptr<const IDiscretization> disc, LevelPair levels,
                 bool shared) : VectorMatrixBase(disc, levels), V(n) {
  if (n < 1) return;
  V[0] = std::make_unique<Vector>(disc, levels);
  for (int i = 1; i < n; ++i) {
    if (shared) V[i] = std::make_unique<Vector>(*V[0]);
    else V[i] = std::make_unique<Vector>(disc, levels);
  }
}

Vectors::Vectors(int n, Scalar b, std::shared_ptr<const IDiscretization> disc, LevelPair levels,
                 bool shared) : VectorMatrixBase(disc, levels), V(n), sharesFlags(shared) {
  if (n < 1) return;
  V[0] = std::make_unique<Vector>(b, disc, levels);
  for (int i = 1; i < n; ++i) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(*V[0]);
    else V[i] = std::make_unique<Vector>(b, disc, levels);
  }
}

Vectors::Vectors(int n, const Vector &u, bool shared) :
    VectorMatrixBase(u), V(n), sharesFlags(shared) {
  for (int i = 0; i < n; ++i) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(u);
    else {
      V[i] = std::make_unique<Vector>(disc, Level());
      *V[i] = u;
    }
  }
}

Vectors::Vectors(int n, const Vectors &u, bool shared) :
    VectorMatrixBase(u), V(n), sharesFlags(shared) {
  for (int i = 0; i < std::min(n, u.size()); ++i) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(u[i]);
    else {
      V[i] = std::make_unique<Vector>(disc, Level());
      *V[i] = u[i];
    }
  }
  for (int i = std::min(n, u.size()); i < n; ++i) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(0.0, u[0]);
    else { V[i] = std::make_unique<Vector>(0.0, disc, Level()); }
  }
}

Vectors::Vectors(const Vectors &u) : VectorMatrixBase(u), V(u.size()), sharesFlags(u.sharesFlags) {
  for (int i = 0; i < V.size(); ++i) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(*(u.V[i]));
    else {
      V[i] = std::make_unique<Vector>(disc, Level());
      *(V[i]) = *(u.V[i]);
    }
  }
}

Vectors::Vectors(Vectors &&u) :
    VectorMatrixBase(u), V(std::move(u.V)), sharesFlags(u.sharesFlags) {}

int computeLength(int length, int startIdx, int endIdx) {
  if (endIdx < 0 || endIdx >= length) endIdx = length - 1;
  return endIdx - startIdx + 1;
}

Vectors::Vectors(const Vectors &u, int startIdx, int endIdx) :
    VectorMatrixBase(u), V(computeLength(u.size(), startIdx, endIdx)), sharesFlags(u.sharesFlags) {
  for (int i = 0; i < V.size(); ++i, ++startIdx) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(*(u.V[startIdx]));
    else {
      V[i] = std::make_unique<Vector>(disc, Level());
      *(V[i]) = *(u.V[startIdx]);
    }
  }
}

Vectors::Vectors(const constAB<Operator, Vectors> &Ov) :
    VectorMatrixBase(Ov.second()), V(Ov.second().size()), sharesFlags(Ov.second().sharesFlags) {
  for (int i = 0; i < V.size(); ++i) {
    if (sharesFlags) V[i] = std::make_unique<Vector>(Ov.second()[i]);
    else {
      V[i] = std::make_unique<Vector>(disc, Level());
      *V[i] = Ov.second()[i];
    }
  }
  Ov.first().multiply(*this, Ov.second());
}

void Vectors::resize(int N) {
  if (N == V.size()) {
    for (int i = 0; i < V.size(); ++i)
      *V[i] = 0.0;
  } else if (N < V.size()) {
    V.erase(V.begin() + N, V.end());
    for (int i = 0; i < V.size(); ++i)
      *V[i] = 0.0;
  } else {
    if (V.size() == 0) V.push_back(std::make_unique<Vector>(0.0, disc, Level()));
    else {
      for (int i = 0; i < V.size(); ++i)
        *V[i] = 0.0;
    }
    for (int i = V.size(); i < N; ++i) {
      if (sharesFlags) V.push_back(std::make_unique<Vector>(0.0, *V[0]));

      else { V.push_back(std::make_unique<Vector>(0.0, disc, Level())); }
    }
  }
}

void Vectors::resizeData(int N) {
  if (N == V.size()) return;
  if (N < V.size()) {
    V.erase(V.begin() + N, V.end());
  } else {
    if (V.size() == 0) V.push_back(std::make_unique<Vector>(0.0, disc, Level()));
    for (int i = V.size(); i < N; ++i)
      V.push_back(std::make_unique<Vector>(0.0, *V[0]));
  }
}

void Vectors::set(const Vectors &u) {
  for (int i = 0; i < min(size(), u.size()); ++i)
    *(V[i]) = u[i];
  for (int i = min(size(), u.size()); i < size(); ++i)
    *(V[i]) = 0.0;
}

void Vectors::erase(int i) { V.erase(V.begin() + i); }

Vectors &Vectors::operator=(const Vector &u) {
#ifdef BUILD_IA
  iadisc = u.GetSharedIADisc();
#endif
  int n = V.size();
  for (int i = 0; i < n; ++i)
    (*V[i]) = u;
  return *this;
}

Vectors &Vectors::operator=(const Vectors &u) {
#ifdef BUILD_IA
  iadisc = u.iadisc;
#endif
  int n = V.size();
  for (int i = 0; i < n; ++i)
    (*V[i]) = *(u.V[i]);
  return *this;
}

Vectors &Vectors::operator=(Scalar b) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) = b;
  return *this;
}

Vectors &Vectors::operator*=(Scalar b) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) *= b;
  return *this;
}

Vectors &Vectors::operator*=(int b) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) *= Scalar(b);
  return *this;
}

Vectors &Vectors::operator/=(Scalar b) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) /= b;
  return *this;
}

Vectors &Vectors::operator+=(const Vector &u) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) += u;
  return *this;
}

Vectors &Vectors::operator+=(const Vectors &u) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) += *(u.V[i]);
  return *this;
}

Vectors &Vectors::operator-=(const Vector &u) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) -= u;
  return *this;
}

Vectors &Vectors::operator-=(const Vectors &u) {
  for (int i = 0; i < V.size(); ++i)
    (*V[i]) -= *(u.V[i]);
  return *this;
}

Vectors &Vectors::operator=(const constAB<Scalar, Vectors> &au) {
  for (int i = 0; i < au.second().size(); ++i)
    *V[i] = constAB<Scalar, Vector>(au.first(), au.second()[i]);
  return *this;
}

Vectors &Vectors::operator+=(const constAB<Scalar, Vectors> &au) {
  for (int i = 0; i < au.second().size(); ++i)
    *V[i] += constAB<Scalar, Vector>(au.first(), au.second()[i]);
  return *this;
}

Vectors &Vectors::operator-=(const constAB<Scalar, Vectors> &au) {
  for (int i = 0; i < au.second().size(); ++i)
    *V[i] -= constAB<Scalar, Vector>(au.first(), au.second()[i]);
  return *this;
}

Vectors &Vectors::operator=(const constAB<int, Vectors> &au) {
  for (int i = 0; i < au.second().size(); ++i)
    *V[i] = constAB<int, Vector>(au.first(), au.second()[i]);
  return *this;
}

Vectors &Vectors::operator=(const constAB<Vectors, Operator> &vO) {
  vO.second().multiply_transpose(*this, vO.first());
  return *this;
}

Vectors &Vectors::operator+=(const constAB<Vectors, Operator> &vO) {
  vO.second().multiply_transpose_plus(*this, vO.first());
  return *this;
}

Vectors &Vectors::operator-=(const constAB<Vectors, Operator> &vO) {
  vO.second().multiply_transpose_minus(*this, vO.first());
  return *this;
}

Vectors &Vectors::operator=(const constAB<Operator, Vectors> &Ov) {
  Ov.first().multiply(*this, Ov.second());
  return *this;
}

Vectors &Vectors::operator+=(const constAB<Operator, Vectors> &Ov) {
  Ov.first().multiply_plus(*this, Ov.second());
  return *this;
}

Vectors &Vectors::operator-=(const constAB<Operator, Vectors> &Ov) {
  Ov.first().multiply_minus(*this, Ov.second());
  return *this;
}

void Vectors::ClearDirichletValues() {
  for (int n = 0; n < size(); ++n)
    V[n]->ClearDirichletValues();
  //  vector<Scalar *> a;
  //  a.resize(V.size());
  //  for (int n = 0; n < V.size(); ++n)
  //    a[n] = (*V[n])();
  //  bool *b = D();
  //  for (int i = 0; i < V[0]->size(); ++i, ++b) {
  //    if (*b)
  //      for (int n = 0; n < V.size(); ++n) *(a[n]) = 0;
  //    for (int n = 0; n < V.size(); ++n) ++(a[n]);
  //  }
}

void Vectors::DirichletConsistent() {
  for (int n = 0; n < size(); ++n)
    V[n]->DirichletConsistent();
}

void Vectors::Accumulate() {
  for (int n = 0; n < size(); ++n)
    V[n]->Accumulate();
}

void Vectors::Collect() {
  for (int n = 0; n < size(); ++n)
    V[n]->Collect();
}

void Vectors::MakeAdditive() {
  for (int n = 0; n < size(); ++n)
    V[n]->MakeAdditive();
}

void Vectors::Average() {
  for (int n = 0; n < size(); ++n)
    V[n]->Average();
}

void Vectors::SetAccumulateFlags(bool flag) {
  for (int n = 0; n < size(); ++n)
    V[n]->SetAccumulateFlag(flag);
}

std::vector<double> Vectors::norm() const {
  std::vector<double> no(V.size());
  for (int i = 0; i < V.size(); ++i)
    no[i] = V[i]->norm();
  return no;
}

void Vectors::Clear() {
  for (int i = 0; i < this->size(); ++i) {
    (*this)[i].Clear();
  }
}

double Vectors::normScalar() const { return std::sqrt((*this) * (*this)); }

void Vectors::ConsistentAddition(const Vectors &u) {
  for (int i = 0; i < u.size(); ++i)
    V[i]->ConsistentAddition(u[i]);
}

double ConsistentScalarProduct(const Vectors &u, const Vectors &v) {
  double s = 0.0;
  for (int i = 0; i < u.size(); ++i) {
    s += ConsistentScalarProduct(u[i], v[i]);
  }
  return s;
}

Saver &Vectors::save(Saver &saver) const {
  saver << size();
  for (int i = 0; i < size(); ++i)
    saver << *(V[i]);
  return saver;
}

Loader &Vectors::load(Loader &loader) {
  int s;
  loader >> s;
  resize(s);
  for (int i = 0; i < s; ++i)
    loader >> *(V[i]);
  return loader;
}

bool Vectors::GetSharesFlags() const { return sharesFlags; }

void Vectors::SetSharesFlags(bool sharesFlags) { Vectors::sharesFlags = sharesFlags; }

Vector operator+(const Vector &u, const Vector &v) {
  if (u.GetAccumulateFlag() != v.GetAccumulateFlag()) {
    Warning("Adding accumulated and not accumulated Vector");
  }
  Vector w = u;
  return w += v;
}

Vector operator-(const Vector &u, const Vector &v) {
  if (u.GetAccumulateFlag() != v.GetAccumulateFlag()) {
    Warning("Subtracting accumulated and not accumulated Vector");
  }
  Vector w = u;
  return w -= v;
}

Scalar operator*(const Vector &u, const Vector &v) {
  // Todo: Is there a case where both Vectors are on different communicators?
  if (u.GetAccumulateFlag() && v.GetAccumulateFlag()) {
    Warning("Both vectors accumulated in scalar product!");
  }
  double value = u.data * v.data;
  value = PPM->SumOnCommSplit(value, u.CommSplit());
  return value;
}

Scalar operator*(const Vectors &u, const Vectors &v) {
  if (u.size() != v.size()) THROW("sizes dont match")
  double s = 0.0;
  for (int i = 0; i < u.size(); ++i) {
    s += u[i] * v[i];
  }
  return s;
}

double norm(const Vector &u) { return u.norm(); }

double norm(const double &u) { return std::abs(u); }

double norm(double &u) { return std::abs(u); }

Vectors operator+(const Vectors &u, const Vectors &v) {
  Vectors w = u;
  return w += v;
}

Vectors operator-(const Vectors &u, const Vectors &v) {
  Vectors w = u;
  return w -= v;
}

Vector HadamardProductVector(const Vector &u, const Vector &v) {
  Vector product(0.0, u);
  for (int i = 0; i < u.size(); ++i) {
    product[i] = u[i] * v[i];
  }
  return product;
}

Vector ComponentSqrt(const Vector &u) {
  Vector root(0.0, u);
  for (int i = 0; i < u.size(); ++i) {
    root[i] = std::sqrt(u[i]);
  }
  return root;
}

Vector ComponentDivide(const Vector &u, const Vector &v) {
  Vector division(0.0, u);
  for (int i = 0; i < u.size(); ++i) {
    if (v[i] == 0)
      Exit("Components of second vector need to be unequal 0") else division[i] = u[i] / v[i];
  }
  return division;
}

Vector &LazyVectors::operator[](int i) {
  if (i < 0 && i >= v.size()) {
    Exit("Trying to access LazyVector with i=" + std::to_string(i) + " but size was "
         + std::to_string(v.size()))
  }
  if (!v[i]) { v[i] = std::make_unique<Vector>(*u); }
  return *v[i];
}

Vector SumVectorOnCommSplit(const Vector &u, int commSplit, int commSplitOfReturn) {
  Vector sumOnCommSplit(0.0, u.GetSharedDisc(), u.Level().WithCommSplit(commSplitOfReturn));
  for (row r = u.rows(); r != u.rows_end(); r++) {
    for (int indexOfDof = 0; indexOfDof < r.NumberOfDofs(); indexOfDof++) {
      sumOnCommSplit(r(), indexOfDof) = PPM->SumOnCommSplit(u(r(), indexOfDof), commSplit);
    }
  }
  return sumOnCommSplit;
}

Vector SumVectorAcrossCommSplit(const Vector &u, int commSplit, int commSplitOfReturn) {
  Vector sumOnCommSplit(0.0, u.GetSharedDisc(), u.Level().WithCommSplit(commSplitOfReturn));
  for (row r = u.rows(); r != u.rows_end(); r++) {
    for (int indexOfDof = 0; indexOfDof < r.NumberOfDofs(); indexOfDof++) {
      double sum = PPM->SumAcrossComm(u(r, indexOfDof), commSplit);
      if (sumOnCommSplit.find_row(r()) != sumOnCommSplit.rows_end()) {
        sumOnCommSplit(r(), indexOfDof) = sum;
      }
    }
  }
  return sumOnCommSplit;
}
