#include "SuperLU.hpp"

#include "Sparse.hpp"
#include "TimeDate.hpp"
#include "Vector.hpp"

using std::set;

void SuperLU::Construct(const Matrix &_A) {
  Date Start;
  Matrix A(_A);
  A.EliminateDirichlet();
  if (!A.identify()) A.Accumulate();
  S = new SparseMatrix(A);
  //	S->CheckDiagonal();
  if (A.identify()) S->pc_mat_convert(A);
  Sol = GetSparseSolver(*S, "SuperLU");
  vout(2) << "   decompose SuperLU " << Date() - Start << endl;
}

void SuperLU::Destruct() {
  if (S) delete S;
  S = 0;
  if (Sol) delete Sol;
  Sol = 0;
}

void SuperLU::multiply(Vector &u, const Vector &b) const {
  Date Start;
  u = b;
  if (PPM->Size(u.CommSplit()) > 1) u.Average();
  else if (u.identify()) S->ShrinkIdentify(u, b);
  Sol->Solve(u());
  u.SetAccumulateFlag(false);
  u.Accumulate();
  if (u.identify())
    Exit("Check if this is needed S->ExpandIdentify(u); instead of AccumulateIdentify();") vout(2)
        << "   solve SuperLU " << Date() - Start << endl;
}

SuperLU_local::SuperLU_local() : S(0), Sol(0), rhs(0), IProc(0) {
  int spec = 0;
  Config::Get("SuperLUVerbose", verbose);
  int TimeLevelLoc = 0;
  Config::Get("TimeLevel", TimeLevelLoc);
  if (TimeLevelLoc == 0) Warning("TimeLevel is deprecated!") Config::Get("SHIFTspecial", spec);
  if (spec == 0) shift_special = false;
  else shift_special = true;
  bool Overlap = false;
  Config::Get("Overlap_Distribution", Overlap);
  if (!Overlap) THROW("SuperLU_local failed: Set Overlap_Distribution = 1")
}

void SuperLU_local::CommunicateMatrix(Matrix &A) {
  Date Start_matrix;
  ExchangeBuffer exBuffer(A.CommSplit());
  Scalar *a = A.GetData()();
  if (!PPM->Master(A.CommSplit())) {
    for (row r = A.rows(); r != A.rows_end(); ++r) {
      int id = r.Id();
      int d = A.Entry(A.Diag(id));
      int n = A.Dof(id);
      exBuffer.Send(0) << r();
      for (short k = 0; k < n; ++k)
        for (short l = 0; l < n; ++l)
          exBuffer.Send(0) << a[d + k * n + l];
      exBuffer.Send(0) << int(r.size());
      for (entry e = r.entries(); e != r.entries_end(); ++e) {
        exBuffer.Send(0) << e();
        int m = 2 * n * A.Dof(e.Id());
        int dd = e.GetEntry();
        for (int j = 0; j < m; ++j)
          exBuffer.Send(0) << a[dd + j];
      }
    }
  }
  exBuffer.Communicate();
  for (short q = 0; q < PPM->Size(A.CommSplit()); ++q) {
    Scalar tmp;
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point x;
      exBuffer.Receive(q) >> x;
      row r = A.find_row(x);
      int id = r.Id();
      int d = A.Entry(A.Diag(id));
      int n = A.Dof(id);
      for (int j = 0; j < n * n; ++j) {
        exBuffer.Receive(q) >> tmp;
        a[d + j] += tmp;
      }
      int s;
      exBuffer.Receive(q) >> s;
      for (int k = 0; k < s; ++k) {
        Point y;
        exBuffer.Receive(q) >> y;
        entry e = r.find_entry(y);
        int dd = e.GetEntry();
        int m = 2 * n * A.Dof(e.Id());
        for (int j = 0; j < m; ++j) {
          exBuffer.Receive(q) >> tmp;
          a[dd + j] += tmp;
        }
      }
    }
  }
  vout(3) << "   SuperLU_local: Communicate Matrix " << Date() - Start_matrix << endl;
}

void SuperLU_local::CreateIProc(const Vector &u) {
  Date Start;
  IProc = new std::unordered_map<Point, set<short>>;
  ExchangeBuffer exBuffer(u.CommSplit());
  for (row r = u.rows(); r != u.rows_end(); ++r) {
    exBuffer.Send(0) << r() << r.n();
    for (int k = 0; k < r.n(); ++k)
      exBuffer.Send(0) << u(r, k);
  }
  exBuffer.Communicate();
  for (short q = 0; q < PPM->Size(u.CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point P;
      int dofs;
      exBuffer.Receive(q) >> P >> dofs;
      auto IProciter = IProc->find(P);
      if (IProciter == IProc->end()) {
        set<short> tmp;
        tmp.insert(q);
        (*IProc)[P] = tmp;
      } else {
        (IProciter->second).insert(q);
      }
      for (int k = 0; k < dofs; ++k) {
        Scalar e;
        exBuffer.Receive(q) >> e;
      }
    }
  vout(3) << "   SuperLU_local: Create IProc " << Date() - Start << endl;
}

void SuperLU_local::CollectResidual(const Vector &b) const {
  Date Start;
  ExchangeBuffer exBuffer(b.CommSplit());
  for (row r = b.rows(); r != b.rows_end(); ++r) {
    exBuffer.Send(0) << r() << r.n();
    for (int k = 0; k < r.n(); ++k)
      exBuffer.Send(0) << b(r, k);
  }

  exBuffer.Communicate();

  (*rhs) = 0.0;

  for (short q = 0; q < PPM->Size(b.CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point P;
      int dofs;
      exBuffer.Receive(q) >> P >> dofs;
      row r = rhs->find_row(P);
      for (int k = 0; k < dofs; ++k) {
        Scalar e;
        exBuffer.Receive(q) >> e;
        (*rhs)(r, k) += e;
      }
    }
  vout(4) << "   SuperLU_local: Collect Residual " << Date() - Start << endl;
}

void SuperLU_local::DistributeSolution(Vector &u) const {
  Date Start;
  u = 0.0;
  ExchangeBuffer exBuffer(u.CommSplit());
  if (PPM->Master(u.CommSplit())) {
    for (row r = rhs->rows(); r != rhs->rows_end(); ++r) {
      auto iprii = IProc->find(r());
      set<short> pset = iprii->second;
      for (auto s : pset) {
        exBuffer.Send(s) << r() << r.n();
        for (int dof = 0; dof < r.n(); ++dof)
          exBuffer.Send(s) << (*rhs)(r, dof);
      }
    }
  }
  exBuffer.Communicate();
  for (short q = 0; q < PPM->Size(u.CommSplit()); ++q)
    while (exBuffer.Receive(q).size() < exBuffer.ReceiveSize(q)) {
      Point P;
      int dofs;
      exBuffer.Receive(q) >> P >> dofs;
      row r = rhs->find_row(P);
      for (int dof = 0; dof < dofs; ++dof) {
        Scalar e;
        exBuffer.Receive(q) >> e;
        u(r, dof) = e;
      }
    }
  vout(4) << "   SuperLU_local: Distribute Solution " << Date() - Start << endl;
}

void SuperLU_local::Construct(const Matrix &_A) {
  Date Start;
  Matrix A(_A);
  Vector u(_A.GetVector());
  rhs = new Vector(u);
  CommunicateMatrix(A);
  if (PPM->Master(u.CommSplit())) {
    Date Start_convert;
    S = new SparseMatrix(A);
    if (A.identify()) S->pc_mat_convert(A);
    vout(3) << "   SuperLU_local: Construct Sparse Matrix " << Date() - Start_convert << endl;
    //             S->CheckDiagonal();
    Date Start_create;
    Sol = GetSparseSolver(*S, "SuperLU");
    vout(2) << "   SuperLU_local: Create Solver " << Date() - Start_create << endl;
  }
  CreateIProc(u);
  vout(2) << "   SuperLU_local: Construct " << Date() - Start << endl;
}

void SuperLU_local::Destruct() {
  if (S) delete S;
  S = 0;
  if (Sol) delete Sol;
  Sol = 0;
  if (rhs) delete rhs;
  rhs = 0;
  if (IProc) delete IProc;
  IProc = 0;
}

void SuperLU_local::multiply(Vector &u, const Vector &b) const {
  Date Start;
  CollectResidual(b);
  if (PPM->Master(u.CommSplit())) {
    Date Start_solve;
    if (u.identify()) {
      Vector tmp(*rhs);
      S->ShrinkIdentify(*rhs, tmp);
    }
    Sol->Solve((*rhs)());
    if (u.identify()) S->ExpandIdentify(*rhs);
    vout(4) << "   SuperLU_local: solve " << Date() - Start_solve << endl;
  }
  DistributeSolution(u);
  vout(3) << "   SuperLU_local: multiply " << Date() - Start << endl;
}
