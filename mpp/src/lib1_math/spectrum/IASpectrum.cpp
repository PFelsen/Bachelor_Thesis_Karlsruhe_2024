#include "IASpectrum.hpp"
#include <IACInterval.hpp>

template<typename IAMATRIX, typename EVALUES, typename EVECTORS>
void EVrealEnclosures(const IAMATRIX &A, const IAMATRIX *B, const EVALUES &l, const EVECTORS &x,
                      IAEigenvalueEnclosures &Lambda) {
  int R = l.size();

  double r0 = 0.0;
  double r1 = 0.0;
  for (int i = 0; i < R; ++i) {
    IAInterval s0(0);
    IAInterval s1(0);
    for (int j = 0; j < R; ++j) {
      std::vector<IAInterval> tmpA_j(R);
      std::vector<IAInterval> tmpB_j(R);
      for (int k = 0; k < R; ++k) {
        RVector col_j = x.col(j);
        tmpA_j[k] = accumulate(A.row(k).asVector(), col_j.asVector());
        if (B) tmpB_j[k] = accumulate(B->row(k).asVector(), col_j.asVector());
        else tmpB_j[k] = col_j[k];
      }
      RVector X_i = x.col(i);
      IAInterval a_ij = accumulate(X_i.asVector(), tmpA_j);
      IAInterval b_ij = accumulate(X_i.asVector(), tmpB_j);

      s0 += abs(a_ij - l[j] * b_ij);
      if (i == j) s1 += abs(b_ij - 1);
      else s1 += abs(b_ij);
    }
    r0 = std::max(r0, sup(s0));
    r1 = std::max(r1, sup(s1));
  }

  if (r1 >= 1) THROW("Eigenvalue verification failed!")

  IAInterval D = IAInterval(-1, 1) * sup(r0 / (1 - IAInterval(r1)));
  //    mout << r0 << "  " << r1 << "  " << D << std::endl;

  Lambda.resize(R);
  for (int i = 0; i < R; ++i)
    Lambda[i] = l[i] + D;

  for (int i = 1; i < R;) {
    int j = i;
    while (sup(Lambda[j - 1]) > inf(Lambda[j])) {
      ++j;
      if (j == R) break;
    }

    if (i == j) {
      i++;
    } else {
      double infValue = inf(Lambda[i - 1]);
      double supValue = sup(Lambda[j - 1]);
      for (int n = i - 1; n < j; ++n) {
        Lambda[n].setInf(infValue);
        Lambda[n].setSup(supValue);
      }
      i = j;
    }
  }
}

void verifiedEVreal(const IASymRMatrix &A, IAEigenvalueEnclosures &Lambda) {
  Eigenvalues l;
  Eigenvectors x;
  EVreal(mid(A), l, x);
  IASymRMatrix *B = nullptr;
  EVrealEnclosures(A, B, l, x, Lambda);
}

void verifiedEVreal(const IASymRMatrix &A, IALowerEigenvalueBounds &lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVreal(A, l);
  l.getLowerBounds(lambda);
}

void verifiedEVreal(const IASymRMatrix &A, IAUpperEigenvalueBounds &Lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVreal(A, l);
  l.getUpperBounds(Lambda);
}

void verifiedEVreal(const IASymRMatrix &A, const IASymRMatrix &B, IAEigenvalueEnclosures &Lambda) {
  Eigenvalues l;
  Eigenvectors x;
  EVreal(mid(A), mid(B), l, x);
  EVrealEnclosures(A, &B, l, x, Lambda);
}

void verifiedEVreal(const IASymRMatrix &A, const IASymRMatrix &B, IALowerEigenvalueBounds &lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVreal(A, B, l);
  l.getLowerBounds(lambda);
}

void verifiedEVreal(const IASymRMatrix &A, const IASymRMatrix &B, IAUpperEigenvalueBounds &Lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVreal(A, B, l);
  l.getUpperBounds(Lambda);
}

template<typename IAMATRIX, typename EVALUES, typename EVECTORS>
void EVcomplexEnclosures(const IAMATRIX &A, const IAMATRIX *B, const EVALUES &l, const EVECTORS &x,
                         IAEigenvalueEnclosures &Lambda) {
  int R = l.size();

  double r0 = 0.0;
  double r1 = 0.0;
  for (int i = 0; i < R; ++i) {
    IAInterval s0;
    IAInterval s1;
    for (int j = 0; j < R; ++j) {
      std::vector<IACInterval> tmpA_j(R);
      std::vector<IACInterval> tmpB_j(R);
      for (int k = 0; k < R; ++k) {
        CVector col_j = x.col(j);
        tmpA_j[k] = accumulate(A.row(k).asVector(), col_j.asVector());
        if (B) tmpB_j[k] = accumulate(B->row(k).asVector(), col_j.asVector());
        else tmpB_j[k] = col_j[k];
      }
      CVector conjX_i = conj(x.col(i));
      IACInterval a_ij = accumulate(conjX_i.asVector(), tmpA_j);
      IACInterval b_ij = accumulate(conjX_i.asVector(), tmpB_j);

      s0 += abs(a_ij - l[j] * b_ij);
      if (i == j) s1 += abs(b_ij - 1);
      else s1 += abs(b_ij);
    }
    r0 = std::max(r0, sup(s0));
    r1 = std::max(r1, sup(s1));
  }

  if (r1 >= 1) THROW("Eigenvalue verification failed!")

  IAInterval D = IAInterval(-1, 1) * sup(r0 / (1 - IAInterval(r1)));

  Lambda.resize(R);
  for (int i = 0; i < R; ++i)
    Lambda[i] = l[i] + D;

  for (int i = 1; i < R;) {
    int j = i;
    while (sup(Lambda[j - 1]) > inf(Lambda[j])) {
      ++j;
      if (j == R) break;
    }

    if (i == j) {
      i++;
    } else {
      double infValue = inf(Lambda[i - 1]);
      double supValue = sup(Lambda[j - 1]);
      for (int n = i - 1; n < j; ++n) {
        Lambda[n].setInf(infValue);
        Lambda[n].setSup(supValue);
      }
      i = j;
    }
  }
}

void verifiedEVcomplex(const IAHermCMatrix &A, IAEigenvalueEnclosures &Lambda) {
  Eigenvalues l;
  CEigenvectors x;
  EVcomplex(mid(A), l, x);
  IAHermCMatrix *B = nullptr;
  EVcomplexEnclosures(A, B, l, x, Lambda);
}

void verifiedEVcomplex(const IAHermCMatrix &A, IALowerEigenvalueBounds &lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVcomplex(A, l);
  l.getLowerBounds(lambda);
}

void verifiedEVcomplex(const IAHermCMatrix &A, IAUpperEigenvalueBounds &Lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVcomplex(A, l);
  l.getUpperBounds(Lambda);
}

void verifiedEVcomplex(const IAHermCMatrix &A, const IAHermCMatrix &B,
                       IAEigenvalueEnclosures &Lambda) {
  Eigenvalues l;
  CEigenvectors x;
  EVcomplex(mid(A), mid(B), l, x);
  EVcomplexEnclosures(A, &B, l, x, Lambda);
}

void verifiedEVcomplex(const IAHermCMatrix &A, const IAHermCMatrix &B,
                       IALowerEigenvalueBounds &lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVcomplex(A, B, l);
  l.getLowerBounds(lambda);
}

void verifiedEVcomplex(const IAHermCMatrix &A, const IAHermCMatrix &B,
                       IAUpperEigenvalueBounds &Lambda) {
  IAEigenvalueEnclosures l;
  verifiedEVcomplex(A, B, l);
  l.getUpperBounds(Lambda);
}
