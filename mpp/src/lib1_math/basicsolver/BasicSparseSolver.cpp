#include "BasicSparseSolver.hpp"

#include <cstdlib>


#undef TRUE
#undef FALSE

#ifdef DCOMPLEX
#define Create_Dense_Matrix zCreate_Dense_Matrix
#define Create_CompCol_Matrix zCreate_CompCol_Matrix
#define Print_CompCol_Matrix zPrint_CompCol_Matrix
#define Print_SuperNode_Matrix zPrint_SuperNode_Matrix
#define gstrs zgstrs
#define gstrf zgstrf
#define SLU_DT SLU_Z

typedef complex<double> doublecomplex;

#else
#define Create_CompCol_Matrix dCreate_CompCol_Matrix
#define Create_Dense_Matrix dCreate_Dense_Matrix
#define Print_CompCol_Matrix dPrint_CompCol_Matrix
#define Print_SuperNode_Matrix dPrint_SuperNode_Matrix
#define gstrs dgstrs
#define gstrf dgstrf
#define SLU_DT SLU_D
#endif

#include "slu_ddefs.h"

class SuperSolver : public SparseSolver {
  SuperMatrix *AA; /* A in SLU_NC format used by the factorization routine.*/
  SuperMatrix AC;  /* Matrix postmultiplied by Pc */
  int *etree;
  trans_t trans;

  SuperMatrix A;
  NCformat *Astore;
  double *a;
  int *asub, *xa;
  int *perm_c;   /* column permutation vector */
  int *perm_r;   /* row permutations from partial pivoting */
  SuperMatrix L; /* factor L */
  SCformat *Lstore;
  SuperMatrix U; /* factor U */
  NCformat *Ustore;

  int ldx, info, m, n, nnz;
  double *xact, *rhs;
  mem_usage_t mem_usage;
  superlu_options_t options;
  SuperLUStat_t stat;
  GlobalLU_t Glu;

  void slu_decompose() {
    int lwork = 0, i;

    /* Set default values for some parameters */
    int panel_size; /* panel size */
    int relax;      /* no of columns in a relaxed snodes */
    int permc_spec;

    trans = NOTRANS;
    if (A.Stype == SLU_NR) {
      NRformat *Astore = (NRformat *)A.Store;
      AA = (SuperMatrix *)SUPERLU_MALLOC(sizeof(SuperMatrix));
      Create_CompCol_Matrix(AA, A.ncol, A.nrow, Astore->nnz, (Scalar *)Astore->nzval,
                            Astore->colind, Astore->rowptr, SLU_NC, A.Dtype, A.Mtype);
      trans = TRANS;
    } else {
      if (A.Stype == SLU_NC) AA = &A;
    }

    permc_spec = options.ColPerm;
    if (permc_spec != MY_PERMC && options.Fact == DOFACT) get_perm_c(permc_spec, AA, perm_c);
    etree = intMalloc(A.ncol);

    sp_preorder(&options, AA, perm_c, etree, &AC);
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    /* Compute the LU factorization of A. */
    gstrf(&options, &AC, relax, panel_size, etree, NULL, lwork, perm_c, perm_r, &L, &U, &Glu, &stat,
          &info);
    if (info != 0) {
      if (info < 0) {
        THROW("Error in superLU: " + std::to_string(-info) + "-th argument had an illegal value")
      } else {
        if (info <= AC.ncol) {
          THROW("Error in superLU: A(" + std::to_string(info) + "," + std::to_string(info)
                + " is exactly zero. The factorization has been completed, but the factor U is "
                  "exactly singular, and division by zero will occur if it is used to solve a "
                  "system of equations.")
        } else {
          THROW(
              "Error in superLU: number of bytes allocated when memory allocation failure occurred")
        }
      }
    }
  }

  void slu_solve(SuperMatrix *B) {
    DNformat *Bstore;
    Bstore = (DNformat *)B->Store;
    if (B->ncol < 0 || Bstore->lda < SUPERLU_MAX(0, A.nrow) || B->Stype != SLU_DN
        || B->Dtype != SLU_DT || B->Mtype != SLU_GE)
      info = -7;
    if (info != 0) { THROW("superlu error in slu_solve") }
    gstrs(trans, &L, &U, perm_c, perm_r, B, &stat, &info);
  }
public:
  SuperSolver(SparseRMatrix &M) {
    set_default_options(&options);
    StatInit(&stat);
    n = m = M.rows();
    perm_r = intMalloc(m);
    perm_c = intMalloc(n);
    Create_CompCol_Matrix(&A, m, n, M.Size(), M.nzval(), M.colptr(), M.rowind(), SLU_NR, SLU_DT,
                          SLU_GE);
    slu_decompose();
  }

  void Solve(Scalar *b, int nrhs = 1) {
    if (nrhs == 0) return;
    SuperMatrix B;

    Create_Dense_Matrix(&B, m, nrhs, b, m, SLU_DN, SLU_DT, SLU_GE);

    slu_solve(&B);
    Destroy_SuperMatrix_Store(&B);
  }

  ~SuperSolver() {
    SUPERLU_FREE(etree);
    Destroy_CompCol_Permuted(&AC);
    if (A.Stype == SLU_NR) {
      Destroy_SuperMatrix_Store(AA);
      SUPERLU_FREE(AA);
    }
    StatFree(&stat);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
  }

  void test_output() {
    char CompColL[] = "CompColL";
    Print_CompCol_Matrix(CompColL, &L);
    char CompColU[] = "CompColU";
    Print_CompCol_Matrix(CompColL, &U);
    char SuperNodeL[] = "SuperNodeL";
    Print_SuperNode_Matrix(SuperNodeL, &L);
  }
};

SparseSolver *GetSparseSolver(SparseRMatrix &S, const string &name) {
  if (name == "SuperLU") return new SuperSolver(S);
  Exit("Sparsesolver " + name + " not compiled. File: " + string(__FILE__))
}
