#ifndef _LAPDEF_H_
#define _LAPDEF_H_

#include "GlobalDefinitions.hpp"

// getrf: LU-decomposition
// getri: inverse of LU-decomposition
// getrs: solving AX=B with A=LU
// gemm : matrix/matrix-operation: C = C+AB
// gemv : matrix/vector-operation: Y = Y+AX
// gees : eigenvalue+Schur-vectors
// gesvd: singular value decomposition A = U*S*V^t
// pptrf: LL^T-decomposition of symmetric packed matrices
// pptrs: solving AX=B with LL^T-decomposition (packed)
// sptrf: LDL^T-decomposition of symmetric packed matrices with pivoting
// sptrs: solving AX=B with LDL^T-decomposition
// sytrf: LDL^T-decomposition of symmetric matrices with pivoting
// potrf: LL^T-decomposition
// trsm: solving AX=B with LL^T-decomposition
// syrk: matrix/matrix-operation with symmetric C: C=C+AB

#ifdef DCOMPLEX

#define GEMM zgemm_
#define GETRS zgetrs_
#define LGEMV zgemv_
#define GETRI zgetri_
#define GETRF zgetrf_
#define GESVD zgesvd_
#define POTRF zpotrf_
#define POTRS zpotrs_
#define POTRI zpotri_
#define SYMM zsymm_
#define TRSM ztrsm_
#define TRMM ztrmm_
#define PPTRF zpptrf_
#define TPSV ztpsv_
#define GEQRF zgeqrf_
#define LACPY zlacpy_
#define ORGQR zungqr_
#define GELQF zgelqf_
#define ORGLQ zunglq_
#define LASET zlaset_
#define GELS zgels_
#define SYRK zsyrk_

#else // of #ifdef DCOMPLEX

#define GEMM dgemm_
#define GETRS dgetrs_
#define LGEMV dgemv_
#define GETRI dgetri_
#define GETRF dgetrf_
#define GESVD dgesvd_
#define POTRF dpotrf_
#define POTRS dpotrs_
#define POTRI dpotri_
#define SYMM dsymm_
#define TRSM dtrsm_
#define TRMM dtrmm_
#define PPTRF dpptrf_
#define TPSV dtpsv_
#define GEQRF dgeqrf_
#define LACPY dlacpy_
#define ORGQR dorgqr_
#define GELQF dgelqf_
#define ORGLQ dorglq_
#define LASET dlaset_
#define GELS dgels_
#define SYRK dsyrk_

#define SYTRI dsytri_
#define SYTRF dsytrf_

#define SPTRF dsptrf_
#define SPTRS dsptrs_

#endif // of #ifdef DCOMPLEX


extern "C" void dsytri_(char *UPLO, int *N, void *A, int *LDA, int *IPIV, void *WORK, int *INFO);
extern "C" void dsytrf_(char *UPLO, int *N, void *A, int *LDA, int *IPIV, void *WORK, void *LWORK,
                        int *INFO);


extern "C" void dgels_(const char *TRANS, const int *M, const int *N, const int *NRHS, void *A,
                       const int *lda, void *B, const int *ldb, void *WORK, const int *LWORK,
                       int *INFO);
extern "C" void zgels_(const char *TRANS, const int *M, const int *N, const int *NRHS, void *A,
                       const int *lda, void *B, const int *ldb, void *WORK, const int *LWORK,
                       int *INFO);

extern "C" void dgetrf_(int *M, int *N, void *A, int *LDA, int *IPIV, int *INFO);
extern "C" void zgetrf_(int *M, int *N, void *A, int *LDA, int *IPIV, int *INFO);

extern "C" void dpptrf_(char *UPLO, int *N, void *A, int *INFO);
extern "C" void zpptrf_(char *UPLO, int *N, void *A, int *INFO);

extern "C" void dgetri_(int *N, void *A, int *LDA, int *IPIV, void *WORK, int *LWORK, int *INFO);
extern "C" void zgetri_(int *N, void *A, int *LDA, int *IPIV, void *WORK, int *LWORK, int *INFO);

extern "C" void dgetrs_(const char *TRANS, const int *N, const int *NRHS, void *A, const int *LDA,
                        int *IPIV, void *B, const int *LDB, int *INFO);
extern "C" void zgetrs_(const char *TRANS, const int *N, const int *NRHS, void *A, const int *LDA,
                        int *IPIV, void *B, const int *LDB, int *INFO);

extern "C" void dpptrs_(char *UPLO, int *N, int *NRHS, void *A, void *B, int *LDB, int *INFO);
extern "C" void zpptrs_(char *UPLO, int *N, int *NRHS, void *A, void *B, int *LDB, int *INFO);

// extern "C" void
// dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K,
//        const Scalar *alpha,
//        const void *A, const int *LDA, const void *B, const int *LDB, const Scalar *beta, void *C,
//        const int *LDC);
extern "C" void zgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N,
                       const int *K, const Scalar *alpha, const void *A, const int *LDA,
                       const void *B, const int *LDB, const Scalar *beta, void *C, const int *LDC);

// extern "C" void dgemv_(char *TRANS, int *M, int *N, const Scalar *alpha,
//                        const void *A, const int *LDA, const void *X, int *INCX,
//                        const Scalar *BETA, void *Y,
//                        int *INCY);
extern "C" void zgemv_(const char *TRANS, const int *M, const int *N, const Scalar *alpha,
                       const void *A, const int *LDA, const void *X, const int *INCX,
                       const Scalar *BETA, void *Y, const int *INCY);

extern "C" void dtpsv_(const char *UPLO, const char *TRANSA, const char *DIAG, const int *N,
                       const void *A, void *X, const int *INCX);
extern "C" void ztpsv_(const char *UPLO, const char *TRANSA, const char *DIAG, const int *N,
                       const void *A, void *X, const int *INCX);

extern "C" void dgees_(char *JOBVS, char *SORT, bool *SELECT, int *N, void *A, int *LDA, int *SDIM,
                       void *WR, void *WI, void *VS, int *LDVS, void *WORK, int *LWORK, bool *BWORK,
                       int *INFO);

extern "C" void dgesvd_(char *JOBU, char *JOBVT, int *M, int *N, void *A, int *LDA, double *S,
                        void *U, int *LDU, void *VT, int *LDVT, void *WORK, int *LWORK, int *INFO);
extern "C" void zgesvd_(char *JOBU, char *JOBVT, int *M, int *N, void *A, int *LDA, double *S,
                        void *U, int *LDU, void *VT, int *LDVT, void *WORK, int *LWORK, int *INFO);


extern "C" void dpotrf_(char *UPLO, int *N, void *A, int *LDA, int *INFO);
extern "C" void zpotrf_(char *UPLO, int *N, void *A, int *LDA, int *INFO);

extern "C" void dpotrs_(char *UPLO, int *N, int *NRHS, void *A, int *LDA, void *B, int *LDB,
                        int *INFO);
extern "C" void zpotrs_(char *UPLO, int *N, int *NRHS, void *A, int *LDA, void *B, int *LDB,
                        int *INFO);

extern "C" void dpotri_(char *UPLO, int *N, void *A, int *LDA, int *INFO);
extern "C" void zpotri_(char *UPLO, int *N, void *A, int *LDA, int *INFO);

extern "C" void dsymm_(char *SIDE, char *UPLO, int *M, int *N, double *alpha, void *A, int *LDA,
                       void *B, int *LDB, double *beta, void *C, int *LDC);
extern "C" void zsymm_(char *SIDE, char *UPLO, int *M, int *N, void *alpha, void *A, int *LDA,
                       void *B, int *LDB, void *beta, void *C, int *LDC);

// extern "C" void
// dtrsm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG, const int *M,
//        const int *N,
//        const double *alpha,
//        const void *A, const int *LDA, void *B, const int *LDB);
extern "C" void ztrsm_(const char *SIDE, const char *UPLO, const char *TRANSA, const char *DIAG,
                       const int *M, const int *N, const double *alpha, const void *A,
                       const int *LDA, void *B, const int *LDB);

extern "C" void dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int *M, int *N,
                       double *alpha, void *A, int *LDA, void *B, int *LDB);
extern "C" void ztrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int *M, int *N,
                       double *alpha, void *A, int *LDA, void *B, int *LDB);

extern "C" void dgeqrf_(int *M, int *N, void *A, int *LDA, void *TAU, void *WORK, int *LWORK,
                        int *INFO);
extern "C" void zgeqrf_(int *M, int *N, void *A, int *LDA, void *TAU, void *WORK, int *LWORK,
                        int *INFO);

extern "C" void dlacpy_(char *UPLO, int *M, int *N, void *A, int *LDA, void *B, int *LDB);
extern "C" void zlacpy_(char *UPLO, int *M, int *N, void *A, int *LDA, void *B, int *LDB);

extern "C" void dorgqr_(int *M, int *N, int *K, void *A, int *LDA, void *TAU, void *WORK,
                        int *LWORK, int *INFO);
extern "C" void zungqr_(int *M, int *N, int *K, void *A, int *LDA, void *TAU, void *WORK,
                        int *LWORK, int *INFO);

extern "C" void dgelqf_(int *M, int *N, void *A, int *LDA, void *TAU, void *WORK, int *LWORK,
                        int *INFO);
extern "C" void zgelqf_(int *M, int *N, void *A, int *LDA, void *TAU, void *WORK, int *LWORK,
                        int *INFO);

extern "C" void dorglq_(int *M, int *N, int *K, void *A, int *LDA, void *TAU, void *WORK,
                        int *LWORK, int *INFO);
extern "C" void zunglq_(int *M, int *N, int *K, void *A, int *LDA, void *TAU, void *WORK,
                        int *LWORK, int *INFO);

extern "C" void dlaset_(char *UPLO, int *M, int *N, void *ALPHA, void *BETA, void *A, int *LDA);
extern "C" void zlaset_(char *UPLO, int *M, int *N, void *ALPHA, void *BETA, void *A, int *LDA);

extern "C" void dsptrf_(const char *UPLO, const int *N, void *A, void *IPIV, int *INFO);

extern "C" void dsptrs_(const char *UPLO, const int *N, const int *NRHS, const void *A,
                        const void *IPIV, void *B, const int *LDB, int *INFO);

extern "C" void zsyrk_(const char *UPLO, const char *TRANS, const int *N, const int *K,
                       const Scalar *ALPHA, const void *A, const int *LDA, const Scalar *BETA,
                       void *C, const int *LDC);
extern "C" void dsyrk_(const char *UPLO, const char *TRANS, const int *N, const int *K,
                       const Scalar *ALPHA, const void *A, const int *LDA, const Scalar *BETA,
                       void *C, const int *LDC);


#endif // of #ifndef _LAPDEF_H_
