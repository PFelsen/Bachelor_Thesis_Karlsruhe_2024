#include "FFT.hpp"
#include <fftw3.h>

void fft_w(int N, fftw_complex *in, fftw_complex *out);

void fft2_w(int N2, int N1, fftw_complex *in, fftw_complex *out);

void ifft_w(int N, fftw_complex *in, fftw_complex *out);

void FFT::RealToRealVector(const RVector &in, RVector &out) {
  CVector in_complex(in);
  CVector out_complex(out);
  ComplexToComplexVector(in_complex, out_complex);
  for (int i = 0; i < in_complex.size(); i++) {
    out[i] = out_complex[i].real();
  }
}

void FFT::ComplexToComplexVector(const CVector &in, CVector &out) {
  int N = in.size();
  fftw_complex in_arr[N], out_arr[N];
  for (int i = 0; i < N; i++) {
    in_arr[i][0] = in[i].real();
    in_arr[i][1] = in[i].imag();
  }
  fft_w(N, in_arr, out_arr);
  out.resize(N);
  for (int i = 0; i < N; i++) {
    out[i] = {out_arr[i][0], out_arr[i][1]};
  }
}

void FFT::InvComplexToComplexVector(const CVector &in, CVector &out) {
  int N = in.size();
  fftw_complex in_arr[N], out_arr[N];
  for (int i = 0; i < N; i++) {
    in_arr[i][0] = in[i].real();
    in_arr[i][1] = in[i].imag();
  }
  ifft_w(N, in_arr, out_arr);
  out.resize(N);
  for (int i = 0; i < N; i++) {
    out[i] = {out_arr[i][0], out_arr[i][1]};
  }
  out /= out.size();
}

void FFT::RealToComplexVector(const RVector &in, CVector &out) {
  int N = in.size();
  fftw_complex in_arr[N], out_arr[N];
  for (int i = 0; i < N; i++) {
    in_arr[i][0] = in[i];
    in_arr[i][1] = 0.0;
  }
  fft_w(N, in_arr, out_arr);
  out.resize(N);
  for (int i = 0; i < N; i++) {
    out[i] = {out_arr[i][0], out_arr[i][1]};
  }
}

void FFT::InvComplexToRealVector(const CVector &in, RVector &out) {
  CVector out_complex(out);
  InvComplexToComplexVector(in, out_complex);
  for (int i = 0; i < out.size(); i++) {
    out[i] = out_complex[i].real();
  }
}

void fft_w(int N, fftw_complex *in, fftw_complex *out) {
  fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();
}

void ifft_w(int N, fftw_complex *in, fftw_complex *out) {
  fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();
}

void FFT::ComplexToComplexMatrix(const CMatrix &in, CMatrix &out) {
  int rows = in.rows();
  int cols = in.cols();
  auto in_arr = new fftw_complex[rows * cols];
  auto out_arr = new fftw_complex[rows * cols];
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      in_arr[i * rows + j][0] = in[i][j].real();
      in_arr[i * rows + j][1] = in[i][j].imag();
    }
  }
  fft2_w(rows, cols, in_arr, out_arr);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      out[i][j] = {out_arr[i * rows + j][0], out_arr[i * rows + j][1]};
    }
  }

  delete[] in_arr;
  delete[] out_arr;
}

void FFT::RealToRealMatrix(const RMatrix &in, RMatrix &out) {
  CMatrix in_complex(in);
  CMatrix out_complex(out);
  ComplexToComplexMatrix(in_complex, out_complex);
  for (int i = 0; i < in_complex.rows(); i++) {
    for (int j = 0; j < in_complex.cols(); j++) {
      out[i][j] = out_complex[i][j].real();
    }
  }
}

void fft2_w(int N2, int N1, fftw_complex *in, fftw_complex *out) {
  fftw_plan plan = fftw_plan_dft_2d(N2, N1, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();
}

void fft3_w(int N3, int N2, int N1, fftw_complex *in, fftw_complex *out) {
  fftw_plan plan = fftw_plan_dft_3d(N3, N2, N1, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_cleanup();
}

void FFT::ComplexToComplexTensor(const CTensor &in, CTensor &out) {
  int FirstDimension = in.FirstDimension();
  int SecondDimension = in.SecondDimension();
  int ThirdDimension = in.ThirdDimension();

  auto in_arr = new fftw_complex[FirstDimension * SecondDimension * ThirdDimension];
  auto out_arr = new fftw_complex[FirstDimension * SecondDimension * ThirdDimension];

  for (int i = 0; i < FirstDimension; i++) {
    for (int j = 0; j < SecondDimension; j++) {
      for (int k = 0; k < ThirdDimension; k++) {
        in_arr[i + j * FirstDimension + k * FirstDimension * SecondDimension][0] =
            in.real()(i, j, k);
        in_arr[i + j * FirstDimension + k * FirstDimension * SecondDimension][1] =
            in.imag()(i, j, k);
      }
    }
  }
  fft3_w(FirstDimension, SecondDimension, ThirdDimension, in_arr, out_arr);
  for (int i = 0; i < FirstDimension; i++) {
    for (int j = 0; j < SecondDimension; j++) {
      for (int k = 0; k < ThirdDimension; k++) {
        out(i, j, k) = {out_arr[i + j * FirstDimension + k * FirstDimension * SecondDimension][0],
                        out_arr[i + j * FirstDimension + k * FirstDimension * SecondDimension][1]};
      }
    }
  }

  delete[] in_arr;
  delete[] out_arr;
};

void FFT::RealToRealTensor(const RTensor &in, RTensor &out) {
  CTensor in_complex(in);
  CTensor out_complex(out);
  ComplexToComplexTensor(in_complex, out_complex);
  for (int i = 0; i < in_complex.FirstDimension(); i++) {
    for (int j = 0; j < in_complex.SecondDimension(); j++) {
      for (int k = 0; k < in_complex.ThirdDimension(); k++) {
        out(i, j, k) = out_complex(i, j, k).real();
      }
    }
  }
}
