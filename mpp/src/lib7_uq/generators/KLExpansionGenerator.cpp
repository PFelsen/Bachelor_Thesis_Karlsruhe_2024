#include "KLExpansionGenerator.hpp"

#include "RVector.hpp"
#include "Random.hpp"
#include "Spectrum.hpp"

// Methode zum Ziehen eines Samples basierend auf der Methode und dem Mesh
void KLExpansionGenerator::drawSample(const SampleID &id) {
  if (method_ == "MC") {
    RVector r_vec_sample = Random::Uniform(id.commSplit, cellCount, -1, 1);
    // RVector r_vec_sample = Random::Normal(id.commSplit, sparseGridGenerator_.Domain().dim);

    for (int i = 0; i < sampleVec.size(); ++i) {
      double X_tk = 0.0;
      for (int j = 0; j < sparseGridGenerator_.Domain().dim; ++j) {
        X_tk += pow(lambda.Eigenvalues.asVector().at(j), 0.5) * lambda.DEigenfunctions.col(j)[i]
                * r_vec_sample[j];
      }
      sampleVec[i] = X_tk;
    }

    sampleWeight = 1;

  } else if (method_ == "SC") {
    sparseGridGenerator_.DrawSample(id);
    RVector r_vec_sample = sparseGridGenerator_.EvalSample();

    for (int i = 0; i < sampleVec.size(); ++i) {
      double X_tk = 0.0;
      for (int j = 0; j < sparseGridGenerator_.Domain().dim; ++j) {
        X_tk += pow(lambda.Eigenvalues.asVector().at(j), 0.5) * lambda.DEigenfunctions.col(j)[i]
                * r_vec_sample[j];
      }
      sampleVec[i] = X_tk; // Todo: Reshape
    }

    sampleWeight = sparseGridGenerator_.EvalWeight();

    vout(2) << "w: " << sampleWeight << endl;

  } else {
    throw std::invalid_argument("Unbekannte Methode: " + method_);
  }
}

//void KLExpansionGenerator::InitGenerator(std::shared_ptr<Meshes> meshes,
//                                         int init) {
//  size = meshes->fine().CellCountGeometry();
//  sampleVec.resize(meshes->fine().CellCountGeometry());
//  lambda = CovarianceSpectrum(*meshes);
//}

// void KLExpansionGenerator::InitGenerator(std::shared_ptr<TimeSeries> ts,
//                                          std::shared_ptr<Meshes> meshes,
//                                          int init) {
//   if (method_ == "MC") {
//     size = ts->MaxStep() + 1;
//     sampleVec.resize(ts->MaxStep() + 1);
//     lambda = CovarianceSpectrum(*ts);
//   } else if (method_ == "SC") {
//     sparseGridGenerator_.InitGenerator(meshes, init);
//   } else {
//     throw std::invalid_argument("Unbekannte Methode: " + method_);
//   }
// }

double KLExpansionGenerator::CovKernel(
    double dist) { // Todo: make non static and make cov_kernel a method of KL expansion
  return pow(KLExpansionGenerator::sigma, 2)
         * exp(
             -pow(abs(dist), KLExpansionGenerator::alpha)
             / pow(KLExpansionGenerator::length, KLExpansionGenerator::alpha)); // Todo: refactor
}

KLExpansionGenerator::Eigenpair KLExpansionGenerator::CovarianceSpectrum(
    Meshes &meshes) {
  RMatrix K(cellCount);
  RMatrix W_root(cellCount);
  RMatrix W_root_inv(cellCount);

  int i = 0;

  // try to generate the complete time series vector or retrieve it somehow

  for (cell c1 = meshes.fine().cells(); c1 != meshes.fine().cells_end(); ++c1) {
    int j = 0;
    for (cell c2 = meshes.fine().cells(); c2 != meshes.fine().cells_end();
         ++c2) {
      K(i, j) = CovKernel(dist(c1.Center(),
               c2.Center()));  // dist(c1.Centre(),c2.Centre()) Can feed in the
      // distance of the cell centers
      if (i == j and i != 0 and i != cellCount - 1) {
        W_root(i, j) = 1 / sqrt((double)(cellCount - 1));
        W_root_inv(i, j) = sqrt((double)(cellCount - 1));
      } else if ((i == j and i == 0) or (i == j and i == cellCount - 1)) {
        W_root(i, j) = 1 / sqrt((double)(2 * (cellCount - 1)));
        W_root_inv(i, j) = sqrt((double)(2 * (cellCount - 1)));
      }
      j += 1;
    }
    i += 1;
  }

  RMatrix A = W_root * K * W_root;

  CEigenvalues lambda;
  CEigenvectors x;
  EVreal(A, lambda, x); // how do we cut off?

  Eigenpair DE;

  DE.Eigenvalues = lambda.real();
  DE.DEigenfunctions = W_root_inv * x.real();

  return DE;
}

KLExpansionGenerator::Eigenpair KLExpansionGenerator::CovarianceSpectrum(
    TimeSeries &ts) {
  // try to generate the complete time series vector or retrieve it somehow
  RVector timeseries(ts.MaxStep() + 1);

  for (int k = 0; k <= timeseries.size(); ++k) {
    timeseries[k] = ts.Time();
    ts.NextTimeStep();
  }

  int n = timeseries.size();
  RMatrix K(n);
  RMatrix W_root(n);
  RMatrix W_root_inv(n);

  for (int i = 0; i < timeseries.size(); ++i) {
    double ti = timeseries[i];
    for (int j = 0; j < timeseries.size(); ++j) {
      double tj = timeseries[j];
      K(i, j) = CovKernel(abs(ti - tj)); // dist(c1.Centre(),c2.Centre()) Can feed
      // in the distance of the cell centers
      if (i == j and i != 0 and i != n - 1) {
        W_root(i, j) = 1 / sqrt((double) (n - 1));
        W_root_inv(i, j) = sqrt((double) (n - 1));
      } else if (i == j and i == 0) {
        W_root(i, j) = 1 / sqrt((double) (2 * (n - 1)));
        W_root_inv(i, j) = sqrt((double) (2 * (n - 1)));
      } else if (i == j and i == n - 1) {
        W_root(i, j) = 1 / sqrt((double) (2 * (n - 1)));
        W_root_inv(i, j) = sqrt((double) (2 * (n - 1)));
      }
    }
    std::cout << endl;
  }

  RMatrix A = W_root * K * W_root;

  CEigenvalues lambda;
  CEigenvectors x;
  EVreal(A, lambda, x);

  Eigenpair DE;

  DE.Eigenvalues = lambda.real();
  DE.DEigenfunctions = W_root_inv * x.real();

  return DE;
}

// Scalar OrnsteinUhlenbeckProcess::eigenfunction(int k, double x) {
//   double lamb_k = eigenvalues[k];
//   double c1 = sqrt((2 - lamb_k) / (lamb_k + 1));
//   double c2 = sqrt(lamb_k / (lamb_k + 1));
//
//   return c1 * std::cos(sqrt(c1 / c2) * x) + c2 * std::sin(sqrt(c1 / c2) * x);
// }
//
// void OrnsteinUhlenbeckProcess::drawSample(const SampleID &id) {
//   //  RVector r_vec_sample =
//   //      Random::Uniform(commSplit, size, -sqrt(3), sqrt(3));
//   //  int k = 0;
//   //  for (cell c = meshes->fine().cells(); c != meshes->fine().cells_end(); ++c) {
//   //    double X_tk = 0.0;
//   //    for (int i = 0; i < eigenvalues.size(); ++i) {
//   //      X_tk += pow(eigenvalues[i], 0.5) * eigenfunction(i, c.Center()[0]) *
//   //              r_vec_sample[i];
//   //    }
//   //    sampleVec[k] = X_tk;
//   //    k += 1;
//   //  }
// }
//
// RVector OrnsteinUhlenbeckProcess::CovarianceSpectrum(int n) {
//   RMatrix K(n + 2);
//   RMatrix W_root(n + 2);
//   RMatrix W_root_inv(n + 2);
//
//   auto t = linspace(0, 1, n + 2);
//
//   for (int i = 0; i < n + 2; ++i) {
//     for (int j = 0; j < n + 2; ++j) {
//       K(i, j) = cov_kernel(t[i] - t[j]);
//
//       if (i == j and i != 0 and i != n + 1) {
//         W_root(i, j) = 1 / sqrt((double) (n + 1));
//         W_root_inv(i, j) = sqrt((double) (n + 1));
//       } else if (i == j and i == 0) {
//         W_root(i, j) = 1 / sqrt((double) (2 * (n + 1)));
//         W_root_inv(i, j) = sqrt((double) (2 * (n + 1)));
//       } else if (i == j and i == n + 1) {
//         W_root(i, j) = 1 / sqrt((double) (2 * (n + 1)));
//         W_root_inv(i, j) = sqrt((double) (2 * (n + 1)));
//       }
//     }
//   }
//
//   RMatrix A = W_root * K * W_root;
//
//   CEigenvalues lambda;
//   EVreal(A, lambda);
//
//   return lambda.real();
// }
//
// void OrnsteinUhlenbeckProcess::InitGenerator(std::shared_ptr<Meshes> meshes, int init) {
//   size = init;
//   sampleVec.resize(init);
//   eigenvalues = CovarianceSpectrum(init);
// }
