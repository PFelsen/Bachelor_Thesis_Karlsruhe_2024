#include <filesystem>

#include "Parallel.hpp"
#include "Quadrature.hpp"

#include "m++.hpp"

// TODO REMOVE THIS. This is not code which should be called at the start of the program
#include "lib4_fem/plot/VtuPlot.hpp"

Mpp *Mpp::mpp = nullptr;

Mpp &Mpp::initialize(int *argc, char **argv) {
  static MppGuard g;
  if (mpp) THROW("Mpp already initialized")
  mpp = new Mpp(argc, argv);
  return *mpp;
}

Mpp &Mpp::instance() {
  if (!mpp) THROW("Mpp not initialized")
  return *mpp;
}

// Todo change signature to
// Mpp::Mpp(int argc, char **argv)
Mpp::Mpp(int *argc, char **argv) {
  ParallelProgrammingModel::Initialize(*argc, argv);
  Config::Initialize(argc, argv);


  // Todo clean this stuff up, rather use filesystem; commented code below found in dead main
  bool clear = false;
  Config::Get("ClearData", clear);
  if (clear && PPM->master()) { int cleared = system("exec rm -rf data/vt*/*"); }
  //  bool clear = false;
  //  Config::Get("ClearData", clear);
  //  if (clear && PPM->master()) {
  //    // The artifacts are stored in the following structure:
  //    // data/
  //    // └── vtu/
  //    //     └── foo.vtu
  //    // └── vtk/
  //    //     └── bar.vtk
  //    // We remove everything within the directories in data/, but keep the
  //    // directory structure itself intact. This is necessary since the plotting
  //    // code does not verify the existence of data/vtu for example, and silently
  //    // fails without it.
  //    for (const auto &entry :
  //        std::filesystem::directory_iterator(Config::GetDataPath())) {
  //      if (entry.is_directory()) {
  //        for (const auto &entry_to_remove :
  //            std::filesystem::directory_iterator(entry)) {
  //          std::filesystem::remove_all(entry_to_remove.path());
  //        }
  //      }
  //    }
  //  }

  // TODO REMOVE THIS. This is not code which should be called at the start of the program
  Config::Get("ParallelPlotting", globalParallelPlotting);
  Config::Get("Compression", globalCompression);
}

Mpp::~Mpp() {
  clearQuadrature();
  ParallelProgrammingModel::Close();
}

Mpp::MppGuard::~MppGuard() {
  if (nullptr != Mpp::mpp) {
    delete Mpp::mpp;
    Mpp::mpp = nullptr;
  }
}
