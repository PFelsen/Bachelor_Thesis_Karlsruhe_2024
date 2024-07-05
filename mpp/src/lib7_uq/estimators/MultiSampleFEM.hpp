#ifndef MPP_MULTISAMPLEFEM_HPP
#define MPP_MULTISAMPLEFEM_HPP

#include "Sample.hpp"
#include "ScalarElement.hpp"

struct MultiSampleFEMConfig {
  int fLevel = 0;

  int cLevel = 0;

  int commSplit = 0;

  std::string protocolName{};

  PDESolverConfig pdeSolverConf;

  std::string quantityOfInterest{};

  MultiSampleFEMConfig(std::string protocolName, std::string quantityOfInterest) :
      protocolName(std::move(protocolName)), quantityOfInterest(std::move(quantityOfInterest)) {
  }

  MultiSampleFEMConfig() {
    if (!Config::IsInitialized()) return;

    Config::Get("Protocol", protocolName);
    Config::Get("Problem", protocolName); // Todo: settle on one name
    Config::Get("Quantity", quantityOfInterest);
  }

  MultiSampleFEMConfig WithFineLevel(int level) {
    fLevel = level;
    return *this;
  }

  MultiSampleFEMConfig WithProtocol(std::string protocol) {
    protocolName = std::move(protocol);
    return *this;
  }

  MultiSampleFEMConfig WithCoarseLevel(int level) {
    cLevel = level;
    return *this;
  }

  MultiSampleFEMConfig WithCommSplit(int communicationSplit) {
    commSplit = communicationSplit;
    return *this;
  }

  MultiSampleFEMConfig WithPDESolverConfig(const PDESolverConfig &conf) {
    pdeSolverConf = conf;
    return *this;
  }
};

class MultiSampleFEM { // Rename to MultiXFEM?
protected:
  int verbose = 0;

  std::string domainName{};

  MultiSampleFEMConfig conf;

  std::shared_ptr<Meshes> meshes;

  virtual SampleSolution protocol(const SampleID &id) = 0;

  virtual SampleSolution protocol(const SampleID &id, const SampleSolution &control) {return emptyProtocol(id); };

  virtual double protocolQOI(const SampleID &id, const SampleSolution &control) {return 0.0; };

  virtual SampleSolution emptyProtocol(const SampleID &id) = 0;

  virtual double selectQoI(const SampleSolution &sample) const {
    return sample.solution.values.at(conf.quantityOfInterest);
  }

public:
  explicit MultiSampleFEM(std::string domainName, MultiSampleFEMConfig conf) :
      domainName(std::move(domainName)), conf(std::move(conf)) {

    Config::Get("MSFEMVerbose", verbose);

    meshes = MeshesCreator(this->domainName)
        .WithPLevel(this->conf.cLevel)
        .WithLevel(this->conf.fLevel)
        .CreateShared();

    meshes->PrintInfo();
  }

  SampleSolution RunProtocol(const SampleID &id) {
    mout.StartBlock("MS-FEM");
    vout(1) << "Run Protocol " << Name() << endl;

    // Todo remove again MatrixGraph, Vectors and Meshes in SLEstimator
    // Todo how to incorporate this in the new design?
    (*meshes)[LevelPair(id.MeshIndex())];

    double start = MPI_Wtime();
    auto sample = protocol(id);
    sample.Q = selectQoI(sample);
    sample.C = MPI_Wtime() - start;

    mout.EndBlock(verbose == 0);
    return sample;
  }

    SampleSolution RunProtocol(const SampleID &id, const SampleSolution &control) {
        mout.StartBlock("MS-FEM");
        vout(1) << "Run Protocol " << Name() << endl;

        // Todo remove again MatrixGraph, Vectors and Meshes in SLEstimator
        (*meshes)[(LevelPair(id.MeshIndex()))];

        double start = MPI_Wtime();
        auto sample = protocol(id, control);
        sample.Q = selectQoI(sample);
        sample.C = MPI_Wtime() - start;

        mout.EndBlock(verbose == 0);
        return sample;
    }

    double RunProtocolQOI(const SampleID &id, const SampleSolution &control) {
      mout.StartBlock("MS-FEM");
      vout(1) << "Run Protocol " << Name() << endl;

      // Todo remove again MatrixGraph, Vectors and Meshes in SLEstimator
      (*meshes)[(LevelPair(id.MeshIndex()))];

      double start = MPI_Wtime();
      auto sample = protocolQOI(id, control);

      mout.EndBlock(verbose == 0);
      return sample;
    }

  void ClearMeshesOnMeshIndex(const LevelPair &meshIndex) { meshes->Erase(meshIndex); }

  SampleSolution EmptyProtocol(const SampleID &id) { return emptyProtocol(id); };

  const Meshes &GetMeshes() const { return *meshes; }

  virtual std::string Name() const = 0;

  virtual ~MultiSampleFEM() = default;
};

MultiSampleFEM *CreateMSFEM(const MultiSampleFEMConfig &conf);

std::unique_ptr<MultiSampleFEM> CreateUniqueMSFEM(const MultiSampleFEMConfig &conf);

std::shared_ptr<MultiSampleFEM> CreateSharedMSFEM(const MultiSampleFEMConfig &conf);

#endif //MPP_MULTISAMPLEFEM_HPP