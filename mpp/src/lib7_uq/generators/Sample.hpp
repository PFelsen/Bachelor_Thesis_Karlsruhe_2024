#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include "PDESolver.hpp"
#include <string>
#include <utility>


struct SampleID {
  std::string name = "";

  int commSplit{};

  int fLevel{};

  int cLevel{};

  int number{};

  bool coarse{};

  explicit SampleID(std::string name = "") : name(std::move(name)) {}

  SampleID(int level, int number, bool coarse, std::string name = "") :
      fLevel(level), cLevel((level > 0) ? (level - 1) : 0),
      number(number), coarse(coarse), name(std::move(name)) {}

  SampleID(int level, int number, bool coarse, int commSplit) :
      fLevel(level), cLevel((level > 0) ? (level - 1) : 0),
      number(number), coarse(coarse), commSplit(commSplit) {}

  SampleID WithCommSplit(int _commSplit) {
    this->commSplit = _commSplit;
    return *this;
  }

  SampleID WithName(const std::string &idName) {
    this->name = idName;
    return *this;
  }


  LevelPair MeshIndex() const {
#ifdef USE_SPACETIME
    if (coarse) return {cLevel, cLevel};
    else return {fLevel, fLevel};
#else
    if (coarse) return {cLevel, -1, 0, commSplit};
    else return {fLevel, -1, 0, commSplit};
#endif
  }

  std::string IdString() const {
    return name + "." + std::to_string(fLevel)
           + "." + std::to_string((int) coarse)
           + "." + std::to_string(number);
  }

  friend std::ostream &operator<<(std::ostream &s, const SampleID &id) {
    return s << id.IdString();
  }
};

struct SampleSolution {
  double Q = 0.0;  // Quantity of interest of sample solution

  double C = 0.0;  // Cost to compute sample solution

  double W = 1.0;  // Weight of the sample solution

  SampleID id;  // Identification of Sample

  Solution solution;

  SampleSolution(std::shared_ptr<const IDiscretization> disc, SampleID id, std::string name = "U") :
      id(std::move(id)), solution(std::move(disc), this->id.MeshIndex()) {
    this->id.name = std::move(name);
  }

  SampleSolution(const Vector &u, SampleID id, std::string name = "U") :
      solution(u), id(std::move(id)) {
    this->id.name = std::move(name);
  }

  SampleSolution(Solution solution, SampleID id, double weight = 1) :
      solution(std::move(solution)), id(std::move(id)), W(weight) {}

  const SampleID &Id() const {
    return id;
  }

  std::string IdString() const {
    return id.IdString();
  }

  friend std::ostream &operator<<(std::ostream &s, const SampleSolution &sol) {
    return s << sol.id << " Q=" << sol.Q << " C=" << sol.C << " W=" << sol.W;
  }
};

#endif //SAMPLE_HPP
