#ifndef STPATHMANAGER_HPP
#define STPATHMANAGER_HPP

#include <vector>
#include "Mesh.hpp"
#include "STAssemble.hpp"

using Path = std::vector<LevelPair>;

class PathStrategy {
public:
  virtual Path createPath(LevelPair from, LevelPair to) const = 0;
  virtual std::string name() const = 0;

  virtual ~PathStrategy() {};
};

class AdaptiveToZeroDirectPathStrategy : public PathStrategy {
  Path createPath(LevelPair from, LevelPair to) const override {
    Path path{from};
    LevelPair current = from;
    while (current.adaptivityLevel > to.adaptivityLevel) {
      current = current.CoarserInAdaptivity();
      path.push_back(current);
    }
    while (current.space > to.space && current.time > to.time) {
      current = current.CoarserInBoth();
      path.push_back(current);
    }
    return path;
  }
  std::string name() const override{
    return "SpaceThenTimePathStrategy";
  }
};


class DirectPathStrategy : public PathStrategy {
    Path createPath(LevelPair from, LevelPair to) const override {
        Path path{from};
        LevelPair current = from;
        while (current.space > to.space && current.time > to.time) {
            current = current.CoarserInBoth();
            path.push_back(current);
        }
        return path;
    }
    std::string name() const override{
        return "SpaceThenTimePathStrategy";
    }
};

class SpaceThenTimePathStrategy : public PathStrategy {
  Path createPath(LevelPair from, LevelPair to) const override {
    Path path{from};
    LevelPair current = from;
    while (current.space > to.space || current.time > to.time) {
      if (current.space > to.space){
        current = current.CoarserInSpace();
      } else {
        current = current.CoarserInTime();
      }
      path.push_back(current);
    }
    return path;
  }
  std::string name() const override{
    return "SpaceThenTimePathStrategy";
  }
};

class TimeThenSpacePathStrategy : public PathStrategy {
  Path createPath(LevelPair from, LevelPair to) const override {
    Path path{from};
    LevelPair current = from;
    while (current.space > to.space || current.time > to.time) {
      if (current.time > to.time){
        current = current.CoarserInTime();
      } else {
        current = current.CoarserInSpace();
      }
      path.push_back(current);
    }
    return path;
  }
  std::string name() const override{
    return "TimeThenSpacePathStrategy";
  }
};

class OnlyTimePathStrategy : public PathStrategy {
  Path createPath(LevelPair from, LevelPair to) const override {
    Path path{from};
    LevelPair current = from;
    while (current.time > to.time) {
      current = current.CoarserInTime();
      path.push_back(current);
    }
    return path;
  }
  std::string name() const override{
    return "OnlyTimePathStrategy";
  }
};

class OnlySpacePathStrategy : public PathStrategy {
  Path createPath(LevelPair from, LevelPair to) const override {
    Path path{from};
    LevelPair current = from;
    while (current.space > to.space) {
      current = current.CoarserInSpace();
      path.push_back(current);
    }
    return path;
  }
  std::string name() const override{
    return "OnlySpacePathStrategy";
  }
};

Path makePathSpaceThenTime(LevelPair levels);

Path makePathTimeThenSpace(LevelPair levels);

Path makeTimePath(LevelPair levels);

Path makeTimePath(LevelPair levels, LevelPair first);

Path makeSpacePath(LevelPair levels);

Path makeSpacePath(LevelPair levels, LevelPair first);

std::unique_ptr<Preconditioner> createPreconditioner(const std::string &pathChoice,
                                                     STAssemble &assemble);

#endif //STPATHMANAGER_HPP
