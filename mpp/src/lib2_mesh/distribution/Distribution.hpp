#ifndef _DISTRIBUTION_H_
#define _DISTRIBUTION_H_

#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Cell.hpp"
#include "Point.hpp"

class Mesh;

class Distribution {
public:
  using Weights = std::unordered_map<Point, short>;
private:
  using MarkedCells = std::vector<std::list<cell>>;

  int verbose = 0;

  int commSplit = 0;

  std::string _distName = "RCB";

  std::vector<Point> distributedCells{};

  std::vector<cell> cells;

  std::vector<int> procsForCells;

  MarkedCells markedCells;

  std::shared_ptr<Weights> weights = std::make_shared<Weights>();

  bool distribute = true;

  std::shared_ptr<const Distribution> spaceDist;

  void init(Mesh &mesh);

  bool clearDistributed = true;

  void emptyData(bool clearAll = true);

  double getMidOfLowestAndHighestCoord(int begin, int end, int coordinate);

  int setProcsForCellsAndReturnMid(int rekDepth, int begin, int end);
public:
  explicit Distribution(int commSplit = 0);

  explicit Distribution(std::string distName, int commSplit = 0);

  Distribution(string distName, std::shared_ptr<const Distribution> spaceDistribution,
               int commSplit = 0);

  Distribution(std::string distName, std::shared_ptr<Weights> weights, int commSplit = 0);

  std::string Name() const;

  int CommSplit() const;

  void Clear() { emptyData(true); }

  ~Distribution() { emptyData(); }

  void DistributeMesh(Mesh &mesh);

  // Todo DistributeMesh should swallow DistributeSTMesh
  // Todo introduce enums for Distributions
  // Todo remove Less_x functions and find occurring bug in RCBGeom
  void DistributeSTMesh(Mesh &mesh);

  void CommunicateIdentifySets(Mesh &mesh);

  void Communicate(Mesh &mesh);

  void RIB(int depth);

  void Stripes_x(int rekDepth, int begin, int end);

  void Stripes_y(int rekDepth, int begin, int end);

  void Stripes_z(int rekDepth, int begin, int end);

  //  void Stripes_t(int rekDepth, int begin, int end);

  void RCB_x(int rekDepth, int begin, int end);

  void RCB_y(int rekDepth, int begin, int end);

  void RCB_z(int rekDepth, int begin, int end);

  //  void RCB_t(int rekDepth, int begin, int end);

  void RCBGeom_x(int rekDepth, int begin, int end);

  void RCBGeom_y(int rekDepth, int begin, int end);

  void RCBGeom_z(int rekDepth, int begin, int end);

  void RCB2D_x(int rekDepth, int begin, int end);

  void RCB2D_y(int rekDepth, int begin, int end);

  // Todo remove the functions below
  void Stripes_Old();

  void RCB_Old(int rekDepth, int begin, int end);

  void RCB_Geom_Old(int rekDepth, int begin, int end, char d);

  // Todo refactor the functions below
  void trcb(int level, int begin, int end, char d);

  void trcb_space(int level, int begin, int end, char d);

  void timeOpt(int tProc, int level, int begin, int end);

  void timercb(int level, int begin, int end, char d);

  void timercb_weighted(int tlevel, int xlevel, int ylevel, int begin, int end, char d, Weights &W,
                        int default_weight = 1);

  void time_stripes(int n_procs, int begin, int end);
};

#endif // of #ifndef _DISTRIBUTION_H_
