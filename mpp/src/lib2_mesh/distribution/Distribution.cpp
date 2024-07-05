#include "Distribution.hpp"

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <valarray>

#include "Mesh.hpp"
#include "Parallel.hpp"
#include "HermCMatrix.hpp"
#include "Spectrum.hpp"

Distribution::Distribution(int commSplit) :
    markedCells(PPM->Size(commSplit)), commSplit(commSplit) {
  Config::Get("DistributionVerbose", verbose);
  Config::Get("Distribution", _distName);
  Config::Get("ClearDistribution", clearDistributed);
}

Distribution::Distribution(std::string distName, int commSplit) :
    markedCells(PPM->Size(commSplit)), _distName(std::move(distName)), commSplit(commSplit) {
  Config::Get("DistributionVerbose", verbose);
  Config::Get("ClearDistribution", clearDistributed);
}

Distribution::Distribution(std::string distName,
                           std::shared_ptr<const Distribution> spaceDistribution, int commSplit) :
    markedCells(PPM->Size(commSplit)), _distName(std::move(distName)),
    spaceDist(std::move(spaceDistribution)), commSplit(commSplit) {
  Config::Get("DistributionVerbose", verbose);
  Config::Get("ClearDistribution", clearDistributed);
}

Distribution::Distribution(std::string distName, std::shared_ptr<Weights> weights, int commSplit) :
    markedCells(PPM->Size(commSplit)), _distName(std::move(distName)), weights(std::move(weights)),
    commSplit(commSplit) {
  Config::Get("DistributionVerbose", verbose);
  Config::Get("ClearDistribution", clearDistributed);
}

void Distribution::DistributeMesh(Mesh &mesh) {
  mout.StartBlock("Distribute");
  int startDist = (int) distributedCells.size();
  init(mesh);
  int rekDepth = int(log2(PPM->Size(commSplit)));

  if (mesh.IsSTMesh()) {
    DistributeSTMesh(mesh);
  } else {
    if (_distName == "Stripes" || _distName == "Stripes_x")
      Stripes_x(rekDepth, startDist, cells.size());
    else if (_distName == "Stripes_y") Stripes_y(rekDepth, startDist, cells.size());
    else if (_distName == "Stripes_z") Stripes_z(rekDepth, startDist, cells.size());

      // Recursive coordinate bisection
    else if (_distName == "RCB" || _distName == "RCBx" || _distName == "RCBNoTime")
      RCB_x(rekDepth, startDist, cells.size());
    else if (_distName == "RCBy") RCB_y(rekDepth, startDist, cells.size());
    else if (_distName == "RCBz") RCB_z(rekDepth, startDist, cells.size());

    else if (_distName == "RCBG" || _distName == "RCBGeom_x" || _distName == "RCBGx")
      RCBGeom_x(rekDepth, startDist, cells.size());
    else if (_distName == "RCBGy" || _distName == "RCBGeom_y")
      RCBGeom_y(rekDepth, startDist, cells.size());
    else if (_distName == "RCBGz" || _distName == "RCBGeom_z")
      RCBGeom_z(rekDepth, startDist, cells.size());

      // Dimension could be considered one time,
      // for now it does the same as RCB.
    else if (_distName == "RCB2D" || _distName == "RCB2Dx")
      RCB2D_x(rekDepth, startDist, cells.size());
    else if (_distName == "RCB2Dy") RCB2D_y(rekDepth, startDist, cells.size());

    else if (_distName == "KAFFPA") THROW("Kaffpa is not included anymore")
      // See https://git.scc.kit.edu/mpp/mpp/-/blob/2.0.0/src/Distribution.C#L1517
      // See https://git.scc.kit.edu/mpp/mpp/-/issues/124
      // See https://github.com/KaHIP/KaHIP

      // Recursive inertia bisection
    else if (_distName == "RIB") {
      RIB(rekDepth);
      // See https://git.scc.kit.edu/mpp/mpp/-/blob/2.0.0/src/Distribution.C#L1587

    }


      // These cases should be removed at some point
    else if (_distName == "Stripes_Old") Stripes_Old();
    else if (_distName == "RCB_Old" || _distName == "RCB2")
      RCB_Old(rekDepth, startDist, cells.size());
    else if (_distName == "RCB_Geom_Old") RCB_Geom_Old(rekDepth, startDist, cells.size(), 'x');

    else if (_distName == "NoDistribution") distribute = false;
    else Exit(_distName + " not implemented")
  }

  if (distribute) {

    for (int i = startDist; i < cells.size(); ++i) {
      distributedCells.emplace_back(cells[i]());
    }

    for (int i = 0; i < cells.size(); ++i)
      markedCells[procsForCells[i]].push_back(cells[i]);

    Communicate(mesh);
  }

  int cellCount = mesh.CellCount();
  pvout(1) << cellCount << " cells on logging process. " << mesh.CellCountGeometry()
           << " cells on all processes." << endl;

  if (PPM->Min(cellCount, commSplit) == 0)
    VerboseWarning("At least one processor has no cells!", 1)

  emptyData(clearDistributed);
  mout.EndBlock(true);
}

void Distribution::CommunicateIdentifySets(Mesh &mesh) {
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    if (!is.master()) continue;
    procset p = mesh.find_procset(is());
    if (p == mesh.procsets_end()) continue;
    for (int j = 0; j < is.size(); ++j) {
      procset q = mesh.find_procset(is[j]);
      if (q == mesh.procsets_end()) continue;
      for (int i = 0; i < q.size(); ++i)
        mesh.procSets.Append(is(), q[i]);
    }
    for (int j = 0; j < is.size(); ++j)
      mesh.procSets.Copy(p, is[j]);
  }
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    if (is.master()) continue;
    procset p = mesh.find_procset(is());
    if (p == mesh.procsets_end()) continue;
    for (int j = 0; j < is.size(); ++j) {
      procset q = mesh.find_procset(is[j]);
      if (q == mesh.procsets_end()) continue;
      for (int i = 0; i < q.size(); ++i)
        mesh.procSets.Append(is(), q[i]);
    }
    for (int j = 0; j < is.size(); ++j)
      mesh.procSets.Copy(p, is[j]);
  }
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    if (!is.master()) continue;
    procset p = mesh.find_procset(is());
    for (int j = 0; j < is.size(); ++j)
      mesh.procSets.Copy(p, is[j]);
  }
  ExchangeBuffer bufferIdentifySets(commSplit);
  for (identifyset is = mesh.identifysets(); is != mesh.identifysets_end(); ++is) {
    procset p = mesh.find_procset(is());
    if (p == mesh.procsets_end()) continue;
    for (int i = 0; i < p.size(); ++i) {
      if (p[i] == PPM->Proc(commSplit)) continue;
      bufferIdentifySets.Send(p[i]) << short(is.size()) << is();
      for (int j = 0; j < is.size(); ++j)
        bufferIdentifySets.Send(p[i]) << is[j];
    }
  }
  mesh.procSets.RemoveSingle();
  for (short q = 0; q < PPM->Size(commSplit); ++q) {
    if (PPM->Proc(commSplit) == q) continue;
    if (markedCells[q].empty()) continue;
    bufferIdentifySets.Send(q) << short(0);
  }
  bufferIdentifySets.Communicate();

  for (short q = 0; q < PPM->Size(commSplit); ++q) {
    if (bufferIdentifySets.Receive(q).Size() == 0) continue;
    short m;
    bufferIdentifySets.Receive(q) >> m;
    while (m) {
      Point x;
      bufferIdentifySets.Receive(q) >> x;
      for (int i = 0; i < m; ++i) {
        Point y;
        bufferIdentifySets.Receive(q) >> y;
        mesh.identifySets.Insert(x, y);
      }
      bufferIdentifySets.Receive(q) >> m;
    }
  }
}

void Distribution::Communicate(Mesh &mesh) {
  vector<std::list<bnd_face>> MarkedBnd(PPM->Size(commSplit));
  for (short q = 0; q < PPM->Size(commSplit); ++q) {
    for (cell &c: markedCells[q]) {
      for (int i = 0; i < c.Corners(); ++i) {
        mesh.procSets.Add(c.Corner(i), q);
      }
      for (int i = 0; i < c.Edges(); ++i) {
        mesh.procSets.Add(c.Edge(i), q);
      }
      for (int i = 0; i < c.Faces(); ++i) {
        Point z = c.Face(i);
        mesh.procSets.Add(z, q);
        bnd_face b = mesh.find_bnd_face(z);
        if (b != mesh.bnd_faces_end()) { MarkedBnd[q].push_back(b); }
      }
    }
  }

  CommunicateIdentifySets(mesh);

  ExchangeBuffer exBuffer(commSplit);
  for (procset p = mesh.procsets(); p != mesh.procsets_end(); ++p)
    for (int i = 0; i < p.size(); ++i) {
      if (p[i] != PPM->Proc(commSplit)) {
        short m = 1;
        face f = mesh.find_face(p());
        if (f != mesh.faces_end()) ++m;
        exBuffer.Send(p[i]) << m << p() << *p;
        if (m > 1) exBuffer.Send(p[i]) << *f;
      }
    }
  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (PPM->Proc(commSplit) == q) continue;
    if (markedCells[q].empty()) continue;
    exBuffer.Send(q) << short(0);
    exBuffer.Send(q) << int(MarkedBnd[q].size());
    for (bnd_face b: MarkedBnd[q]) {
      exBuffer.Send(q) << b;
    }
    for (cell &c: markedCells[q]) {
      exBuffer.Send(q) << c;
#ifdef USE_DATAMESH
      for (int i = 0; i < c.Corners(); i++) {
        exBuffer.Send(q) << mesh.find_vertex(c.Corner(i)).GetData();
      }
#endif
      mesh.removeCell(c);
    }
  }

  exBuffer.Communicate();
  for (short q = 0; q < PPM->Size(commSplit); ++q) {
    if (exBuffer.Receive(q).Size() == 0) continue;
    short m;
    exBuffer.Receive(q) >> m;
    while (m) {
      Point z;
      ProcSet P(mesh.CommSplit());
      exBuffer.Receive(q) >> z >> P;
      mesh.procSets.Insert(z, P);
      if (m > 1) {
        Face F;
        exBuffer.Receive(q) >> F;
        mesh.insertFace(z, F);
      }
      exBuffer.Receive(q) >> m;
    }
    int n;
    exBuffer.Receive(q) >> n;
    for (int i = 0; i < n; ++i) {
      exBuffer.Receive(q) >> mesh.bnd_ref();
    }
    while (exBuffer.Receive(q).size() < exBuffer.Receive(q).Size()) {
      exBuffer.Receive(q) >> mesh;
    }
  }
  exBuffer.ClearBuffers();
  mesh.procSets.Clean();
  for (int q = 0; q < PPM->Size(commSplit); ++q)
    mesh.procSets.Add(Infty, q);
  mesh.finishParallel();
}

int Distribution::setProcsForCellsAndReturnMid(int rekDepth, int begin, int end) {
  int mid = begin + (end - begin) / 2;
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1); // 2^(rekDepth - 1)
  return mid;
}

bool Less_x(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P[0] < Q[0] - GeometricTolerance) return true;
  if (P[0] > Q[0] + GeometricTolerance) return false;
  if (P[1] < Q[1] - GeometricTolerance) return true;
  if (P[1] > Q[1] + GeometricTolerance) return false;
  if (P[2] < Q[2] - GeometricTolerance) return true;
  return false;
}

bool Less_y(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P[1] < Q[1] - GeometricTolerance) return true;
  if (P[1] > Q[1] + GeometricTolerance) return false;
  if (P[2] < Q[2] - GeometricTolerance) return true;
  if (P[2] > Q[2] + GeometricTolerance) return false;
  if (P[0] < Q[0] - GeometricTolerance) return true;
  return false;
}

bool Less_z(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P[2] < Q[2] - GeometricTolerance) return true;
  if (P[2] > Q[2] + GeometricTolerance) return false;
  if (P[0] < Q[0] - GeometricTolerance) return true;
  if (P[0] > Q[0] + GeometricTolerance) return false;
  if (P[1] < Q[1] - GeometricTolerance) return true;
  return false;
}

bool Less_t(const cell &c0, const cell &c1) { return (c0() < c1()); }

bool TLess_x(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P[0] < Q[0] - GeometricTolerance) return true;
  if (P[0] > Q[0] + GeometricTolerance) return false;
  if (P[1] < Q[1] - GeometricTolerance) return true;
  if (P[1] > Q[1] + GeometricTolerance) return false;
  if (P.t() < Q.t() - GeometricTolerance) return true;
  if (P.t() > Q.t() + GeometricTolerance) return false;
  if (P[2] < Q[2] - GeometricTolerance) return true;
  if (P[2] > Q[2] + GeometricTolerance) return false;
  return false;
}

bool TLess_y(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P[1] < Q[1] - GeometricTolerance) return true;
  if (P[1] > Q[1] + GeometricTolerance) return false;
  if (P[0] < Q[0] - GeometricTolerance) return true;
  if (P[0] > Q[0] + GeometricTolerance) return false;
  if (P.t() < Q.t() - GeometricTolerance) return true;
  if (P.t() > Q.t() + GeometricTolerance) return false;
  if (P[2] < Q[2] - GeometricTolerance) return true;
  if (P[2] > Q[2] + GeometricTolerance) return false;
  return false;
}

bool TLess_z(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P[2] < Q[2] - GeometricTolerance) return true;
  if (P[2] > Q[2] + GeometricTolerance) return false;
  if (P[0] < Q[0] - GeometricTolerance) return true;
  if (P[0] > Q[0] + GeometricTolerance) return false;
  if (P[1] < Q[1] - GeometricTolerance) return true;
  if (P[1] > Q[1] + GeometricTolerance) return false;
  if (P.t() < Q.t() - GeometricTolerance) return true;
  if (P.t() > Q.t() + GeometricTolerance) return false;
  return false;
}

bool TLess_t(const cell &c0, const cell &c1) {
  const Point &P = c0();
  const Point &Q = c1();
  if (P.t() < Q.t() - GeometricTolerance) return true;
  if (P.t() > Q.t() + GeometricTolerance) return false;
  if (P[0] < Q[0] - GeometricTolerance) return true;
  if (P[0] > Q[0] + GeometricTolerance) return false;
  if (P[1] < Q[1] - GeometricTolerance) return true;
  if (P[1] > Q[1] + GeometricTolerance) return false;
  if (P[2] < Q[2] - GeometricTolerance) return true;
  if (P[2] > Q[2] + GeometricTolerance) return false;
  return false;
}

std::function<bool(const cell &, const cell &)> GetLess(char var) {
  switch (var) {
    case 'x':
      return Less_x;
    case 'y':
      return Less_y;
    case 'z':
      return Less_z;
    case 't':
      return Less_t;
    default:
      THROW("Dimension " + std::to_string(var) + " not known")
  }
}

void Distribution::RCB_x(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_x);
  int mid = setProcsForCellsAndReturnMid(rekDepth, begin, end);
  RCB_y(rekDepth - 1, begin, mid);
  RCB_y(rekDepth - 1, mid, end);
}

void Distribution::RCB_y(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_y);
  int mid = setProcsForCellsAndReturnMid(rekDepth, begin, end);
  RCB_z(rekDepth - 1, begin, mid);
  RCB_z(rekDepth - 1, mid, end);
}

void Distribution::RCB_z(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_z);
  int mid = setProcsForCellsAndReturnMid(rekDepth, begin, end);
  RCB_x(rekDepth - 1, begin, mid);
  RCB_x(rekDepth - 1, mid, end);
}

// void Distribution::RCB_t(int rekDepth, int begin, int end) {
//   if (rekDepth == 0) return;
//   sort(cells.begin() + begin, cells.begin() + end, Less_t);
//   int mid = setProcsForCellsAndReturnMid(rekDepth, begin, end);
//   RCB_t(rekDepth - 1, begin, mid);
//   RCB_t(rekDepth - 1, mid, end);
// }

void Distribution::Stripes_x(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_x);
  int mid = begin + (end - begin) / 2;
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCB_x(rekDepth - 1, begin, mid);
  RCB_x(rekDepth - 1, mid, end);
}

void Distribution::Stripes_y(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_y);
  int mid = begin + (end - begin) / 2;
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCB_y(rekDepth - 1, begin, mid);
  RCB_y(rekDepth - 1, mid, end);
}

void Distribution::Stripes_z(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_z);
  int mid = begin + (end - begin) / 2;
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCB_z(rekDepth - 1, begin, mid);
  RCB_z(rekDepth - 1, mid, end);
}

// void Distribution::Stripes_t(int rekDepth, int begin, int end) {
//   if (rekDepth == 0) return;
//   sort(cells.begin() + begin, cells.begin() + end, Less_t);
//   int mid = begin + (end - begin) / 2;
//   for (int i = mid; i < end; ++i) procsForCells[i] += 1 << (rekDepth - 1);
//   RCB_t(rekDepth - 1, begin, mid);
//   RCB_t(rekDepth - 1, mid, end);
// }


void Distribution::RCB2D_x(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_x);
  int mid = begin + (end - begin) / 2;
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCB2D_y(rekDepth - 1, begin, mid);
  RCB2D_y(rekDepth - 1, mid, end);
}

void Distribution::RCB2D_y(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_y);
  int mid = begin + (end - begin) / 2;
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCB2D_x(rekDepth - 1, begin, mid);
  RCB2D_x(rekDepth - 1, mid, end);
}

double Distribution::getMidOfLowestAndHighestCoord(int begin, int end, int coordinate) {
  double low = infty;
  double high = -infty;
  for (int i = begin; i < end; ++i) {
    if (low > cells[i]()[coordinate]) low = cells[i]()[coordinate];
    if (high < cells[i]()[coordinate]) high = cells[i]()[coordinate];
  }
  return (low + high) / 2.0;
}

void Distribution::RCBGeom_x(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_x);
  double midlh = getMidOfLowestAndHighestCoord(begin, end, 0);
  int mid = end;
  for (int i = begin; i < end; ++i) {
    bool con = false;
    if (midlh > cells[i]()[0]) con = true;
    if (con) continue;
    mid = i;
    break;
  }
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCBGeom_y(rekDepth - 1, begin, mid);
  RCBGeom_y(rekDepth - 1, mid, end);
}

void Distribution::RCBGeom_y(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_y);
  auto midlh = getMidOfLowestAndHighestCoord(begin, end, 1);
  int mid = end;
  for (int i = begin; i < end; ++i) {
    bool con = false;
    if (midlh > cells[i]()[1]) con = true;
    if (con) continue;
    mid = i;
    break;
  }
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCBGeom_x(rekDepth - 1, begin, mid);
  RCBGeom_x(rekDepth - 1, mid, end);
  //    RCBGeom_z(rekDepth - 1, begin, mid); <- Was used in old implementation
  //    RCBGeom_z(rekDepth - 1, mid, end);
}

void Distribution::RCBGeom_z(int rekDepth, int begin, int end) {
  if (rekDepth == 0) return;
  sort(cells.begin() + begin, cells.begin() + end, Less_z);
  auto midlh = getMidOfLowestAndHighestCoord(begin, end, 2);
  int mid = end;
  for (int i = begin; i < end; ++i) {
    bool con = false;
    if (midlh > cells[i]()[2]) con = true;
    if (con) continue;
    mid = i;
    break;
  }
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (rekDepth - 1);
  RCBGeom_x(rekDepth - 1, begin, mid);
  RCBGeom_x(rekDepth - 1, mid, end);
}

std::string Distribution::Name() const { return _distName; }

void Distribution::init(Mesh &mesh) {
  // Ensure correct numeration of previously distributed cells
  for (auto p: distributedCells) {
    cells.emplace_back(mesh.find_cell(p));
  }

  // Add new cells from mesh
  for (cell c = mesh.cells(); c != mesh.cells_end(); ++c) {
    if (std::none_of(distributedCells.begin(), distributedCells.end(),
                     [&c](const Point &p) { return c() == p; })) {
      cells.emplace_back(c);
      procsForCells.push_back(0);
      (*weights)[c()] = 1;
    }
  }
}

void Distribution::emptyData(bool clearAll) {
  cells.clear();
  markedCells.clear();
  if (clearAll) {
    markedCells.clear();
    procsForCells.clear();
    distributedCells.clear();
  }
  markedCells.resize(PPM->Size(commSplit));
}

// STARTING RIB
namespace RIB_DIST {

  inline bool common_face(const cell &c, const cell &d) {
    for (int i = 0; i < c.Faces(); ++i)
      for (int j = 0; j < d.Faces(); ++j)
        if (c.Face(i) == d.Face(j)) return true;
    return false;
  }


  inline void fiedler(const vector<cell> &C, std::valarray<double> &F) {
    int N = int(C.size());
    SymRMatrix G(0.0, N);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < i; ++j) {
        if (common_face(C[i], C[j])) {
          // Symmetric entry is inserted by default for SymRMatrix
          G(-1.0, i, j);
        }
      }
    }
    for (int i = 0; i < N; ++i) {
      int sum = 0;
      for (int j = 0; j < N; ++j) {
        if (G(i, j) == -1.0) ++sum;
      }
      G(sum, i, i);
    }
    Eigenvector lambda(0.0, N);
    Eigenvectors E(N);

    // maybe just calc the second EV
    EVreal(G, lambda, E);

    for (int i = 0; i < N; ++i) {
      F[i] = E[i][1];
    }
  }


  inline int pow_int(int b, int e) {
    int p = 1;
    while (e > 1) {
      p *= b;
      --e;
    }
    return p;
  }


  using Bundle = std::tuple<cell, int, double>;

  inline bool bundle_compare(const Bundle &a, const Bundle &b) {
    return std::get<2>(a) <= std::get<2>(b);
  }

  void rib(int level, vector<int> &D, size_t begin, size_t end, vector<cell> &C) {
    size_t N = end - begin;
    vector<cell> tmp(N);
    for (int i = 0; i < N; ++i)
      tmp[i] = C[begin + i];

    std::valarray<double> F(N);
    fiedler(tmp, F);
    vector<Bundle> B(N);
    for (size_t i = 0; i < N; ++i)
      B[i] = std::make_tuple(C[begin + i], D[begin + i], F[i]);
    sort(B.begin(), B.end(), bundle_compare);
    for (size_t i = 0; i < N; ++i) {
      C[begin + i] = std::get<0>(B[i]);
      D[begin + i] = std::get<1>(B[i]);
      F[i] = std::get<2>(B[i]);
    }
    if (level == 0) return;
    else {
      size_t mid = begin + (end - begin) / 2;
      rib(level - 1, D, begin, mid, C);
      for (size_t i = mid; i < end; ++i)
        D[i] += pow_int(2, level);
      rib(level - 1, D, mid, end, C);
    }
  }
}

void Distribution::RIB(int depth) {
  std::fill(procsForCells.begin(), procsForCells.end(), 0);
  Date StartF;
  if (PPM->master()) RIB_DIST::rib(depth, procsForCells, 0, cells.size(), cells);
  vout(2) << "find fiedler vectors " << Date() - StartF << "\n";
}
// ENDING RIB


void Distribution::Stripes_Old() {
  sort(cells.begin(), cells.end(), Less_x);
  int m = (int(cells.size()) + PPM->Size(commSplit) - 1) / PPM->Size(commSplit);
  int r = int(cells.size()) - (m - 1) * PPM->Size(commSplit);
  int k = 0;
  int q = 0;
  for (int i = 0; i < cells.size(); ++i) {
    if (k >= m) {
      k = 0;
      ++q;
      if (q == r) --m;
    }
    ++k;
    procsForCells[i] = q;
  }
}

void Distribution::RCB_Old(int rekDepth, int begin, int end) {
  double dist_x = 0;
  double dist_y = 0;
  double dist_z = 0;
  for (int i = begin; i < end; ++i) {
    for (int j = begin; j < end; ++j) {
      if (dist_x < abs(cells[i]()[0] - cells[j]()[0])) dist_x = abs(cells[i]()[0] - cells[j]()[0]);
      if (dist_y < abs(cells[i]()[1] - cells[j]()[1])) dist_y = abs(cells[i]()[1] - cells[j]()[1]);
      if (dist_z < abs(cells[i]()[2] - cells[j]()[2])) dist_z = abs(cells[i]()[2] - cells[j]()[2]);
    }
  }
  if (dist_x < dist_y)
    if (dist_z < dist_y) sort(cells.begin() + begin, cells.begin() + end, LessCell_y);
    else sort(cells.begin() + begin, cells.begin() + end, LessCell_z);
  else if (dist_z < dist_x) sort(cells.begin() + begin, cells.begin() + end, LessCell_x);
  else sort(cells.begin() + begin, cells.begin() + end, LessCell_z);

  if (rekDepth == 0) {
    return;
  } else {
    int mid = begin + (end - begin) / 2;
    for (int i = mid; i < end; ++i)
      procsForCells[i] += 1 << (rekDepth - 1);
    RCB_Old(rekDepth - 1, begin, mid);
    RCB_Old(rekDepth - 1, mid, end);
  }
}

void Distribution::RCB_Geom_Old(int rekDepth, int begin, int end, char d) {
  switch (d) {
    case 'x':
      sort(cells.begin() + begin, cells.begin() + end, Less_x);
      break;
    case 'y':
      sort(cells.begin() + begin, cells.begin() + end, Less_y);
      break;
    case 'z':
      sort(cells.begin() + begin, cells.begin() + end, Less_z);
      break;
    default:
      break;
  }
  if (rekDepth == 0) {
    return;
  } else {
    int mid;
    double low = infty;
    double high = -infty;
    for (int i = begin; i < end; ++i) {
      switch (d) {
        case 'x':
          if (low > cells[i]()[0]) low = cells[i]()[0];
          if (high < cells[i]()[0]) high = cells[i]()[0];
          break;
        case 'y':
          if (low > cells[i]()[1]) low = cells[i]()[1];
          if (high < cells[i]()[1]) high = cells[i]()[1];
          break;
        case 'z':
          if (low > cells[i]()[2]) low = cells[i]()[2];
          if (high < cells[i]()[2]) high = cells[i]()[2];
          break;
        default:
          break;
      }
    }
    double midlh = (high + low) / 2.;
    mid = end;
    for (int i = begin; i < end; ++i) {
      bool con = false;
      switch (d) {
        case 'x':
          if (midlh > cells[i]()[0]) con = true;
          break;
        case 'y':
          if (midlh > cells[i]()[1]) con = true;
          break;
        case 'z':
          if (midlh > cells[i]()[2]) con = true;
          break;
        default:
          break;
      }
      if (con) continue;
      mid = i;
      break;
    }
    for (int i = mid; i < end; ++i)
      procsForCells[i] += 1 << (rekDepth - 1);
    switch (d) {
      case 'x':
        RCB_Geom_Old(rekDepth - 1, begin, mid, 'y');
        RCB_Geom_Old(rekDepth - 1, mid, end, 'y');
        break;
      case 'y':
      case 'z':
        RCB_Geom_Old(rekDepth - 1, begin, mid, 'x');
        RCB_Geom_Old(rekDepth - 1, mid, end, 'x');
        break;
      default:
        break;
    }
  }
}

void Distribution::trcb(int level, int begin, int end, char d) {
  int N = end - begin;
  switch (d) {
    case 'x':
    case 'y':
    case 'z':
      sort(cells.begin() + begin, cells.begin() + end, TLess_t);
    case 'a':
    case 'c':
    case 'e':
      sort(cells.begin() + begin, cells.begin() + end, TLess_x);
      break;
    case 'b':
    case 'd':
    case 'f':
      sort(cells.begin() + begin, cells.begin() + end, TLess_y);
      break;
    default:
      break;
  }
  if (level == 0) {
    return;
  } else {
    int mid = begin + (end - begin) / 2;
    for (int i = mid; i < end; ++i)
      procsForCells[i] += 1 << (level - 1);
    switch (d) {
      case 'x':
        trcb(level - 1, begin, mid, 'y');
        trcb(level - 1, mid, end, 'y');
        break;
      case 'y':
        trcb(level - 1, begin, mid, 'z');
        trcb(level - 1, mid, end, 'z');
        break;
      case 'z':
        trcb(level - 1, begin, mid, 'a');
        trcb(level - 1, mid, end, 'a');
        break;
      case 'a':
        trcb(level - 1, begin, mid, 'b');
        trcb(level - 1, mid, end, 'b');
        break;
      case 'b':
        trcb(level - 1, begin, mid, 'c');
        trcb(level - 1, mid, end, 'c');
        break;
      case 'c':
        trcb(level - 1, begin, mid, 'd');
        trcb(level - 1, mid, end, 'd');
        break;
      case 'd':
        trcb(level - 1, begin, mid, 'e');
        trcb(level - 1, mid, end, 'e');
        break;
      case 'e':
        trcb(level - 1, begin, mid, 'f');
        trcb(level - 1, mid, end, 'f');
        break;
      case 'f':
        trcb(level - 1, begin, mid, 'x');
        trcb(level - 1, mid, end, 'x');
        break;
      default:
        break;
    }
  }
}

void Distribution::trcb_space(int level, int begin, int end, char d) {
  int N = end - begin;
  switch (d) {
    case 'x':
      sort(cells.begin() + begin, cells.begin() + end, TLess_x);
      break;
    case 'y':
      sort(cells.begin() + begin, cells.begin() + end, TLess_y);
      break;
    default:
      break;
  }
  if (level == 0) {
    return;
  } else {
    int mid = begin + (end - begin) / 2;
    for (int i = mid; i < end; ++i)
      procsForCells[i] += 1 << (level - 1);
    switch (d) {
      case 'x':
        trcb_space(level - 1, begin, mid, 'y');
        trcb_space(level - 1, mid, end, 'y');
        break;
      case 'y':
        trcb_space(level - 1, begin, mid, 'x');
        trcb_space(level - 1, mid, end, 'x');
        break;
      default:
        break;
    }
  }
}

void Distribution::timeOpt(int tProc, int level, int begin, int end) {
  int N = end - begin;
  sort(cells.begin() + begin, cells.begin() + end, TLess_t);

  if (PPM->Proc(commSplit) == 0) {
    vector<int> P(tProc);
    for (int i = 0; i < tProc; ++i)
      P[i] = ((end - 1) / tProc) * (i + 1);
    P[tProc - 1] = end - 1;

    for (int i = 0; i < tProc - 1; ++i) {
      int old_mid = P[i];
      int r_mid, l_mid;
      Point a = cells[P[i]]().CopyWithT(0.0);
      Point b = cells[P[i] - 1]().CopyWithT(0.0);

      while (cells[P[i]]().t() == cells[P[i] - 1]().t() || a == b) {
        ++P[i];
        a = cells[P[i]]().CopyWithT(0.0);
        b = cells[P[i] - 1]().CopyWithT(0.0);
      }
      r_mid = P[i];

      P[i] = old_mid;
      a = cells[P[i]]().CopyWithT(0.0);
      b = cells[P[i] - 1]().CopyWithT(0.0);

      while (cells[P[i]]().t() == cells[P[i] - 1]().t() || a == b) {
        --P[i];
        a = cells[P[i]]().CopyWithT(0.0);
        b = cells[P[i] + 1]().CopyWithT(0.0);
        if (P[i] == 0) break;
      }
      l_mid = P[i];

      if (old_mid - l_mid < r_mid - old_mid) P[i] = l_mid - 1;
      else P[i] = r_mid - 1;
    }

    if (level > 0) trcb_space(level, 0, P[0], 'x');

    int spaceProc = int(std::pow(2.0, double(level)) + 1e-10);

    if (tProc > 1)
      for (int i = 1; i < tProc; ++i) {
        for (int j = P[i - 1] + 1; j <= P[i]; ++j)
          procsForCells[j] += i * spaceProc;
        if (level > 0) trcb_space(level, P[i - 1] + 1, P[i], 'x');
      }
  }
}

void Distribution::timercb(int level, int begin, int end, char d) {
  if (cells.size() == 0) return;
  int N = end - begin;
  sort(cells.begin() + begin, cells.begin() + end, GetLess(d));
  if (level == 0) return;
  int mid = begin + (end - begin) / 2;
  if (PPM->Proc(commSplit) == 0) {
    int old_mid = mid;
    int r_mid, l_mid;
    Point a = cells[mid]().CopyWithT(0.0);
    Point b = cells[mid - 1]().CopyWithT(0.0);
    while ((mid + 1 < cells.size()) && (cells[mid]().t() == cells[mid - 1]().t() || a == b)) {
      ++mid;
      a = cells[mid]().CopyWithT(0.0);
      b = cells[mid - 1]().CopyWithT(0.0);
    }
    if (mid == end) --mid;
    // if (mid == cells.size()) --mid;
    r_mid = mid;
    mid = old_mid;
    a = cells[mid]().CopyWithT(0.0);
    b = cells[mid - 1]().CopyWithT(0.0);
    while ((mid > 0) && (cells[mid]().t() == cells[mid - 1]().t() || a == b)) {
      --mid;
      a = cells[mid]().CopyWithT(0.0);
      b = cells[mid + 1]().CopyWithT(0.0);
    }
    if (mid == begin) ++mid;
    // if (mid == 0) ++mid;
    l_mid = mid;
    if (old_mid - l_mid < r_mid - old_mid) mid = l_mid;
    else mid = r_mid;
  }
  for (int i = mid; i < end; ++i)
    procsForCells[i] += 1 << (level - 1);
  if (cells[0].dim() == 1) {
    switch (d) {
      case 'x':
        timercb(level - 1, begin, mid, 't');
        timercb(level - 1, mid, end, 't');
        break;
      case 't':
        timercb(level - 1, begin, mid, 'x');
        timercb(level - 1, mid, end, 'x');
        break;
      default:
        break;
    }
  } else if (cells[0].dim() == 2) {
    switch (d) {
      case 'x':
        timercb(level - 1, begin, mid, 'y');
        timercb(level - 1, mid, end, 'y');
        break;
      case 'y':
        timercb(level - 1, begin, mid, 't');
        timercb(level - 1, mid, end, 't');
        break;
      case 't':
        timercb(level - 1, begin, mid, 'x');
        timercb(level - 1, mid, end, 'x');
        break;
      default:
        break;
    }
  } else {
    switch (d) {
      case 'x':
        timercb(level - 1, begin, mid, 'y');
        timercb(level - 1, mid, end, 'y');
        break;
      case 'y':
        timercb(level - 1, begin, mid, 'z');
        timercb(level - 1, mid, end, 'z');
        break;
      case 'z':
        timercb(level - 1, begin, mid, 't');
        timercb(level - 1, mid, end, 't');
        break;
      case 't':
        timercb(level - 1, begin, mid, 'x');
        timercb(level - 1, mid, end, 'x');
        break;
      default:
        break;
    }
  }
}

void Distribution::timercb_weighted(int tlevel, int xlevel, int ylevel, int begin, int end, char d,
                                    Weights &W, int default_weight) {
  int N = end - begin;
  int level = tlevel + xlevel + ylevel;
  sort(cells.begin() + begin, cells.begin() + end, GetLess(d));
  if (PPM->Proc(commSplit) == 0) {
    int W_sum = 0;
    for (int i = begin; i < end; ++i) {
      if (W.find(cells[i]()) != W.end()) W_sum += W[cells[i]()];
      else W_sum += default_weight;
    }

    int W_mid = 0;
    int mid = begin;
    while (W_mid < W_sum / 2) {
      if (W.find(cells[mid]()) != W.end()) W_mid += W[cells[mid]()];
      else W_mid += default_weight;
      ++mid;
    }

    int old_mid = mid;
    int r_mid, l_mid;
    Point a = cells[mid]().CopyWithT(0.0);
    Point b = cells[mid - 1]().CopyWithT(0.0);

    while (cells[mid]().t() == cells[mid - 1]().t() || a == b) {
      ++mid;
      if (mid >= end) break;
      a = cells[mid]().CopyWithT(0.0);
      b = cells[mid - 1]().CopyWithT(0.0);
    }
    r_mid = mid;

    mid = old_mid;
    a = cells[mid]().CopyWithT(0.0);
    b = cells[mid - 1]().CopyWithT(0.0);

    while (cells[mid]().t() == cells[mid - 1]().t() || a == b) {
      --mid;
      if (mid < begin || mid == 0) break;
      a = cells[mid]().CopyWithT(0.0);
      b = cells[mid + 1]().CopyWithT(0.0);
    }
    l_mid = mid;

    if (old_mid - l_mid < r_mid - old_mid) mid = l_mid;
    else mid = r_mid;

    if (mid == begin || mid == end) { mid = old_mid; }
    for (int i = mid; i < end; ++i)
      procsForCells[i] += 1 << (level);

    if (tlevel > 0) {
      timercb_weighted(tlevel - 1, xlevel, ylevel, begin, mid, 't', W, default_weight);
      timercb_weighted(tlevel - 1, xlevel, ylevel, mid, end, 't', W, default_weight);
      return;
    }
    if (xlevel > 0) {
      timercb_weighted(tlevel, xlevel - 1, ylevel, begin, mid, 'x', W, default_weight);
      timercb_weighted(tlevel, xlevel - 1, ylevel, mid, end, 'x', W, default_weight);
      return;
    }
    if (ylevel > 0) {
      timercb_weighted(tlevel, xlevel, ylevel - 1, begin, mid, 'y', W, default_weight);
      timercb_weighted(tlevel, xlevel, ylevel - 1, mid, end, 'y', W, default_weight);
      return;
    }
  }
}

void Distribution::time_stripes(int n_procs, int begin, int end) {
  int N = end - begin;
  sort(cells.begin() + begin, cells.begin() + end, TLess_t);

  if (PPM->Proc(commSplit) == 0) {
    vector<int> P(n_procs);
    for (int i = 0; i < n_procs; ++i)
      P[i] = ((end - 1) / n_procs) * (i + 1);
    P[n_procs - 1] = end - 1;

    for (int i = 0; i < n_procs - 1; ++i) {
      int old_mid = P[i];
      int r_mid, l_mid;
      Point a = cells[P[i]]().CopyWithT(0.0);
      Point b = cells[P[i] - 1]().CopyWithT(0.0);

      while (cells[P[i]]().t() == cells[P[i] - 1]().t() || a == b) {
        ++P[i];
        a = cells[P[i]]().CopyWithT(0.0);
        b = cells[P[i] - 1]().CopyWithT(0.0);
      }
      r_mid = P[i];

      P[i] = old_mid;
      a = cells[P[i]]().CopyWithT(0.0);
      b = cells[P[i] - 1]().CopyWithT(0.0);

      while (cells[P[i]]().t() == cells[P[i] - 1]().t() || a == b) {
        --P[i];
        a = cells[P[i]]().CopyWithT(0.0);
        b = cells[P[i] + 1]().CopyWithT(0.0);
      }
      l_mid = P[i];

      if (old_mid - l_mid < r_mid - old_mid) P[i] = l_mid - 1;
      else P[i] = r_mid - 1;
    }

    for (int j = 0; j <= P[0]; ++j)
      procsForCells[j] = 0;

    for (int i = 1; i < n_procs; ++i)
      for (int j = P[i - 1] + 1; j <= P[i]; ++j)
        procsForCells[j] = i;
  }
}

void Distribution::DistributeSTMesh(Mesh &mesh) {
  if (PPM->Size(commSplit) == 1) return;
  if (_distName == "NoDistribution") return;
  const int N = mesh.CellCount();
  int P = PPM->Size(commSplit);
  int rekDepth = int(log2(PPM->Size(commSplit)));
  Weights &weights = *this->weights;

  int default_weight = 1;
  Config::Get("lb_default_weight", default_weight);
  if (_distName == "LikeSpace") {
    if (spaceDist != nullptr) {
      std::unordered_map<Point, int> spacePointToProc;
      for (size_t i = 0; i < spaceDist->distributedCells.size(); i++) {
        spacePointToProc[spaceDist->distributedCells[i]] = spaceDist->procsForCells[i];
      }
      sort(cells.begin(), cells.end(), Less_t);
      for (int i = 0; i < cells.size(); ++i) {
        Point spacePoint = cells[i].SpaceCell()();
        if (spacePointToProc.find(spacePoint) == spacePointToProc.end()) {
          THROW("SpacePoint not found!");
        }
        procsForCells[i] = spacePointToProc[spacePoint];
      }
    } else {
      THROW("Distribution 'LikeSpace' needs spatial distribution to be set!")
    }
  } else if (_distName == "Stripes") {
    sort(cells.begin(), cells.end(), Less_t);
    int m = (N + PPM->Size(commSplit) - 1) / PPM->Size(commSplit);
    int r = N - (m - 1) * PPM->Size(commSplit);
    int k = 0;
    int q = 0;
    for (int i = 0; i < cells.size(); ++i) {
      if (k >= m) {
        k = 0;
        ++q;
        if (q == r) --m;
      }
      ++k;
      markedCells[q].push_back(cells[i]);
    }
  } else if (_distName == "st-opt") {
    int S = mesh.steps() - 1;
    int T_Proc = int(double(P) / std::pow(2.0, double(rekDepth)) + 1e-10);
    T_Proc = std::min(T_Proc, S);
    timeOpt(T_Proc, rekDepth, 0, cells.size());
  } else if (_distName == "RCB") {
    // RCB_x(rekDepth, 0, cells.size());
    trcb(rekDepth, 0, cells.size(), 'x');
  } else if (_distName == "RCBNoTime") {
    trcb(rekDepth, 0, cells.size(), 'a');
  } else if (_distName == "deformed") {
    timercb(rekDepth, 0, cells.size(), 't');
  } else if (_distName == "deformed_optimized") {
    int S = mesh.steps() - 1;
    int Slice_Maximum = int(log(double(S)) / log(2.0) + 1e-10);
    int X_Maximum = int(log(double(sqrt(N / S))) / log(2.0) + 1e-10);
    int Y_Maximum = X_Maximum;

    if (Slice_Maximum == 0 && X_Maximum == 0 && Y_Maximum == 0) {
      THROW("Not enough cells to distribute. Level: " + mesh.Level().str())
    }

    int tL = 0;
    int yL = 0;
    int xL = 0;
    int L0 = rekDepth;
    for (int i = 0; i < rekDepth; ++i) {
      if (L0 == 0) break;
      if (tL < Slice_Maximum) {
        ++tL;
        --L0;
      }
      if (L0 == 0) break;
      if (xL < X_Maximum) {
        ++xL;
        --L0;
      }
      if (L0 == 0) break;
      if (yL < Y_Maximum) {
        ++yL;
        --L0;
      }
    }

    vout(1) << "Distribution:" << endl;
    vout(1) << " " << tL << " refinements in Time" << endl;
    vout(1) << " " << xL << " refinements in x" << endl;
    vout(1) << " " << yL << " refinements in y" << endl;

    if (tL > 0) timercb_weighted(tL - 1, xL, yL, 0, cells.size(), 't', weights, default_weight);
    else if (xL > 0)
      timercb_weighted(tL, xL - 1, yL, 0, cells.size(), 'x', weights, default_weight);
    else if (yL > 0)
      timercb_weighted(tL, xL, yL - 1, 0, cells.size(), 'y', weights, default_weight);
  } else if (_distName == "optimized_doubleslit") {
    int S = mesh.steps() - 2;
    int Slice_Maximum = int(log(double(S)) / log(2.0) + 1e-10) - 1;
    int X_Maximum = int(log(double(sqrt(N / (2 * S)))) / log(2.0) + 1e-10);
    int Y_Maximum = 1 + int(log(double(sqrt(N / (2 * S)))) / log(2.0) + 1e-10);
    int tL = 0;
    int yL = 0;
    int xL = 0;
    int L0 = rekDepth;
    if (tL < Slice_Maximum) {
      ++tL;
      --L0;
    }
    if (tL < Slice_Maximum) {
      ++tL;
      --L0;
    }
    if (yL < Y_Maximum) {
      ++yL;
      --L0;
    }
    for (int i = 0; i < rekDepth; ++i) {
      if (L0 == 0) break;
      if (tL < Slice_Maximum) {
        ++tL;
        --L0;
      }
      if (L0 == 0) break;
      if (xL < X_Maximum) {
        ++xL;
        --L0;
      }
      if (L0 == 0) break;
      if (yL < Y_Maximum) {
        ++yL;
        --L0;
      }
    }

    vout(1) << "Distribution:" << endl;
    vout(1) << " " << tL << " refinements in Time" << endl;
    vout(1) << " " << xL << " refinements in x" << endl;
    vout(1) << " " << yL << " refinements in y" << endl;

    for (int i = 0; i < N; ++i) {
      if (cells[i]()[0] < 0) weights[cells[i]()] = 0;
    }

    if (tL > 0) timercb_weighted(tL - 1, xL, yL, 0, cells.size(), 't', weights, default_weight);
    else if (xL > 0)
      timercb_weighted(tL, xL - 1, yL, 0, cells.size(), 'x', weights, default_weight);
    else if (yL > 0)
      timercb_weighted(tL, xL, yL - 1, 0, cells.size(), 'y', weights, default_weight);
  } else if (_distName == "time_stripes") {
    time_stripes(P, 0, cells.size());
  } else if (_distName == "x_stripes") {
    timercb_weighted(0, rekDepth - 1, 0, 0, cells.size(), 'x', weights);
  } else if (_distName == "y_stripes") {
    timercb_weighted(0, 0, rekDepth - 1, 0, cells.size(), 'y', weights);
  }
}

int Distribution::CommSplit() const { return commSplit; }
