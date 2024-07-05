#include "Mesh.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Buffer.hpp"
#include "Cell.hpp"
#include "Celltype.hpp"
#include "CoarseGeometry.hpp"
#include "Config.hpp"
#include "Edge.hpp"
#include "ExchangeBuffer.hpp"
#include "Face.hpp"
#include "ICell.hpp"
#include "Identify.hpp"
#include "LevelPair.hpp"
#include "MeshInfo.hpp"
#include "MeshPart.hpp"
#include "MeshSettings.hpp"
#include "Parallel.hpp"
#include "Point.hpp"
#include "ProcSet.hpp"
#include "Rule.hpp"

namespace {

int find_faceid_on_bndface(const cell &c, bnd_face b) {
  for (int i = 0; i < c.SpaceCell().Faces(); ++i) {
    if (c.Face(i) != b()) continue;
    return i;
  }
  THROW("faceid for bndface not found.");
}

TimeSteps refineTimeSteps(const TimeSteps &base) {
  TimeSteps refinedTimeSteps(2 * base.size() - 1);

  for (int k = 0, l = 0; k < base.size(); ++k, ++l) {
    refinedTimeSteps[l] = base[k];
    ++l;
    if (l >= refinedTimeSteps.size()) break;
    refinedTimeSteps[l] = 0.5 * (base[k] + base[k + 1]);
  }
  return refinedTimeSteps;
}

} // namespace

bool Mesh::isIncluded(int subdomain) const {
  return (settings.include.empty()
          || std::find(settings.include.begin(), settings.include.end(), subdomain)
                 != settings.include.end());
}

void Mesh::removeFace(const Point &F, const face &f, const Point &C) {
  if (find_cell(C) != cells_end()) return;
  if (C == Infty) procSets.RemoveIfOn(F, PPM->Proc(level.commSplit));
  meshFaces.Remove(F);
}

void Mesh::removeFace(const cell &c, const Point &F) {
  face f = find_face(F);
  if (f == faces_end()) return;
  if (find_cell(f.Left()) == c) removeFace(F, f, f.Right());
  else if (find_cell(f.Right()) == c) removeFace(F, f, f.Left());

  bnd_face bf = find_bnd_face(F);
  if (bf != bnd_faces_end()) removeBoundaryFace(F);
}

void Mesh::finishParallel() { Parallel = PPM->Or(procsets() != procsets_end(), level.commSplit); }

void Mesh::finish() {
  Parallel = PPM->Or(procsets() != procsets_end(), level.commSplit);
  const int dimension = dim();
  for (cell c2 = cells(); c2 != cells_end(); c2++) {
    if (c2.dim() != dimension)
      Exit("No good Mesh_ " + to_string(c2.dim()) + " - " + to_string(dimension))
  }
  for (bnd_face b = bnd_faces(); b != bnd_faces_end(); ++b) {
    face f = find_face(b());
    if (f == faces_end()) { THROW("face not found") }
    cell c = find_cell_on_bndface(f);
    int faceid = find_faceid_on_bndface(c, b);
    identifySets.Identify(b(), b.Part());
    for (int j = 0; j < c.SpaceCell().FaceCorners(faceid); ++j)
      identifySets.Identify(c.SpaceCell().FaceCorner(faceid, j), b.Part());
    for (int j = 0; j < c.SpaceCell().FaceEdges(faceid); ++j)
      identifySets.Identify(c.SpaceCell().FaceEdge(faceid, j), b.Part());
  }
  identifyBnd = PPM->Or(identifysets() != identifysets_end(), level.commSplit);
  findMinMaxMeshWidth();
  if (!Parallel) return;
  ExchangeBuffer E(level.commSplit);
  for (identifyset is = identifysets(); is != identifysets_end(); ++is) {
    procset p = find_procset(is());
    if (p == procsets_end()) continue;
    for (int i = 0; i < p.size(); ++i) {
      if (p[i] == PPM->Proc(level.commSplit)) continue;
      E.Send(p[i]) << short(is.size()) << is();
      for (int j = 0; j < is.size(); ++j)
        E.Send(p[i]) << is[j];
    }
  }
  for (short q = 0; q < PPM->Size(level.commSplit); ++q)
    if (E.Send(q).size()) E.Send(q) << short(0);
  E.Communicate();
  for (short q = 0; q < PPM->Size(level.commSplit); ++q) {
    if (E.Receive(q).Size() == 0) continue;
    short m;
    E.Receive(q) >> m;
    while (m) {
      Point x;
      E.Receive(q) >> x;
      for (int i = 0; i < m; ++i) {
        Point y;
        E.Receive(q) >> y;
        identifySets.Insert(x, y);
      }
      E.Receive(q) >> m;
    }
  }
  E.ClearBuffers();
}

bool Mesh::onBnd(const Point &x) const { return (find_bnd_face(x) != bnd_faces_end()); }

bool Mesh::onBoundary(const Cell &c, int f) const {
#ifdef USE_SPACETIME
  face ff = find_face(c.Face(f));
  return ((ff.Right() == Infty) || (ff.Left() == Infty));
#else
  return onBnd(c.Face(f));
#endif
}

void Mesh::findMinMaxMeshWidth() {
  double h_min = infty;
  double h_max = 0;
  for (cell c = cells(); c != cells_end(); ++c)
    if (this->dim() != 1) {
      for (int i = 0; i < c.SpaceCell().Edges(); ++i) {
        double h = dist(c.SpaceCell().EdgeCorner(i, 0), c.SpaceCell().EdgeCorner(i, 1));
        if (h < h_min) h_min = h;
        if (h > h_max) h_max = h;
      }
    } else {
      double h = dist(c.SpaceCell().Corner(0), c.SpaceCell().Corner(1));
      if (h < h_min) h_min = h;
      if (h > h_max) h_max = h;
    }
  minMeshWidth = PPM->Min(h_min, level.commSplit);
  maxMeshWidth = PPM->Max(h_max, level.commSplit);
}

double Mesh::MaxMeshWidth() const { return maxMeshWidth; }

double Mesh::MinMeshWidth() const { return minMeshWidth; }

Mesh::~Mesh() {
  while (CellCount())
    removeCell(cells());
  procSets.clear();
  meshBndFaces.clear();
}

BFParts::BFParts(const Mesh &M, const Cell &c) : n(c.Faces()), onbnd(false) {
  for (int i = 0; i < n; ++i) {
    bnd[i] = M.BoundaryFacePart(c.Face(i));
    if (bnd[i] != -1) onbnd = true;
  }
}

std::ostream &operator<<(std::ostream &s, const BFParts &BF) {
  for (int i = 0; i < BF.size(); ++i)
    s << BF[i] << " ";
  return s << endl;
}

bool Mesh::refineEdge(const procset &p, const Mesh &base) {
  edge e = base.find_edge(p());
  if (e != base.edges_end()) {
    procSets.Copy(p, 0.5 * (e.Left() + e()));
    procSets.Copy(p, 0.5 * (e.Right() + e()));
    return true;
  }
  return false;
}

void Mesh::refineFace(const procset &p, const Mesh &base) {
  face f = base.find_face(p());
  if (f == base.faces_end()) return;
  cell c = base.find_cell(f.Left());
  if (c == base.cells_end()) c = base.find_cell(f.Right());
  if (c == base.cells_end()) return;
  int i = c.facecorner(p());
  if (i == -1) return;
  Point P = f->first;
  Point E0 = 0.5 * (c.FaceCorner(i, 0) + c.FaceCorner(i, 1));
  Point E1 = 0.5 * (c.FaceCorner(i, 1) + c.FaceCorner(i, 2));
  if (c.FaceCorners(i) == 4) {
    Point E2 = 0.5 * (c.FaceCorner(i, 2) + c.FaceCorner(i, 3));
    Point E3 = 0.5 * (c.FaceCorner(i, 3) + c.FaceCorner(i, 0));
    procSets.Copy(p, 0.5 * (E0 + P));
    procSets.Copy(p, 0.5 * (E1 + P));
    procSets.Copy(p, 0.5 * (E2 + P));
    procSets.Copy(p, 0.5 * (E3 + P));
    procSets.Copy(p, 0.25 * (c.FaceCorner(i, 0) + E0 + P + E3));
    procSets.Copy(p, 0.25 * (c.FaceCorner(i, 1) + E1 + P + E0));
    procSets.Copy(p, 0.25 * (c.FaceCorner(i, 2) + E2 + P + E1));
    procSets.Copy(p, 0.25 * (c.FaceCorner(i, 3) + E3 + P + E2));
  } else if (c.FaceCorners(i) == 3) {
    Point E2 = 0.5 * (c.FaceCorner(i, 2) + c.FaceCorner(i, 0));
    procSets.Copy(p, 0.5 * (E0 + E1));
    procSets.Copy(p, 0.5 * (E1 + E2));
    procSets.Copy(p, 0.5 * (E0 + E2));
    procSets.Copy(p, (1 / 3.0) * (c.FaceCorner(i, 0) + E0 + E2));
    procSets.Copy(p, (1 / 3.0) * (c.FaceCorner(i, 1) + E0 + E1));
    procSets.Copy(p, (1 / 3.0) * (c.FaceCorner(i, 2) + E1 + E2));
    procSets.Copy(p, (1 / 3.0) * (E0 + E1 + E2));
  }
}

Mesh::Mesh(const MeshSettings &_settings, int commSplit) :
    settings(_settings), level(LevelPair{0, -1, 0, commSplit}),
    refineInSpaceWrapper(std::mem_fn(
        _settings.barycentricRefinement ? &Mesh::refineInSpaceBaryCentric : &Mesh::refineInSpace)),
    Parallel(false) {
  Config::Get("MeshVerbose", verbose);
  procSets.SetCommSplit(commSplit);
  fill();
  finish();
}

Mesh::Mesh(const Mesh &base, bool refine, bool refineInTime, bool mapGeometry) :
    verbose(base.verbose), settings(base.settings),
    level(LevelPair(refine ? (refineInTime ? base.level.NextInTime() : base.level.NextInSpace()) :
                             base.level)),
    timesteps(base.timesteps), Parallel(base.Parallel),
    refineInSpaceWrapper(base.refineInSpaceWrapper), refineInTimeWrapper(base.refineInTimeWrapper) {
  procSets.SetCommSplit(level.commSplit);

  if (mapGeometry) {
    if (refine) {
      const Mesh refinedMesh(base, true, refineInTime);
      this->mapGeometry(refinedMesh);
    } else {
      this->mapGeometry(base);
    }
    return;
  }

  if (refine) {
    if (refineInTime) {
      std::invoke(refineInTimeWrapper, *this, base);
    } else {
      std::invoke(refineInSpaceWrapper, *this, base);
    }
    return;
  }

  // Copy from base mesh if no map or refinement was requested
  // Copy vertices, edges, faces and cells
  std::for_each(base.cells(), base.cells_end(), [&](const auto &c) {
    cell C = insertCell(c.Type(), c.Subdomain(), c.AsVector());
#ifdef USE_DATAMESH
    C.SetData(c.GetData());
    vector<DataContainer> d(C.Corners());
    for (int k = 0; k < C.Corners(); ++k) {
      find_vertex(C.Corner(k)).SetData(base.find_vertex(c.Corner(k)).GetData());
    }
#endif
  });
  // Copy remaining members
  meshBndFaces = BoundaryFaces(base.meshBndFaces);
  procSets = ProcSets(base.procSets);
  identifySets = IdentifySets(base.identifySets);
  timeMidPointValues = std::vector<double>(base.timeMidPointValues);
  minMeshWidth = base.minMeshWidth;
  maxMeshWidth = base.maxMeshWidth;
}

void Mesh::fill() {
  const CoarseGeometry &_cGeo = *settings.coarseGeometry;
  if (!PPM->Master(level.commSplit)) return;

#ifdef USE_DATAMESH
  int cnt = 0;
#endif
  for (const auto &Cell_Id : _cGeo.GetCellIds()) {
    if (isIncluded(Cell_Id.Subdomain())) {
      vector<Point> corners = _cGeo.Corners(Cell_Id);
      auto newCell = insertCell(Cell_Id.Type(), Cell_Id.Subdomain(), corners);
#ifdef USE_DATAMESH
      for (int i = 0; i < corners.size(); ++i)
        find_vertex(corners[i]).SetData(_cGeo.GetVertexDataList()[Cell_Id.CornerIndices()[i]]);
      newCell.SetData(_cGeo.GetCellDataList()[cnt++]);
#endif
    }
  }

  for (const auto &Face_Id : _cGeo.GetFaceIds()) {
    Point faceCenter = _cGeo.FaceCenter(Face_Id);
    if (find_face(faceCenter) != faces_end()) insertBoundaryFace(faceCenter, Face_Id.BndType());
  }
  int bndFaceCount = 0;
  for (face f = faces(); f != faces_end(); ++f) {
    if (f.Right() == Infty) {
      bndFaceCount++;
      insertBoundaryFace(f(), 0);
    }
  }

  if (settings.sdSetter) {
    for (cell c = cells(); c != cells_end(); ++c)
      c->second->subdomain = settings.sdSetter(*c);
  }
}

void Mesh::fillST(const Mesh &base) {
  // Switch to ST refinement functions for subsequent calls
  refineInSpaceWrapper = std::mem_fn(&Mesh::refineSTMeshInSpace);
  refineInTimeWrapper = std::mem_fn(&Mesh::refineSTMeshInTime);

  timesteps = settings.coarseGeometry->timeSteps;

  for (cell c = base.cells(); c != base.cells_end(); ++c) {
    insertCellsWithTime(c);
  }

  for (face f = faces(); f != faces_end(); ++f) {
    Point pf = f().CopyWithT(0.0);
    bnd_face bf = base.find_bnd_face(pf);
    if (bf != base.bnd_faces_end()) { insertBoundaryFace(f(), bf.Part()); }
  }
  finish();
}

int find_other_proc(const procset &ps) {
  if (ps.size() != 2) THROW("Face-procset should only be on 2 processes.");
  int p1 = ps->second[0];
  int p2 = ps->second[1];
  return p1 == PPM->Proc(ps->second.CommSplit()) ? p2 : p1;
}

void Mesh::refineInSpaceInternal(const Mesh &base, const Cells &baseCells,
                                 const Mesh::CellRefinementFunction &cellRefinement) {
  for (procset p = base.procsets(); p != base.procsets_end(); ++p) {
    procSets.Copy(p);
    if (!refineEdge(p, base)) refineFace(p, base);
  }
  std::vector<Point> z;
  std::for_each(cell(baseCells.Begin()), cell(baseCells.End()), [&](const auto &c) {
    BFParts bnd(base, c);
    const std::vector<Rule> &rules = std::invoke(cellRefinement, c, z);
    for (const Rule &rule : rules) {
      cell new_cell = insertCell(rule.type(), c.Subdomain(), rule.getCorners(z), c.min(), c.max());
      if (bnd.onBnd()) {
        for (int j = 0; j < new_cell.SpaceCell().Faces(); ++j) {
          int faceid = rule.face[j];
          if (faceid != -1 && bnd[faceid] != -1) {
            insertBoundaryFace(new_cell.Face(j), bnd[faceid]);
          }
        }
      }
    }
  });
  finish();
  ExchangeBuffer buffer(base.level.commSplit);
  std::for_each(cells(), cells_end(), [&](const auto &fine_child) {
    BFParts bnd(*this, fine_child);
    for (int fid = 0; fid < fine_child.Faces(); fid++) {
      face fine_face = find_face(fine_child.Face(fid));
      if (fine_face.Right() == Infty && bnd[fid] == -1) {
        procset ps = find_procset(fine_face());
        if (ps == procsets_end()) continue;
        int other_process = find_other_proc(ps);
        buffer.Send(other_process) << fine_face() << fine_face.Left();
      }
    }
  });
  buffer.Communicate();
  for (int q = 0; q < PPM->Size(base.level.commSplit); ++q) {
    while (buffer.Receive(q).size() < buffer.ReceiveSize(q)) {
      Point fine_face, fine_face_right;
      buffer.Receive(q) >> fine_face >> fine_face_right;
      insertFace(fine_face, fine_face_right);
    }
  }
}

void Mesh::refineInSpace(const Mesh &base) {
  const CellRefinementFunction cellRefinement =
      [](const Cell &c, std::vector<Point> &z) -> const std::vector<Rule> & {
    return c.SpaceCell().Refine(z);
  };

  refineInSpaceInternal(base, base.meshCells, cellRefinement);
}

void Mesh::refineInSpaceBaryCentric(const Mesh &base) {
  Cells baseCells; // TODO: use different container - must be the same as
                   // decltype(meshCells)
  baseCells.reserve(base.CellCount() * 2);

  std::vector<Point> z;
  std::for_each(base.cells(), base.cells_end(), [&](const auto &cell) {
    const std::vector<Rule> &rules = cell.Refine(z);
    for (const auto &rule : rules) {
      Cell *insertedCell = CreateCell(rule.type(), cell.Subdomain(), rule.getCorners(z));
      baseCells.emplace((*insertedCell)(), insertedCell);
    }
  });

  const double shiftWidth = base.minMeshWidth / 10;
  const CellRefinementFunction baryCentricCellRefinement =
      [shiftedCenter = base.settings.shiftedCenter,
       shiftWidth](const Cell &c, std::vector<Point> &z) -> const std::vector<Rule> & {
    return c.RefineBarycentric(z, shiftedCenter ? Point::Random(-shiftWidth, shiftWidth) : Origin);
  };

  refineInSpaceInternal(base, baseCells, baryCentricCellRefinement);
}

void Mesh::refineSTInternal(const Mesh &base, const STCellsRefinementFunction &cellsRefinement) {
  for (procset p = base.procsets(); p != base.procsets_end(); ++p) {
    if (base.find_vertex(p()) != base.vertices_end()) { procSets.Copy(p); }
  }

  std::invoke(cellsRefinement, base);

  std::for_each(cells(), cells_end(), [&](const auto &new_cell) {
    for (int ei = 0; ei < new_cell.Edges(); ++ei) {
      procset e1p = find_procset(new_cell.EdgeCorner(ei, 0));
      procset e2p = find_procset(new_cell.EdgeCorner(ei, 1));
      if (e1p != procsets_end() && e2p != procsets_end()) {
        for (short q : *e1p) {
          if (e2p.in(q)) { procSets.Add(new_cell.Edge(ei), q); }
        }
      }
    }
  });

  for (procset ps = procsets(); ps != procsets_end(); ps++) {
    if (find_vertex(ps()) != vertices_end()) continue;
    if (find_face(ps()) != faces_end()) continue;
    if (find_edge(ps()) != edges_end()) continue;
    if (ps() == Infty) continue;
    procSets.erase(ps());
  }

  finish();
}

void Mesh::mapGeometry(const Mesh &base) {
  const CoarseGeometry &coarseGeometry = *settings.coarseGeometry;

  std::for_each(base.cells(), base.cells_end(), [&](const auto &c) {
    std::vector<Point> corners(c.Corners());
    for (int i = 0; i < corners.size(); ++i)
      corners[i] = coarseGeometry(c[i]);
    auto newCell = insertCell(c.Type(), c.Subdomain(), corners);
    BFParts bnd(base, c);
    for (int f = 0; f < c.Faces(); ++f) {
      bnd_face f_bnd = base.find_bnd_face(c.Face(f));
      if (f_bnd == base.bnd_faces_end()) { continue; }
      insertBoundaryFace(newCell.Face(f), f_bnd.Part());
    }
  });
  finish();
}

void Mesh::refineSTMeshInSpace(const Mesh &base) {
  const STCellsRefinementFunction cellsSpaceRefinement = [&](const Mesh &base) {
    std::vector<Point> z;
    std::for_each(base.cells(), base.cells_end(), [&](const auto &c) {
      BFParts bnd(base, c);
      const std::vector<Rule> &rules = c.SpaceCell().Refine(z);
      for (const Rule &rule : rules) {
        cell new_cell = insertCell(SpaceTimeCellType(rule.type()), c.Subdomain(),
                                   rule.getCorners(z), c.min(), c.max());
        if (bnd.onBnd()) {
          for (int j = 0; j < new_cell.SpaceCell().Faces(); ++j) {
            int faceid = rule.face[j];
            if (faceid != -1 && bnd[faceid] != -1) {
              insertBoundaryFace(new_cell.Face(j), bnd[faceid]);
            }
          }
        }
        for (int j = 0; j < new_cell.Faces(); ++j) {
          int faceid = j < new_cell.SpaceCell().Faces() ? rule.face[j] : j;
          if (faceid < 0) continue;
          procset pf = base.find_procset(c.Face(faceid));
          if (pf == base.procsets_end()) continue;
          procSets.Add(new_cell.Face(j), pf);
        }
        for (int ci = 0; ci < new_cell.Corners(); ci++) {
          procset pci = base.find_procset(new_cell.Corner(ci));
          if (pci != base.procsets_end()) { procSets.Add(new_cell.Corner(ci), pci); }
        }
      }
    });
  };

  refineSTInternal(base, cellsSpaceRefinement);
}

void Mesh::refineSTMeshInTime(const Mesh &base) {
  timesteps = refineTimeSteps(base.timesteps);

  const STCellsRefinementFunction cellsTimeRefinement = [&](const Mesh &base) {
    std::vector<Point> z;
    std::for_each(base.cells(), base.cells_end(), [&](const auto &c) {
      BFParts bnd(base, c);
      const std::vector<Rule> &rules = c.TimeCell().Refine(z);
      for (const Rule &rule : rules) {
        std::vector<Point> x;
        rule(z, x);
        cell new_cell = insertCell(SpaceTimeCellType(c.SpaceCell().Type()), c.Subdomain(),
                                   c.SpaceCell().AsVector(), x[0][0], x[1][0]);
        if (bnd.onBnd()) {
          for (int j = 0; j < new_cell.SpaceCell().Faces(); ++j) {
            if (bnd[j] != -1) { insertBoundaryFace(new_cell.Face(j), bnd[j]); }
          }
        }
        for (int j = 0; j < new_cell.Faces(); ++j) {
          Point fj = new_cell.Face(j);
          if (c().t() == fj.t()) continue;
          procset pf = base.find_procset(c.Face(j));
          if (pf == base.procsets_end()) continue;
          procSets.Add(fj, pf);
        }
        for (int ci = 0; ci < new_cell.Corners(); ci++) {
          procset pci = base.find_procset(new_cell.Corner(ci));
          if (pci == base.procsets_end()) continue;
          if (find_procset(pci()) != procsets_end()) { procSets.Add(new_cell.Corner(ci), pci); }
        }
      }
    });
  };

  refineSTInternal(base, cellsTimeRefinement);
}

void Mesh::insertVertex(const Point &p) {
  auto check = find_vertex(p);
  if (check == vertices_end()) {
    meshVertices.Insert(p, new Vertex(0));
  } else check->second->inc();
}

bool Mesh::removeVertex(const Point &p) {
  auto v = find_vertex(p);
  if (v->second->dec() > 0) return false;
  meshVertices.Remove(p);
  return true;
}

void Mesh::insertEdge(const Point &left, const Point &right) {
  Point middle = 0.5 * (left + right);
  auto e = meshEdges.Find(middle);
  if (e == meshEdges.End()) meshEdges.Insert(middle, new Edge(left, right));
  else e->second->inc();
}

bool Mesh::removeEdge(const Point &edgeCenter) {
  auto e = meshEdges.Find(edgeCenter);
  if (e->second->dec() > 0) return false;
  meshEdges.Remove(edgeCenter);
  return true;
}

void Mesh::insertFace(const Point &faceCenter, const Point &cellCenter) {
  auto f = meshFaces.Find(faceCenter);
  if (f == meshFaces.End()) {
    meshFaces.Insert(faceCenter, Face(cellCenter));
  } else if (f->second.Right() == Infty) {
    if (f->second.Left() != cellCenter) {
      meshFaces.Replace(faceCenter, Face(f->second.Left(), cellCenter));
    }
  }
}

void Mesh::insertFace(const Point &faceCenter, const Face &faceValue) {
  meshFaces.Insert(faceCenter, Face(faceValue));
}

void Mesh::insertBoundaryFace(const Point &z, int part) { meshBndFaces.Insert(z, part); }

void Mesh::removeBoundaryFace(const Point &x) { meshBndFaces.Remove(x); }

cell Mesh::insertCell(CELLTYPE tp, short subdomain, const std::vector<Point> &x, double a,
                      double b) {
  auto *insertedCell = CreateCell(tp, subdomain, x, a, b);
  cell c = meshCells.Insert((*insertedCell)(), insertedCell);
  for (int i = 0; i < c.Corners(); ++i) {
    insertVertex(c.Corner(i));
  }
  for (int i = 0; i < c.Edges(); ++i) {
    insertEdge(c.EdgeCorner(i, 0), c.EdgeCorner(i, 1));
  }
  for (int i = 0; i < c.Faces(); ++i) {
    insertFace(c.Face(i), c());
  }

  return c;
}

void Mesh::insertCellsWithTime(cell &c) {
  static int truncate = -1;
  Point Pmid;
  double vMax = 2;
  double rhoMin = 0.25;
  if (truncate == -1 || truncate > 0) {
    Config::Get("ProblemMid", Pmid);
    Config::Get("MaxKappa", vMax);
    Config::Get("RhoMin", rhoMin);
    Config::Get("truncateSTMesh", truncate);
    if (truncate == -1) { truncate = 0; }
  }

  int numberofsubtimesteps = 2;
  // Config::Get("subtimesamples", numberofsubtimesteps);

  const Point center = c.Center();
  for (int k = 0; k < timesteps.size() - 1; ++k) {
    vector<double> subtimesteps(numberofsubtimesteps);
    for (int a = 0; a < numberofsubtimesteps; ++a)
      subtimesteps[a] = t(k) + double(a) / double(numberofsubtimesteps - 1) * (t(k + 1) - t(k));

    if (truncate == 1) {
      if (dist(center, Pmid) > vMax / 2.25 * (t(k + 1) + 0.5)) continue;
      double x = center[0];
      if (x < 4.75 - vMax * (4 - t(k))) continue;
      if (x > 7.25 + vMax * (4 - t(k))) continue;

      if (x < 4.75 - vMax / 2 * (0.5 + 4 - t(k))) continue;
      if (x > 7.25 + vMax / 2 * (0.5 + 4 - t(k))) continue;
      double y = center[1];
      if (y > 0.4 + vMax / 2 * (1 + 4 - t(k))) continue;
      if (y > 0.4 + vMax * (4 - t(k))) continue;
    }

    if (truncate == 2) {
      Point pMid{0.5, 0.75, 0.0};
      vMax = sqrt(1 / rhoMin);
      if (dist(center, pMid) > vMax * (t(k + 1) + 0.0625)) continue;
      //      double x = center[0];
      //      if (x < 4.75 - vMax * (4 - t(k))) continue;
      //      if (x > 7.25 + vMax * (4 - t(k))) continue;
      //
      //      if (x < 4.75 - vMax / 2 * (0.5 + 4 - t(k))) continue;
      //      if (x > 7.25 + vMax / 2 * (0.5 + 4 - t(k))) continue;
      //      double y = center[1];
      //      if (y > 0.4 + vMax / 2 * (1 + 4 - t(k))) continue;
      //      if (y > 0.4 + vMax * (4 - t(k))) continue;
    }

    if (truncate == 3) {
      double x = center[0];
      double y = center[1];
      double t = timesteps[k];
      if (x - 1.5 > t) continue;
    }

    if (truncate == 4) {
      double x = center[0];
      double y = center[1];
      double t = timesteps[k];
      if (y < 0.5 - t) continue;
      if (y > 1.5 - t) continue;
    }

    for (int a = 0; a < numberofsubtimesteps - 1; ++a) {
      cell tc = insertCell(SpaceTimeCellType(c.Type()), 0, c.AsVector(), subtimesteps[a],
                           subtimesteps[a + 1]);
    }
  }
}

void Mesh::removeCell(cell c) {
  for (int i = 0; i < c.Corners(); ++i) {
    if (removeVertex(c.Corner(i))) { procSets.RemoveIfOn(c.Corner(i), PPM->Proc(level.commSplit)); }
  }
  for (int i = 0; i < c.Edges(); ++i) {
    if (removeEdge(c.Edge(i))) { procSets.RemoveIfOn(c.Edge(i), PPM->Proc(level.commSplit)); }
  }
  for (int i = 0; i < c.Faces(); ++i) {
    removeFace(c, c.Face(i));
  }
  meshCells.Remove(c());
}

vertex Mesh::vertices() const { return {meshVertices.Begin()}; }

vertex Mesh::vertices_end() const { return {meshVertices.End()}; }

vertex Mesh::find_vertex(const Point &p) const { return {meshVertices.Find(p)}; }

cell Mesh::cells() const { return {meshCells.Begin()}; }

cell Mesh::cells_end() const { return {meshCells.End()}; }

cell Mesh::find_cell(const Point &center) const { return {meshCells.Find(center)}; }

cell Mesh::find_cell_on_bndface(const face f) const {
  cell c = find_cell(f.Left());
  if (c == cells_end()) c = find_cell(f.Right());
  if (c == cells_end()) THROW("Cell on Bnd not found.");
  return c;
}

int Mesh::find_neighbour_face_id(const Point &f_c, const cell &cf) const {
  for (int f1 = 0; f1 < cf.Faces(); ++f1) {
    if (f_c == cf.Face(f1)) return f1;
  }
  THROW("Mesh::find_neighbour_face_id Face not found!")
}

edge Mesh::edges() const { return {meshEdges.Begin()}; }

edge Mesh::edges_end() const { return {meshEdges.End()}; }

edge Mesh::find_edge(const Point &center) const { return {meshEdges.Find(center)}; }

face Mesh::faces() const { return {meshFaces.Begin()}; }

face Mesh::faces_end() const { return {meshFaces.End()}; }

face Mesh::find_face(const Point &center) const { return {meshFaces.Find(center)}; }

bnd_face Mesh::bnd_faces() const { return {meshBndFaces.Begin()}; }

bnd_face Mesh::bnd_faces_end() const { return {meshBndFaces.End()}; }

bnd_face Mesh::find_bnd_face(const Point &z) const { return {meshBndFaces.Find(z)}; }

procset Mesh::procsets() const { return {procSets.Begin()}; }

procset Mesh::procsets_end() const { return {procSets.End()}; }

procset Mesh::find_procset(const Point &z) const { return {procSets.Find(z)}; }

bool Mesh::master(const Point &z) const { return procSets.master(z); }

identifyset Mesh::identifysets() const { return {identifySets.Begin()}; }

identifyset Mesh::identifysets_end() const { return {identifySets.End()}; }

identifyset Mesh::find_identifyset(const Point &z) const { return {identifySets.Find(z)}; }

const std::string &Mesh::DistributionName() const { return settings.distributionName; }

std::string Mesh::meshWidthStr() const {
  std::string strMeshWidth =
      "[" + std::to_string(minMeshWidth) + ", " + std::to_string(maxMeshWidth) + "]";
  return strMeshWidth;
}

std::string TimeStepsString(const std::vector<double> &timesteps) {
  std::string str = "[ ";
  for (auto ts : timesteps)
    str += to_string(ts) + " ";
  str += "]";
  return str;
}

void Mesh::PrintBndIdsPerProc() const {
  std::unordered_map<int, size_t> bndValueCount;
  int cnt = 0;
  for (bnd_face f = bnd_faces(); f != bnd_faces_end(); f++) {
    cnt += 1;
    auto bndValueIter = bndValueCount.find(f.Part());
    if (bndValueIter == bndValueCount.end()) { bndValueCount[f.Part()] = 0; }
    bndValueCount[f.Part()] = bndValueCount[f.Part()] + 1;
  }
  std::stringstream ss;
  ss << "[";
  for (auto &[part, count] : bndValueCount) {
    ss << part << ":" << count << ", ";
  }
  ss << "]";

  pout << "cnt:" << cnt << " bnd_faces: " << ss.str() << endl;
}

void Mesh::PrintInfo() const {
  ExchangeBuffer exBuffer(level.commSplit);
  exBuffer.Send(0) << MeshInfoOnProc(*this);
  exBuffer.Communicate();
  MeshInfoOnProcs infos;
  if (PPM->Proc(level.commSplit) == 0) {
    for (int p = 0; p < PPM->Size(level.commSplit); p++) {
      MeshInfoOnProc info;
      exBuffer.Receive(p) >> info;
      infos.push_back(info);
    }
  }

  mout.PrintInfo("Mesh", verbose, PrintInfoEntry("Mesh Name", Name(), 1),
                 PrintInfoEntry("Distribution", settings.distributionName, 1),
                 PrintInfoEntry("CommSplit", level.commSplit, 1),
                 PrintInfoEntry("Level", level.str(), 1),
                 PrintInfoEntry("pLevel", settings.distributeLevel, 1),
                 PrintInfoEntry("Dimension", dim(), 1),
                 PrintInfoEntry("Cells", CellCountGeometry(), 1),
                 PrintInfoEntry("Vertices", VertexCountGeometry(), 1),
                 PrintInfoEntry("Edges", EdgeCountGeometry(), 1),
                 PrintInfoEntry("Faces", FaceCountGeometry(), 1),
                 PrintInfoEntry("BND Faces", infos.BNDFacesCountString(), 2),
                 PrintInfoEntry("ProcSets", ProcSetsCountGeometryWithoutInfty(), 1),
                 PrintInfoEntry("Mesh width [min, max]", meshWidthStr(), 1),
                 PrintInfoEntry("Time steps", timesteps.size(), 2),
                 PrintInfoEntry("Time steps", TimeStepsString(timesteps), 3),
                 PrintInfoEntry("Min cells on proc", PPM->Min(CellCount(), level.commSplit), 2),
                 PrintInfoEntry("Max cells on proc", PPM->Max(CellCount(), level.commSplit), 2),
                 PrintInfoEntry("Processes", infos.ProcString(), 2),
                 PrintInfoEntry("Cell distribution", infos.CellCountString(), 2),
                 PrintInfoEntry("Vertex distribution", infos.VertexCountString(), 2),
                 PrintInfoEntry("Edges distribution", infos.EdgesCountString(), 2),
                 PrintInfoEntry("Faces distribution", infos.FacesCountString(), 2),
                 PrintInfoEntry("ProcSets distribution", infos.ProcSetsCountString(), 2));
}

int Mesh::CellCount() const { return (int)meshCells.size(); }

int Mesh::CellCountGeometry() const {
  return PPM->SumOnCommSplit(meshCells.size(), level.commSplit);
}

int Mesh::SpaceCellCountGeometry(int slice_index) const {
  if (IsSTMesh()) {
    int count = 0;
    for (auto [p, c] : meshCells) {
      if (abs(p.t() - slice(slice_index)) < GeometricTolerance) count++;
    }
    return PPM->SumOnCommSplit(count, level.commSplit);
  } else {
    return CellCountGeometry();
  }
}

int Mesh::FaceCount() const { return (int)meshFaces.size(); }

int Mesh::FaceCountGeometry() const {
  int count = 0;
  for (face f = this->faces(); f != this->faces_end(); ++f)
    if (this->procSets.master(f())) count++;
  return PPM->SumOnCommSplit(count, level.commSplit);
}

int Mesh::EdgeCount() const { return (int)meshEdges.size(); }

int Mesh::EdgeCountGeometry() const {
  int count = 0;
  for (edge e = this->edges(); e != this->edges_end(); ++e)
    if (this->procSets.master(e())) count++;
  return PPM->SumOnCommSplit(count, level.commSplit);
}

int Mesh::VertexCount() const { return (int)meshVertices.size(); }

int Mesh::VertexCountGeometry() const {
  int count = 0;
  for (vertex v = this->vertices(); v != this->vertices_end(); ++v)
    if (this->procSets.master(v())) count++;
  return PPM->SumOnCommSplit(count, level.commSplit);
}

int Mesh::SpaceVertexCountGeometry() const {
  if (IsSTMesh()) return VertexCountGeometry() / steps();
  else return VertexCountGeometry();
}

int Mesh::ProcSetsCount() const { return (int)procSets.size(); }

int Mesh::ProcSetsCountWithoutInfty() const {
  if (find_procset(Infty) == procsets_end()) {
    return (int)procSets.size();
  } else {
    return (int)procSets.size() - 1;
  }
}

int Mesh::ProcSetsCountGeometry() const {
  int count = 0;
  for (procset p = this->procsets(); p != this->procsets_end(); ++p)
    if (this->procSets.master(p())) count++;
  return PPM->SumOnCommSplit(count, level.commSplit);
}

int Mesh::ProcSetsCountGeometryWithoutInfty() const {
  int count = 0;
  for (procset p = this->procsets(); p != this->procsets_end(); ++p)
    if (this->procSets.master(p()) && p() != Infty) count++;
  return PPM->SumOnCommSplit(count, level.commSplit);
}

std::pair<double, double> Mesh::MeshWidth() const { return {minMeshWidth, maxMeshWidth}; }

bool Mesh::IsSTMesh() const { return (!timesteps.empty()); }

int Mesh::BoundaryFaceCount() const { return (int)meshBndFaces.size(); }

int Mesh::CommSplit() const { return level.commSplit; }

int Mesh::GetBNDFaceCount(int i) const {
  int count{0};
  for (bnd_face f = bnd_faces(); f != bnd_faces_end(); ++f) {
    count += f.Part() == i;
  }
  return count;
}

bool Mesh::onBoundary(const Cell &c) const {
  for (int i = 0; i < c.Faces(); ++i) {
    int part = BoundaryFacePart(c.Face(i));
    if (part != -1) return true;
  }
  return false;
}

double Mesh::slice(int k) const {
  if (IsSTMesh()) {
    if (timeMidPointValues.empty()) { InitTValues(); }
    if (k >= 0 && k < timeMidPointValues.size()) {
      return timeMidPointValues[k];
    } else if (k < 0 && k >= -timeMidPointValues.size()) {
      return timeMidPointValues[timeMidPointValues.size() + k];
    }
    THROW("No Slice with index " + std::to_string(k));
  }
  return 0.0;
}

void Mesh::InitTValues() const {
  std::set<double> times;
  for (cell c = cells(); c != cells_end(); c++) {
    times.insert(c().t());
  }
  ExchangeBuffer buffer(CommSplit());
  for (int q = 0; q < PPM->Size(CommSplit()); q++) {
    if (q != PPM->Proc(CommSplit())) { buffer.Send(q) << times; }
  }
  buffer.Communicate();
  for (int q = 0; q < PPM->Size(CommSplit()); q++) {
    if (q != PPM->Proc(CommSplit())) { buffer.Receive(q) >> times; }
  }
  timeMidPointValues.clear();
  std::copy(times.begin(), times.end(), std::back_inserter(timeMidPointValues));
  std::sort(timeMidPointValues.begin(), timeMidPointValues.end());
}

Buffer &operator>>(Buffer &b, Mesh &M) {
  if (!M.IsSTMesh()) {
    int tp;
    short n, sd;

    b >> tp >> sd >> n;

    vector<Point> x(n);
    for (int i = 0; i < n; ++i)
      b >> x[i];
    cell c = M.insertCell(CELLTYPE(tp), sd, x);
#ifdef USE_DATAMESH
    DataContainer data;
    b >> data;
    c.SetData(data);

    for (int i = 0; i < n; ++i) {
      b >> data;
      M.find_vertex(c.Corner(i)).SetData(data);
    }
#endif
    return b;
  } else {
#ifdef USE_SPACETIME
    int tp;
    short n, sd;
    double min, max;
    b >> tp >> sd >> n >> min >> max;
    vector<Point> x(n);
    for (int i = 0; i < n; ++i)
      b >> x[i];
    M.insertCell(SpaceTimeCellType(static_cast<CELLTYPE>(tp)), sd, x, min, max);
    return b;
#endif
    return b;
  }
}

const ProcSets &Mesh::GetProcSets() const { return procSets; }

std::string Mesh::Name() const { return settings.coarseGeometry->Name(); }

const LevelPair &Mesh::Level() const { return level; }

int Mesh::BoundaryFacePart(const Point &bndFaceCenter) const {
  auto b = meshBndFaces.Find(bndFaceCenter);
  if (b == meshBndFaces.End()) return -1;
  return b->second;
}

double Mesh::GetEndTime() const { return timesteps[timesteps.size() - 1]; }

vector<double> Mesh::GetTimeMidpointValues() const { return timeMidPointValues; }

const TimeSteps &Mesh::GetTimesteps() const { return timesteps; }

int Mesh::slices() const { return (int)timesteps.size() - 1; }

int Mesh::steps() const { return (int)timesteps.size(); }

bool Mesh::identify() const { return identifyBnd; }

bool Mesh::parallel() const { return Parallel; }

const std::vector<int> &Mesh::IncludedSubdomains() const { return settings.include; }

double Mesh::t(int k) const {
  if (IsSTMesh()) return timesteps[k];
  else return 0.0;
}
