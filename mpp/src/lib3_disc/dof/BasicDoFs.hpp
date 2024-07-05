#ifndef BASICDOFS_HPP
#define BASICDOFS_HPP

#include "MixedDoF.hpp"

class CellDoF : public IDoF {
  short m;
public:
  explicit CellDoF(short M = 1) : m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override { return 1; }

  std::vector<Point> GetNodalPoints(const Cell &c) const override {
    return std::vector<Point>{c()};
  }

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>{m};
  }

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override { return 0; }

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override {
    THROW("Should never reach this function")
  }

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override { return 0; }

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId,
                             int k) const override{THROW("Should never reach this function")}

  string Name() const override {
    return "CellDoF";
  }
};

class VertexDoF : public IDoF {
  short m;
public:
  explicit VertexDoF(short M = 1, int n = 0, bool bnd = false) : IDoF(n, bnd), m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override { return c.Corners(); }

  std::vector<Point> GetNodalPoints(const Cell &c) const override {
    std::vector<Point> sp(NumberOfNodalPoints(c));
    for (int i = 0; i < sp.size(); ++i)
      sp[i] = c[i];
    return sp;
  }

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>(NumberOfNodalPoints(c), m);
  }

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override {
    return c.FaceCorners(faceId);
  }

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override {
    return c.facecorner(faceId, k);
  }

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override {
    return c.EdgeCorners(edgeId);
  }

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override {
    return c.edgecorner(edgeId, k);
  }

  string Name() const override { return "VertexDoF"; }
};

class EdgeDoF : public IDoF {
  short m;
public:
  explicit EdgeDoF(short M = 1) : m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override { return c.Edges(); }

  std::vector<Point> GetNodalPoints(const Cell &c) const override {
    std::vector<Point> sp(NumberOfNodalPoints(c));
    for (int i = 0; i < c.Edges(); ++i)
      sp[i] = c.Edge(i);
    return sp;
  }

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>(NumberOfNodalPoints(c), m);
  }

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override {
    return c.FaceEdges(faceId);
  }

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override {
    return c.faceedge(faceId, k);
  }

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override { return 1; }

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId, int k) const override { return edgeId; }

  string Name() const override { return "EdgeDoF"; }
};

class FaceDoF : public IDoF {
  short m;
public:
  explicit FaceDoF(short M = 1, int n = 0, bool bnd = false) : IDoF(n, bnd), m(M) {}

  short NumberOfComponents() const override { return m; }

  short NumberOfNodalPoints(const Cell &c) const override { return c.Faces(); }

  std::vector<Point> GetNodalPoints(const Cell &c) const override {
    std::vector<Point> sp(NumberOfNodalPoints(c));
    for (int i = 0; i < c.Faces(); ++i)
      sp[i] = c.Face(i);
    return sp;
  }

  std::vector<short> DoFSizesAtNodalPoints(const Cell &c) const override {
    return std::vector<short>(NumberOfNodalPoints(c), m);
  }

  short NumberOfNodalPointsOnFace(const Cell &c, int faceId) const override { return 1; }

  short IdOfNodalPointOnFace(const Cell &c, int faceId, int k) const override { return faceId; }

  short NumberOfNodalPointsOnEdge(const Cell &c, int edgeId) const override { return 0; }

  short IdOfNodalPointOnEdge(const Cell &c, int edgeId,
                             int k) const override{THROW("Should never reach this function")}

  string Name() const override {
    return "FaceDoF";
  }
};

class VertexCellDoF : public MixedDoF {
  std::vector<std::unique_ptr<IDoF>> createDoFs(short M, int n, bool bnd) {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<VertexDoF>(M, n, bnd);
    dofs[1] = std::make_unique<CellDoF>(M);
    return dofs;
  }
public:
  explicit VertexCellDoF(short M = 1, int n = 0, bool bnd = false) :
      MixedDoF(std::move(createDoFs(M, n, bnd)), n, bnd) {}

  std::string Name() const override { return "VertexCellDoF"; }
};

class VertexEdgeDoF : public MixedDoF {
  std::vector<std::unique_ptr<IDoF>> createDoFs(short M, int n, bool bnd) {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<VertexDoF>(M, n, bnd);
    dofs[1] = std::make_unique<EdgeDoF>(M);
    return dofs;
  }
public:
  explicit VertexEdgeDoF(short M = 1, int n = 0, bool bnd = false) :
      MixedDoF(std::move(createDoFs(M, n, bnd)), n, bnd) {}

  std::string Name() const override { return "VertexEdgeDoF"; }
};

class VertexEdgeCellDoF : public MixedDoF {
  std::vector<std::unique_ptr<IDoF>> createDoFs(short M, int n, bool bnd) {
    std::vector<std::unique_ptr<IDoF>> dofs(3);
    dofs[0] = std::make_unique<VertexDoF>(M, n, bnd);
    dofs[1] = std::make_unique<EdgeDoF>(M);
    dofs[2] = std::make_unique<CellDoF>(M);
    return dofs;
  }
public:
  explicit VertexEdgeCellDoF(short M = 1, int n = 0, bool bnd = false) :
      MixedDoF(std::move(createDoFs(M, n, bnd)), n, bnd) {}

  std::string Name() const override { return "VertexEdgeCellDoF"; }
};

class VertexEdgeFaceCellDoF : public MixedDoF {
  std::vector<std::unique_ptr<IDoF>> createDoFs(short M, int n, bool bnd) {
    std::vector<std::unique_ptr<IDoF>> dofs(4);
    dofs[0] = std::make_unique<VertexDoF>(M, n, bnd);
    dofs[1] = std::make_unique<EdgeDoF>(M);
    dofs[2] = std::make_unique<FaceDoF>(M);
    dofs[3] = std::make_unique<CellDoF>(M);
    return dofs;
  }
public:
  explicit VertexEdgeFaceCellDoF(short M = 1, int n = 0, bool bnd = false) :
      MixedDoF(std::move(createDoFs(M, n, bnd)), n, bnd) {}

  std::string Name() const override { return "VertexEdgeFaceCellDoF"; }
};

class FaceCellDoF : public MixedDoF {
  std::vector<std::unique_ptr<IDoF>> createDoFs(short M, int n, bool bnd) {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<FaceDoF>(M, n, bnd);
    dofs[1] = std::make_unique<CellDoF>(M);
    return dofs;
  }
public:
  explicit FaceCellDoF(short M = 1, int n = 0, bool bnd = false) :
      MixedDoF(std::move(createDoFs(M, n, bnd)), n, bnd) {}

  std::string Name() const override { return "FaceCellDoF"; }
};

class VertexFaceDoF : public MixedDoF {
  std::vector<std::unique_ptr<IDoF>> createDoFs(int _mvertex, int _mface) {
    std::vector<std::unique_ptr<IDoF>> dofs(2);
    dofs[0] = std::make_unique<VertexDoF>(_mvertex);
    dofs[1] = std::make_unique<FaceDoF>(_mface);
    return dofs;
  }
public:
  explicit VertexFaceDoF(int _mvertex = 1, int _mface = 1) :
      MixedDoF(std::move(createDoFs(_mvertex, _mface))) {}

  std::string Name() const override { return "VertexFaceDoF"; }
};

#endif // BASICDOFS_HPP
