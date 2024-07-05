#ifndef MESHINFO_HPP
#define MESHINFO_HPP

#include "Mesh.hpp"

struct MeshInfoOnProc {
  int proc = 0;
  int cellCount = 0;
  int vertexCount = 0;
  int facesCount = 0;
  int edgesCount = 0;
  int procSetsCount = 0;
  int bndFaceCount = 0;

  MeshInfoOnProc() {}

  MeshInfoOnProc(const Mesh &mesh) :
      proc(PPM->Proc(mesh.CommSplit())), cellCount(mesh.CellCount()),
      vertexCount(mesh.VertexCount()), facesCount(mesh.FaceCount()), edgesCount(mesh.EdgeCount()),
      procSetsCount(mesh.ProcSetsCountWithoutInfty()), bndFaceCount(mesh.BoundaryFaceCount()) {}
};

struct MeshInfoOnProcs : public std::vector<MeshInfoOnProc> {
  std::string ProcString();

  std::string CellCountString();

  std::string VertexCountString();

  std::string EdgesCountString();

  std::string FacesCountString();

  std::string BNDFacesCountString();

  std::string ProcSetsCountString();
};

Buffer &operator<<(Buffer &buffer, const MeshInfoOnProc &info);

Buffer &operator>>(Buffer &buffer, MeshInfoOnProc &info);

#endif // MESHINFO_HPP
