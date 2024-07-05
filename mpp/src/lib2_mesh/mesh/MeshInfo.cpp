#include "MeshInfo.hpp"

std::string MeshInfoOnProcs::ProcString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.proc) + " ";
  str += "]";
  return str;
}

std::string MeshInfoOnProcs::CellCountString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.cellCount) + " ";
  str += "]";
  return str;
}

std::string MeshInfoOnProcs::VertexCountString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.vertexCount) + " ";
  str += "]";
  return str;
}

std::string MeshInfoOnProcs::EdgesCountString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.edgesCount) + " ";
  str += "]";
  return str;
}

std::string MeshInfoOnProcs::FacesCountString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.facesCount) + " ";
  str += "]";
  return str;
}

std::string MeshInfoOnProcs::BNDFacesCountString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.bndFaceCount) + " ";
  str += "]";
  return str;
}

std::string MeshInfoOnProcs::ProcSetsCountString() {
  std::string str = "[ ";
  for (auto info : *this)
    str += to_string(info.procSetsCount) + " ";
  str += "]";
  return str;
}

Buffer &operator<<(Buffer &buffer, const MeshInfoOnProc &info) {
  buffer << info.proc << info.cellCount << info.vertexCount << info.facesCount << info.edgesCount
         << info.procSetsCount << info.bndFaceCount;
  return buffer;
}

Buffer &operator>>(Buffer &buffer, MeshInfoOnProc &info) {
  buffer >> info.proc >> info.cellCount >> info.vertexCount >> info.facesCount >> info.edgesCount
      >> info.procSetsCount >> info.bndFaceCount;
  return buffer;
}
