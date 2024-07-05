#ifndef ELEMENTPOOL_HPP
#define ELEMENTPOOL_HPP

#include <queue>
#include "VectorMatrixBase.hpp"

template<typename ElementType>
class ElementPool {
  bool poolElements = false;
  std::unordered_map<Point, std::unique_ptr<ElementType>> elementMap{};
  std::queue<std::unique_ptr<ElementType>> elementQueue;
public:
  ElementPool() { Config::Get("ElementPool", poolElements); }

  ElementPool(bool poolElements) : poolElements(poolElements) {}

  void Initialize(const VectorMatrixBase &u) {
    Initialize(u, [](const VectorMatrixBase &u, const Cell &c) {
      return std::make_unique<ElementType>(u, c);
    });
  }

  void
  Initialize(const VectorMatrixBase &u,
             std::function<std::unique_ptr<ElementType>(const VectorMatrixBase &, const Cell &)>
                 elementGenerator) {
    if (poolElements) {
      for (cell c = u.cells(); c != u.cells_end(); ++c) {
        elementMap.try_emplace(c(), elementGenerator(u, *c));
      }
      for (cell c = u.overlap(); c != u.overlap_end(); ++c) {
        elementMap.try_emplace(c(), elementGenerator(u, *c));
      }
    }
  }

  const ElementType &Get(const Point &p) const { return *(elementMap.find(p)->second); }

  const ElementType &Get(
      const VectorMatrixBase &u, const Cell &c,
      std::function<std::unique_ptr<ElementType>(const VectorMatrixBase &, const Cell &)>
          elementGenerator = [](const VectorMatrixBase &u, const Cell &c) {
            return std::make_unique<ElementType>(u, c);
          }) {
    if (poolElements) {
      return Get(c());
    } else {
      if (elementQueue.size() > 1) { elementQueue.pop(); }
      elementQueue.emplace(elementGenerator(u, c));
      return *(elementQueue.back());
    }
  }

  bool Empty() const { return elementMap.empty(); }

  size_t Size() const { return elementMap.size(); }
};

template<typename FaceElementType>
class FaceElementPool {
  bool poolElements = false;
  std::unordered_map<Point, std::unordered_map<int, std::unique_ptr<FaceElementType>>>
      faceElementMap{};
  std::queue<std::unique_ptr<FaceElementType>> faceElementQueue;
public:
  FaceElementPool() { Config::Get("ElementPool", poolElements); }

  FaceElementPool(bool poolElements) : poolElements(poolElements) {}

  void Initialize(const VectorMatrixBase &u) {
    Initialize(u,
               [](const VectorMatrixBase &u, const Cell &c, int f) {
                 return std::make_unique<FaceElementType>(u, c, f);
               },
               {-1});
  }

  void Initialize(
      const VectorMatrixBase &u,
      std::function<std::unique_ptr<FaceElementType>(const VectorMatrixBase &, const Cell &, int f)>
          elementGenerator) {
    Initialize(u, elementGenerator, {-1});
  }

  void Initialize(
      const VectorMatrixBase &u,
      std::function<std::unique_ptr<FaceElementType>(const VectorMatrixBase &, const Cell &, int f)>
          elementGenerator,
      const std::vector<int> &bndTypes) {
    if (poolElements) {
      if (std::find(bndTypes.begin(), bndTypes.end(), -1) != bndTypes.end()) {
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
          for (int f = 0; f < c.Faces(); ++f) {
            faceElementMap[c()][f] = elementGenerator(u, *c, f);
          }
        }
        for (cell c = u.overlap(); c != u.overlap_end(); ++c) {
          for (int f = 0; f < c.Faces(); ++f) {
            faceElementMap[c()][f] = elementGenerator(u, *c, f);
          }
        }
      } else {
        for (cell c = u.cells(); c != u.cells_end(); ++c) {
          for (int f = 0; f < c.Faces(); ++f) {
            auto bndFace = u.find_bnd_face(c.Face(f));
            if (bndFace != u.bnd_faces_end()) {
              if (std::find(bndTypes.begin(), bndTypes.end(), bndFace.Part()) != bndTypes.end()) {
                faceElementMap[c()][f] = elementGenerator(u, *c, f);
              }
            }
          }
        }
        for (cell c = u.overlap(); c != u.overlap_end(); ++c) {
          for (int f = 0; f < c.Faces(); ++f) {
            auto bndFace = u.find_bnd_face(c.Face(f));
            if (bndFace != u.bnd_faces_end()) {
              if (std::find(bndTypes.begin(), bndTypes.end(), bndFace.Part()) != bndTypes.end()) {
                faceElementMap[c()][f] = elementGenerator(u, *c, f);
              }
            }
          }
        }
      }
    }
  }

  const FaceElementType &Get(const Point &p, int f) const {
    return *((faceElementMap.find(p)->second).find(f)->second);
  }

  const FaceElementType &Get(
      const VectorMatrixBase &u, const Cell &c, int f,
      std::function<std::unique_ptr<FaceElementType>(const VectorMatrixBase &, const Cell &, int f)>
          elementGenerator = [](const VectorMatrixBase &u, const Cell &c, int f) {
            return std::make_unique<FaceElementType>(u, c, f);
          }) {
    if (poolElements) {
      return Get(c(), f);
    } else {
      if (faceElementQueue.size() > 1) { faceElementQueue.pop(); }
      faceElementQueue.emplace(elementGenerator(u, c, f));
      return *(faceElementQueue.back());
    }
  }

  bool Empty() const { return faceElementMap.empty(); }

  size_t Size() const { return faceElementMap.size(); }
};

#endif // ELEMENTPOOL_HPP
