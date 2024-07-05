#ifndef SHAPEITERATOR_HPP
#define SHAPEITERATOR_HPP

template<typename DATA>
struct ShapeIterator {
protected:
  DATA data;
public:
  ShapeIterator(DATA &&data) : data(std::move(data)) {}

  ShapeIterator &operator++() {
    data.Incr();
    return *this;
  }

  bool operator==(const ShapeIterator &other) const { return data.id == other.data.id; }

  bool operator!=(const ShapeIterator &other) const { return data.id != other.data.id; }

  DATA &operator*() { return data; }
};

struct ShapeId {
  int id;

  ShapeId(int id) : id(id) {}
};

template<typename T>
struct IteratorProvider {
private:
  T &t;
  int n;
public:
  IteratorProvider(T &t, int n) : t(t), n(n) {}

  auto begin() { return t.begin(n); }

  auto end() { return t.end(n); }
};

#endif // SHAPEITERATOR_HPP
