#include "Buffer.hpp"

Buffer::Buffer(size_t m) : n(m), max_n(0) {
  if (m) b = new char[m];
  else b = 0;
  p = b;
}

Buffer::~Buffer() { Destruct(); }

void Buffer::Destruct() {
  if (b) {
    delete[] b;
    b = 0;
    max_n = n;
  }
  p = 0;
  n = 0;
}

size_t Buffer::size() const { return size_t(p - b); }

size_t Buffer::Size() const { return n; }

void Buffer::rewind() {
  if (max_n > 0 && !b) resize(max_n);
  p = b;
}

char *Buffer::operator()() { return b; }

void Buffer::resize(size_t m) {
  int s = size();
  if (s && b) {
    char *tmp = new char[m];
    memcpy(tmp, b, s);
    delete[] b;
    p = tmp + s;
    b = tmp;
  } else {
    b = new char[m];
    p = b;
  }
  n = m;
}
