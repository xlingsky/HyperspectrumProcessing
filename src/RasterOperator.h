#ifndef RASTER_OPERATOR_H
#define RASTER_OPERATOR_H

#include <vector>

namespace xlingsky {

namespace raster {

class Operator {
 public:
  Operator() {}
  virtual ~Operator() {}
  virtual bool operator()(void* data, int size[3], int space[3],
                          int prior[3]) = 0;
};

class ComboOperator : public Operator {
 public:
  virtual ~ComboOperator() {
    for (auto& op : _ops) delete op;
    _ops.clear();
  }
  void Add(Operator* op) { _ops.push_back(op); }
  bool operator()(void* data, int size[3], int space[3],
                  int prior[3]) override {
    for (auto& op : _ops)
      if (!op->operator()(data, size, space, prior)) return false;
    return true;
  }
  size_t size() const { return _ops.size(); }

 protected:
  std::vector<Operator*> _ops;
};


class FrameIterator : public Operator {
 public:
  bool operator()(void* data, int size[3], int space[3],
                  int prior[3]) override {
    bool ret = true;
    for (int r = 0; r < size[prior[2]]; ++r)
      if (this->operator()(r, (char*)data + (size_t)space[prior[2]] * r,
                           size[prior[0]], size[prior[1]])) {
      }
    return ret;
  }
  virtual bool operator()(int r, void* data, int cols, int rows) = 0;
};


};

};

#endif