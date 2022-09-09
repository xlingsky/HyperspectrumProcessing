#ifndef RASTER_OPERATOR_H
#define RASTER_OPERATOR_H

#include <vector>
#include <numeric>

namespace xlingsky {

namespace raster {

class Operator {
 public:
  Operator() {}
  virtual ~Operator() {}
  virtual bool operator()(void* data, int imoff[3], int size[3], int space[3],
                          int prior[3]) = 0;
};

class ComboOperator : public Operator {
 public:
  virtual ~ComboOperator() {
    for (auto& op : _ops) delete op;
    _ops.clear();
  }
  void Add(Operator* op) { _ops.push_back(op); }
  bool operator()(void* data, int imoff[3], int size[3], int space[3],
                  int prior[3]) override {
    for (auto& op : _ops)
      if (!op->operator()(data, imoff, size, space, prior)) return false;
    return true;
  }
  size_t size() const { return _ops.size(); }

 protected:
  std::vector<Operator*> _ops;
};

class OperatorEx : public Operator {
 public:
  OperatorEx() : _nodata(std::numeric_limits<double>::min()) {}
  void SetNoDataValue(double d) { _nodata = d; }
  double GetNoDataValue() const { return _nodata; }
  template<typename T>
  bool IsNoData(T& d) const {
    return (T)_nodata == d;
  }
  template<typename T>
  int FindValid(T* data, int n, int& st, int& ed){
    st = 0;
    ed = n-1;
    while(IsNoData(data[st]) && st<n) ++st;
    while(st<=ed && IsNoData(data[ed])) --ed;
    return ed-st+1;
  }
 private:
  double _nodata;
};

class FrameIterator : public OperatorEx {
 public:
  bool operator()(void* data, int imoff[3], int size[3], int space[3],
                  int prior[3]) override {
    bool ret = true;
    for (int r = 0; r < size[prior[2]]; ++r)
      if (this->operator()(r+imoff[prior[2]], imoff[prior[0]], imoff[prior[1]], (char*)data + (size_t)space[prior[2]] * r,
                           size[prior[0]], size[prior[1]])) {
      }
    return ret;
  }
  virtual bool operator()(int b, int xoff, int yoff, void* data, int cols, int rows) = 0;
};


};

};

#endif
