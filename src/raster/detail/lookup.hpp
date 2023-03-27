#ifndef XLINGSKY_RASTER_DETAIL_LOOKUP_HPP
#define XLINGSKY_RASTER_DETAIL_LOOKUP_HPP

#include <vector>
#include <limits>
#include <cmath>
#include "raster/histogram.hpp"

namespace xlingsky{
namespace raster{
namespace detail{

template<typename _SrcType, typename _DstType>
class Lookup {
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  Lookup(){}
  virtual ~Lookup(){}
  virtual DstType operator[](SrcType x) const = 0;
  void operator()(void* p){
    SrcType* data = (SrcType*)p;
    *data = (SrcType)this->operator[](*data);
  }
};

template<typename _SrcType, typename _DstType>
class LookupLinear : public Lookup<_SrcType, _DstType>{
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
 protected:
  DstType _a;
  DstType _b;
 public:
  LookupLinear(DstType a, DstType b) : _a(a), _b(b) {}
  virtual ~LookupLinear() {}
  DstType operator[](SrcType x) const override{
    return _a*x+_b;
  }
};

template<typename _SrcType, typename _DstType>
class LookupTable : public Lookup<_SrcType, _DstType>{
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
 protected:
  std::vector<DstType> _table;
 public:
  LookupTable(){}
  template<class Iterator>
  LookupTable(Iterator first, Iterator last) : _table(first, last) {}
  DstType operator[](SrcType x) const override{
    // int i;
    // if(x<0) i=0;
    // else if(x>_table.size()-1) i = (int)_table.size()-1;
    // else i = (int)x;
    return _table[(int)x];
  }
  template<class Iterator>
  void set_table(Iterator first, Iterator last){
    _table.assign(first, last);
  }
};

template<typename _SrcType, typename _DstType>
class LookupMap : public LookupTable<_SrcType, _DstType>{
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  typedef Range<double> Domain;
  typedef LookupTable<_SrcType, _DstType> Base;
 protected:
  Domain _domain;
 public:
  LookupMap() {}
  void set_domain(double start , double end , double step ){
    _domain.reset(start, end, step);
  }
  DstType operator[](SrcType x) const override{
    return Base::operator[](_domain.bin(x));
  }
};

};
};
};

#endif
