#ifndef XLINGSKY_RASTER_DETAIL_LOOKUP_HPP
#define XLINGSKY_RASTER_DETAIL_LOOKUP_HPP

#include <vector>
#include <limits>
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
    int i;
    if(x<0) i=0;
    else if(x>_table.size()-1) i = (int)_table.size()-1;
    else i = (int)x;
    return _table[i];
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

template<typename _SrcType, typename _DstType>
class LookupCreator{
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  typedef xlingsky::raster::Histogram<SrcType> HistType;
  typedef typename HistType::index_type IndexType;
 protected:
  SrcType _src_minimum;
  SrcType _src_unbounded_maximum;
  DstType _dst_minimum;
  DstType _dst_maximum;

  SrcType _src_step;
  float _cut_ratio_upper;
  float _cut_ratio_lower;
 public:
  LookupCreator(SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax)
      : _src_minimum(src_min), _src_unbounded_maximum(src_umax),
        _dst_minimum(dst_min), _dst_maximum(dst_umax),
        _src_step(2), _cut_ratio_upper(-1), _cut_ratio_lower(-1) {
  }
  void set_src_step(SrcType step) {
    _src_step = step;
  }
  void set_cut_ratio(float lower, float upper){
    _cut_ratio_lower = lower;
    _cut_ratio_upper = upper;
  }
  LookupLinear<SrcType, DstType>* Create() const {
    return Create(_src_minimum, _src_unbounded_maximum, _dst_minimum, _dst_maximum);
  }
  LookupLinear<SrcType, DstType>* Create(const SrcType* data, IndexType cols, IndexType rows, IndexType pixelspace, IndexType linespace) const{
    HistType hist(_src_minimum, _src_unbounded_maximum, _src_step);
    hist.adds(data, cols, rows, pixelspace, linespace);
    if(_cut_ratio_lower>=0 && _cut_ratio_upper>=0) hist.cut(_cut_ratio_lower, _cut_ratio_upper);
    return Create(hist.value(0), hist.value(hist.range().bins()-1), _dst_minimum, _dst_maximum);
  }
  LookupMap<SrcType, DstType>* Create(const SrcType* data, IndexType cols, IndexType rows, IndexType pixelspace, IndexType linespace, float clipratio) const {
    HistType hist(_src_minimum, _src_unbounded_maximum, _src_step);
    hist.adds(data, cols, rows, pixelspace, linespace);
    if(_cut_ratio_lower>=0 && _cut_ratio_upper>=0) hist.cut(_cut_ratio_lower, _cut_ratio_upper);
    size_t cliplimit = 0;
    if(clipratio>0){
      cliplimit = (size_t)(clipratio*cols*rows/hist.range().bins());
      if(cliplimit<1) cliplimit = 1;
    }else cliplimit = (std::numeric_limits<size_t>::max)();
    hist.redistribution(cliplimit);

    std::vector<DstType> table;
    {
      auto accum = hist.accumulation();
      LookupLinear<SrcType, DstType>* p = Create(accum.front(), accum.back(), _dst_minimum, _dst_maximum);
      table.resize(accum.size());
      for(int i=0; i<accum.size(); ++i)
        table[i] = p->operator[](accum[i]);
      delete p;
    }
    LookupMap<SrcType, DstType>* p = new LookupMap<SrcType, DstType>();
    p->set_domain(hist.range().bounded_minimum(), hist.range().unbounded_maximun(), hist.range().step());
    p->set_table(table.begin(), table.end());
    return p;
  }
 private:
  static LookupLinear<SrcType, DstType>* Create(SrcType src_min, SrcType src_max, DstType dst_min, DstType dst_max){
    DstType a = 0;
    if(src_max-src_min>std::numeric_limits<SrcType>::epsilon())
      a = (dst_max-dst_min)/(src_max-src_min);
    return new LookupLinear<SrcType, DstType>(a,dst_min-a*src_min);
  }
};

};
};
};

#endif
