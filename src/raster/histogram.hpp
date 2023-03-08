#ifndef XLINGSKY_RASTER_HISTOGRAM_HPP
#define XLINGSKY_RASTER_HISTOGRAM_HPP

//may use boost histogram instead in future

#include <deque>
#include <vector>
#include <type_traits>

#include "raster/transform.hpp"

namespace xlingsky {
namespace raster {

template <typename T>
struct Range {
  typedef T value_type;
  value_type   	_bounded_minimum;
  value_type   	_unbounded_maximun;
  value_type 		_step;
  int _bins;
  Range(value_type start = 0, value_type end = 256, value_type step = 1){
    reset(start, end, step);
  }
  void reset(value_type start, value_type end, value_type step){
    _bounded_minimum = start;
    _unbounded_maximun = end;
    _step = step;
    _bins = int((end-start)/step);
  }
  int bin(value_type v) const {
    return int((v-_bounded_minimum)/_step);
  }
  int bins() const { return _bins; }
  value_type bounded_minimum(int bin = 0) const { return _bounded_minimum + bin*_step; };
  value_type unbounded_maximun() const { return _unbounded_maximun; };
  value_type unbounded_maximun(int bin) const { return bounded_minimum(bin+1); };
  value_type step() const { return _step; }
  value_type random(int bin) const{
    const static float factor = _step==1?0:0.5;//0.5;
    return bounded_minimum(bin)+(value_type)(factor*_step);
  }
};

template <typename T>
class Histogram{
 public:
  typedef T value_type;
  typedef const T* reference_type;
  typedef Range<T> range_type;
  typedef int index_type;
  typedef Histogram<T> self;
 protected:
  std::vector<index_type> _count;
  range_type _range;
  index_type _valid_count;
  index_type _under_minimum_count;
  index_type _beyond_maximum_count;
 public:
  Histogram(value_type start = 0, value_type end = 256, value_type step = 1)
      : _range(start, end, step), _valid_count(0), _under_minimum_count(0), _beyond_maximum_count(0) {
    _count.resize(_range.bins(), 0);
 }
  virtual void clear() {
    for(auto& c : _count)
      c = 0;
    _valid_count = _under_minimum_count = _beyond_maximum_count = 0;
  }
  void reset(value_type start, value_type end, value_type step) {
    resize(start, end, step);
    clear();
  }
  void cut( float cut_lower, float cut_upper){
    int cut_lower_count = (int)(_valid_count*cut_lower);
    int cut_upper_count = (int)(_valid_count*cut_upper);

    int st = 0, ed = (int)_count.size()-1;
    while( st<_count.size() && (cut_lower_count-=_count[st])>=0) ++st;
    while(ed>=0 && (cut_upper_count-=_count[ed])>=0) --ed;
    cut(st, ed);
  }
  void redistribution(index_type cliplimit){
    int clipped = 0;
    for (auto& c : _count) {
      if (c > cliplimit) {
        clipped += c - cliplimit;
        c = cliplimit;
      }
    }

    // redistribute clipped pixels
    int redistBatch = clipped / _count.size();
    int residual = clipped - redistBatch * _count.size();

    for (auto& c : _count) c += redistBatch;

    if (residual != 0) {
      int residualStep = MAX(_count.size() / residual, 1);
      for (int i = 0; i < _count.size() && residual > 0;
           i += residualStep, residual--)
        _count[i]++;
    }
    return;

    size_t excess = 0;
    for(const auto& c : _count) {
      long bin_excess = (long)c-(long)cliplimit;
      if(bin_excess>0) excess += bin_excess;
    }
    if(excess==0) return;

    size_t bin_increment = excess/_count.size();
    size_t upper = cliplimit-bin_increment; /* Bins larger than upper set to cliplimit */

    for(auto& c : _count){
      if(c>cliplimit) c = cliplimit;
      else if(c > upper){
        excess -= (cliplimit-c);
        c = cliplimit;
      }else{
        excess -= bin_increment;
        c += bin_increment;
      }
    }

    while(excess){ /* Redistribute remaining excess  */
      for(int i=0; excess && i<_count.size(); ++i){
        index_type step = _count.size()/excess;
        if(step<1) step = 1;
        for(int j=i; j < _count.size() && excess; j += step ){
          if(_count[j]<cliplimit){
            ++_count[j]; --excess;
          }
        }
      }
    }

  }
  void add(reference_type orig){
    if(*orig<_range.bounded_minimum()) ++_under_minimum_count;
    else if(*orig>=_range.unbounded_maximun()) ++_beyond_maximum_count;
    else{
      add(_range.bin(*orig), orig);
      ++_valid_count;
    }
  }
  void del(reference_type orig){
    if(*orig<_range.bounded_minimum()) --_under_minimum_count;
    else if(*orig>=_range.unbounded_maximun()) --_beyond_maximum_count;
    else{
      del(_range.bin(*orig), orig);
      --_valid_count;
    }
  }
#define BATCH_FUNCTION(name)                                            \
  void name##s(reference_type src, index_type cols, index_type rows, index_type pixelspace, index_type linespace) { \
    if(cols<=0 || rows<=0 ) return;                                     \
    struct _##name{                                                     \
      self* _p;                                                         \
      _##name(self* p) : _p(p) {}                                       \
      void operator()(const void* data){                                \
        _p->name((reference_type)data);                                 \
      }                                                                 \
    } op(this);                                                         \
    xlingsky::raster::transform(src, cols, rows, pixelspace*sizeof(value_type), linespace*sizeof(value_type), op); \
  }

  BATCH_FUNCTION(add);
  BATCH_FUNCTION(del);

  int bin(value_type& v) const { return _range.bin(v); }
  const range_type& range() const { return _range; }
  index_type count(int bin) const { return _count[bin]; }
  value_type median(value_type _default) const{
    index_type i = median(_valid_count);
    if(i==(index_type)-1) return _default;
    return random_elem(i);
  }
  value_type value(int bin) const {
    return random_elem(bin);
  }
  index_type under_minimum_count() const { return _under_minimum_count; }
  index_type beyond_maximum_count() const { return _beyond_maximum_count; }
  index_type valid_count() const { return _valid_count; }

  std::vector<size_t> accumulation() const {
    std::vector<size_t> accum;
    if(_count.size()>0){
      accum.resize(_count.size());
      accum[0] = _count[0];
      for(int i=1; i<_count.size(); ++i){
        accum[i] = accum[i-1]+_count[i];
      }
    }
    return accum;
  }
 protected:
  virtual void resize(value_type start, value_type end, value_type step) {
    _range.reset(start, end, step);
    _count.resize(_range.bins());
  }
  virtual value_type random_elem (int bin) const
  {
    return _range.random(bin);
  }
  index_type median(int count) const{
    if(!count) return (index_type)-1;
    count = (count+1)/2;
    index_type i=0;
    while((count-=_count[i])>0) ++i;
    return i;
  }
  virtual void add(int bin, reference_type){
    ++_count[bin];
  }
  virtual void del(int bin, reference_type){
    --_count[bin];
  }
  virtual void cut(int st,int ed){
    if (st > ed) return;
    int count = 0;
    for(int i=0; i<st; ++i) count += _count[i];
    _under_minimum_count += count;
    _valid_count -= count;
    count = 0;
    for(int i=(int)_count.size()-1; i>ed; --i) count += _count[i];
    _beyond_maximum_count += count;
    _valid_count -= count;

    if(st==0){
      if(ed!=_count.size()-1){
        _count.resize(ed+1);
        _range._unbounded_maximun = _range.unbounded_maximun(ed);
        _range._bins = (int)_count.size();
      }
    }else{
      for(int i=st; i<=ed; ++i)
        _count[i-st] = _count[i];
      if(ed!=_count.size()-1) _range._unbounded_maximun = _range.unbounded_maximun(ed);
      _range._bounded_minimum = _range.bounded_minimum(st);
      _count.resize(ed-st+1);
      _range._bins = (int)_count.size();
    }
  }
};

template <typename T, typename S = const T*>
class RefHistogram : public Histogram<T>
{
 public:
  typedef Histogram<T> base;
  typedef typename base::value_type value_type;
  typedef typename base::reference_type reference_type;
  typedef typename base::index_type index_type;
  typedef S store_type;
 protected:
  std::vector< std::deque<store_type> > _origs;
 public:
  RefHistogram(value_type start = 0, value_type end = 256, value_type step = 1)
      : base(start, end, step){
    _origs.resize(base::range().bins());
  }
  void clear() override{
    for(auto& o : _origs)
      o.clear();
    base::clear();
  }
protected:
  void add(int bin, reference_type orig) override{
    if(std::is_same<store_type, value_type>::value)
      _origs[bin].push_back(*orig);
    else
      _origs[bin].push_back(orig);
    base::add(bin, orig);
  }
  void del(int bin, reference_type orig) override{
    _origs[bin].pop_front();
    base::del(bin, orig);
  }
  value_type random_elem (int bin) const override
  {
    auto& o = _origs[bin];
    if(std::is_same<store_type, value_type>::value)
      return o[rand()%o.size()];
    else
      return *(o[rand()%o.size()]);
  }
  void resize(value_type start, value_type end, value_type step) override{
    base::resize(start, end, step);
    _origs.resize(base::range().bins());
  }
};

};  // namespace raster
};  // namespace xlingsky

#endif
