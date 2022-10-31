#ifndef XLINGSKY_HISTOGRAM_HPP
#define XLINGSKY_HISTOGRAM_HPP

#include <deque>
#include <vector>
#include <type_traits>

namespace xlingsky {
namespace raster {
namespace detail {

template <typename T>
struct Range {
  typedef T value_type;
  value_type   	_bounded_minimum;
  value_type   	_unbounded_maximun;
  value_type 		_step;
  int _buckets;
  Range(value_type start = 0, value_type end = 256, value_type step = 1){
    reset(start, end, step);
  }
  void reset(value_type start, value_type end, value_type step){
    _bounded_minimum = start;
    _unbounded_maximun = end;
    _step = step;
    _buckets = int((end-start)/step);
  }
  int bucket(const value_type &v) const {
    return int((v-_bounded_minimum)/_step);
  }
  int buckets() const { return _buckets; }
  value_type bounded_minimum(int bucket = 0) const { return _bounded_minimum + bucket*_step; };
  value_type unbounded_maximun() const { return _unbounded_maximun; };
  value_type unbounded_maximun(int bucket) const { return bounded_minimum(bucket+1); };
  value_type step() const { return _step; }
};

template <typename T>
class Histogram{
 public:
  typedef T value_type;
  typedef const T* reference_type;
  typedef Range<T> Buckets;
  typedef int index_type;
  typedef Histogram<T> self;
 protected:
  std::vector<index_type> _count;
  Buckets _buckets;
  index_type _valid_count;
  index_type _under_minimum_count;
  index_type _beyond_maximum_count;
 public:
  Histogram(value_type start = 0, value_type end = 256, value_type step = 1)
      : _buckets(start, end, step), _valid_count(0), _under_minimum_count(0), _beyond_maximum_count(0) {
    _count.resize(_buckets.buckets(), 0);
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
  void add(reference_type orig){
    if(*orig<_buckets.bounded_minimum()) ++_under_minimum_count;
    else if(*orig>=_buckets.unbounded_maximun()) ++_beyond_maximum_count;
    else{
      add(_buckets.bucket(*orig), orig);
      ++_valid_count;
    }
  }
  void del(reference_type orig){
    if(*orig<_buckets.bounded_minimum()) --_under_minimum_count;
    else if(*orig>=_buckets.unbounded_maximun()) --_beyond_maximum_count;
    else{
      del(_buckets.bucket(*orig), orig);
      --_valid_count;
    }
  }
  int bucket(value_type& v) const { return _buckets.bucket(v); }
  int buckets() const { return _buckets.buckets(); }
  index_type count(int bucket) const { return _count[bucket]; }
  value_type median(value_type _default) const{
    index_type i = median(_valid_count);
    if(i==(index_type)-1) return _default;
    return random_elem(i);
  }
  index_type under_minimum_count() const { return _under_minimum_count; }
  index_type beyond_maximum_count() const { return _beyond_maximum_count; }
 protected:
  virtual void resize(value_type start, value_type end, value_type step) {
    _buckets.reset(start, end, step);
    _count.resize(_buckets.buckets());
  }
  virtual value_type random_elem (int bucket) const
  {
    const static float factor = _buckets.step()==1?0:0.5;//0.5;
    return _buckets.bounded_minimum(bucket)+(value_type)(factor*_buckets.step());
  }
  index_type median(int count) const{
    if(!count) return (index_type)-1;
    count = (count+1)/2;
    index_type i=0;
    while((count-=_count[i])>0) ++i;
    return i;
  }
  virtual void add(int bucket, reference_type){
    ++_count[bucket];
  }
  virtual void del(int bucket, reference_type){
    --_count[bucket];
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
      : base(start, end, step), _origs(base::buckets()){}
  void clear() override{
    for(auto& o : _origs)
      o.clear();
    base::clear();
  }
protected:
  void add(int bucket, reference_type orig) override{
    if(std::is_same<store_type, value_type>::value)
      _origs[bucket].push_back(*orig);
    else
      _origs[bucket].push_back(orig);
    base::add(bucket, orig);
  }
  void del(int bucket, reference_type orig) override{
    _origs[bucket].pop_front();
    base::del(bucket, orig);
  }
  value_type random_elem (int bucket) const override
  {
    auto& o = _origs[bucket];
    if(std::is_same<store_type, value_type>::value)
      return o[rand()%o.size()];
    else
      return *(o[rand()%o.size()]);
  }
  void resize(value_type start, value_type end, value_type step) override{
    base::resize(start, end, step);
    _origs.resize(base::buckets());
  }
};

};
};  // namespace raster
};  // namespace xlingsky

#endif
