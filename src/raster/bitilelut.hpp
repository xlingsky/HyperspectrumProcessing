#ifndef XLINGSKY_RASTER_BITILELUT_HPP
#define XLINGSKY_RASTER_BITILELUT_HPP

#include <algorithm>

#include "raster/operator.h"
#include "raster/detail/lookup.hpp"
#include "util/TileManager.hpp"

namespace xlingsky {
namespace raster {

namespace lookup{
template <typename _SrcType, typename _DstType>
class Creator {
public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  typedef xlingsky::raster::detail::Lookup<SrcType, DstType>* LookupPtr;
  typedef xlingsky::raster::Histogram<SrcType> HistType;
  typedef typename HistType::index_type IndexType;
 protected:
  SrcType _src_minimum;
  SrcType _src_unbounded_maximum;
  DstType _dst_minimum;
  DstType _dst_maximum;
 public:
  Creator(SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax)
      : _src_minimum(src_min), _src_unbounded_maximum(src_umax),
        _dst_minimum(dst_min), _dst_maximum(dst_umax) {};
  virtual ~Creator() {}
  
  SrcType src_minimum() const { return _src_minimum; }
  SrcType src_unbounded_maximum() const { return _src_unbounded_maximum;}
  DstType dst_minimum() const { return _dst_minimum; }
  DstType dst_maximum() const { return _dst_maximum; }
  
  void Destroy(LookupPtr p) {
    if(p!=nullptr) delete p;
  }
  bool Check(SrcType v){
    return _src_minimum<=v && v < _src_unbounded_maximum;
  }
  virtual LookupPtr Create(const SrcType* data, IndexType cols, IndexType rows, IndexType pixelspace, IndexType linespace){
    return Create(_src_minimum, _src_unbounded_maximum, _dst_minimum, _dst_maximum);
  }
 public:
  static detail::LookupLinear<SrcType, DstType>* Create(SrcType src_min, SrcType src_max, DstType dst_min, DstType dst_max){
    DstType a = 0;
    if(src_max-src_min>std::numeric_limits<SrcType>::epsilon())
      a = (dst_max-dst_min)/(src_max-src_min);
    return new detail::LookupLinear<SrcType, DstType>(a,dst_min-a*src_min);
  }
};

template <typename _SrcType, typename _DstType>
class LinearCreator : public Creator<_SrcType, _DstType>{
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  typedef Creator<_SrcType, _DstType> Base;
  typedef typename Base::LookupPtr LookupPtr;
  typedef typename Base::IndexType IndexType;
  typedef typename Base::HistType HistType;
 protected:
  SrcType _src_step;
  float _cut_ratio_upper;
  float _cut_ratio_lower;
 public:
  LinearCreator(SrcType src_step, SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax)
      : Base(src_min, src_umax, dst_min, dst_umax), _src_step(src_step) {}

  void set_cut_ratio(float lower, float upper){
    _cut_ratio_lower = lower;
    _cut_ratio_upper = upper;
  }
  LookupPtr Create(const SrcType* data, IndexType cols, IndexType rows, IndexType pixelspace, IndexType linespace) override{
    HistType hist = ComputeHist(data, cols, rows, pixelspace, linespace);
    if(hist.valid_count()==0) return nullptr;
    return Base::Create(hist.value(0), hist.value(hist.range().bins()-1), Base::_dst_minimum, Base::_dst_maximum);
  }
 protected:
  virtual HistType ComputeHist(const SrcType* data, IndexType cols, IndexType rows, IndexType pixelspace, IndexType linespace){
    HistType hist(Base::_src_minimum, Base::_src_unbounded_maximum, _src_step);
    hist.adds(data, cols, rows, pixelspace, linespace);
    if(_cut_ratio_lower>=0 && _cut_ratio_upper>=0) hist.cut(_cut_ratio_lower, _cut_ratio_upper);
    return hist;
  }
};

template <typename _SrcType, typename _DstType>
class ClaheCreator : public LinearCreator<_SrcType, _DstType>{
 protected:
  float _clip_ratio;
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  typedef LinearCreator<_SrcType, _DstType> Base;
  typedef typename Base::LookupPtr LookupPtr;
  typedef typename Base::IndexType IndexType;
  typedef typename Base::HistType HistType;
  ClaheCreator(float clipratio, SrcType src_step, SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax)
      : _clip_ratio(clipratio), Base(src_step, src_min, src_umax, dst_min, dst_umax) {}
  LookupPtr Create(const SrcType* data, IndexType cols, IndexType rows, IndexType pixelspace, IndexType linespace) override{
    HistType hist = ComputeHist(data, cols, rows, pixelspace, linespace);
    if(hist.valid_count()==0) return nullptr;

    std::vector<DstType> table;
    {
      auto accum = hist.accumulation();
      detail::LookupLinear<SrcType, DstType>* p = Creator<SrcType, DstType>::Create(accum.front(), accum.back(), Base::_dst_minimum, Base::_dst_maximum);
      table.resize(accum.size());
      for(int i=0; i<accum.size(); ++i)
        table[i] = p->operator[](accum[i]);
      Base::Destroy(p);
    }
    detail::LookupMap<SrcType, DstType>* p = new detail::LookupMap<SrcType, DstType>();
    p->set_domain(hist.range().bounded_minimum(), hist.range().unbounded_maximun(), hist.range().step());
    p->set_table(table.begin(), table.end());
    return p;
  }
 protected:
   HistType ComputeHist(const SrcType *data, IndexType cols,
                                IndexType rows, IndexType pixelspace,
                                IndexType linespace) override{
     HistType hist = Base::ComputeHist(data, cols, rows, pixelspace, linespace);
     size_t cliplimit = 0;
     if(_clip_ratio>0){
       cliplimit = (size_t)(_clip_ratio*cols*rows/hist.range().bins());
       if(cliplimit<1) cliplimit = 1;
     }else cliplimit = (std::numeric_limits<size_t>::max)();
     hist.redistribution(cliplimit);
     return hist;
   }
};
template <typename _SrcType, typename _DstType>
class WallisCreator : public Creator<_SrcType, _DstType> {
 public:
  typedef _SrcType SrcType;
  typedef _DstType DstType;
  typedef Creator<_SrcType, _DstType> Base;
  typedef typename Base::LookupPtr LookupPtr;
  typedef typename Base::IndexType IndexType;
 protected:
  SrcType _src_step;
  float _dst_mean;
  float _dst_std;
  float _factor_contrast;//0-1
  float _factor_bright;//0-1
 public:
  WallisCreator(SrcType src_step, SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax)
      : Base(src_min, src_umax, dst_min, dst_umax), _src_step(src_step), _factor_contrast(0.8), _factor_bright(0.9){
    _dst_mean = (Base::_dst_maximum+Base::_dst_minimum)/2;
    _dst_std = (Base::_dst_maximum-Base::_dst_minimum)/2;
  }
  void Setup(float dst_mean, float dst_std, float c, float b){
    _dst_mean = dst_mean;
    _dst_std = dst_std;
    _factor_contrast = c;
    _factor_bright = b;
  }
  LookupPtr Create(const SrcType *data, IndexType cols, IndexType rows,
                           IndexType pixelspace, IndexType linespace) override{
    typedef Range<SrcType> RangeType;
    RangeType range(Base::_src_minimum, Base::_src_unbounded_maximum, _src_step);
    struct _op{
      RangeType& _range;
      float _mean;
      float _std;
      SrcType _min;
      SrcType _max;
      IndexType _count;
      _op(RangeType& range) : _range(range), _mean(0), _std(0), _count(0) {
        _min = _range.unbounded_maximun();
        _max = _range.bounded_minimum();
      }
      void operator()(const void* data){
        const SrcType* pdata = (const SrcType*)data;
        if(*pdata < _range.bounded_minimum()) return;
        if(*pdata >= _range.unbounded_maximun()) return;
        _mean += (float)*pdata;
        _std += (float)*pdata**pdata;
        if(*pdata<_min) _min = *pdata;
        else if(*pdata>_max) _max = *pdata;
        ++_count;
      }
      void process(){
        if(_count==0) return;
        _mean /= _count;
        _std = _std/_count-_mean*_mean;
        if(_std<0) _std = 0;
        else _std = std::sqrt(_std);
        if(_min<_max) _range.reset(_min, _max, _range.step());
      }
    } op(range);
    xlingsky::raster::transform(data, cols, rows, pixelspace*sizeof(SrcType), linespace*sizeof(SrcType), op);
    op.process();

    if(op._count==0) return nullptr;
    float r1 = (_factor_contrast*_dst_std)/(_factor_contrast*op._std+_dst_std/_factor_contrast);
    float r0 = _factor_bright*_dst_mean+(1-_factor_bright-r1)*op._mean;

    std::vector<DstType> table;
    {
      table.resize(range.bins());
      for(int i=0; i<range.bins(); ++i){
        float v = r1*range.random(i)+r0;
        if(v<Base::_dst_minimum) table[i] = Base::_dst_minimum;
        else if(v>Base::_dst_maximum) table[i] = Base::_dst_maximum;
        else table[i] = (DstType)v;
      }
    }
    detail::LookupMap<SrcType, DstType>* p = new detail::LookupMap<SrcType, DstType>();
    p->set_domain(range.bounded_minimum(), range.unbounded_maximun(), range.step());
    p->set_table(table.begin(), table.end());
    return p;
  }
};

template <class Config, typename _SrcType, typename _DstType>
Creator<_SrcType, _DstType>* Create(const Config& config){
  float src_min = config. template get<float>("src_min", 7);
  float src_umax = config. template get<float>("src_umax", 4000);
  float dst_min = config. template get<float>("dst_min", 0);
  float dst_max = config. template get<float>("dst_max", 255);
  std::string m = config.template get<std::string>("mode", "WALLIS|MINMAX|CUT");
  std::transform(m.begin(), m.end(), m.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  int mode = 0;
  if (m.find("minmax") != std::string::npos)
    mode |= 0x01;
  if (m.find("cut") != std::string::npos)
    mode |= 0x02;
  if (m.find("hist") != std::string::npos)
    mode |= 0x04;
  if (m.find("wallis") != std::string::npos)
    mode |= 0x08;

  if ((mode & 0x02) || (mode & 0x04) || (mode & 0x08)) {
    float src_step = config. template get<float>("src_step", 1);
    if (mode & 0x08) {
      float dst_mean = config. template get<float>("dst_mean", 127);
      float dst_std = config. template get<float>("dst_std", 70); // 40-70
      float c = config. template get<float>("c", 0.8);            // 0-1
      float b = config. template get<float>("b", 0.9);            // 0-1
      auto *l = new WallisCreator<_SrcType, _DstType>( src_step, src_min, src_umax, dst_min, dst_max);
      l->Setup(dst_mean, dst_std, c, b);
      return l;
    } else {
      float cut_ratio_lower = config. template get<float>("cut_lower", 0.002);
      float cut_ratio_upper = config. template get<float>("cut_upper", 0.002);
      if (mode & 0x02) {
        auto *l = new LinearCreator<_SrcType, _DstType>( src_step, src_min, src_umax, dst_min, dst_max);
        l->set_cut_ratio(cut_ratio_lower, cut_ratio_upper);
        return l;
      } else {
        float hist_clip = config. template get<float>("hist_clip", 10);
        auto *l = new ClaheCreator<_SrcType, _DstType>( hist_clip, src_step, src_min, src_umax, dst_min, dst_max);
        l->set_cut_ratio(cut_ratio_lower, cut_ratio_upper);
        return l;
      }
    }
  } else {
    return new Creator<_SrcType, _DstType>(src_min, src_umax, dst_min, dst_max);
  }
}

};

class BiTileLut : public FrameIterator {
public:
  enum Mode{
    NONE = 0x00,
    GLOBAL = 0x01
  };
  typedef float SrcType;
  typedef SrcType DstType;
  using SegContainer = xlingsky::TileManager::SegContainer;
  typedef xlingsky::raster::detail::Lookup<SrcType, DstType>* LookupPtr;
  typedef xlingsky::raster::Histogram<SrcType>::index_type IndexType;
  typedef xlingsky::raster::lookup::Creator<SrcType, DstType>* LutCreatorPtr;
 protected:
  unsigned int _tile_col;
  unsigned int _tile_row;
  int _mode;
  unsigned int _resample_col_step;
  unsigned int _resample_row_step;

  LutCreatorPtr _lookup_creator;
  std::vector< LookupPtr > _tile_lookup;
  SegContainer _tile_colseg;
  SegContainer _tile_rowseg;
 protected:
  void Clear(){
    for(auto& p : _tile_lookup){
      _lookup_creator->Destroy(p);
      p = nullptr;
    }
    _tile_lookup.clear();
  }
  void Interpolate(SrcType* data, IndexType pixelspace, IndexType linespace, LookupPtr plu, LookupPtr pru, LookupPtr plb, LookupPtr prb, IndexType cols, IndexType rows){
    unsigned int num = cols*rows;
    for(int coef_y=0, coef_invy=rows; coef_y < rows; ++coef_y, --coef_invy, data+=linespace){
      SrcType* pr = data;
      for(int coef_x=0, coef_invx=cols; coef_x < cols; ++coef_x, --coef_invx, pr += pixelspace){
        if (_lookup_creator->Check(*pr)) {
          *pr =
              (DstType)((coef_invy *
                             (coef_invx * (*plu)[*pr] + coef_x * (*pru)[*pr]) +
                         coef_y *
                             (coef_invx * (*plb)[*pr] + coef_x * (*prb)[*pr])) /
                        num);
        }
      }
    }
  }
 public:
  BiTileLut(unsigned int tile_col, unsigned int tile_row, int mode = NONE) : _lookup_creator(nullptr), _tile_col(tile_col), _tile_row(tile_row), _mode(mode) {
    _resample_col_step = tile_col/16;
    if(_resample_col_step<1) _resample_col_step = 1;
    _resample_row_step = tile_row/16;
    if(_resample_row_step<1) _resample_row_step = 1;
  }
  virtual ~BiTileLut(){
    Clear();
    if (_lookup_creator) {
      delete _lookup_creator;
      _lookup_creator = nullptr;
    }
  }
  // template<class Config>
  // bool load_config(const Config& config){
  //     int resample_col_step = v.second.get<float>("resample_col_step", 3);
  //     int resample_row_step = v.second.get<float>("resample_row_step", 3);
  //     unsigned int tile_cols = v.second.get<unsigned int>("tile_cols", 0);
  //     unsigned int tile_rows = v.second.get<unsigned int>("tile_rows", 0);
  //     unsigned int grid_x = v.second.get<unsigned int>("grid_x", 8);
  //     unsigned int grid_y = v.second.get<unsigned int>("grid_y", 8);
  //     if (tile_cols == 0) {
  //       if (grid_x == 0)
  //         tile_cols = 512;
  //       else {
  //         tile_cols = src_size[store_prior[0]] / grid_x;
  //         if (tile_cols < 1)
  //           tile_cols = 1;
  //       }
  //     }
  //     if (tile_rows == 0) {
  //       if (grid_y == 0)
  //         tile_rows = 512;
  //       else {
  //         tile_rows = src_size[store_prior[1]] / grid_y;
  //         if (tile_rows < 1)
  //           tile_rows = 1;
  //       }
  //     }
  // }

  void set_lut_creator(LutCreatorPtr lookup_creator){
    if (_lookup_creator) {
      delete _lookup_creator;
      _lookup_creator = nullptr;
    }
    _lookup_creator = lookup_creator;
  }
  void set_resample_interval(int col, int row){
    _resample_col_step = col;
    _resample_row_step = row;
  }
  bool operator()(int , int , int , void* data, int cols, int rows) override {

    if(_tile_lookup.size()==0){
      _tile_colseg = TileManager::Tiling(cols, _tile_col, 0, _tile_col);
      _tile_rowseg = TileManager::Tiling(rows, _tile_row, 0, _tile_row);
      _tile_lookup.resize(_tile_colseg.size()*_tile_rowseg.size());

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
      for (int r = 0; r < _tile_rowseg.size(); ++r) {
        auto& tr = _tile_rowseg[r];
        auto p = _tile_lookup.data()+r*_tile_colseg.size();
        for (int c = 0; c < _tile_colseg.size(); ++c) {
          auto &tc = _tile_colseg[c];
          p[c] = _lookup_creator->Create(
              (const SrcType *)data + tr.first * cols + tc.first,
              std::max<IndexType>(1, tc.second/_resample_col_step), std::max<IndexType>(1, tr.second/_resample_row_step), _resample_col_step, cols*_resample_row_step);
        }
      }
    }

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for(int r=0; r<=_tile_rowseg.size(); ++r){
      SrcType* pdata = (SrcType*)data;
      unsigned int sub_y, upper_y, bottom_y;
      if(r==0){
        sub_y = _tile_rowseg[r].second >> 1;
        upper_y = bottom_y = 0;
      }else if(r==_tile_rowseg.size()){
        sub_y = (_tile_rowseg[r-1].second+1) >> 1;
        upper_y = bottom_y = r-1;
        pdata += (_tile_rowseg[r-1].first+(_tile_rowseg[r-1].second>>1))*cols;
      }else{
        sub_y = ((_tile_rowseg[r-1].second+1)>>1)+(_tile_rowseg[r].second>> 1);
        upper_y = r-1; bottom_y = r;
        pdata += (_tile_rowseg[r-1].first+(_tile_rowseg[r-1].second>>1))*cols;
      }
      upper_y *= _tile_colseg.size();
      bottom_y *= _tile_colseg.size();
      for(int c=0; c<=_tile_colseg.size(); ++c){
        unsigned int sub_x, left_x, right_x;
        if(c==0){
          sub_x = _tile_colseg[c].second >> 1;
          left_x = right_x = 0;
        }else if(c==_tile_colseg.size()){
          sub_x = (_tile_colseg[c-1].second+1) >> 1;
          left_x = right_x = c-1;
        }else{
          sub_x = ((_tile_colseg[c-1].second+1)>>1) +(_tile_colseg[c].second>>1);
          left_x = c-1; right_x = c;
        }
        LookupPtr p[4] = {
            _tile_lookup[upper_y + left_x], _tile_lookup[upper_y + right_x],
            _tile_lookup[bottom_y + left_x], _tile_lookup[bottom_y + right_x]};
        for(int i=1; i<4; ++i){
          if(p[i]==nullptr){
            p[i] = p[i-1];
          }else{
            for(int j=i-1; j>=0; --j){
              if(p[j]!=nullptr) break;
              p[j] = p[i];
            }
          }
        }
        if(p[0]) Interpolate(pdata, 1, cols, p[0], p[1], p[2], p[3], sub_x, sub_y);
        pdata += sub_x;
      }
    }

    if(!(_mode&GLOBAL)){
      Clear();
    }
    return true;
  }
};

typedef lookup::WallisCreator<BiTileLut::SrcType,BiTileLut::DstType> LutWallis;
typedef lookup::ClaheCreator<BiTileLut::SrcType,BiTileLut::DstType>  LutClahe;
typedef lookup::LinearCreator<BiTileLut::SrcType,BiTileLut::DstType> LutLinear;
typedef lookup::Creator<BiTileLut::SrcType,BiTileLut::DstType> LutCreator;

template<class Config>
LutCreator* LutCreate(const Config& config) {
  return lookup::Create<Config, BiTileLut::SrcType, BiTileLut::DstType>(config);
}

};  // namespace raster
};  // namespace xlingsky

#endif
