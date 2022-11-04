#ifndef XLINGSKY_RASTER_ENHANCEMENT_H
#define XLINGSKY_RASTER_ENHANCEMENT_H

#include "raster/operator.h"
#include <opencv2/imgproc.hpp>
#include "raster/detail/despeckle.hpp"
#include "raster/detail/lookup.hpp"

namespace xlingsky{
namespace raster{
namespace enhancement {

class Destripe : public FrameIterator{
 public:
  typedef float DataType;
  typedef double HistType;
 private:
  int _tile_width;
  HistType _hist_mean;
  DataType _v_min;
  DataType _v_max;
 public:
  Destripe(int tile_size) : _tile_width(tile_size), _hist_mean(std::numeric_limits<HistType>::min()), 
      _v_min(0), _v_max(std::numeric_limits<DataType>::max()) {}
  virtual ~Destripe(){}
  bool operator()(int , int , int , void* data, int cols, int rows) override {
  std::vector<HistType> hist(rows, 0), corr(rows);
    DataType *pdata = (DataType *)data;

  /*
   * collect "histogram" data.
   */
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
  for (int r=0; r < rows; ++r) {
    DataType *t = pdata + r*cols;
    for (int c=0; c < cols; c++) {
      hist[r] += (HistType)t[c];
    }
   }

  /*
   * average out histogram
   */
  {
     int extend = _tile_width / 2;
     HistType *h = &hist[0] - extend;
     HistType *c = &corr[0] - extend;
     HistType sum = 0;
     int cnt = 0;

     for (int x = -extend; x < rows; ++x) {
       if (x + extend < rows) {
         sum += h[extend];
         cnt++;
       }

       if (x - extend >= 0) {
         sum -= h[-extend];
         cnt--;
       }

       if (x >= 0) {
         if (*h) {
             //*c = ((sum / cnt - *h) << 10) / *h;
             *c = (sum / cnt - *h) / *h;
         }
         else
           *c = std::numeric_limits<HistType>::max();
       }

       ++h;
       ++c;
     }
  }

  /*
   * remove stripes.
   */
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
  for (int r=0; r < rows; ++r) {
    DataType *src = pdata + r * cols;
    DataType *dst = src;
    if (_hist_mean==std::numeric_limits<HistType>::min()) {
      for (int c = 0; c < cols; ++c) {
        HistType v = src[c] + (src[c] * corr[r]);//src[c]+(src[c]*corr[r]>>10)
        if (v < _v_min)
          dst[c] = _v_min;
        else if (v > _v_max)
          dst[c] = _v_max;
        else
          dst[c] = (DataType)v;
      }
    } else {
      for (int c = 0; c < cols; ++c) {
        HistType v = _hist_mean + (src[c] * corr[r]);//_hist_mean+(src[c]*corr[r]>>10)
        if (v < _v_min)
          dst[c] = _v_min;
        else if (v > _v_max)
          dst[c] = _v_max;
        else
          dst[c] = (DataType)v;
      }
    }
  }

  return true;
  }
};

class Despeckle : public FrameIterator {
 public:
  typedef float DataType;
  typedef xlingsky::raster::enhancement::detail::DespeckleHistogram<DataType> HistType;
  enum FilterType{
    ADAPTIVE = 0x01,
    RECURSIVE = 0x02
  };
 protected:
  int _radius;
  int _filter_type;
  DataType _minimum;
  DataType _unbounded_maximum;
  DataType _step;
 public:
  Despeckle(int radius = 3, DataType start = 7, DataType end = 4000, DataType step = 1, int filter_type = ADAPTIVE)
      : _radius(radius), _filter_type(filter_type), _minimum(start), _unbounded_maximum(end), _step(step) {}
  virtual ~Despeckle(){}
  bool operator()(int , int , int , void* data, int cols, int rows) override {
    int adapt_radius = _radius;
    int pixelspace = 1;
    int linespace = cols;
    std::vector<DataType> buffer;
    DataType* pdata = (DataType*)data;
    DataType* src = pdata;
    DataType* dst = nullptr;
    if(_filter_type & RECURSIVE){
      dst = src;
    }else{
      buffer.resize((size_t)cols*rows);
      dst = &buffer[0];
    }
    for (int y = 0; y < rows; y++) {
      int x = 0;
      int ymin = MAX (0, y - adapt_radius);
      int ymax = MIN (rows - 1, y + adapt_radius);
      int xmin = MAX (0, x - adapt_radius);
      int xmax = MIN (cols - 1, x + adapt_radius);
      HistType hist(_minimum, _unbounded_maximum, _step);
      hist.set_rect(xmin, ymin, xmax, ymax);
      hist.adds(src, pixelspace, linespace, xmin, ymin, xmax, ymax);

      for (x = 0; x < cols; x++)
      {
        ymin = MAX (0, y - adapt_radius); /* update ymin, ymax when adapt_radius changed (FILTER_ADAPTIVE) */
        ymax = MIN (rows- 1, y + adapt_radius);
        xmin = MAX (0, x - adapt_radius);
        xmax = MIN (cols - 1, x + adapt_radius);

        hist.update(src, pixelspace, linespace, xmin, ymin, xmax, ymax);

        int pos = x*pixelspace + (y * linespace);
        auto pixel = hist.median (src[pos]);

        if (_filter_type & RECURSIVE)
        {
          hist.del (src+pos);
          hist.add (&pixel);
        }

        dst[pos] = pixel;
        /*
         * Check the histogram and adjust the diameter accordingly...
         */
        if (_filter_type & ADAPTIVE)
        {
          if (hist.under_minimum_count() >= adapt_radius || hist.beyond_maximum_count() >= adapt_radius)
          {
            if (adapt_radius < _radius)
              adapt_radius++;
          }
          else if (adapt_radius > 1)
          {
            adapt_radius--;
          }
        }
      }
    }
    if(dst!=pdata)
      memcpy(pdata, dst, sizeof(DataType)*cols*rows);
    return true;
  }
};

class Wallis: public FrameIterator {
 public:
  Wallis() {}
  virtual ~Wallis(){}
  bool operator()(int , int , int , void* data, int cols,
      int rows) override {
    return true;
  }
};

class Render : public FrameIterator {
 public:
  typedef float SrcType;
  typedef SrcType DstType;
  typedef xlingsky::raster::Histogram<SrcType> HistType;
  enum Mode{
    MINMAX = 0x01,
    CLIP = 0x02,
    HIST_EQU = 0x04,
    GLOBAL = 0x08
  };
 protected:
  SrcType _src_minimum;
  SrcType _src_unbounded_maximum;
  DstType _dst_minimum;
  DstType _dst_maximum;
  int _mode;

  SrcType _src_step;
  int _hist_col_step;
  int _hist_row_step;
  float _cut_ratio_upper;
  float _cut_ratio_lower;

  xlingsky::raster::detail::Lookup<SrcType, DstType>* _lookup;

 public:
  Render(SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax, int mode = MINMAX|CLIP)
      :  _src_minimum(src_min), _src_unbounded_maximum(src_umax),
         _dst_minimum(dst_min), _dst_maximum(dst_umax),
         _mode(mode), _lookup(nullptr),
         _src_step(1), _hist_col_step(3), _hist_row_step(3), _cut_ratio_upper(0.002), _cut_ratio_lower(0.002) {}
  virtual ~Render(){
    if(_lookup) delete _lookup;
  }
  void set_src_step(SrcType step) {
    _src_step = step;
  }
  void set_hist_interval(int col, int row){
    _hist_col_step = col;
    _hist_row_step = row;
  }
  void set_clip_ratio(float lower, float upper){
    _cut_ratio_lower = lower;
    _cut_ratio_upper = upper;
  }
  bool operator()(int , int , int , void* data, int cols,
      int rows) override {
    if(_lookup==nullptr){
      if((_mode&CLIP) || (_mode&HIST_EQU)){
        HistType hist(_src_minimum, _src_unbounded_maximum, _src_step);
        hist.adds((const SrcType*)data, MAX(1, cols/_hist_col_step), MAX(1, rows/_hist_row_step), _hist_col_step, cols*_hist_row_step);
        if(_mode&CLIP) hist.cut(_cut_ratio_lower, _cut_ratio_upper);
        if(_mode&HIST_EQU){
          std::vector<DstType> table;
          {
            auto accum = hist.accumulation();
            auto minimum = accum.front();
            auto maximum = accum.back();
            DstType a = 0;
            if(maximum-minimum>std::numeric_limits<SrcType>::epsilon())
              a = (_dst_maximum-_dst_minimum)/(maximum-minimum);
            xlingsky::raster::detail::LookupLinear<SrcType, DstType> p(a, (DstType)(_dst_minimum-a*minimum));
            table.resize(accum.size());
            for(int i=0; i<accum.size(); ++i)
              table[i] = p[accum[i]];
          }
          xlingsky::raster::detail::LookupMap<SrcType, DstType>* p = new xlingsky::raster::detail::LookupMap<SrcType, DstType>();
          p->set_buckets(_src_minimum, _src_unbounded_maximum, _src_step);
          p->set_table(table.begin(), table.end());
          _lookup = p;
        }else{
          SrcType minimum = hist.value(0);
          SrcType maximum = hist.value(hist.range().buckets()-1);
          DstType a = 0;
          if(maximum-minimum>std::numeric_limits<SrcType>::epsilon())
            a = (_dst_maximum-_dst_minimum)/(maximum-minimum);
          _lookup = new xlingsky::raster::detail::LookupLinear<SrcType, DstType>(a,_dst_minimum-a*minimum);
        }
      }else{
        assert((_src_unbounded_maximum-_src_minimum)>std::numeric_limits<SrcType>::epsilon());
        auto a = (_dst_maximum-_dst_minimum)/(_src_unbounded_maximum-_src_minimum);
        _lookup = new xlingsky::raster::detail::LookupLinear<SrcType, DstType>(a,_dst_minimum-a*_src_minimum);
      }
    }

    transform(data, cols, rows, sizeof(SrcType), cols*sizeof(SrcType), *_lookup);

    if(!(_mode&GLOBAL)){
      delete _lookup;
      _lookup = nullptr;
    }
    return true;
  }
};

class Clahe : public FrameIterator {
 public:
  Clahe() {}
  virtual ~Clahe(){}
  bool operator()(int b, int xoff, int yoff, void* data, int cols,
      int rows) override {
    return true;
  }
};

};
}; // namespace raster
}; // namespace xlingsky

#endif
