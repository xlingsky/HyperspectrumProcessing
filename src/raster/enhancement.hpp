#ifndef XLINGSKY_RASTER_ENHANCEMENT_H
#define XLINGSKY_RASTER_ENHANCEMENT_H

#include "raster/operator.h"
#include <opencv2/imgproc.hpp>
#include "raster/detail/despeckle.hpp"
#include "raster/detail/lookup.hpp"
#include "TileManager.hpp"

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

class Render : public FrameIterator {
 public:
  typedef float SrcType;
  typedef SrcType DstType;
  enum Mode{
    MINMAX = 0x01,
    CLIP = 0x02,
    HIST_EQU = 0x04,
    GLOBAL = 0x08
  };
 protected:
  int _mode;

  int _hist_col_step;
  int _hist_row_step;
  float _hist_clip_ratio;

  xlingsky::raster::detail::LookupCreator<SrcType, DstType> _lookup_creator;
  xlingsky::raster::detail::Lookup<SrcType, DstType>* _lookup;

 public:
  Render(SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax, int mode = MINMAX|CLIP)
      :  _lookup_creator(src_min, src_umax, dst_min, dst_umax),
         _mode(mode), _lookup(nullptr),
         _hist_col_step(3), _hist_row_step(3), _hist_clip_ratio(80){}
  virtual ~Render(){
    if(_lookup) delete _lookup;
  }
  void set_src_step(SrcType step) {
    _lookup_creator.set_src_step(step);
  }
  void set_cut_ratio(float lower, float upper){
    _lookup_creator.set_cut_ratio(lower, upper);
  }
  void set_hist_interval(int col, int row){
    _hist_col_step = col;
    _hist_row_step = row;
  }
  void set_hist_clip_ratio(float ratio){
    _hist_clip_ratio = ratio;
  }
  bool operator()(int , int , int , void* data, int cols,
      int rows) override {
    if(_lookup==nullptr){
      if(_mode&HIST_EQU){
        _lookup = _lookup_creator.Create((const SrcType*)data, MAX(1, cols/_hist_col_step), MAX(1, rows/_hist_row_step), _hist_col_step, cols*_hist_row_step, _hist_clip_ratio);
      }else if(_mode&CLIP){
        _lookup = _lookup_creator.Create((const SrcType*)data, MAX(1, cols/_hist_col_step), MAX(1, rows/_hist_row_step), _hist_col_step, cols*_hist_row_step);
      }else{
        _lookup = _lookup_creator.Create();
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
  enum Mode{
    MINMAX = 0x01,
    WALLIS = 0x02,
    CLAHE = 0x04,
    GLOBAL = 0x08
  };
 protected:
  typedef float SrcType;
  typedef SrcType DstType;
  using SegContainer = xlingsky::TileManager::SegContainer;
  typedef xlingsky::raster::detail::Lookup<SrcType, DstType>* LookupPtr;
  unsigned int _tile_col;
  unsigned int _tile_row;
  int _hist_col_step;
  int _hist_row_step;
  float _hist_clip_ratio;
  xlingsky::raster::detail::LookupCreator<SrcType, DstType> _lookup_creator;
  std::vector< LookupPtr > _tile_lookup;
  SegContainer _tile_colseg;
  SegContainer _tile_rowseg;

  float _dst_mean;
  float _dst_std;
  float _c;
  float _b;

  int _mode;
 protected:
  void ClearLookup(){
    for(auto& p : _tile_lookup){
      delete p;
      p = nullptr;
    }
    _tile_lookup.clear();
  }
  void Interpolate(SrcType* data, int pixelspace, int linespace, LookupPtr plu, LookupPtr pru, LookupPtr plb, LookupPtr prb, unsigned int cols, unsigned int rows){
    unsigned int num = cols*rows;
    for(int coef_y=0, coef_invy=rows; coef_y < rows; ++coef_y, --coef_invy, data+=linespace){
      SrcType* pr = data;
      for(int coef_x=0, coef_invx=cols; coef_x < cols; ++coef_x, --coef_invx, pr += pixelspace){
        if(*pr < _lookup_creator.src_minimum() ) continue;
        if(*pr >= _lookup_creator.src_unbounded_maximum()) continue;
        *pr = (SrcType)((coef_invy*(coef_invx*(*plu)[*pr]+coef_x*(*pru)[*pr])+coef_y*(coef_invx*(*plb)[*pr]+coef_x*(*prb)[*pr]))/num);
      }
    }
  }
 public:
  Clahe(
      SrcType src_min, SrcType src_umax, DstType dst_min, DstType dst_umax,
      unsigned int tile_col, unsigned int tile_row, float clipratio, int mode)
      : _lookup_creator(src_min, src_umax, dst_min, dst_umax), _tile_col(tile_col), _tile_row(tile_row), _hist_clip_ratio(clipratio), _mode(mode), _dst_std(0), _dst_mean(0), _c(-1), _b(-1) {
    _hist_col_step = tile_col/16;
    if(_hist_col_step<1) _hist_col_step = 1;
    _hist_row_step = tile_row/16;
    if(_hist_row_step<1) _hist_row_step = 1;
  }
  virtual ~Clahe(){
    ClearLookup();
  }
  void set_src_step(SrcType step) {
    _lookup_creator.set_src_step(step);
  }
  void set_cut_ratio(float lower, float upper){
    _lookup_creator.set_cut_ratio(lower, upper);
  }
  void set_hist_interval(int col, int row){
    _hist_col_step = col;
    _hist_row_step = row;
  }
  void set_hist_clip_ratio(float ratio){
    _hist_clip_ratio = ratio;
  }
  void set_wallis_pars(float dst_mean, float dst_std, float c, float b){
    _dst_mean = dst_mean;
    _dst_std = dst_std;
    _c = c;
    _b = b;
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
        if(_mode&CLAHE){
          for (int c = 0; c < _tile_colseg.size(); ++c) {
            auto& tc = _tile_colseg[c];
            p[c] = _lookup_creator.Create(
                (const SrcType*)data + tr.first*cols+tc.first,
                MAX(1, tc.second/_hist_col_step), MAX(1, tr.second/_hist_row_step),
                _hist_col_step, cols*_hist_row_step, _hist_clip_ratio);
          }
        }else if(_mode&WALLIS){
          for (int c = 0; c < _tile_colseg.size(); ++c) {
            auto& tc = _tile_colseg[c];
            p[c] = _lookup_creator.Create(
                (const SrcType*)data + tr.first*cols+tc.first,
                MAX(1, tc.second/_hist_col_step), MAX(1, tr.second/_hist_row_step),
                _hist_col_step, cols*_hist_row_step, _dst_mean, _dst_std, _c, _b);
          }
        }else {
          for (int c = 0; c < _tile_colseg.size(); ++c) {
            auto& tc = _tile_colseg[c];
            p[c] = _lookup_creator.Create(
                (const SrcType*)data + tr.first*cols+tc.first,
                MAX(1, tc.second/_hist_col_step), MAX(1, tr.second/_hist_row_step),
                _hist_col_step, cols*_hist_row_step);
          }
        }
      }
    }

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for(int r=0; r<=_tile_rowseg.size(); ++r){//
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
        auto plu = _tile_lookup[upper_y*_tile_colseg.size()+left_x];
        auto pru = _tile_lookup[upper_y*_tile_colseg.size()+right_x];
        auto plb = _tile_lookup[bottom_y*_tile_colseg.size()+left_x];
        auto prb = _tile_lookup[bottom_y*_tile_colseg.size()+right_x];
        Interpolate(pdata, 1, cols, plu, pru, plb, prb, sub_x, sub_y);
        pdata += sub_x;
      }
    }

    if(!(_mode&GLOBAL)){
      ClearLookup();
    }
    return true;
  }
};

};
}; // namespace raster
}; // namespace xlingsky

#endif
