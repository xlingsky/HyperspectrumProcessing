#ifndef XLINGSKY_RASTER_ENHANCEMENT_H
#define XLINGSKY_RASTER_ENHANCEMENT_H

#include <opencv2/imgproc.hpp>

#include "raster/operator.h"
#include "raster/detail/lookup.hpp"
#include "raster/bitilelut.hpp"

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

class Render : public FrameIterator {
 public:
  typedef float SrcType;
  typedef SrcType DstType;
  typedef xlingsky::raster::lookup::Creator<SrcType, DstType>* LutCreatorPtr;
  enum Mode{
    NONE = 0x00,
    GLOBAL = 0x08
  };
 protected:
  int _mode;
  unsigned int _resample_col_step;
  unsigned int _resample_row_step;

  LutCreatorPtr _lookup_creator;
  xlingsky::raster::detail::Lookup<SrcType, DstType>* _lookup;

 public:
  Render(int mode = NONE)
      :  _lookup_creator(nullptr),
         _mode(mode), _lookup(nullptr)
         {}
  virtual ~Render(){
    if(_lookup)
      _lookup_creator->Destroy(_lookup);
    if(_lookup_creator)
      delete _lookup_creator;
  }
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
  bool operator()(int , int , int , void* data, int cols,
      int rows) override {
    if(_lookup==nullptr){
      _lookup = _lookup_creator->Create((const SrcType*)data, MAX(1, cols/_resample_col_step), MAX(1, rows/_resample_row_step), _resample_col_step, cols*_resample_row_step);
    }

    transform(data, cols, rows, sizeof(SrcType), cols*sizeof(SrcType), *_lookup);

    if(!(_mode&GLOBAL)){
      _lookup_creator->Destroy(_lookup);
      _lookup = nullptr;
    }
    return true;
  }
};

};
}; // namespace raster
}; // namespace xlingsky

#endif
