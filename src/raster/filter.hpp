#ifndef XLINGSKY_RASTER_FILTER_HPP
#define XLINGSKY_RASTER_FILTER_HPP

#include <opencv2/core/core.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/imgproc.hpp>

#include "raster/operator.h"
#include "raster/detail/despeckle.hpp"
namespace xlingsky{
namespace raster{
namespace filter{

class MedianBlur : public FrameIterator {
 private:
  int _ksize;

 public:
  typedef float DataType;
  MedianBlur() {}
  MedianBlur(int ksize) : _ksize(ksize) {}
  bool operator()(int , int , int , void* data, int cols, int rows) override {
    cv::Mat m(rows, cols, cv::DataType<DataType>::type, data);
    cv::medianBlur(m, m, _ksize);
    return true;
  }
  template<class Config>
  bool load_config(const Config& config){
    _ksize = config.template get<int>("ksize", 3);
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
        if (!hist.is_valid()) continue;

        int pos = x*pixelspace + (y * linespace);

        if (std::isnan(src[pos]) || std::isinf(src[pos])){
          dst[pos] = src[pos];
          continue;
        }

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
  template<class Config>
  bool load_config(const Config& config){
    _radius = config.template get<int>("ksize", 3);
    _minimum = config.template get<float>("hist_min", 7);
    _unbounded_maximum = config.template get<float>("hist_max", 4000);
    _step = config.template get<float>("hist_step", 1);
    return true;
  }
};

class LinearFilter : public FrameIterator {
 private:
  int _ksize;
  int _stband;
  cv::Mat _kernel;
  cv::Mat _kernel_ex;

 public:
  typedef float DataType;
  LinearFilter() {}
  LinearFilter(int ksize, int stband, int dim = -1) {
    reset(ksize, stband, dim);
  }
  template<class Config>
  bool load_config(Config& config){
      int band = config.template get<int>("band", 0);
      int ksize = config.template get<int>("ksize", 3);
      int dim = config.template get<int>("dim", 1);

      reset(ksize, band, dim);

      return true;
  }

  void reset(int ksize, int stband, int dim) {
    _ksize = ksize;
    _stband = stband;
    if(dim==1)
      _kernel = cv::getGaussianKernel(_ksize, 0, cv::DataType<DataType>::type);
    else{
      _kernel_ex = cv::getGaussianKernel(_ksize, 0, cv::DataType<DataType>::type);
      _kernel = _kernel_ex.t();
      if(dim==0)
        _kernel_ex.release();
    }
  }

  bool operator()(int , int, int, void* data, int cols, int rows) override {
    cv::Mat o(rows, cols, cv::DataType<DataType>::type, data);
    cv::Mat m = o(cv::Rect(0, _stband, cols, rows - _stband));
    cv::filter2D(m, m, m.depth(), _kernel);
    if(!_kernel_ex.empty())
      cv::filter2D(m, m, m.depth(), _kernel_ex);
    return true;
  }
};

};
};
};

#endif
