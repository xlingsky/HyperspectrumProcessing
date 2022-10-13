#ifndef RASTER_FILTER_HPP
#define RASTER_FILTER_HPP

#include "raster/RasterOperator.h"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>

namespace xlingsky{
namespace raster{
namespace filter{

class MedianBlur : public FrameIterator {
 private:
  int _ksize;

 public:
  typedef float DataType;
  MedianBlur(int ksize) : _ksize(ksize) {}
  bool operator()(int , int , int , void* data, int cols, int rows) override {
    cv::Mat m(rows, cols, cv::DataType<DataType>::type, data);
    cv::medianBlur(m, m, _ksize);
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
  LinearFilter(int ksize, int stband, int dim = -1) : _ksize(ksize), _stband(stband) {
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
