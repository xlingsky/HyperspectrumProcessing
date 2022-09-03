#ifndef RASTER_BLUR_HPP
#define RASTER_BLUR_HPP

#include "RasterOperator.h"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>

namespace xlingsky{
namespace raster{
namespace blur{

class MedianBlur : public FrameIterator {
 private:
  int _ksize;

 public:
  typedef float DataType;
  MedianBlur(int ksize) : _ksize(ksize) {}
  bool operator()(int r, void* data, int cols, int rows) override {
    cv::Mat m(rows, cols, cv::DataType<DataType>::type, data);
    cv::medianBlur(m, m, _ksize);
    return true;
  }
};

class GaussianBlur : public FrameIterator {
 private:
  int _ksize;
  int _stband;

 public:
  typedef float DataType;
  GaussianBlur(int ksize, int stband) : _ksize(ksize), _stband(stband) {}
  bool operator()(int r, void* data, int cols, int rows) override {
    cv::Mat o(rows, cols, cv::DataType<DataType>::type, data);
    cv::Mat m = o(cv::Rect(0, _stband, cols, rows - _stband));
    cv::GaussianBlur(m, m, cv::Size(_ksize, _ksize), 0);
    return true;
  }
};


};
};
};

#endif
