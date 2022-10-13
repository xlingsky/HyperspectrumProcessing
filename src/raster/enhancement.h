#ifndef RASTER_ENHANCEMENT_H
#define RASTER_ENHANCEMENT_H

#include "RasterOperator.h"
#include <opencv2/imgproc.hpp>

namespace xlingsky{
namespace raster{
namespace enhancement {

class Destripe : public FrameIterator{
 public:
  typedef float DataType;
  Destripe() {}
  virtual ~Destripe(){}
  bool operator()(int b, int xoff, int yoff, void* data, int cols, int rows) override;
};

class Despike : public FrameIterator {
 public:
};

class Wallis: public FrameIterator {
 public:
};

class HistEqualization : public FrameIterator {
 public:
};

class Clahe : public FrameIterator {
 public:
};

};
}; // namespace raster
}; // namespace xlingsky

#endif
