#ifndef RASTER_ENHANCEMENT_H
#define RASTER_ENHANCEMENT_H

#include "RasterOperator.h"

namespace xlingsky{
namespace raster{
namespace enhancement {

class Destripe : public FrameIterator{
 public:
  typedef float DataType;
  Destripe() {}
  virtual ~Destripe(){}
  bool operator()(int r, void* data, int cols, int rows) override;
};

};
}; // namespace raster
}; // namespace xlingsky

#endif
