#ifndef XLINGSKY_RASTER_DETAIL_INPAINT_HPP
#define XLINGSKY_RASTER_DETAIL_INPAINT_HPP

#include <opencv2/core/core.hpp>
#include <opencv2/photo.hpp>

namespace xlingsky{
namespace raster{
namespace detail{

class InpaintOp{
 private:
  double _radius;
  int _flags;
  cv::Mat& _src;
  cv::Mat& _dst;
  cv::Mat& _mask;
 public:
  InpaintOp(cv::Mat& src, cv::Mat& dst, cv::Mat& mask) : _src(src), _dst(dst), _mask(mask), _radius(3), _flags(cv::INPAINT_TELEA){
  }
  bool operator()(cv::Rect& src_tile, cv::Rect& aoi, cv::Rect dst_tile) {
    cv::Mat tile;//(src_tile.height, src_tile.width, _src.type() );
    cv::inpaint(_src(src_tile), _mask(src_tile), tile, _radius, _flags);
    tile(aoi).copyTo(_dst(dst_tile));
    return true;
  }
};

};
};
};

#endif
