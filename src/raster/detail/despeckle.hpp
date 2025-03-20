#ifndef XLINGSKY_RASTER_ENHANCEMENT_DETAIL_DESPECKLE_HPP
#define XLINGSKY_RASTER_ENHANCEMENT_DETAIL_DESPECKLE_HPP

#include "raster/histogram.hpp"

namespace xlingsky{
namespace raster{
namespace enhancement {
namespace detail{

template<typename T>
class DespeckleHistogram : public xlingsky::raster::Histogram<T>{
 public:
  typedef xlingsky::raster::Histogram<T> base;
  typedef typename base::value_type value_type;
  typedef typename base::reference_type reference_type;
  typedef typename base::index_type index_type;
 protected:
  int       _xmin;
  int       _ymin;
  int       _xmax;
  int       _ymax; /* Source rect */
 public:
  DespeckleHistogram(value_type start = 0, value_type end = 256, value_type step = 1)
      : base(start, end, step){}
  void set_rect(int xmin, int ymin, int xmax, int ymax){
    _xmin = xmin;
    _ymin = ymin;
    _xmax = xmax;
    _ymax = ymax;
  }
  bool is_valid() const { return base::valid_count() > 0; }
  void adds(reference_type src, int pixelspace, int linespace, int xmin,
            int ymin, int xmax, int ymax) {
    if (xmin > xmax)
      return;

    base::adds(src+xmin*pixelspace+ymin*linespace, xmax-xmin+1, ymax-ymin+1, pixelspace, linespace);
  }
  void dels(reference_type src, int pixelspace, int linespace, int xmin,
            int ymin, int xmax, int ymax) {
    if (xmin > xmax)
      return;

    base::dels(src+xmin*pixelspace+ymin*linespace, xmax-xmin+1, ymax-ymin+1, pixelspace, linespace);
  }
  void update(reference_type src, int pixelspace, int linespace, int xmin,
              int ymin, int xmax, int ymax){
    /* assuming that radious of the box can change no more than one
       pixel in each call */
    /* assuming that box is moving either right or down */

    dels (src, pixelspace, linespace, _xmin, _ymin, xmin - 1, _ymax);
    dels (src, pixelspace, linespace, xmin, _ymin, xmax, ymin - 1);
    dels (src, pixelspace, linespace, xmin, ymax + 1, xmax, _ymax);

    adds (src, pixelspace, linespace, _xmax + 1, ymin, xmax, ymax);
    adds (src, pixelspace, linespace, xmin, ymin, _xmax, _ymin - 1);
    adds (src, pixelspace, linespace, _xmin, _ymax + 1, _xmax, ymax);

    _xmin = xmin;
    _ymin = ymin;
    _xmax = xmax;
    _ymax = ymax;
  }
};

};
};
}; // namespace raster
}; // namespace xlingsky

#endif
