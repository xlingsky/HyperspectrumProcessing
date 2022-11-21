#ifndef XLINGSKY_RASTER_BITILELUT_HPP
#define XLINGSKY_RASTER_BITILELUT_HPP

#include "raster/operator.h"
#include "raster/detail/lookup.hpp"
#include "TileManager.hpp"

namespace xlingsky {
namespace raster {

class BiTileLut : public FrameIterator {
 public:
 protected:
  typedef float SrcType;
  typedef SrcType DstType;
  using SegContainer = xlingsky::TileManager::SegContainer;
  typedef xlingsky::raster::detail::Lookup<SrcType, DstType>* LookupPtr;

  unsigned int _tile_col;
  unsigned int _tile_row;
  int _mode;

  xlingsky::raster::detail::LookupCreator<SrcType, DstType>* _lookup_creator;
  std::vector< LookupPtr > _tile_lookup;
  SegContainer _tile_colseg;
  SegContainer _tile_rowseg;
 protected:
  void Clear(){
    for(auto& p : _tile_lookup){
      delete p;
      p = nullptr;
    }
    _tile_lookup.clear();
    if (_lookup_creator) {
      delete _lookup_creator;
      _lookup_creator = nullptr;
    }
  }
  void Interpolate(SrcType* data, int pixelspace, int linespace, LookupPtr plu, LookupPtr pru, LookupPtr plb, LookupPtr prb, unsigned int cols, unsigned int rows){
    unsigned int num = cols*rows;
    for(int coef_y=0, coef_invy=rows; coef_y < rows; ++coef_y, --coef_invy, data+=linespace){
      SrcType* pr = data;
      for(int coef_x=0, coef_invx=cols; coef_x < cols; ++coef_x, --coef_invx, pr += pixelspace){
        if (_lookup_creator->Check(*pr)) {
          *pr =
              (DstType)((coef_invy *
                             (coef_invx * (*plu)[*pr] + coef_x * (*pru)[*pr]) +
                         coef_y *
                             (coef_invx * (*plb)[*pr] + coef_x * (*prb)[*pr])) /
                        num);
        }
      }
    }
  }
 public:
  BiTileLut(unsigned int tile_col, unsigned int tile_row, int mode) : _lookup_creator(nullptr), _tile_col(tile_col), _tile_row(tile_row), _mode(mode) {}
  virtual 
};

};  // namespace raster
};  // namespace xlingsky

#endif