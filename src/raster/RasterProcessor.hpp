#ifndef RASTER_PROCESSOR_HPP
#define RASTER_PROCESSOR_HPP

#include "RasterOperator.h"
#include "TileManager.hpp"

#include <gdal_priv.h>
#include <numeric>
#include "gdal_traits.hpp"

#ifdef _LOGGING
#include <glog/logging.h>
#define _LOG_LEVEL_RASTERPROCESSOR 3
#endif

namespace xlingsky {

namespace raster {

//all indices starting from 0
struct Patch {
  GDALDataset* dataset;
  int win[6];
  Patch() : dataset(nullptr) { memset(win, 0, sizeof(win)); }
  void SetDataset(GDALDataset* d) {
    dataset = d;
    if (dataset)
      SetWin(0, dataset->GetRasterXSize(), 0, dataset->GetRasterYSize(), 0,
             dataset->GetRasterCount());
  }
  void SetWin(int xoff, int xsize, int yoff, int ysize, int bandoff,
              int bandsize) {
    win[0] = xoff;
    win[1] = yoff;
    win[2] = bandoff;
    win[3] = xsize;
    win[4] = ysize;
    win[5] = bandsize;
  }
};

class Processor {
 public:
  using DimSeg = xlingsky::TileManager::Seg;
  using IndexType = xlingsky::TileManager::IndexType;

 protected:
  Patch _src;
  Patch _dst;
  GDALDataType _datatype;
  std::vector<char> _buffer;
  std::vector<int> _bandlist;
  int _storeorder[3];
  DimSeg _seg[3];
  Operator* _op;
  int _dim_flushcache;

 public:
  Processor(Operator* op, GDALDataType datatype = GDT_Byte, int flushcache = -1)
      : _op(op), _datatype(datatype), _dim_flushcache(flushcache) {
    _storeorder[0] = 0;
    _storeorder[1] = 2;
    _storeorder[2] = 1;
  }
  void SetStoreOrder(int* order) {
    memcpy(_storeorder, order, sizeof(_storeorder));
  }
  void ReserveBufferSize(size_t size) {
    _buffer.resize(size * GetDataTypeSize());
  }
  void SetDataType(GDALDataType datatype) { _datatype = datatype; }
  int GetDataTypeSize() const { return GDALGetDataTypeSizeBytes(_datatype); }
  void SetSource(GDALDataset* src) { _src.SetDataset(src); }
  void SetDestination(GDALDataset* dst) { _dst.SetDataset(dst); }
  Patch& source() { return _src; }
  Patch& destination() { return _dst; }
  bool Begin(DimSeg& seg, int dim) {
    _seg[dim] = seg;
    if (dim == 2) {
      _bandlist.resize(seg.second);
      std::iota(_bandlist.begin(), _bandlist.end(), _src.win[2] + seg.first + 1);
    }
    return true;
  }
  bool End(DimSeg& seg, int dim) { 
    if(dim==_dim_flushcache) {
      if(_src.dataset) _src.dataset->FlushCache();
      if(_dst.dataset) _dst.dataset->FlushCache();
    }
    return true; 
  }
  bool Apply() {
    int store_space[3];
    store_space[_storeorder[0]] = GetDataTypeSize();
    store_space[_storeorder[1]] =
        GetDataTypeSize() * _seg[_storeorder[0]].second;
    store_space[_storeorder[2]] =
        store_space[_storeorder[1]] * _seg[_storeorder[1]].second;

#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RASTERPROCESSOR) << "RasterIO reading...";
#endif
    if (_src.dataset &&
        _src.dataset->RasterIO(
            GF_Read, _src.win[0] + _seg[0].first, _src.win[1] + _seg[1].first,
            _seg[0].second, _seg[1].second, &_buffer[0], _seg[0].second,
            _seg[1].second, _datatype, _seg[2].second, &_bandlist[0],
            store_space[0], store_space[1], store_space[2]) == CE_None) {
      IndexType store_size[3] = {_seg[0].second, _seg[1].second, _seg[2].second};
      IndexType imoff[3] = {_src.win[0] + _seg[0].first, _src.win[1] + _seg[1].first, _src.win[2]+_seg[2].first};
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RASTERPROCESSOR) << "Operators starting...";
#endif
      if (_op->operator()(&_buffer[0], imoff, store_size, store_space, _storeorder)) {
        if (_dst.dataset) {
          int* dst_bandlist = nullptr;
          if (store_size[2] != _seg[2].second) {
            dst_bandlist = new int[store_size[2]];
            std::iota(dst_bandlist, dst_bandlist + store_size[2],
                      _dst.win[2] + _seg[2].first + 1);
          }
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RASTERPROCESSOR) << "RasterIO writing..." ;
#endif
          if (_dst.dataset->RasterIO(
                  GF_Write, _dst.win[0] + _seg[0].first,
                  _dst.win[1] + _seg[1].first, store_size[0], store_size[1],
                  &_buffer[0], store_size[0], store_size[1], _datatype,
                  store_size[2], dst_bandlist ? dst_bandlist : &_bandlist[0],
                  store_space[0], store_space[1], store_space[2]) == CE_None) {
          }
          if (dst_bandlist) delete[] dst_bandlist;
        }
      }
    }
    return true;
  }
};

};  // namespace raster

};  // namespace xlingsky

#endif
