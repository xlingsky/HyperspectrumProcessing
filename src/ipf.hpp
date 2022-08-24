#ifndef IMAGE_PROCESS_FRAMEWORK_HPP
#define IMAGE_PROCESS_FRAMEWORK_HPP

#include <gdal_priv.h>
#include "gdal_traits.hpp"

#include <opencv2/core/core.hpp>
#include <numeric>

#include "TileManager.hpp"

namespace ipf{

class BufferOperator {
 public:
  BufferOperator() {}
  virtual ~BufferOperator() {}
  virtual bool operator()(void* data, int size[3], int space[3], int prior[3]) = 0;
};

struct RasterPatch{
  GDALDataset* dataset;
  int win[6];
  RasterPatch() : dataset(nullptr){
      memset(win, 0, sizeof(win));
  }
  void SetDataset(GDALDataset* d){
    dataset = d;
    if(dataset) SetWin(0, dataset->GetRasterXSize(), 0, dataset->GetRasterYSize(), 1, dataset->GetRasterCount());
  }
  void SetWin(int xoff, int xsize, int yoff, int ysize, int bandoff, int bandsize){
    win[0] = xoff;
    win[1] = yoff;
    win[2] = bandoff;
    win[3] = xsize;
    win[4] = ysize;
    win[5] = bandsize;
  }
};

class RasterOperator {
public:
    using DimSeg = xlingsky::TileManager::Seg;
protected:
  RasterPatch _src;
  RasterPatch _dst;
  GDALDataType _datatype;
  std::vector<char> _buffer;
  std::vector<int> _bandlist;
  int _storeorder[3];
  DimSeg _seg[3];
  BufferOperator* _op;
public:
    RasterOperator(BufferOperator* op, GDALDataType datatype = GDT_Byte) : _op(op), _datatype(datatype) {
        _storeorder[0] = 0;
        _storeorder[1] = 2;
        _storeorder[2] = 1;
    }
    void SetStoreOrder(int* order) {
        memcpy(_storeorder, order, sizeof(_storeorder));
    }
    void ReserveBufferSize(size_t size) {
        _buffer.resize(size*GetDataTypeSize());
    }
    int GetDataTypeSize() const {
        return GDALGetDataTypeSizeBytes(_datatype);
    }
  void SetSource(GDALDataset* src){
    _src.SetDataset(src);
  }
  void SetDestination(GDALDataset* dst){
    _dst.SetDataset(dst);
  }
  RasterPatch& source() { return _src; }
  RasterPatch& destination() { return _dst; }
  bool Begin(DimSeg& seg, int dim) {
      _seg[dim] = seg;
      if (dim == 2) {
          _bandlist.resize(seg.second);
		  std::iota(_bandlist.begin(), _bandlist.end(), _src.win[2] + seg.first);
      }
      return true;
  }
  bool End(DimSeg& seg, int dim) {
      return true;
  }
  bool Apply() {
	  int store_space[3];
	  store_space[_storeorder[0]] = GetDataTypeSize();
	  store_space[_storeorder[1]] = GetDataTypeSize()*_seg[_storeorder[0]].second;
	  store_space[_storeorder[2]] = store_space[_storeorder[1]] * _seg[_storeorder[1]].second;
    if(_src.dataset&&_src.dataset->RasterIO(GF_Read,
        _src.win[0]+_seg[0].first, _src.win[1]+_seg[1].first, _seg[0].second, _seg[1].second,
        &_buffer[0], _seg[0].second, _seg[1].second, _datatype, _seg[2].second, &_bandlist[0],
        store_space[0], store_space[1], store_space[2]) == CE_None) {
		  int store_size[3] = { _seg[0].second, _seg[1].second, _seg[2].second };
		  if (_op->operator()(&_buffer[0], store_size, store_space, _storeorder)) {
			  if (_dst.dataset ) {
				  int* dst_bandlist = nullptr;
                  if (store_size[2] != _seg[2].second) {
                      dst_bandlist = new int[store_size[2]];
					  std::iota(dst_bandlist, dst_bandlist+store_size[2], _dst.win[2] + _seg[2].first);
                  }
				  if (_dst.dataset->RasterIO(GF_Write,
					  _dst.win[0] + _seg[0].first, _dst.win[1] + _seg[1].first, store_size[0], store_size[1],
					  &_buffer[0], store_size[0], store_size[1], _datatype, store_size[2], dst_bandlist?dst_bandlist:&_bandlist[0],
					  store_space[0], store_space[1], store_space[2]) == CE_None) {
				  }
				  if (dst_bandlist) delete[] dst_bandlist;
			  }
		  }
	  }
    return true;
  }
};

template<class Operator>
bool DimProcessing(xlingsky::TileManager& manager, Operator* op, int dim) {
    for (int i = 0; i < manager.Size(dim); ++i) {
        auto seg = manager.Segment(dim, i);
        op->Begin(seg, dim);
        if (dim == 0) op->Apply();
        else DimProcessing(manager, op, dim-1);
        op->End(seg, dim);
    }
    return true;
}

template<class Operator, int dims = 3>
bool TileProcessing(int win_size[], int buffer_size[], Operator* op) {
    xlingsky::TileManager manager;
    using Seg = xlingsky::TileManager::Seg;
    for(int i=0; i<dims; ++i)
      manager.AppendDimension(win_size[i], buffer_size[i]);
    DimProcessing(manager, op, dims-1);
    return true;
}

template<class TileOp>
void Tile(int width, int height, int tilesize, int margin, TileOp& op){
    int tile_x_num = width/tilesize;
    int tile_y_num = height/tilesize;
    if(tile_x_num<1) tile_x_num = 1;
    if(tile_y_num<1) tile_y_num = 1;

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for(int ty=0; ty<tile_y_num; ++ty){
      int sty = ty*tilesize;
      int h = tilesize;
      if(sty+2*tilesize>height) h = height-sty;
      cv::Rect rc1,rc2;
      if(sty>margin) {
        rc1.y = sty-margin;
        rc2.y = margin;
      }else{
        rc1.y = 0;
        rc2.y = 0;
      }
      if(sty+h+margin>height){
        rc1.height = height-rc1.y;
      }else rc1.height = sty+h+margin-rc1.y;
      rc2.height = h;
      for(int tx=0; tx<tile_x_num; ++tx){
        int stx = tx*tilesize;
        int w = tilesize;
        if(stx+2*tilesize>width) w = width-stx;
        cv::Rect rc0(stx,sty,w,h);
        if(stx>margin) {
          rc1.x = stx-margin;
          rc2.x = margin;
        }else{
          rc1.x = 0;
          rc2.x = 0;
        }
        rc2.width = w;
        if(stx+w+margin>width) rc1.width = width-rc1.x;
        else rc1.width = stx+w+margin-rc1.x;
        op(rc1, rc2, rc0);
      }
    }
}

};

#endif
