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
    win[1] = xsize;
    win[2] = yoff;
    win[3] = ysize;
    win[4] = bandoff;
    win[5] = bandsize;
  }
};

class Framework {
 public:
  Framework() : _datatype(GDT_Byte), _datatype_size(sizeof(unsigned char)) {
  }
  virtual ~Framework() {}
  void SetDataType(GDALDataType datatype) { _datatype = datatype; }
  void SetDataTypeSize(int size) { _datatype_size = size; }
  int GetDataTypeSize() const { return _datatype_size; }
  GDALDataType GetDataType() const { return _datatype; }

  void SetSource(GDALDataset* src){
    _src.SetDataset(src);
  }
  void SetDestination(GDALDataset* dst){
    _dst.SetDataset(dst);
  }

  bool Apply(int buffer_size[3], int store_prior[3], BufferOperator* op) {
    xlingsky::TileManager manager;
    using Seg = xlingsky::TileManager::Seg;
    for(int i=0; i<3; ++i)
      manager.AppendDimension(_src.win[2*i+1], buffer_size[i]);
    std::vector<char> buffer((size_t)buffer_size[0]*buffer_size[1]*buffer_size[2]*_datatype_size);
    int store_space[3];
    Seg segs[3];

    for(int b=0; b<manager.Size(2); ++b){
      segs[2] = manager.Segment(2, b);
      std::vector<int> bandlist(segs[2].second);
      std::iota(bandlist.begin(), bandlist.end(), _src.win[4]+segs[2].first);
      for(int r=0; r<manager.Size(1); ++r){
        segs[1] = manager.Segment(1, r);
        for(int c=0; c<manager.Size(0); ++c){
          segs[0] = manager.Segment(0, c);
          store_space[store_prior[0]] = _datatype_size;
          store_space[store_prior[1]] = _datatype_size*segs[store_prior[0]].second;
          store_space[store_prior[2]] = store_space[store_prior[1]]*segs[store_prior[1]].second;
          if(Preprocessing(&buffer[0], segs[0].first, segs[1].first, segs[0].second, segs[1].second, segs[2].second, &bandlist[0], store_space[0], store_space[1], store_space[2])){
            int store_size[3] = {segs[0].second, segs[1].second, segs[2].second};
            if(op->operator()(&buffer[0], store_size, store_space, store_prior)){
              Postprocessing(&buffer[0], segs[0].first, segs[1].first, segs[0].second, segs[1].second, segs[2].second, &bandlist[0], store_space[0], store_space[1], store_space[2]);
            }
          }
        }
      }
    }
    return true;
  }
  virtual bool Preprocessing(void* data, int xoff, int yoff, int xsize, int ysize, int bandcount, int* bandlist, int xspace, int yspace, int bandspace){
    if(_src.dataset&&_src.dataset->RasterIO(GF_Read, _src.win[0]+xoff, _src.win[2]+yoff, xsize, ysize, data, xsize, ysize, GetDataType(), bandcount, bandlist, xspace, yspace, bandspace) != CE_None)
      return false;
    return true;
  }
  virtual bool Postprocessing(void* data, int xoff, int yoff, int xsize, int ysize, int bandcount, int* bandlist, int xspace, int yspace, int bandspace){
    if(_dst.dataset&&_dst.dataset->RasterIO(GF_Write, _dst.win[0]+xoff, _dst.win[2]+yoff, xsize, ysize, data, xsize, ysize, GetDataType(), bandcount, bandlist, xspace, yspace, bandspace) != CE_None)
      return false;
    return true;
  }
 protected:
  int _datatype_size;
  GDALDataType _datatype;
  RasterPatch _src;
  RasterPatch _dst;
};

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
