#ifndef IMAGE_PROCESS_FRAMEWORK_HPP
#define IMAGE_PROCESS_FRAMEWORK_HPP

#include <gdal_priv.h>
#include "gdal_traits.hpp"

#include <opencv2/core/core.hpp>

namespace ipf{

template<class RowOp>
void RowMajor(GDALDataset* src, GDALDataset* dst, RowOp& op){
  using DataType = typename RowOp::DataType;
  int buffer_cols = src->GetRasterXSize()>dst->GetRasterXSize()?src->GetRasterXSize():dst->GetRasterXSize();
  int buffer_rows = src->GetRasterCount()>dst->GetRasterCount()?src->GetRasterCount():dst->GetRasterCount();

  std::vector<DataType> data((size_t)buffer_cols*buffer_rows);

  int pixelspace = sizeof(DataType);
  int bandspace = sizeof(DataType)*buffer_cols;
  int linespace = 0;

  for(int r=0; r<src->GetRasterYSize(); ++r){
    if(src->RasterIO(GF_Read, 0, r, src->GetRasterXSize(), 1, &data[0], src->GetRasterXSize(), 1, gdal::DataType<DataType>::type(), buffer_rows, nullptr, pixelspace, linespace, bandspace)){}
    if( op(r, &data[0], buffer_cols, buffer_rows) && dst->RasterIO(GF_Write, 0, r, dst->GetRasterXSize(), 1, &data[0], dst->GetRasterXSize(), 1, gdal::DataType<DataType>::type(), buffer_rows, nullptr, pixelspace, linespace, bandspace)){}
  }
}

template<class BandOp, bool ColMajor = false>
void BandMajor(GDALDataset* src, GDALDataset* dst, BandOp& op){
  using DataType = typename BandOp::DataType;
  int buffer_cols = src->GetRasterXSize()>dst->GetRasterXSize()?src->GetRasterXSize():dst->GetRasterXSize();
  int buffer_rows = src->GetRasterYSize()>dst->GetRasterYSize()?src->GetRasterYSize():dst->GetRasterYSize();

  std::vector<DataType> data((size_t)buffer_cols*buffer_rows);

  int pixelspace = sizeof(DataType);
  int linespace = sizeof(DataType)*buffer_cols;

  if(ColMajor){
    std::swap(buffer_cols, buffer_rows);
    pixelspace = sizeof(DataType)*buffer_cols;
    linespace = sizeof(DataType);
  }

  for(int b=1; b<=src->GetRasterCount(); ++b){
    if(src->GetRasterBand(b)->RasterIO(GF_Read, 0, 0, src->GetRasterXSize(), src->GetRasterYSize(), &data[0], src->GetRasterXSize(), src->GetRasterYSize(), gdal::DataType<DataType>::type(), pixelspace, linespace)){}
    if( op(b, &data[0], buffer_cols, buffer_rows) && dst->GetRasterBand(b)->RasterIO(GF_Write, 0, 0, dst->GetRasterXSize(), dst->GetRasterYSize(), &data[0], dst->GetRasterXSize(), dst->GetRasterYSize(), gdal::DataType<DataType>::type(), pixelspace, linespace)){}
  }
}

template<class BandOp, bool ColMajor = false>
void BandMajor(GDALDataset* src, BandOp& op){
  using DataType = typename BandOp::DataType;
  int buffer_cols = src->GetRasterXSize();
  int buffer_rows = src->GetRasterYSize();

  std::vector<DataType> data((size_t)buffer_cols*buffer_rows);

  int pixelspace = sizeof(DataType);
  int linespace = sizeof(DataType)*buffer_cols;

  if(ColMajor){
    std::swap(buffer_cols, buffer_rows);
    pixelspace = sizeof(DataType)*buffer_cols;
    linespace = sizeof(DataType);
  }

  for(int b=1; b<=src->GetRasterCount(); ++b){
    if(src->GetRasterBand(b)->RasterIO(GF_Read, 0, 0, src->GetRasterXSize(), src->GetRasterYSize(), &data[0], src->GetRasterXSize(), src->GetRasterYSize(), gdal::DataType<DataType>::type(), pixelspace, linespace)){}
    if( op(b, &data[0], buffer_cols, buffer_rows) ){}
  }
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
