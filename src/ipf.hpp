#ifndef IMAGE_PROCESS_FRAMEWORK_HPP
#define IMAGE_PROCESS_FRAMEWORK_HPP

#include <gdal_priv.h>
#include "gdal_traits.hpp"

#include <opencv2/core/core.hpp>
#include <numeric>

namespace ipf{

    class BufferOperator {
    public:
        BufferOperator() {}
        virtual ~BufferOperator() {}
        virtual bool operator()(int r, void* data, int cols, int rows) = 0;
    };

    class Framework {
    public:
        Framework() : _datatype(GDT_Byte), _datatype_size(sizeof(unsigned char)), _src(nullptr), _max_buffer_size(1024*1024*1024) {
        }
        virtual ~Framework() {}
        void SetDataType(GDALDataType datatype) { _datatype = datatype; }
        void SetDataTypeSize(int size) { _datatype_size = size; }
        int GetDataTypeSize() const { return _datatype_size; }
        GDALDataType GetDataType() const { return _datatype; }
        size_t GetMaxBufferSize() const { return _max_buffer_size; }

        virtual bool Apply(BufferOperator* op) = 0;
        bool Apply(int srcwin[6], int bufferwin[3], BufferOperator* op) {
            l;
        }
    protected:
        int _datatype_size;
        GDALDataType _datatype;
        size_t _max_buffer_size;
        GDALDataset* _src;
    };

    class RowMajor : public Framework {
    public:
        RowMajor(GDALDataset* src, GDALDataset* dst) : _src(src), _dst(dst) {}
        virtual ~RowMajor() {}
        bool Apply(BufferOperator* op) override {
            int buffer_cols = HasSrcWinX() ?  GetSrcWin()[1]: _src->GetRasterXSize();
            int buffer_rows = _src->GetRasterCount();
            std::vector<int> bandlist;
            if (HasSrcWinZ()) {
                buffer_rows = GetSrcWin()[5];
                bandlist.resize(buffer_rows);
                std::iota(bandlist.begin(), bandlist.end(), GetSrcWin()[4]);
            }

            size_t buffer_size = (size_t)buffer_cols * buffer_rows;
            int buffer_count = GetMaxBufferSize() / buffer_size;

            std::vector<char> data();

            int pixelspace = GetDataTypeSize();
            int bandspace = GetDataTypeSize() * buffer_cols;
            int linespace = bandspace * buffer_rows;

            for (int r = 0; r < _src->GetRasterYSize(); ++r) {
                if (_src->RasterIO(GF_Read, 0, r, _src->GetRasterXSize(), 1, &data[0], _src->GetRasterXSize(), 1, _datatype, buffer_rows, nullptr, pixelspace, linespace, bandspace)) {}
                if ((*op)(r, &data[0], buffer_cols, buffer_rows) && _dst->RasterIO(GF_Write, 0, r, _dst->GetRasterXSize(), 1, &data[0], _dst->GetRasterXSize(), 1, _datatype, buffer_rows, nullptr, pixelspace, linespace, bandspace)) {}
            }
        }
    protected:
        GDALDataset* _src;
        GDALDataset* _dst;
    };



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
