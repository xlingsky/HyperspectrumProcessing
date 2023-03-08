#ifndef XLINGSKY_RASTER_COMMON_HPP
#define XLINGSKY_RASTER_COMMON_HPP

#include "raster/operator.h"
#include "raster/gdalex.hpp"
#include "raster/gdal_traits.hpp"

#include <fstream>
#include <boost/algorithm/string.hpp>

#ifdef _LOGGING
#include <glog/logging.h>
#define _LOG_LEVEL_COMMON 1
#endif

namespace xlingsky{

namespace raster{

namespace common{

class Sort : public FrameIterator {
 public:
  typedef float DataType;
  Sort(){}
  bool operator()(int , int , int, void* data, int cols, int rows) override{
    DataType* pdata = (DataType*)data;
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for(int r=0; r<rows; ++r){
      auto t = pdata+r*cols;
      std::sort( t, t + cols);
    }
    return true;
  }
};

template<typename T>
bool load(const char* filepath, int cols, int rows, T* data){
    if (IsRasterDataset(filepath)) {
      GDALDataset* dataset = (GDALDataset*)GDALOpen(filepath, GA_ReadOnly);
      if(dataset->GetRasterXSize() != cols || dataset->GetRasterYSize() != rows){
        GDALClose(dataset);
        return false;
      }
      if( dataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, cols, rows, data, cols, rows, gdal::DataType<T>::type(), 0, 0) ){
      }
      GDALClose(dataset);
      return true;
    }
    std::ifstream file(filepath);
    if(!file.is_open()) return false;
    std::string line;
    int r=0;
    auto token = boost::is_any_of(" ,;|\t");
    while (std::getline(file, line)) {
      std::vector<std::string> vs;
      boost::trim_if(line, token);
      boost::algorithm::split(vs, line, token, boost::token_compress_on);
      if(vs.size() == cols){
        T* pm = data+r*cols;
        for(int c=0; c<cols; ++c)
          pm[c] = (T)std::stof(vs[c]);
      }else if(vs.size()==2){
        int tc = std::stoi(vs[1]), tr = std::stoi(vs[0]);
        if ( tr > 0 && tr <= rows && tc > 0 && tc <= cols)
          data[(tr - 1)*cols+(tc-1)] = 1;
        else return false;
      }else return false;
      ++r;
    }
    return true;
}

class Extend : public FrameIterator{
 protected:
  typedef float DataType;
  DataType* _data;
  bool _own_data;
  void reset(){
    if(_own_data && _data)
    {
      delete[] _data;
      _data = nullptr;
      _own_data = false;
    }
  }
 public:
  Extend( DataType* data = nullptr) : _data(data), _own_data(false) {
  }
  ~Extend(){
    reset();
  }
  bool load(const char* filepath, int cols, int rows) {
    DataType* d = new DataType[(size_t)cols*rows];
    if (!xlingsky::raster::common::load(filepath, cols, rows, d)){
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_COMMON) << "[Extend] size not match " << filepath;
#endif
      delete[] d;
      return false;
    }
    reset();
    _data = d;
    _own_data = true;
    return true;
  }
  bool operator()(int , int , int, void* data, int cols, int rows) override{
    DataType* pdata = (DataType*)data;
    memcpy(pdata, _data, sizeof(DataType)*cols*rows);
    return true;
  }
};

};
};
};

#endif
