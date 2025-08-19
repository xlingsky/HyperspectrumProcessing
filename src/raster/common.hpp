#ifndef XLINGSKY_RASTER_COMMON_HPP
#define XLINGSKY_RASTER_COMMON_HPP

#include "raster/operator.h"
#include "util/gdal_traits.hpp"
#include "util/gdalex.hpp"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <gdal.h>

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

class Minus : public Operator {
  private:
  GDALDataset* _dataset;
  int _interval;
  bool _owned;
  public:
  typedef float DataType;
  Minus() : _dataset(nullptr), _interval(30), _owned(false) {}
  Minus(GDALDataset* dataset, int interval) : _dataset(dataset), _interval(interval), _owned(false) {}
  ~Minus(){
    if ( _owned && _dataset) {
      GDALClose(_dataset);
      _dataset = nullptr;
    }
  }
  template<class Config>
  bool load_config(const Config& config) {
    _interval = config.template get<int>("interval", 30);
    std::string filepath = config.template get<std::string>("background", "");
    if (!filepath.empty()) {
      set_dataset(filepath.c_str());
      return good();
    }
    return true;
  }
  bool set_dataset(const char* filepath) {
    _dataset = (GDALDataset *)GDALOpen(filepath, GA_ReadOnly);
    _owned = true;
    return good();
  }
  bool good() const { return _dataset != nullptr; }
  bool operator()(void* data, int imoff[3], int size[3], int space[3], int prior[3]) override{
    int blknum = (int)std::ceil(size[prior[2]] / (float)_interval);
    const unsigned bytes = (unsigned)sizeof(DataType);
    unsigned stride[3] = {space[0] / bytes, space[1] / bytes, space[2] / bytes};

    int imoff_new[3];
    imoff_new[prior[0]] = imoff[prior[0]];
    imoff_new[prior[1]] = imoff[prior[1]];
    imoff_new[prior[2]] = imoff[prior[2]] / _interval;
    int size_new[3];
    size_new[prior[0]] = size[prior[0]];
    size_new[prior[1]] = size[prior[1]];
    size_new[prior[2]] = blknum;
    std::vector<int> bandlist(size[2]);
    std::iota(bandlist.begin(), bandlist.end(),
              imoff[2]+1);
    std::vector<DataType> bgdata(size[prior[0]] * size[prior[1]] * blknum);
    if (_dataset->RasterIO(GF_Read, imoff_new[0], imoff_new[1], size_new[0],
                           size_new[1], bgdata.data(), size_new[0], size_new[1],
                           gdal::DataType<DataType>::type(), size_new[2], &bandlist[0],
                           space[0], space[1], space[2]) != CE_None) {
      return false;
    }
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < blknum; ++i) {
      int st = i * _interval;
      int sz = i == blknum - 1 ? size[prior[2]] - st : _interval;
      DataType *pdata = (DataType *)data + st * stride[prior[2]];
      DataType *bgptr =
          bgdata.data() + (size_t)i * size[prior[0]] * size[prior[1]];
      for (int j = 0; j < size[prior[0]] * size[prior[1]]; ++j) {
        for (int k = 0; k < sz; ++k)
          pdata[k * stride[prior[2]]] -= *bgptr;
        ++pdata;
        ++bgptr;
      }
    }
    return true;
  }
};

};
};
};

#endif
