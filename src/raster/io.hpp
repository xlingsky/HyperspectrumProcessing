#ifndef XLINGSKY_RASTER_IO_HPP
#define XLINGSKY_RASTER_IO_HPP

#include "operator.h"
#include "util/gdalex.hpp"
#include <gdal.h>
#include <gdal_priv.h>

namespace xlingsky {
namespace raster {

class RasterIO : public Operator {
protected:
  GDALDataset* _dataset;
  GDALDataType _datatype;
  bool _owned;
 public:
   RasterIO(GDALDataset *dataset, GDALDataType datatype, bool owned = true)
       : _dataset(dataset), _datatype(datatype), _owned(owned) {}
   virtual ~RasterIO() {
    close();
   }
   int width() const { return _dataset->GetRasterXSize(); }
   int height() const { return _dataset->GetRasterYSize(); }
   int channel() const { return _dataset->GetRasterCount(); }
   bool good() const { return _dataset != nullptr; }
   GDALDataset* ptr() const { return _dataset; }
   virtual void close() {
    if (_owned && _dataset) {
      GDALClose(_dataset);
      _dataset = nullptr;
    }
   }
};

class ImageReader : public RasterIO {
 public:
  ImageReader(GDALDataset* dataset, GDALDataType datatype, bool owned = true) : RasterIO(dataset, datatype, owned) {}
  ImageReader(const char* filepath, GDALDataType type = GDT_Unknown) : RasterIO(nullptr, type, false) {
    open(filepath, type);
  }
  virtual ~ImageReader() {}
  bool open(const char* filepath, GDALDataType type = GDT_Unknown) {
    GDALDataset* dataset = (GDALDataset*)GDALOpen(filepath, GA_ReadOnly);
    if (!dataset) return false;
    _dataset = dataset;
    _datatype = (type==GDT_Unknown? _dataset->GetRasterBand(1)->GetRasterDataType(): type);
    _owned = true;
    return true;
  }
  bool operator()(void *data, int imoff[3], int size[3], int space[3],
                          int prior[3]) override {
    space[prior[0]] = GDALGetDataTypeSizeBytes(_datatype); 
    space[prior[1]] = space[prior[0]] * size[prior[0]];
    space[prior[2]] = space[prior[1]] * size[prior[1]];
    std::vector<int> bandlist(size[2]);
    std::iota(bandlist.begin(), bandlist.end(),
              imoff[2]+1);
    return _dataset->RasterIO(
               GF_Read, imoff[0], imoff[1], size[0], size[1], data, size[0],
               size[1], _datatype, size[2],
               &bandlist[0], space[0], space[1], space[2]) == CE_None;
  }
};

class ImageWriter : public RasterIO {
 public:
  ImageWriter(GDALDataset* dataset, GDALDataType datatype) : RasterIO(dataset, datatype) {}
  ImageWriter(const char* filepath, int width, int height, int channel, GDALDataType type) : RasterIO(nullptr, type, false) {
    open(filepath, width, height, channel, type);
  }
  virtual ~ImageWriter() {}
  bool open(const char* filepath, int width, int height, int channel, GDALDataType type) {
    GDALDataset* dataset = GDALCreate(filepath, width, height, channel, type);
    if (!dataset) return false;
    _dataset = dataset;
    _datatype = type;
    _owned = true;
    return true;
  }
  bool operator()(void *data, int imoff[3], int size[3], int space[3],
                          int prior[3]) override{
    std::vector<int> bandlist(size[2]);
    std::iota(bandlist.begin(), bandlist.end(),
              imoff[2]+1);
    return _dataset->RasterIO(
               GF_Write, imoff[0], imoff[1], size[0], size[1], data, size[0],
               size[1], _datatype, size[2],
               &bandlist[0], space[0], space[1], space[2]) == CE_None;
  }
};

};
};

#endif