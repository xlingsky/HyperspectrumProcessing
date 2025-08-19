#ifndef XLINGSKY_RASTER_RADIOMETRIC_HPP
#define XLINGSKY_RASTER_RADIOMETRIC_HPP

#include <assert.h>
#include <vector>
#include <sstream>

#include "util/TileManager.hpp"
#include "util/gdalex.hpp"
#include "util/gdal_traits.hpp"

#include "raster/operator.h"
#include "raster/detail/inpaint.hpp"
#include "raster/common.hpp"

//#define DEBUG

#if defined(DEBUG) || defined(_DEBUG)
#include <boost/dll/runtime_symbol_info.hpp>
#include "raster/gdalex.hpp"
#include "raster/gdal_traits.hpp"
#endif

#ifdef _LOGGING
#include <glog/logging.h>
#define _LOG_LEVEL_RADIOMETRIC 1
#endif
// #include <sstream>
//#define DEFECTIVEPIXEL_COUNTING

namespace xlingsky {

namespace raster {

namespace radiometric {

template <typename T>
bool save(const char* filepath, T* data, size_t count, int newline) {
  std::ostringstream out;
  size_t i = 0;
  while (i < count) {
    out << data[i] << "\t";
    ++i;
    if (i % newline == 0) {
      out << std::endl;
    }
  }
  FILE* fp = fopen(filepath, "w");
  if (fp == nullptr) {
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "CANNOT create " << filepath;
#endif
    return false;
  }
  fputs(out.str().c_str(), fp);
  fclose(fp);
  return true;
}

class PixelCorrection : public FrameIterator {
 public:
  typedef float DataType;
 protected:
  int _cols;
  int _rows;
  int _threshold;
  DataType _dst_minimum;
  DataType _dst_maximum;
#ifdef DEFECTIVEPIXEL_COUNTING
  std::vector<std::vector<int> > _list;
#endif

 public:
  PixelCorrection(int cols, int rows, DataType dst_min = 0, DataType dst_max = 255)
      : _cols(cols), _rows(rows), _dst_minimum(dst_min), _dst_maximum(dst_max),
#ifdef DEFECTIVEPIXEL_COUNTING
    _list(rows),
#endif
    _threshold(10) {}
  void set_destination_range(DataType min, DataType max) {
    _dst_minimum = min;
    _dst_maximum = max;
  }
  virtual DataType correct(int, DataType, int) = 0;

  bool operator()(int b, int xoff, int yoff, void* data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
// #ifdef _USE_OPENMP
// #if _USE_OPENMP > 4
// #pragma omp declare reduction (merge : std::vector<int> :
// omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) #pragma omp
// parallel for reduction(merge: _list) #else #pragma omp parallel for
// reduction(+: cnt) #endif #endif
#ifndef DEFECTIVEPIXEL_COUNTING
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
#endif
    for (int r=0; r<rows; ++r){
      size_t i0 = (size_t)(r+yoff)*_cols+xoff;
      size_t i1 = (size_t)r*cols;
      for(int c=0; c<cols; ++c){
        size_t i = i1+c;
        if(IsNoData(pdata[i])) continue;
        DataType d = correct(b, pdata[i], i0+c);
        if (d < _dst_minimum) {
          pdata[i] = _dst_minimum;
#ifdef DEFECTIVEPIXEL_COUNTING
          if (d < -_threshold) _list[b].push_back(i);
#endif
        } else if(d>_dst_maximum){
          d = _dst_maximum;
        }else pdata[i] = d;
      }
    }
    return true;
  }
  int cols() const { return _cols; }
  int rows() const { return _rows; }
    
#ifdef DEFECTIVEPIXEL_COUNTING
  int max_bad_pixel_num() const {
    int cnt = 0;
    for (auto it = _list.begin(); it != _list.end(); ++it)
      if (cnt < it->size()) cnt = it->size();
    return cnt;
  }
#endif
};

class DefectivePixelCorrection : public Operator {
 public:
  typedef float DataType;
  typedef void* Mask;
  DefectivePixelCorrection() {}
  ~DefectivePixelCorrection(){}

  virtual Mask init(int imoff[3], int size[3], int prior[3]) = 0;
  virtual cv::Mat begin(int b, Mask mask) = 0;
  virtual void end(int b, Mask mask) {}
  virtual void exit(Mask mask) = 0;

  template<class Config>
  static DefectivePixelCorrection* Create(const Config& config, int src_win[3], int store_prior[3]);

  bool operator()(void* data, int imoff[3], int size[3], int space[3],
                  int prior[3]) override {

    Mask mask = init(imoff, size, prior);
    if(mask == nullptr) return false;

    int rows = size[prior[1]];
    int cols = size[prior[0]];

    xlingsky::TileManager manager;
    const int margin = 10;
    int tilesize = rows;
    if(cols<tilesize) tilesize = cols;
    if(500<tilesize) tilesize = 500;
    manager.AppendDimension(cols, tilesize, margin, margin);
    manager.AppendDimension(rows, tilesize, margin, margin);

    for (int b = 0; b < size[prior[2]]; ++b){
      cv::Mat m = begin(b, mask);
      cv::Mat d( rows, cols, cv::DataType<DataType>::type, (char*)data+(size_t)space[prior[2]]*b);
      detail::InpaintOp op(d, d, m);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
      for (int r = 0; r < manager.Size(1); ++r) {
        auto segr = manager.Segment(1, r);
        int h = (r+1==manager.Size(1)?segr.second:tilesize);
        for (int c = 0; c < manager.Size(0); ++c) {
          auto segc = manager.Segment(0, c);
          int w = (c+1==manager.Size(0)?segc.second:tilesize);
          cv::Rect rc0(segc.first, segr.first, w, h),
              rc1(segc.first, segr.first, segc.second, segr.second), rc2(0,0,w,h);
          op(rc1, rc2, rc0);
        }
      }
      end(b, mask);
    }
    exit(mask);

    return true;
  }
};

class DefectivePixelCorrectionV1 : public DefectivePixelCorrection{
protected:
  typedef float DataType;
  struct MM{
    cv::Mat m0;
    int w;
    int h;
    bool flag;
  };

  std::vector<unsigned char> _data;
  int _width;
  int _height;
  int _band_dim;
 public:
  DefectivePixelCorrectionV1() {}
  ~DefectivePixelCorrectionV1(){}
  bool load(const char* filepath, int cols, int rows, int band_dim){
    _data.resize((size_t)cols*rows);
    if (xlingsky::raster::common::load(filepath, cols, rows, &_data[0])) {
      _width = cols;
      _height = rows;
      _band_dim = band_dim;
      return true;
    }
    return false;
  }

  Mask init(int imoff[3], int size[3], int prior[3]) override {
    MM *mm = new MM;
    cv::Mat t = cv::Mat( _height, _width, cv::DataType<unsigned char>::type,
                    &_data[0]);
    mm->w = size[prior[0]];
    mm->h = size[prior[1]];
    if(prior[2] == _band_dim) {
      mm->flag = true;
      mm->m0 = t(cv::Rect(imoff[prior[0]], imoff[prior[1]], size[prior[0]], size[prior[1]]));
    }
    else{
      mm->flag = false;
      mm->m0 = t(cv::Rect(imoff[prior[0]], imoff[prior[2]], size[prior[0]], size[prior[2]]));
    }
    return mm;
  }
  cv::Mat begin(int b, Mask mask) override{
    MM* mm = (MM*)mask;
    if(mm->flag){
      return mm->m0;
    }else{
      cv::Mat m = cv::Mat(mm->h, mm->w, mm->m0.type(), mm->m0.ptr(b), 0);
      m.step[0] = 0;
      return m;
    }
  }
  void exit(Mask mask) override{
    MM* mm = (MM*)mask;
    delete mm;
  }
};

class DefectivePixelCorrectionV2 : public DefectivePixelCorrection {
 protected:
  struct MM{
    unsigned char* data;
    int space;
    int cols;
    int rows;
  };
  GDALDataset* _mask;

 public:
  typedef float DataType;
  DefectivePixelCorrectionV2() : _mask(nullptr){
  }
  ~DefectivePixelCorrectionV2(){
    if(_mask) GDALClose(_mask);
  }
  bool load(const char* maskpath, int cols, int rows, int bands){
    if(!IsRasterDataset(maskpath)) return false;
    GDALDataset* d = (GDALDataset*)GDALOpen(maskpath, GA_ReadOnly);
    if(d==nullptr){
      return false;
    }
    if(d->GetRasterXSize()<cols || d->GetRasterYSize()<rows || d->GetRasterCount()<bands){
      return false;
    }
    if(_mask) GDALClose(_mask);
    _mask = d;
    return true;
  }
  Mask init(int imoff[3], int size[3], int prior[3]) override{
    MM* mm = new MM;
    mm->data = new unsigned char[(size_t)size[0]*size[1]*size[2]];
    int mmspace[3];
    mmspace[prior[0]] = sizeof(unsigned char);
    mmspace[prior[1]] = mmspace[prior[0]]*size[prior[0]];
    mmspace[prior[2]] = mmspace[prior[1]]*size[prior[1]];
    std::vector<int> bandlist;
    bandlist.resize(size[2]);
    std::iota(bandlist.begin(), bandlist.end(), imoff[2]+1);
    if(_mask->RasterIO(GF_Read, imoff[0], imoff[1], size[0], size[1], &mm->data[0], size[0], size[1], gdal::DataType<unsigned char>::type(), size[2], &bandlist[0], mmspace[0], mmspace[1], mmspace[2])!=CE_None){
      delete[] mm->data;
      delete mm;
      return nullptr;
    }
    mm->space = mmspace[prior[2]];
    mm->cols = size[prior[0]];
    mm->rows = size[prior[1]];
    return mm;
  }
  cv::Mat begin(int b, Mask mask) override{
    MM* mm = (MM*)mask;
    return cv::Mat( mm->rows, mm->cols, cv::DataType<unsigned char>::type,
                    mm->data+b*mm->space);
  }
  void exit(Mask mask) override{
    MM* mm = (MM*)mask;
    if(mm->data) delete[] mm->data;
    delete mm;
  }
};

template<class Config>
DefectivePixelCorrection* DefectivePixelCorrection::Create(const Config& config, int src_size[3], int store_prior[3]){
  xlingsky::raster::radiometric::DefectivePixelCorrectionV2 *p1 =
      new xlingsky::raster::radiometric::DefectivePixelCorrectionV2;
  if (p1->load(config.filepath("file").string().c_str(), src_size[0], src_size[1], src_size[2])) {
    return p1;
  } 
  delete p1;
  xlingsky::raster::radiometric::DefectivePixelCorrectionV1 *p2 =
      new xlingsky::raster::radiometric::DefectivePixelCorrectionV1;
  int pt[3];
  {
    std::string prior = config.template get<std::string>("dim_prior", "");
    std::vector<std::string> result;
    boost::split(result, prior, boost::is_any_of(","));
    if (result.size() == 3) {
      pt[0] = std::stoi(result[0]);
      pt[1] = std::stoi(result[1]);
      pt[2] = std::stoi(result[2]);
    } else {
      memcpy(pt, store_prior, sizeof(int) * 3);
    }
  }
  if (!p2->load( config.filepath("file").string().c_str(), src_size[pt[0]],
                src_size[pt[1]], pt[2])) {
    delete p2;
    return nullptr;
  }
  return p2;
}

class MeanStdCalculator : public FrameIterator {
 private:
  std::vector<double> _mean;
  std::vector<double> _std;
  int _width;
  char _filepath[512];
  char _xmlpath[512];

  int _dim_order[3];

  float _cut_ratio_upper;
  float _cut_ratio_lower;

 public:
  typedef float DataType;
  MeanStdCalculator(int w, float cut_ratio_lower = 0, float cut_ratio_upper = 0) : _width(w), _cut_ratio_lower(cut_ratio_lower), _cut_ratio_upper(cut_ratio_upper) {
    _xmlpath[0] = 0;
  }
  ~MeanStdCalculator() { 
      if (save(_filepath)) {
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "Means was saved to " << _filepath;
#endif
      }
      if(_xmlpath[0] && save2xml(_xmlpath)){
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "xml was saved to " << _xmlpath;
#endif
      }
  }
  std::pair<double, double> compute(DataType* data, int n, float cut_ratio_lower, float cut_ratio_upper) {
    std::pair<double, double> ret = std::make_pair(0, 0);

    DataType* end = RemoveNoData(data, data+n);
    n = end-data;
    int st, ed;
    st = int(n*cut_ratio_lower);
    ed = n - int(n*cut_ratio_upper);
    n = ed-st;

    if(n<=0) return ret;

    std::sort(data, end);

    int i = st;
    while(i<ed){
      ret.first += data[i];
      ret.second += (double)data[i] * data[i];
      ++i;
    }
    ret.first /= n;
    ret.second = std::sqrt(ret.second / n - ret.first * ret.first);
    return ret;
  }
  void SetFilePath(const char* filepath) {
    strcpy(_filepath, filepath);
  }
  void SetXmlPath(const char* xmlpath){
    if(xmlpath) strcpy(_xmlpath, xmlpath);
    else _xmlpath[0] = 0;
  }
  bool save(const char* filepath) {
    if (::xlingsky::raster::radiometric::save(filepath, _mean.data(), _mean.size(), _width)) {
      char path[512];
      strcpy(path, filepath);
      strcpy(strrchr(path, '.'), "_std.txt");
      if (::xlingsky::raster::radiometric::save(path, _std.data(), _std.size(), _width))
        return true;
    }

    return false;
  }
  bool operator()(int b, int , int , void* data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
    std::vector<double> mean(rows), std(rows);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for (int r = 0; r < rows; ++r) {
      auto ret = compute(pdata + r * cols, cols, _cut_ratio_lower, _cut_ratio_upper);
      mean[r] = ret.first;
      std[r] = ret.second;
    }
    _mean.insert(_mean.end(), mean.begin(), mean.end());
    _std.insert(_std.end(), std.begin(), std.end());
    return true;
  }

  void SetDimOrder(int order[3]) { memcpy(_dim_order, order, 3*sizeof(int)); }
  bool save2xml(const char* xmlpath){
    std::vector<double> a, b;
    for(int i = 0;i<_mean.size(); i+=_width){
      double mt = 0, st = 0;
      int cnt = 0;
      std::vector<bool> flag(_width);
      for(int j=0; j < _width; ++j){
        if(fabs(_std[i+j])<=std::numeric_limits<float>::epsilon()){
          flag[j] = false;
          continue;
        }
        mt += _mean[i+j];
        st += _std[i+j];
        ++cnt;
        flag[j] = true;
      }
      if(cnt>0){
        mt /= cnt;
        st /= cnt;
      }
      for(int j=0; j < _width; ++j){
        if(flag[j]){
          a.push_back(st/_std[i+j]);
          b.push_back(mt-_mean[i+j]*a.back());
        }else{
          a.push_back(1);
          b.push_back(0);
        }
      }
    }

    char path[512], patha[512], pathb[512];
    strcpy(path, xmlpath);
    char* ps = strrchr(path, '.');
    strcpy(ps, "_a.txt");
    if (!::xlingsky::raster::radiometric::save(path, a.data(), a.size(), _width))
      return false;
    strcpy(patha, path);
    strcpy(ps, "_b.txt");
    if (!::xlingsky::raster::radiometric::save(path, b.data(), b.size(), _width))
      return false;
    strcpy(pathb, path);

    FILE* fp = fopen(xmlpath, "w");
    if(fp){
      fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      fprintf(fp, "<HSP>\n");
      fprintf(fp, "\t<dim_prior>%d,%d,%d</dim_prior>\n", _dim_order[1], _dim_order[2], _dim_order[0]);
      fprintf(fp, "\t<task name=\"uniform\">\n");
      fprintf(fp, "\t\t<a>%s</a>\n", patha);
      fprintf(fp, "\t\t<b>%s</b>\n", pathb);
      fprintf(fp, "\t</task>\n");
      // if (_bp_path[0] && _hi_path[0]) {
      //   fprintf(fp, "\t<task name=\"dpc\">\n");
      //   fprintf(fp, "\t\t<file>%s</file>\n", _bp_path);
      //   fprintf(fp, "\t</task>\n");
      // }
      fprintf(fp, "</HSP>\n");
      fclose(fp);
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "NUC xml was saved to " << xmlpath;
#endif
      return true;
    }
    return false;
  }
};

class MedianCalculator : public FrameIterator {
 public:
  typedef float DataType;

 private:
  std::vector<DataType> _median;
  int _width;
  char _filepath[512];

 public:
  MedianCalculator(int w) : _width(w) {}
  ~MedianCalculator() { 
      if (save(_filepath)) {
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "Medians was saved to " << _filepath;
#endif
    }
  }
  void SetFilePath(const char* filepath) { strcpy(_filepath, filepath); }
  DataType compute(DataType* data, int n) {
    DataType* end = RemoveNoData(data, data+n);
    std::sort(data, end);
    return data[(end-data) >> 1];
  }
  bool save(const char* filepath) {
    return ::xlingsky::raster::radiometric::save(filepath, _median.data(), _median.size(),
                               _width);
  }
  bool operator()(int , int , int, void* data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
    std::vector<double> median(rows);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for (int r = 0; r < rows; ++r) {
      median[r] = compute(pdata + r * cols, cols);
    }
    _median.insert(_median.end(), median.begin(), median.end());
    return true;
  }
};

};  // namespace radiometric
};  // namespace raster

};  // namespace xlingsky

#endif
