#ifndef XLINGSKY_RASTER_RADIOMETRIC_HPP
#define XLINGSKY_RASTER_RADIOMETRIC_HPP

#include <assert.h>
#include <vector>
#include <sstream>

#include "util/TileManager.hpp"
#include "util/InterpolatorAdaptor.hpp"
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

template <class TileOp>
void Tile(int width, int height, int tilesize, int margin, TileOp& op) {
  int tile_x_num = width / tilesize;
  int tile_y_num = height / tilesize;
  if (tile_x_num < 1) tile_x_num = 1;
  if (tile_y_num < 1) tile_y_num = 1;

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
  for (int ty = 0; ty < tile_y_num; ++ty) {
    int sty = ty * tilesize;
    int h = tilesize;
    if (sty + 2 * tilesize > height) h = height - sty;
    cv::Rect rc1, rc2;
    if (sty > margin) {
      rc1.y = sty - margin;
      rc2.y = margin;
    } else {
      rc1.y = 0;
      rc2.y = 0;
    }
    if (sty + h + margin > height) {
      rc1.height = height - rc1.y;
    } else
      rc1.height = sty + h + margin - rc1.y;
    rc2.height = h;
    for (int tx = 0; tx < tile_x_num; ++tx) {
      int stx = tx * tilesize;
      int w = tilesize;
      if (stx + 2 * tilesize > width) w = width - stx;
      cv::Rect rc0(stx, sty, w, h);
      if (stx > margin) {
        rc1.x = stx - margin;
        rc2.x = margin;
      } else {
        rc1.x = 0;
        rc2.x = 0;
      }
      rc2.width = w;
      if (stx + w + margin > width)
        rc1.width = width - rc1.x;
      else
        rc1.width = stx + w + margin - rc1.x;
      op(rc1, rc2, rc0);
    }
  }
}

template <typename T>
bool load(const char* filepath, size_t count, T* data) {
  if(IsRasterDataset(filepath) ){
    GDALDataset* dataset = (GDALDataset*)GDALOpen(filepath, GA_ReadOnly);
    if((size_t)dataset->GetRasterXSize()*dataset->GetRasterYSize() < count){
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "Factors NOT enough " << filepath;
#endif
      GDALClose(dataset);
    return false;
    }
    int cols = dataset->GetRasterXSize();
    int rows = count/cols;
    if( dataset->RasterIO(GF_Read, 0, 0, cols, rows, data, cols, rows, gdal::DataType<T>::type(), 1, nullptr, 0, 0, 0) ){
    }
    GDALClose(dataset);
    return true;
  }
  FILE* fp = fopen(filepath, "r");
  if (fp == nullptr) {
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "CANNOT open " << filepath;
#endif
    return false;
  }
  size_t i;
  for (i = 0; i < count; ++i) {
    double t;
    if (fscanf(fp, "%lf%*c", &t) != 1) break;
    data[i] = (T)t;
  }

  fclose(fp);

#ifdef _LOGGING
  if (i!=count) {
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "ERROR format " << filepath;
  }
#endif

  return i == count;
}

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

// dv is a bandsXcols  dark level composation value matrix stored in row-major
class DarkBackgroundCorrection : public PixelCorrection {
 public:
  typedef PixelCorrection::DataType DataType;

 protected:
  DataType* _data;

 public:
  DarkBackgroundCorrection(int cols, int rows, DataType dst_min, DataType dst_max)
      : _data(new DataType[(size_t)cols * rows]),
        PixelCorrection(cols, rows, dst_min, dst_max) {}
  ~DarkBackgroundCorrection() {
    if (_data) delete[] _data;
  }
  bool load(const char* filepath) {
    return ::xlingsky::raster::radiometric::load(filepath, (size_t)cols() * rows(), _data);
  }
  DataType correct(int, DataType d, int i) override { return d - _data[i]; }
};

class DarkBackgroundLinear : public PixelCorrection{
public:
  typedef PixelCorrection::DataType DataType;
protected:
    DataType* _a;
    DataType* _b;
    size_t* _index;
    int _bands;
public:
  DarkBackgroundLinear(int cols, int rows, int bands, DataType dst_min, DataType dst_max) : PixelCorrection(cols, rows, dst_min, dst_max), _bands(bands){
        size_t sz = (size_t)cols * rows;
        _a = new DataType[sz];
        _b = new DataType[sz];
        _index = new size_t[bands];
    }
    ~DarkBackgroundLinear(){
        if(_a) delete[] _a;
        if(_b) delete[] _b;
        if(_index) delete[] _index;
    }
    bool load(const char* a, const char* b, const char* index){
        if (::xlingsky::raster::radiometric::load(a, (size_t)cols() * rows(), _a) &&
            ::xlingsky::raster::radiometric::load(b, (size_t)cols() * rows(), _b)){
            FILE* fp = fopen(index, "r");
            if(fp){
                size_t cnt = 0;
                char strline[1000000];
                while(fgets(strline, 1000000, fp) && cnt<_bands){
                    size_t d;
                    if(sscanf(strline, "%ld",&d)==1){
                        _index[cnt++] = d;
                    }
                }
                fclose(fp);
                return cnt==_bands;
            }
        }
        return false;
    }
  DataType correct(int b, DataType d, int i) override {
      return d-(_a[i] * _index[b] + _b[i]);
  }
};

class NonUniformCorrection : public PixelCorrection {
 public:
  typedef PixelCorrection::DataType DataType;

 protected:
  DataType* _a;
  DataType* _b;

 public:
  NonUniformCorrection(int cols, int rows) : PixelCorrection(cols, rows) {
    size_t sz = (size_t)cols * rows;
    _a = new DataType[sz];
    _b = new DataType[sz];
  }
  ~NonUniformCorrection() {
    if (_a) delete[] _a;
    if (_b) delete[] _b;
  }
  template <class Config>
  bool load_config(const Config& config) {
    std::string a, b;
    try {
      a = config.filepath("a").string();
      b = config.filepath("b").string();
    } catch (const boost::property_tree::ptree_error &e) {
      throw std::runtime_error("NonUniformCorrection: " +
                               std::string(e.what()));
      return false;
    }
    if (a.empty() || b.empty()){
      throw std::runtime_error("NonUniformCorrection: a or b is empty");
      return false;
    }
    float dst_min =
        config.template get<float>("dst_min", (std::numeric_limits<float>::min)());
    float dst_max =
        config.template get<float>("dst_max", (std::numeric_limits<float>::max)());
    set_destination_range(dst_min, dst_max);
    return load(a.c_str(), b.c_str());
  }
  bool load(const char* a, const char* b) {
    if (::xlingsky::raster::radiometric::load(a, (size_t)cols() * rows(), _a) &&
        ::xlingsky::raster::radiometric::load(b, (size_t)cols() * rows(), _b))
      return true;
    return false;
  }
  DataType correct(int, DataType d, int i) override { return _a[i] * d + _b[i]; }
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
  // bool load(const char* filepath, int cols, int rows, int bands) {
  //   int c[] = {cols, cols, rows};
  //   int r[] = {rows, bands, bands};
  //   int b[] = {2, 1, 0};
  //   for(int i=0; i<3; ++i){
  //     _data.resize((size_t)c[i]*r[i]);
  //     if(xlingsky::raster::common::load(filepath, c[i], r[i], &_data[0])){
  //       _width = c[i];
  //       _height = r[i];
  //       _band_dim = b[i];
  //       return true;
  //     }else if(xlingsky::raster::common::load(filepath, r[i], c[i], &_data[0])){
  //       _width = r[i];
  //       _height = c[i];
  //       _band_dim = b[i];
  //       return true;
  //     }
  //   }
  //   return false;
  // }
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

class NucCalculator : public FrameIterator{
 public:
  typedef float DataType;
  enum Mode{
    SCALE,
    OFFSET
  };
  enum InterpType{
    BARYCENTRIC,
    PCHIP,
    MAKIMA
  };
 private:
  float _cut_ratio_bright;
  float _cut_ratio_dark;
  int _sample_maxnum_bright;
  int _sample_maxnum_dark;
  int _sample_minnum_bright;
  int _sample_minnum_dark;
  int _line_tile_size;
  int _line_tile_overlap;
  float _value_ratio_threshold_bright;
  float _value_ratio_threshold_dark;
  int _mode;
  InterpType _interp_type;

  std::vector<double> _a;
  std::vector<double> _b;
  std::vector<double> _hi;
  std::vector<double> _lo;
  std::vector<int> _defective_pixels;
  int _width;

  char _a_path[512];
  char _b_path[512];
  char _bp_path[512];
  char _xml_path[512];
  char _hi_path[512];
  char _lo_path[512];

  int _dim_order[3];
 public:
  NucCalculator(int w, float cut_ratio_dark = 0.03, float cut_ratio_bright = 0.1, float ratio_threshold_dark = 0.03, float ratio_threshold_bright = 0.03, int sample_num_dark = 40, int sample_num_bright = 40, int line_tile_size = 30, int line_tile_overlap = 15, InterpType type = PCHIP, int mode = SCALE) : _width(w), _cut_ratio_dark(cut_ratio_dark), _cut_ratio_bright(cut_ratio_bright), _sample_maxnum_dark(sample_num_dark), _sample_maxnum_bright(sample_num_bright), _sample_minnum_dark(3), _sample_minnum_bright(3), _value_ratio_threshold_dark(ratio_threshold_dark), _value_ratio_threshold_bright(ratio_threshold_bright), _mode(mode), _line_tile_size(line_tile_size), _line_tile_overlap(line_tile_overlap), _interp_type(type) {
    _a_path[0] = _b_path[0] = _bp_path[0] = _xml_path[0] = 0;
    _hi_path[0] = _lo_path[0] = 0;
    _dim_order[0] = 0;
    _dim_order[1] = 1;
    _dim_order[2] = 2;
  }
  ~NucCalculator() {
    if(_a_path[0]) {
      if(::xlingsky::raster::radiometric::save(_a_path, _a.data(), _a.size(), _width)){
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "Linear factor was saved to " << _a_path;
#endif
      }
    }
    if(_b_path[0] && ::xlingsky::raster::radiometric::save(_b_path, _b.data(), _b.size(), _width)){
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "Offset factor was saved to " << _b_path;
#endif
    }
    if(_bp_path[0] && ::xlingsky::raster::radiometric::save(_bp_path, _defective_pixels.data(), _defective_pixels.size(), _width)) {
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "bad pixel list was saved to " << _bp_path;
#endif
    }
    if(_xml_path[0] && _a_path[0] && _b_path[0]){
      FILE* fp = fopen(_xml_path, "w");
      if(fp){
        fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fp, "<HSP>\n");
        fprintf(fp, "\t<dim_prior>%d,%d,%d</dim_prior>\n", _dim_order[1], _dim_order[2], _dim_order[0]);
        fprintf(fp, "\t<task name=\"uniform\">\n");
        fprintf(fp, "\t\t<a>%s</a>\n", _a_path);
        fprintf(fp, "\t\t<b>%s</b>\n", _b_path);
        fprintf(fp, "\t</task>\n");
        if(_bp_path[0] && _hi_path[0]){
          fprintf(fp, "\t<task name=\"dpc\">\n");
          fprintf(fp, "\t\t<file>%s</file>\n", _bp_path);
          fprintf(fp, "\t</task>\n");
        }
        fprintf(fp, "</HSP>\n");
        fclose(fp);
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "NUC xml was saved to " << _xml_path;
#endif
      }
      if(_hi_path[0] && ::xlingsky::raster::radiometric::save(_hi_path, _hi.data(), _hi.size(), _width)){
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "bright dn was saved to " << _hi_path;
#endif
      }
      if(_lo_path[0] && ::xlingsky::raster::radiometric::save(_lo_path, _lo.data(), _lo.size(), _width)){
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_RADIOMETRIC) << "dark dn was saved to " << _lo_path;
#endif
      }
    }
  }
  void SetDimOrder(int order[3]) { memcpy(_dim_order, order, 3*sizeof(int)); }
  void SetFilePath(const char *apath, const char *bpath, const char *bppath, const char* xmlpath, const char* hi_path, const char* lo_path) {
    if(apath) strcpy(_a_path, apath);
    if(bpath) strcpy(_b_path, bpath);
    if(bppath) strcpy(_bp_path, bppath);
    if(xmlpath) strcpy(_xml_path, xmlpath);
    if(hi_path) strcpy(_hi_path, hi_path);
    if(lo_path) strcpy(_lo_path, lo_path);
  }
  template<class Config>
  bool load_config(const Config& config){
      _cut_ratio_dark = config. template get<float>("cut_dark", 0.03);
      _cut_ratio_bright = config. template get<float>("cut_bright", 0.1);
      _value_ratio_threshold_dark = config. template get<float>("threshold_dark", 0.03);
      _value_ratio_threshold_bright = config. template get<float>("threshold_bright", 0.03);
      _sample_maxnum_dark = config. template get<int>("sample_dark", 40);
      _sample_maxnum_bright = config. template get<int>("sample_bright", 40);
      _line_tile_size = config. template get<int>("tile_size", (std::numeric_limits<int>::max)());
      _line_tile_overlap = config. template get<int>("tile_overlap", -1);
      bool preferred_a = config. template get<bool>("a", true);
      _mode = preferred_a?xlingsky::raster::radiometric::NucCalculator::SCALE:xlingsky::raster::radiometric::NucCalculator::OFFSET;
      std::string t = config. template get<std::string>("interp", "pchip");
      if(t=="makima")
        _interp_type = xlingsky::raster::radiometric::NucCalculator::MAKIMA;
      else if(t=="barycentric")
        _interp_type = xlingsky::raster::radiometric::NucCalculator::BARYCENTRIC;
      else
        _interp_type = xlingsky::raster::radiometric::NucCalculator::PCHIP;
      if(_line_tile_overlap<0) _line_tile_overlap = _line_tile_size/2;
      return true;
  }
  bool operator()(int b, int, int, void *data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
    float factor_bright = 1-_value_ratio_threshold_bright;
    float factor_dark = 1+_value_ratio_threshold_dark;
    struct BW {
      double v_dark;
      double v_bright;
      int cnt_dark;
      int cnt_bright;
      int stat;
      BW() : cnt_dark(0), cnt_bright(0), v_dark(0), v_bright(0), stat(0) {}
    };
    std::vector<BW> bws(rows);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for (int r = 0; r < rows; ++r) {
      DataType* t = pdata+r*cols;
      DataType* end = RemoveNoData(t, t+cols);
      std::sort( t, end);
      int n = end-t;
      int st = int(n*_cut_ratio_dark);
      int ed = n - 1 - int(n * _cut_ratio_bright);
      if (ed < st) {
        bws[r].stat = rows+1;
        continue;
      }
      int cnt_dark = 1, cnt_bright = 1;
      double v_dark = (double)t[st], v_bright = (double)t[ed];
      while (st + cnt_dark < ed && cnt_dark < _sample_maxnum_dark) {
        double v = v_dark/cnt_dark;
        v *= factor_dark;
        if(t[st+cnt_dark]>v) break;
        v_dark += t[st+cnt_dark];
        ++cnt_dark;
      }
      while (ed-cnt_bright>st && cnt_bright < _sample_maxnum_bright) {
        double v = v_bright/cnt_bright;
        v *= factor_bright;
        if(t[ed-cnt_bright] < v ) break;
        v_bright += t[ed-cnt_bright];
        ++cnt_bright;
      }
      bws[r].v_dark = v_dark;
      bws[r].v_bright = v_bright;
      bws[r].cnt_dark = cnt_dark;
      bws[r].cnt_bright = cnt_bright;
      bws[r].stat = (cnt_dark >= _sample_minnum_dark &&
                     cnt_bright >= _sample_minnum_bright &&
                     (v_dark / cnt_dark) * factor_dark <
                     (v_bright / cnt_bright) * factor_bright)?0:1;
    }

// #if defined(DEBUG) || defined(_DEBUG)
//       {
//           boost::filesystem::path dirpath(boost::dll::program_location().parent_path());
//           boost::filesystem::path filepath = dirpath;
//           char name[128];
//           sprintf(name, "b%d.tif", b);
//           filepath /= name;
//           GDALDataset* dataset = GDALCreate(filepath.string().c_str(), cols, rows, 1, GDT_UInt16);
//           if(dataset->RasterIO(GF_Write, 0, 0, cols, rows, data, cols, rows, gdal::DataType<DataType>::type(), 1, nullptr, 0, 0, 0)){}
//           GDALClose(dataset);
//       }
// #endif

    xlingsky::TileManager manager;
    manager.AppendDimension(rows, _line_tile_size>rows?rows:_line_tile_size, _line_tile_overlap>rows?rows:_line_tile_overlap, _line_tile_overlap);

    std::vector<BW> tile_sum(manager.Size(0));
    std::vector< std::pair<int, float> > valid_tiles;
    for(int t=0; t<manager.Size(0); ++t){
      auto seg = manager.Segment(0, t);
      auto& sum = tile_sum[t];
      int tilesize = (t+1==manager.Size(0)?seg.second:_line_tile_size);
      double v_dark[2] = {0,0}, v_bright[2] = {0,0};
      int cnt_dark[2]={0,0}, cnt_bright[2]={0,0};
      for(int r=seg.first; r<seg.first+tilesize; ++r){
        switch(bws[r].stat){
          case 0:{
            v_dark[0] += bws[r].v_dark;
            v_bright[0] += bws[r].v_bright;
            cnt_dark[0] += bws[r].cnt_dark;
            cnt_bright[0] += bws[r].cnt_bright;
          }break;
          case 1:{
            v_dark[1] += bws[r].v_dark+bws[r].v_bright;
            cnt_dark[1] += bws[r].cnt_dark+bws[r].cnt_bright;
          }break;
          default:
            break;
        }
      }
      if(cnt_bright[0]>0){
        valid_tiles.emplace_back(t, seg.first+(double)tilesize/2.0);
        sum.stat = 0;
        sum.cnt_dark = cnt_dark[0];
        sum.cnt_bright = cnt_bright[0];
        sum.v_dark = v_dark[0]/cnt_dark[0];
        sum.v_bright = v_bright[0]/cnt_bright[0];
      }else if(cnt_dark[1]>0){
        sum.stat = 1;
        sum.cnt_dark = cnt_dark[1];
        sum.v_dark = v_dark[1]/cnt_dark[1];
      }else{
        sum.stat = 2;
      }
    }
    if(valid_tiles.size()==0){
      for (int t = 0; t < manager.Size(0); ++t) {
        auto seg = manager.Segment(0, t);
        auto &sum = tile_sum[t];
        int tilesize = (t+1==manager.Size(0)?seg.second:_line_tile_size);
        if(sum.stat==1){
          for (int r = seg.first; r < seg.first + tilesize; ++r) {
            double v = bws[r].v_dark+bws[r].v_bright;
            int cnt = bws[r].cnt_dark+bws[r].cnt_bright;
            v /= cnt;
            switch(_mode){
              case OFFSET:{
                _a.push_back(1);
                _b.push_back(v-sum.v_dark);
                _defective_pixels.push_back(0);
              }break;
              default:{
                _a.push_back(v/sum.v_dark);
                _b.push_back(0);
                _defective_pixels.push_back(0);
              };
            }
          }
        }else{
          for (int r = 0; r < tilesize; ++r) {
            _a.push_back(1);
            _b.push_back(0);
            _defective_pixels.push_back(1);
          }
        }
      }
      return true;
    }

    std::vector<double> dn_high(rows, 0), dn_low(rows, 0);
    if (valid_tiles.size() == 1) {
      for (int r = 0; r < rows; ++r) {
        dn_high[r] = tile_sum[valid_tiles.front().first].v_bright;
        dn_low[r] = tile_sum[valid_tiles.front().first].v_dark;
      }
    }else {
      xlingsky::InterpolatorAdaptor* interp = nullptr;
      std::vector<double> x1, x2, y;
      x1.reserve(valid_tiles.size()+2);
      x1.push_back(0);
      for (auto& t : valid_tiles) x1.push_back(t.second);
      x1.push_back(rows);
      x2 = x1;

      y.reserve(valid_tiles.size()+2);
      y.push_back(tile_sum[valid_tiles.front().first].v_bright);
      for (auto& t : valid_tiles) y.push_back(tile_sum[t.first].v_bright);
      y.push_back(tile_sum[valid_tiles.back().first].v_bright);
      switch (_interp_type) { 
          case PCHIP:
            interp = new xlingsky::pchip(std::move(x1), std::move(y));
            break;
          case MAKIMA:
            interp = new xlingsky::makima(std::move(x1), std::move(y));
            break;
          case BARYCENTRIC:
          default:
            interp =
                new xlingsky::barycentric_rational(std::move(x1), std::move(y));
            break;
      }
      for (int r = 0; r < rows; ++r) dn_high[r] = interp->operator()(r);
      delete interp;

      y.reserve(valid_tiles.size()+2);
      y.push_back(tile_sum[valid_tiles.front().first].v_dark);
      for (auto& t : valid_tiles) y.push_back(tile_sum[t.first].v_dark);
      y.push_back(tile_sum[valid_tiles.back().first].v_dark);
      switch (_interp_type) { 
          case PCHIP:
            interp = new xlingsky::pchip(std::move(x2), std::move(y));
            break;
          case MAKIMA:
            interp = new xlingsky::makima(std::move(x2), std::move(y));
            break;
          case BARYCENTRIC:
          default:
            interp =
                new xlingsky::barycentric_rational(std::move(x2), std::move(y));
            break;
      }
      for (int r = 0; r < rows; ++r) dn_low[r] = interp->operator()(r);
      delete interp;
    }

    if(_hi_path[0] || _lo_path[0]){
      for(int r=0; r<rows; ++r){
        _hi.push_back(dn_high[r]);
        _lo.push_back(dn_low[r]);
      }
    }

    for (int r = 0; r < rows; ++r) {
      if (bws[r].stat < 2) {
        auto su = bws[r].v_bright / bws[r].cnt_bright;
        auto sl = bws[r].v_dark / bws[r].cnt_dark;
        if (su - sl > std::numeric_limits<double>::epsilon()) {
          auto tu = dn_high[r];
          auto tl = dn_low[r];
          _a.push_back((tu - tl) / (su - sl));
          _b.push_back(tl - _a.back() * sl);
          _defective_pixels.push_back(0);
          continue;
        }
      }
      _a.push_back(1);
      _b.push_back(0);
      _defective_pixels.push_back(1);
    }
    /*
    {
      const double sigma = (double)rows*rows/2;
      for (int t = 0; t < tile_sum.size(); ++t) {
        auto& sum = tile_sum[t];
        if (sum.stat == 1) {
          auto seg = manager.Segment(0, t);
          int tilesize = (t+1==manager.Size(0)?seg.second:_line_tile_size);
          double x = seg.first + (double)tilesize / 2;
          double w = 0, hi = 0, lo = 0;
          for (int i = 0; i < valid_tiles.size(); ++i) {
            double wt = exp(-pow(x-valid_tiles[i].second,2.0)/sigma);
            w += wt;
            hi += wt*tile_sum[valid_tiles[i].first].v_bright;
            lo += wt*tile_sum[valid_tiles[i].first].v_dark;
          }
          hi /= w;
          lo /= w;
          if (sum.v_dark < (hi + lo) / 2) {
            sum.v_bright = hi;
          } else {
            sum.v_bright = sum.v_dark;
            sum.v_dark = lo;
          }
          sum.cnt_dark = 1;
          sum.cnt_bright = 1;
        }
      }
    }
    std::vector<double> dn_high(rows,0), dn_low(rows,0), dn_w(rows, 0);
    for (int t = 0; t < manager.Size(0); ++t) {
      auto& sum = tile_sum[t];
      if(sum.stat>1) continue;
      auto seg = manager.Segment(0, t);
      int i=0;
      while(i<seg.second/2){
        double w = i+1;
        int id = seg.first+i;
        dn_high[id] += w*sum.v_bright;
        dn_low[id] += w*sum.v_dark;
        dn_w[id] += w;
        ++i;
      }
      while(i<seg.second){
        double w = seg.second-i;
        int id = seg.first+i;
        dn_high[id] += w*sum.v_bright;
        dn_low[id] += w*sum.v_dark;
        dn_w[id] += w;
        ++i;
      }
    }

    if(_hi_path[0] || _lo_path[0]){
      for(int r=0; r<rows; ++r){
        if(dn_w[r]>std::numeric_limits<double>::epsilon()){
          _hi.push_back(dn_high[r]/dn_w[r]);
          _lo.push_back(dn_low[r]/dn_w[r]);
        }else{
          _hi.push_back(-1);
          _lo.push_back(-1);
        }
      }
    }

    for(int r=0; r<rows; ++r){
      if(bws[r].stat<2 && dn_w[r] > std::numeric_limits<double>::epsilon()){
        auto su = bws[r].v_bright/bws[r].cnt_bright;
        auto sl = bws[r].v_dark/bws[r].cnt_dark;
        if(su-sl>std::numeric_limits<double>::epsilon()){
          auto tu = dn_high[r]/dn_w[r];
          auto tl = dn_low[r]/dn_w[r];
          _a.push_back((tu-tl)/(su-sl));
          _b.push_back(tl-_a.back()*sl);
          _DefectivePixels.push_back(0);
          continue;
        }
      }
      _a.push_back(1);
      _b.push_back(0);
      _DefectivePixels.push_back(1);
    }*/

    return true;
  }
};

};  // namespace radiometric
};  // namespace raster

};  // namespace xlingsky

#endif
