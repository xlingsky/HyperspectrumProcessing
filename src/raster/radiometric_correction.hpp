#ifndef RADIOMETRIC_CORRECTION_HPP
#define RADIOMETRIC_CORRECTION_HPP

#include <assert.h>
#include <vector>
#include <sstream>
#include <string>
#include <boost/filesystem.hpp>

#include "RasterOperator.h"
#include "inpaint.hpp"
#include "TileManager.hpp"
#include "InterpolatorAdaptor.hpp"

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
//#define BADPIXEL_COUNTING

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
    if (fscanf(fp, "%lf", &t) != 1) break;
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
#ifdef BADPIXEL_COUNTING
  std::vector<std::vector<int> > _list;
#endif

 public:
  PixelCorrection(int cols, int rows)
      : _cols(cols), _rows(rows),
#ifdef BADPIXEL_COUNTING
    _list(rows),
#endif
    _threshold(10) {}
  virtual DataType correct(int, DataType, int) = 0;

  bool operator()(int b, int xoff, int yoff, void* data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
// #ifdef _USE_OPENMP
// #if _USE_OPENMP > 4
// #pragma omp declare reduction (merge : std::vector<int> :
// omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) #pragma omp
// parallel for reduction(merge: _list) #else #pragma omp parallel for
// reduction(+: cnt) #endif #endif
#ifndef BADPIXEL_COUNTING
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
        if (d < 0) {
          pdata[i] = 0;
#ifdef BADPIXEL_COUNTING
          if (d < -_threshold) _list[b].push_back(i);
#endif
        } else
          pdata[i] = d;
      }
    }
    return true;
  }
  int cols() const { return _cols; }
  int rows() const { return _rows; }
    
#ifdef BADPIXEL_COUNTING
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
  DarkBackgroundCorrection(int cols, int rows)
      : _data(new DataType[(size_t)cols * rows]),
        PixelCorrection(cols, rows) {}
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
    DarkBackgroundLinear(int cols, int rows, int bands) : PixelCorrection(cols, rows), _bands(bands){
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
                char strline[512];
                while(fgets(strline, 512, fp) && cnt<_bands){
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
  NonUniformCorrection(int cols, int rows)
      : PixelCorrection(cols, rows) {
    size_t sz = (size_t)cols * rows;
    _a = new DataType[sz];
    _b = new DataType[sz];
  }
  ~NonUniformCorrection() {
    if (_a) delete[] _a;
    if (_b) delete[] _b;
  }
  bool load(const char* a, const char* b) {
    if (::xlingsky::raster::radiometric::load(a, (size_t)cols() * rows(), _a) &&
        ::xlingsky::raster::radiometric::load(b, (size_t)cols() * rows(), _b))
      return true;
    return false;
  }
  DataType correct(int, DataType d, int i) override { return _a[i] * d + _b[i]; }
};

class BadPixelCorrection : public FrameIterator {
 protected:
  cv::Mat _mask;

 public:
  typedef float DataType;
  BadPixelCorrection(int cols, int rows)
      : _mask(cv::Mat::zeros(rows, cols, CV_8U)) {}
  bool load(const char* filepath) {
    FILE* fp = fopen(filepath, "r");
    if (fp == nullptr) return false;
    char strline[512];
    while (fgets(strline, 512, fp)) {
      int r, c;
      if (sscanf(strline, "%d%d", &r, &c) == 2) {
        if (r > 0 && r <= _mask.rows && c > 0 && c <= _mask.cols)
          _mask.at<unsigned char>(r - 1, c - 1) = 1;
      }
    }
    return true;
  }
  bool operator()(int b, int xoff, int yoff, void* data, int cols, int rows) override {
    cv::Mat mask = _mask(cv::Rect(xoff, yoff, cols, rows));

    cv::Mat m(rows, cols, cv::DataType<DataType>::type, data);
    InpaintOp op(m, m, mask);
    xlingsky::TileManager manager;
    const int margin = 10;
    const int tilesize = rows;
    manager.AppendDimension(cols, tilesize, margin, margin);
    manager.AppendDimension(rows, tilesize, margin, margin);

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
    return true;
  }
};

class MeanStdCalculator : public FrameIterator {
 private:
  std::vector<double> _mean;
  std::vector<double> _std;
  int _width;
  char _filepath[512];
  float _cut_ratio_upper;
  float _cut_ratio_lower;

 public:
  typedef float DataType;
  MeanStdCalculator(int w, float cut_ratio_lower = 0, float cut_ratio_upper = 0) : _width(w), _cut_ratio_lower(cut_ratio_lower), _cut_ratio_upper(cut_ratio_upper) {}
  ~MeanStdCalculator() { save(_filepath); }
  std::pair<double, double> compute(DataType* data, int n, float cut_ratio_lower, float cut_ratio_upper) {
    std::pair<double, double> ret = std::make_pair(0, 0);
    int cnt = 0;
    std::sort(data, data+n);
    int st, ed;
    n = FindValid(data, n, st, ed);
    st += int(n*cut_ratio_lower);
    ed -= int(n*cut_ratio_upper);

    int i = st;
    while(i<=ed){
      ret.first += data[i];
      ret.second += data[i] * data[i];
    }
    n = ed-st+1;
    if(n>0){
      ret.first /= n;
      ret.second = std::sqrt(ret.second / n - ret.first * ret.first);
    }
    return ret;
  }
  void SetFilePath(const char* filepath) { strcpy(_filepath, filepath); }
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
  ~MedianCalculator() { save(_filepath); }
  void SetFilePath(const char* filepath) { strcpy(_filepath, filepath); }
  DataType compute(DataType* data, int n) {
    std::sort(data, data + n);
    int st, ed;
    FindValid(data, n, st, ed);
    return data[(st+ed) >> 1];
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
    BSPLINE_CUBIC,
    BSPLINE_QUINTIC,
    BSPLINE_QUADRATIC,
    BARYCENTRIC
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

  std::vector<double> _a;
  std::vector<double> _b;
  std::vector<double> _hi;
  std::vector<double> _lo;
  std::vector<int> _badpixels;
  int _width;

  char _a_path[512];
  char _b_path[512];
  char _bp_path[512];
  char _xml_path[512];
  char _hi_path[512];
  char _lo_path[512];
 public:
  NucCalculator(int w, float cut_ratio_dark = 0.03, float cut_ratio_bright = 0.1, float ratio_threshold_dark = 0.03, float ratio_threshold_bright = 0.03, int sample_num_dark = 40, int sample_num_bright = 40, int line_tile_size = 30, int line_tile_overlap = 15, int mode = SCALE) : _width(w), _cut_ratio_dark(cut_ratio_dark), _cut_ratio_bright(cut_ratio_bright), _sample_maxnum_dark(sample_num_dark), _sample_maxnum_bright(sample_num_bright), _sample_minnum_dark(3), _sample_minnum_bright(3), _value_ratio_threshold_dark(ratio_threshold_dark), _value_ratio_threshold_bright(ratio_threshold_bright), _mode(mode), _line_tile_size(line_tile_size), _line_tile_overlap(line_tile_overlap) {
    _a_path[0] = _b_path[0] = _bp_path[0] = _xml_path[0] = 0;
    _hi_path[0] = _lo_path[0] = 0;
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
    if(_bp_path[0] && ::xlingsky::raster::radiometric::save(_bp_path, _badpixels.data(), _badpixels.size(), _width)) {
#ifdef _LOGGING
      VLOG(_LOG_LEVEL_RADIOMETRIC) << "bad pixel list was saved to " << _bp_path;
#endif
    }
    if(_xml_path[0] && _a_path[0] && _b_path[0]){
      FILE* fp = fopen(_xml_path, "w");
      if(fp){
        fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fp, "<HSP>\n");
        fprintf(fp, "\t<dim_prior>0,2,1</dim_prior>\n");
        fprintf(fp, "\t<task name=\"uniform\">\n");
        fprintf(fp, "\t\t<a>%s</a>\n", _a_path);
        fprintf(fp, "\t\t<b>%s</b>\n", _b_path);
        fprintf(fp, "\t</task>\n");
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
  void SetFilePath(const char *apath, const char *bpath, const char *bppath, const char* xmlpath, const char* hi_path, const char* lo_path) {
    if(apath) strcpy(_a_path, apath);
    if(bpath) strcpy(_b_path, bpath);
    if(bppath) strcpy(_bp_path, bppath);
    if(xmlpath) strcpy(_xml_path, xmlpath);
    if(hi_path) strcpy(_hi_path, hi_path);
    if(lo_path) strcpy(_lo_path, lo_path);
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
      std::sort( t, t+cols);
      int st, ed;
      int n = FindValid(t, cols, st, ed);
      st += int(n*_cut_ratio_dark);
      ed -= int(n*_cut_ratio_bright);
      if (ed - st < 0) {
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
                _badpixels.push_back(0);
              }break;
              default:{
                _a.push_back(v/sum.v_dark);
                _b.push_back(0);
                _badpixels.push_back(0);
              };
            }
          }
        }else{
          for (int r = 0; r < tilesize; ++r) {
            _a.push_back(1);
            _b.push_back(0);
            _badpixels.push_back(1);
          }
        }
      }
      return true;
    }
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
    // {
    //   xlingsky::InterpolatorAdaptor* interp = nullptr;
    //   switch
    // }
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
          _badpixels.push_back(0);
          continue;
        }
      }
      _a.push_back(1);
      _b.push_back(0);
      _badpixels.push_back(1);
    }

    return true;
  }
};

};  // namespace radiometric
};  // namespace raster

};  // namespace xlingsky

#endif
