#ifndef RADIOMETRIC_CORRECTION_HPP
#define RADIOMETRIC_CORRECTION_HPP

#include <assert.h>
#include <vector>
#include <fstream>

#include <opencv2/core/core.hpp>

#include "inpaint.hpp"
#include "ipf.hpp"

// #include <sstream>
//#define BADPIXEL_COUNTING

namespace radiometric{

template<typename T>
bool load(const char* filepath, size_t count, T* data){
  FILE* fp = fopen(filepath, "r");
  if(fp==nullptr) return false;
  size_t i;
  for ( i = 0; i < count; ++i) {
      double t;
      if (fscanf(fp, "%lf", &t) != 1) break;
      data[i] = (T)t;
  }
  
  fclose(fp);

  return i==count;
}

template <typename T> bool save(const char *filepath, T *data, size_t count, int newline) {
  std::ofstream fp(filepath);
  if(!fp.is_open()) return false;
  size_t i=0;
  while(i<count){
    fp << data[i] << "\t";
    ++i;
    if(i%newline==0){
      fp << std::endl;
    }
  }
  fp.close();
  return true;
}

using BufferOperator = ipf::BufferOperator;

class ComboOperator : public BufferOperator {
public:
    virtual ~ComboOperator() {
        for (auto& op : _ops)
            delete op;
        _ops.clear();
    }
    void Add(BufferOperator* op) {
        _ops.push_back(op);
    }
    bool operator()(void* data, int size[3], int space[3], int prior[3]) override {
        for (auto& op : _ops)
            if (!op->operator()(data,size,space,prior)) return false;
        return true;
    }
protected:
    std::vector<BufferOperator*> _ops;
};

class FrameIterator : public BufferOperator{
 public:
  bool operator()(void* data, int size[3], int space[3], int prior[3]) override{
      bool ret = true;
    for(int r=0; r<size[prior[2]]; ++r)
        if (this->operator()(r, (char*)data + (size_t)space[prior[2]]*r, size[prior[0]], size[prior[1]])) {   }
    return ret;
  }
  virtual bool operator()(int r, void* data, int cols, int rows) = 0;
};

class PixelCorrection : public FrameIterator{
 protected:
  int _cols;
  int _rows;
  int _threshold;
  std::vector<std::vector<int> > _list;
 public:
  typedef float DataType;
  PixelCorrection(int cols, int rows, int bands) : _cols(cols), _rows(bands), _list(rows), _threshold(10) {}
  virtual DataType correct(DataType , int ) = 0;

  bool operator()(int r, void* data, int cols, int rows) override{
    assert(cols==_cols);
    assert(rows==_rows);
    size_t sz = (size_t)cols*rows;
    DataType* pdata = (DataType*)data;
// #ifdef _USE_OPENMP
// #if _USE_OPENMP > 4
// #pragma omp declare reduction (merge : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
// #pragma omp parallel for reduction(merge: _list)
// #else
// #pragma omp parallel for reduction(+: cnt)
// #endif
// #endif
#ifndef BADPIXEL_COUNTING
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int i=0; i<sz; ++i){
      DataType d = correct(pdata[i], i);
      if(d<0){
        pdata[i] = 0;
#ifdef BADPIXEL_COUNTING
        if(d<-_threshold) _list[r].push_back(i);
#endif
      }else pdata[i] = d;
    }
    return true;
  }
  int cols() const { return _cols; }
  int rows() const { return _rows; }
  int max_bad_pixel_num() const {
    int cnt = 0;
    for(auto it=_list.begin(); it!=_list.end(); ++it)
      if(cnt < it->size()) cnt = it->size();
    return cnt;
  }
};

// dv is a bandsXcols  dark level composation value matrix stored in row-major
class DarkLevelCorrection : public PixelCorrection{
 public:
  typedef PixelCorrection::DataType DataType;
 protected:
  DataType* _data;
 public:
  DarkLevelCorrection(int cols, int rows, int bands) : _data(new DataType[(size_t)cols*bands]), PixelCorrection(cols, rows, bands) {}
  ~DarkLevelCorrection(){
    if(_data) delete[] _data;
  }
  bool load(const char* filepath){
    return ::radiometric::load(filepath, (size_t)cols()*rows(), _data);
  }
  DataType correct(DataType d, int i) override{
    return d-_data[i];
  }
};

class NonUniformCorrection : public PixelCorrection{
 public:
  typedef PixelCorrection::DataType DataType;
 protected:
  DataType* _a;
  DataType* _b;
 public:
  NonUniformCorrection(int cols, int rows, int bands) : PixelCorrection(cols, rows, bands) {
    size_t sz = (size_t)cols*bands;
    _a = new DataType[sz];
    _b = new DataType[sz];
  }
  ~NonUniformCorrection(){
    if(_a) delete[] _a;
    if(_b) delete[] _b;
  }
  bool load(const char* a, const char* b){
    if(::radiometric::load(a, (size_t)cols()*rows(), _a) && ::radiometric::load(b, (size_t)cols()*rows(), _b))
      return true;
    return false;
  }
  DataType correct(DataType d, int i) override{
    return _a[i]*d+_b[i];
  }
};

class BadPixelCorrection : public FrameIterator{
 protected:
  cv::Mat _mask;
 public:
  typedef float DataType;
  BadPixelCorrection(int cols, int , int bands) : _mask(cv::Mat::zeros(bands, cols, CV_8U)) {
  }
  bool load(const char* filepath){
    FILE* fp = fopen(filepath, "r");
    if(fp==nullptr) return false;
    char strline[512];
    while(fgets(strline, 512, fp)){
      int r, c;
      if(sscanf(strline, "%d%d", &r, &c)==2){
        if( r>0 && r<=_mask.rows && c>0 && c<=_mask.cols )
          _mask.at<unsigned char>(r-1, c-1) = 1;
      }
    }
    return true;
  }
  bool operator()(int r, void* data, int cols, int rows) override{
    assert(cols==_mask.cols);
    assert(rows==_mask.rows);

    cv::Mat m(rows, cols, cv::DataType<DataType>::type, data);
    InpaintOp op( m, m, _mask);
    ipf::Tile(cols, rows, rows, 10, op);

    return true;
  }
};

class MedianBlur : public FrameIterator{
 private:
  int _ksize;
 public:
  typedef float DataType;
  MedianBlur() : _ksize(5) {}
  bool operator()(int r, void *data, int cols, int rows) override {
    cv::Mat m(rows, cols, cv::DataType<DataType>::type, data);
    cv::medianBlur(m, m, _ksize);
    return true;
  }
};

class GaussianBlur : public FrameIterator {
private:
    int _ksize;
    int _stband;
public:
    typedef float DataType;
    GaussianBlur(int ksize, int stband) : _ksize(ksize), _stband(stband) {}
    bool operator()(int r, void* data, int cols, int rows) override {
        cv::Mat o(rows, cols, cv::DataType<DataType>::type, data);
        cv::Mat m = o(cv::Rect(0, _stband, cols, rows - _stband));
        cv::GaussianBlur(m, m, cv::Size(_ksize,_ksize), 0);
        return true;
    }
};

class MeanStdCalculator : public FrameIterator{
 private:
  std::vector<double> _mean;
  std::vector<double> _std;
  int _width;
 public:
  typedef float DataType;
  MeanStdCalculator(int w) : _width(w) {}
  std::pair<double, double> compute(DataType* data, int n){
    std::pair<double, double> ret = std::make_pair(0, 0);
    for(int i=0; i<n; ++i){
      ret.first += data[i];
      ret.second += data[i]*data[i];
    }
    ret.first /= n;
    ret.second = std::sqrt(ret.second/n-ret.first*ret.first);
    return ret;
  }
  bool save(const char* filepath){

    if(::radiometric::save(filepath, _mean.data(), _mean.size(), _width)){
      char path[512];
      strcpy(path, filepath);
      strcpy(strrchr(path, '.'), "_std.txt");
      if(::radiometric::save(path, _std.data(), _std.size(), _width))
        return true;
    }

    return false;
  }
  bool operator()(int r, void *data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
    std::vector<double> mean(rows), std(rows);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for(int r=0; r<rows; ++r){
      auto ret = compute(pdata+r*cols, cols);
      mean[r] = ret.first;
      std[r] = ret.second;
    }
    _mean.insert(_mean.end(), mean.begin(), mean.end());
    _std.insert(_std.end(), std.begin(), std.end());
    return true;
  }
};

class MedianCalculator : public FrameIterator{
 public:
  typedef float DataType;
 private:
  std::vector<DataType> _median;
  int _width;
 public:
  MedianCalculator(int w) : _width(w) {}
  DataType compute(DataType* data, int n){
    std::sort(data, data+n);
    return data[n>>1];
  }
  bool save(const char* filepath){
    return ::radiometric::save(filepath, _median.data(), _median.size(), _width);
  }
  bool operator()(int r, void *data, int cols, int rows) override {
    DataType* pdata = (DataType*)data;
    std::vector<double> median(rows);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
    for(int r=0; r<rows; ++r){
      median[r]= compute(pdata+r*cols, cols);
    }
    _median.insert(_median.end(), median.begin(), median.end());
    return true;
  }
};

};

#endif
