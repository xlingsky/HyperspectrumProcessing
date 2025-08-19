#ifndef XLINGSKY_ANOMALYDETECTION_HPP
#define XLINGSKY_ANOMALYDETECTION_HPP

#include <Eigen/Dense>
#include <cstddef>

#include "raster/operator.h"

namespace xlingsky {
namespace raster {

template<class Exporter>
class BackgroundExtraction : public Operator {
    public:
    typedef float DataType;
    static DataType average(DataType* data, int n, int stride) {
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += data[i * stride];
        }
        return (DataType)(sum / n);
    }
    static DataType median(DataType* data, int n, int stride) {
        std::vector<DataType> v(n);
        for (int i = 0; i < n; i++) {
            v[i] = data[i * stride];
        }
        std::sort(v.begin(), v.end());
        return v[n >> 1];
    }
    private:
    int _ksize;
    DataType (*_bg)(DataType*, int, int);
    Exporter* _exporter;
    public:
    BackgroundExtraction() : _exporter(nullptr), _ksize(30) {}
    void set_exporter(Exporter* exporter) {
        _exporter = exporter;
    }
    int get_ksize() const {
      return _ksize;
    }
    template<class Config>
    bool load_config(const Config& config) {
        _ksize = config.template get<int>("ksize", 30);
        std::string method = config.template get<std::string>("method", "average");
        if (method == "average")
            _bg = BackgroundExtraction::average; 
        else
            _bg = BackgroundExtraction::median;
        return true;
    }
    int get_blknum(int size) const {
        return (int)std::ceil(size / (float)_ksize);
    }

    bool operator()(void *data, int imoff[3], int size[3], int space[3],
                    int prior[3]) override {
      int blknum = get_blknum(size[prior[2]]);
      const unsigned bytes = (unsigned)sizeof(DataType);
      unsigned stride[3] = {space[0] / bytes, space[1] / bytes, space[2] / bytes};
      std::vector<DataType> bgdata(size[prior[0]] * size[prior[1]] * blknum);
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
      for (int i=0; i<blknum; ++i) {
        int st = i * _ksize;
        int sz = i==blknum-1 ? size[prior[2]]-st : _ksize;
        DataType* pdata = (DataType*)data + st * stride[prior[2]];
        DataType* bgptr = bgdata.data() + (size_t)i * size[prior[0]] * size[prior[1]];
        for (int j=0; j<size[prior[0]]*size[prior[1]]; ++j) {
          *bgptr = (*_bg)(pdata, sz, stride[prior[2]]);
          ++pdata;
          ++bgptr;
        }
      }
      int imoff_new[3];
      imoff_new[prior[0]] = imoff[prior[0]];
      imoff_new[prior[1]] = imoff[prior[1]];
      imoff_new[prior[2]] = imoff[prior[2]]/_ksize;
      int size_new[3];
      size_new[prior[0]] = size[prior[0]];
      size_new[prior[1]] = size[prior[1]];
      size_new[prior[2]] = blknum;
      return _exporter->operator()(bgdata.data(), imoff_new, size_new, space, prior);
    }
};

class RXDetector : public Operator {
  private:
  int _ksize;
  public:
    typedef float DataType;
    RXDetector() {}
    template<class Config>
    bool load_config(const Config& config) {
      _ksize = config.template get<int>("ksize", 30);
      return true;
    }
    int get_blknum(int size) const {
        return (int)std::ceil(size / (float)_ksize);
    }
    Eigen::Matrix<DataType, Eigen::Dynamic, 1> run(DataType* data, size_t cols, size_t rows) {
      // Implement RX detection logic here
      Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>>
          mat(data, cols, rows);
      Eigen::Matrix<DataType, 1, Eigen::Dynamic> mean = mat.colwise().mean();
      Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> centered = mat.rowwise() - mean;
      Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic> cov = (centered.adjoint()*centered) / double(mat.rows() - 1);
      if (cov.determinant() == 0) {
        cov += Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>::Identity(cov.rows(), cov.rows()) * 1e-6;
      }
      auto inv_cov = cov.inverse();
      return (centered * inv_cov).cwiseProduct(centered).rowwise().sum();
    }
    bool operator()(void* data, int imoff[3], int size[3], int space[3],
                    int prior[3]) override {
      int blknum = get_blknum(size[prior[2]]);
      size_t stride = (size_t)size[prior[0]] * size[prior[1]];
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
      for (int i=0; i<blknum; ++i) {
        int st = i * _ksize;
        int sz = i==blknum-1 ? size[prior[2]]-st : _ksize;
        DataType* pdata = (DataType*)data + st * stride;
        Eigen::Map<Eigen::Matrix<DataType, Eigen::Dynamic, 1>> rx(
            (DataType *)data + i * stride,
            stride);
        rx.noalias() = run(pdata, stride, sz);
      }
      size[prior[2]] = blknum;
      return true;
    }
};

};
};

#endif