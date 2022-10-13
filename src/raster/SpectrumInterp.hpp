#ifndef SPECTRUM_INTERP_HPP
#define SPECTRUM_INTERP_HPP

#include "raster/RasterOperator.h"
#include "InterpolatorAdaptor.hpp"

namespace xlingsky {

namespace raster {
namespace spectrum {

class Interpolator : public Operator {
 public:
  typedef float DataType;
  enum InterpType {
    BSPLINE_CUBIC,
    BSPLINE_QUINTIC,
    BSPLINE_QUADRATIC,
    PCHIP,
    MAKIMA,
    BARYCENTRIC
  };

 protected:
  std::vector<double> _wl_old;
  std::vector<double> _wl_new;
  std::vector<DataType> _temp;
  InterpType _type;
  int _num_lines;
  int _num_samples_old;
  int _num_samples_new;

 public:
  Interpolator(int num_lines, std::vector<double>& wl_old, int num_samples_old, std::vector<double>& wl_new, int num_samples_new)
      : _type(PCHIP) {
    SetWaveLength(num_lines, wl_old, num_samples_old,  wl_new, num_samples_new);
  }
  virtual ~Interpolator() {}
  void SetInterpType(InterpType type) { _type = type; }
  bool SetWaveLength(int num_lines, std::vector<double>& wl_old, int num_samples_old, std::vector<double>& wl_new, int num_samples_new) {
    if(wl_old.size() >= (size_t)num_samples_old*num_lines){
      _wl_old.insert(_wl_old.end(), wl_old.begin(), wl_old.begin()+(size_t)num_samples_old*num_lines);
    }else if(wl_old.size()>=(size_t)num_samples_old){
      for(int i=0; i<num_lines; ++i)
        _wl_old.insert(_wl_old.end(), wl_old.begin(), wl_old.begin()+num_samples_old);
    }else return false;

    if(wl_new.size() >= (size_t)num_samples_new*num_lines){
      _wl_new.insert(_wl_new.end(), wl_new.begin(), wl_new.begin()+(size_t)num_samples_new*num_lines);
    }else if(wl_new.size()>=(size_t)num_samples_new){
      for(int i=0; i<num_lines; ++i)
        _wl_new.insert(_wl_new.end(), wl_new.begin(), wl_new.begin()+num_samples_new);
    }else return false;

    _num_lines = num_lines;
    _num_samples_new = num_samples_new;
    _num_samples_old = num_samples_old;
    return true;
  }
  bool operator()(void* data, int imoff[3], int size[3], int space[3],
                  int prior[3]) override {
    int dim_order[3] = {2,0,1};
    if(false){
      dim_order[0] = dim_order[1] = dim_order[2] = -1;
      for(int i=0; i<3; ++i){
        if(dim_order[1]<0 && _num_lines ==  size[i]) dim_order[1] = i;
        else if(dim_order[0]<0 && _num_samples_old == size[i]) dim_order[0] = i;
        else dim_order[2] = i;
      }
    }
    bool ret = true;
    _temp.resize((size_t)size[0] * size[1] * size[2]);
    memcpy(_temp.data(), data, sizeof(DataType) * _temp.size());
    for (int i = 0; i < 3; ++i) space[i] /= sizeof(DataType);
    int num_samples_old = _num_samples_old;
    int num_samples_new = _num_samples_new;
    int size_new[3];
    size_new[dim_order[0]] = num_samples_new;
    size_new[dim_order[1]] = size[dim_order[1]];
    size_new[dim_order[2]] = size[dim_order[2]];
    int space_new[3];
    {
      int i=0;
      while(i<3) {
        space_new[prior[i]] = space[prior[i]];
        if(prior[i++] == dim_order[0]) { break; }
      }
      while(i<3){
        space_new[prior[i]] = space[prior[i]] / num_samples_old*num_samples_new;
        ++i;
      }
    }
    for (int r = 0; r < size[dim_order[2]]; ++r) {
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
      for (int c = 0; c < size[dim_order[1]]; ++c) {
        std::vector<double> y(num_samples_old), x(num_samples_old);
        DataType* pdata = _temp.data()+c * space[dim_order[1]] + r * space[dim_order[2]];
        double * pwl = _wl_old.data()+c*num_samples_old;
        for (int b = 0; b < size[dim_order[0]]; ++b){
          y[b] = pdata[b * space[dim_order[0]]];
          x[b] = pwl[b];
        }
        xlingsky::InterpolatorAdaptor* interp = nullptr;
        switch (_type) {

#if 1
          case BSPLINE_CUBIC:
            interp = new xlingsky::cardinal_cubic_b_spline(
                y.data(), y.size(), x[0], x[1] - x[0]);
            break;
          case BSPLINE_QUINTIC:
            interp = new xlingsky::cardinal_quintic_b_spline(
                y.data(), y.size(), x[0], x[1] - x[0]);
            break;
          case PCHIP:
            interp = new xlingsky::pchip(std::move(x), std::move(y));
            break;
          case MAKIMA:
            interp = new xlingsky::makima(std::move(x), std::move(y));
            break;
#endif
          case BSPLINE_QUADRATIC:
            interp = new xlingsky::cardinal_quadratic_b_spline(
                y.data(), y.size(), x[0], x[1] - x[0]);
            break;
          case BARYCENTRIC:
          default:
            interp =
                new xlingsky::barycentric_rational(std::move(x), std::move(y));
            break;
        }
        pdata = (DataType*)data + c * space_new[dim_order[1]] + r * space_new[dim_order[2]];
        pwl = _wl_new.data()+c*num_samples_new;
        for (int b = 0; b < size_new[dim_order[0]]; ++b)
          *(pdata + b * space_new[dim_order[0]]) = (DataType)interp->operator()(pwl[b]);

        if (interp) delete interp;
      }
    }
    memcpy(size, size_new, sizeof(int) * 3);
    memcpy(space, space_new, sizeof(int) * 3);
    for (int i = 0; i < 3; ++i) space[i] *= sizeof(DataType);
    return ret;
  }
};

};  // namespace spectrum
};  // namespace raster

};  // namespace xlingsky

#endif
