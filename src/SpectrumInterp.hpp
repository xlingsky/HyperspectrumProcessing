#ifndef SPECTRUM_INTERP_HPP
#define SPECTRUM_INTERP_HPP

#include "RasterOperator.h"
#include "InterpolatorAdaptor.hpp"

namespace xlingsky {

namespace raster {
namespace spectrum {

class Interpolator : public Operator {
 public:
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
  std::vector<float> _temp;
  InterpType _type;

 public:
  Interpolator(std::vector<double>&& wl_old, std::vector<double>&& wl_new)
      : _wl_old(std::move(wl_old)), _wl_new(std::move(wl_new)), _type(PCHIP) {}
  virtual ~Interpolator() {}
  void SetInterpType(InterpType type) { _type = type; }
  void SetWaveLength(std::vector<double>& wl_old, std::vector<double>& wl_new) {
    _wl_old = wl_old;
    _wl_new = wl_new;
  }
  bool operator()(void* data, int size[3], int space[3],
                  int prior[3]) override {
    bool ret = true;
    _temp.resize((size_t)size[0] * size[1] * size[2]);
    memcpy(_temp.data(), data, sizeof(float) * _temp.size());
    for (int i = 0; i < 3; ++i) space[i] /= sizeof(float);
    int size_new[3] = {size[0], size[1], (int)_wl_new.size()};
    int space_new[3];
    if (prior[0] == 2) {
      space_new[0] = space[0] / _wl_old.size() * _wl_new.size();
      space_new[1] = space[1] / _wl_old.size() * _wl_new.size();
      space_new[2] = space[2];
    } else if (prior[1] == 2) {
      space_new[prior[0]] = space[prior[0]];
      space_new[prior[1]] = space[prior[1]];
      space_new[prior[2]] = space[prior[2]] / _wl_old.size() * _wl_new.size();
    } else {
      memcpy(space_new, space, sizeof(int) * 3);
    }
    for (int r = 0; r < size[1]; ++r) {
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
      for (int c = 0; c < size[0]; ++c) {
        std::vector<double> x(size[2]), wl_old(_wl_old);
        size_t id = c * space[0] + r * space[1];
        for (int b = 0; b < size[2]; ++b) x[b] = _temp[id + b * space[2]];
        xlingsky::InterpolatorAdaptor* interp = nullptr;
        switch (_type) {
          case BSPLINE_CUBIC:
            interp = new xlingsky::cardinal_cubic_b_spline(
                x.data(), x.size(), _wl_old[0], _wl_old[1] - _wl_old[0]);
            break;
          case BSPLINE_QUADRATIC:
            interp = new xlingsky::cardinal_quadratic_b_spline(
                x.data(), x.size(), _wl_old[0], _wl_old[1] - _wl_old[0]);
            break;
          case BSPLINE_QUINTIC:
            interp = new xlingsky::cardinal_quintic_b_spline(
                x.data(), x.size(), _wl_old[0], _wl_old[1] - _wl_old[0]);
            break;
          case PCHIP:
            interp = new xlingsky::pchip(std::move(wl_old), std::move(x));
            break;
          case MAKIMA:
            interp = new xlingsky::makima(std::move(wl_old), std::move(x));
            break;
          case BARYCENTRIC:
            interp =
                new xlingsky::barycentric_rational(std::move(wl_old), std::move(x));
            break;
        }
        float* pdata = (float*)data + c * space_new[0] + r * space_new[1];
        for (int b = 0; b < _wl_new.size(); ++b)
          *(pdata + b * space_new[2]) = (float)interp->operator()(_wl_new[b]);

        if (interp) delete interp;
      }
    }
    memcpy(size, size_new, sizeof(int) * 3);
    memcpy(space, space_new, sizeof(int) * 3);
    for (int i = 0; i < 3; ++i) space[i] *= sizeof(float);
    return ret;
  }
};

};  // namespace spectrum
};  // namespace raster

};  // namespace xlingsky

#endif