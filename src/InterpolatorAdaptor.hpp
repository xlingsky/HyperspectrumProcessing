#ifndef INTERPOLATOR_ADAPTOR_HPP
#define INTERPOLATOR_ADAPTOR_HPP

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>
#include <boost/math/interpolators/cardinal_quadratic_b_spline.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/pchip.hpp>

namespace xlingsky {

typedef std::vector<double> InterpValueContainer;

class InterpolatorAdaptor {
 public:
  typedef InterpValueContainer::value_type value_type;
  virtual ~InterpolatorAdaptor() {}
  virtual value_type operator()(value_type x) const = 0;
};

#define BEGIN_INTERPOLATOR_DERIVED_CLASS(name, datatype) \
class name : public InterpolatorAdaptor, public boost::math::interpolators::name<datatype> {\
public:\
  using Real = InterpolatorAdaptor::value_type;\
  using Base = boost::math::interpolators::name<datatype>;\
  Real operator()(Real x) const override{\
    return Base::operator()(x);\
  }\

#define END_INTERPOLATOR_DERIVED_CLASS(name, datatype) \
}\

#define INTERPOLATOR_SPLINE(name) \
BEGIN_INTERPOLATOR_DERIVED_CLASS(name, InterpolatorAdaptor::value_type)\
  name(const Real* const f, size_t length, Real left_endpoint, Real step_size) : Base(f, length, left_endpoint, step_size) {} \
END_INTERPOLATOR_DERIVED_CLASS(name, InterpolatorAdaptor::value_type)\

#define INTERPOLATOR_WITHOUT_DERIVATIVES(name)\
BEGIN_INTERPOLATOR_DERIVED_CLASS(name, InterpValueContainer)\
using Container = InterpValueContainer;\
  name(Container&& abscissas, Container&& ordinates) : Base(std::move(abscissas),std::move(ordinates)) {}\
END_INTERPOLATOR_DERIVED_CLASS(name, InterpValueContainer)\

#define INTERPOLATOR_WITH_DERIVATIVES(name)\
BEGIN_INTERPOLATOR_DERIVED_CLASS(name, InterpValueContainer)\
using Container = InterpValueContainer;\
  name(Container&& abscissas, Container&& ordinates, Container&& derivatives) : Base(std::move(abscissas),std::move(ordinates),std::move(derivatives)) {}\
END_INTERPOLATOR_DERIVED_CLASS(name, InterpValueContainer)\

INTERPOLATOR_SPLINE(cardinal_cubic_b_spline);
INTERPOLATOR_SPLINE(cardinal_quintic_b_spline);
INTERPOLATOR_SPLINE(cardinal_quadratic_b_spline);
INTERPOLATOR_WITHOUT_DERIVATIVES(pchip);
INTERPOLATOR_WITHOUT_DERIVATIVES(makima);
INTERPOLATOR_WITH_DERIVATIVES(cubic_hermite);

BEGIN_INTERPOLATOR_DERIVED_CLASS(barycentric_rational, InterpolatorAdaptor::value_type)
barycentric_rational(InterpValueContainer&& x, InterpValueContainer&& y, size_t order = 3) : Base(std::move(x), std::move(y), order) {}
END_INTERPOLATOR_DERIVED_CLASS(barycentric_rational, InterpolatorAdaptor::value_type);

};

#endif