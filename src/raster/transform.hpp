#ifndef XLINGSKY_RASTER_TRANSFORM_HPP
#define XLINGSKY_RASTER_TRANSFORM_HPP

#include "transform_inl.hpp"
#define TRANSFORM_CONST
#include "transform_inl.hpp"
#undef TRANSFORM_CONST

#ifdef _USE_OPENMP
#define TRANSFORM_OMP
#include "transform_inl.hpp"

#define TRANSFORM_CONST
#include "transform_inl.hpp"
#undef TRANSFORM_CONST
#undef TRANSFORM_OMP
#endif

#endif
