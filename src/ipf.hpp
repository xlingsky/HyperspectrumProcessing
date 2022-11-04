#ifndef XLINGSKY_IMAGE_PROCESS_FRAMEWORK_HPP
#define XLINGSKY_IMAGE_PROCESS_FRAMEWORK_HPP

#include "TileManager.hpp"
#ifdef _LOGGING
#include <glog/logging.h>
#define _LOG_LEVEL_IPF 1
#endif

namespace xlingsky {

namespace ipf{

template<class Operator>
bool DimProcessing(xlingsky::TileManager& manager, Operator* op, int dim) {
    for (int i = 0; i < manager.Size(dim); ++i) {
        auto seg = manager.Segment(dim, i);
#ifdef _LOGGING
        VLOG(_LOG_LEVEL_IPF + 1)
            << "DIM[" << dim << "] SEG[" << i+1 << "/" << manager.Size(dim) << "]";
        VLOG(_LOG_LEVEL_IPF + 2) << "SEG start at " << seg.first << " with "
                                 << seg.second << " pixels";
#endif
        op->Begin(seg, dim);
        if (dim == 0) op->Apply();
        else DimProcessing(manager, op, dim-1);
        op->End(seg, dim);
    }
    return true;
}

template<class Operator, int dims = 3>
bool TileProcessing(int win_size[], int buffer_size[], Operator* op) {
    xlingsky::TileManager manager;
    using Seg = xlingsky::TileManager::Seg;
    for(int i=0; i<dims; ++i)
      manager.AppendDimension(win_size[i], buffer_size[i]);
#ifdef _LOGGING
    VLOG(_LOG_LEVEL_IPF) << "Total number of tiles is " << manager.Size();
#endif
    DimProcessing(manager, op, dims-1);
    return true;
}

};

};
#endif
