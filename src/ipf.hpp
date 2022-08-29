#ifndef IMAGE_PROCESS_FRAMEWORK_HPP
#define IMAGE_PROCESS_FRAMEWORK_HPP

#include "TileManager.hpp"

namespace xlingsky {

namespace ipf{

template<class Operator>
bool DimProcessing(xlingsky::TileManager& manager, Operator* op, int dim) {
    for (int i = 0; i < manager.Size(dim); ++i) {
        auto seg = manager.Segment(dim, i);
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
    DimProcessing(manager, op, dims-1);
    return true;
}

};

};
#endif
