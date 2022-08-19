#ifndef TILE_MANAGER_HPP
#define TILE_MANAGER_HPP

#include <vector>

namespace xlingsky{

class TileManager{
 public:
  using IndexType = int;
  using Seg = std::pair<IndexType, IndexType>;
  using SegContainer = std::vector<Seg>;
  using TileInfo = std::vector<Seg>;
  static SegContainer Tiling(IndexType length, IndexType individual_length, IndexType overlap_length, IndexType length_to_merge = 0){
    IndexType count = length/individual_length;
    if(count<1) count = 1;
    if((double)length-(double)count*individual_length > (double)length_to_merge) ++count;
    SegContainer segs(count);
    IndexType tile_size = individual_length+overlap_length;
    IndexType i=0;
    for(; i<count-1; ++i){
      IndexType st = i*individual_length;
      IndexType l = length-st;
      segs[i] = std::make_pair(st, l<tile_size?l:tile_size);
    }
    IndexType st = i*individual_length;
    segs[i] = std::make_pair(st, length-st);
    return segs;
  }
  TileManager(){}
  virtual ~TileManager(){}
  void AppendDimension(IndexType length, IndexType individual_length, IndexType overlap_length = 0, IndexType length_to_merge = 0){
    _info.emplace_back(Tiling(length, individual_length, overlap_length, length_to_merge));
    if(_accumulation.size()==0) _accumulation.push_back(_info.back().size());
    else _accumulation.push_back(_accumulation.back()*_info.back().size());
  }
  IndexType Dims() const { return _info.size(); }
  IndexType Size() const { return _accumulation.back(); }
  IndexType Size(int dim) const { return _info[dim].size(); }
  TileInfo Tile(IndexType i) const {
    std::vector<IndexType> indices(Dims());
    for (IndexType d = Dims()-1; d>0; --d) {
      indices[d] = i/_accumulation[d-1];
      i -= indices[d]*_accumulation[d-1];
    }
    indices[0] = i;
    TileInfo t(Dims());
    for(IndexType d=0; d<Dims(); ++d)
      t[d] = Segment(d, indices[d]);
    return t;
  }
  Seg Segment(IndexType dim, IndexType i) const {
    return _info[dim][i];
  }

protected:
  std::vector<SegContainer> _info;
  std::vector<IndexType> _accumulation;
};

};

#endif
