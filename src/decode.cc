#include "decode.h"

#include <stdio.h>
#include <iostream>

#include "gdalex.hpp"

#include "bigfileio.h"
// Segments parse(const Keyword::DataType* data, size_t count, const Keywords& keywords){
//   Segments segments;
//   std::vector<Keyword::Index> ks(keywords.size(), 0);
//   for(size_t i=0; i<count; ++i){
//     for(Segment::Index j=0; j<keywords.size(); ++j){
//       auto& key = keywords[j];
//       auto& k = ks[j];
//       while(k!=(Keyword::Index)-1 && data[i]!=key.key[k]){
//         k = key.lps[k];
//       }
//       ++k;
//       if(k==key.key.size()){
//         segments.emplace_back(j, data+i+1);
//         k=0;
//       }
//     }
//   }
//   if(segments.size()<1) return segments;
//   auto it1 = segments.begin();
//   auto it0 = it1++;
//   while(it1!=segments.end()){
//     segments.modify(it0, [&](Segment& s) {
//       s.count = it1->data-s.data-keywords[it1->id].key.size();
//     });
//   }
//   segments.modify(it0, [&](Segment& s) {
//     s.count = data+count-s.data-keywords[it1->id].key.size();
//   });
//   return segments;
// }

// typedef std::pair<Segment::DataType*, Segment::Size> SegInfo;

// static bool operator<(const SegInfo& l, const SegInfo& r){
//   return l.first < r.first;
// }
// static bool operator==(const SegInfo& l, const SegInfo& r){
//   return l.first == r.first;
// }

// typedef std::set<SegInfo> SegInfoContainer;

// SegInfoContainer find(Segments& segments, Segment::Index keyid) {
//   std::set<SegInfo> ret;
//   auto& segs = segments.get<1>();
//   auto range = segs.equal_range(keyid);
//   while(range.first!=range.second){
//     ret.insert( std::make_pair(range.first->data, range.first->count) );
//     ++range.first;
//   }
//   return ret;
// }

namespace HSP{

template<typename DstIterator>
size_t cvt_copy(const void* src, size_t count, DstIterator dst){
  const unsigned short* data = (const unsigned short*)src;
  size_t num = count>>1;
  for(size_t i=0; i<num; ++i, ++dst, ++data){
    *dst = *data;
  }
  return num;
}

class UnzipBase{
 public:
  virtual size_t apply(const void*, size_t) = 0;
};

bool Aux2Pos(const char* auxpath, const char* pospath){
  FILE* aux = fopen(auxpath, "rb");
  if(aux==nullptr) return false;
  FILE* pos = fopen(pospath, "w");
  if(pos==nullptr) { fclose(aux); return false; }

  while(1){
    int frame, count;
    if(fread(&frame, sizeof(int), 1, aux)!=1) break;
    if(fread(&count, sizeof(int), 1, aux)!=1) break;
    if(count==0) continue;
    std::vector<char> buffer(count);
    if(fread(buffer.data(), sizeof(char), count, aux)!=count) break;

    fprintf(pos, "%d\t", frame);
    std::vector<char> data(count*2.0/3);
    HSP::cvt_12to16( buffer.data(), count, data.begin());
    fprintf(pos, "%s\n", (char*)data.data()+1);
  }

  fclose(aux);
  fclose(pos);

  return true;
}

template<class DstIterator>
class UnzipT : public UnzipBase{
 protected:
  DstIterator _it;
 public:
  void SetDst(DstIterator& it) { _it = it; }
};

template<class DstIterator>
class Unzip12bits : public UnzipT<DstIterator>{
 public:
  typedef UnzipT<DstIterator> Base;
  static float ratio() { return 1.5; }
  size_t apply(const void* src, size_t count) override{
    return cvt_12to16(src, count, Base::_it);
  }
};

template<class DstIterator>
class Unzip16bits : public UnzipT<DstIterator>{
 public:
  typedef UnzipT<DstIterator> Base;
  static float ratio() { return 2; }
  size_t apply(const void* src, size_t count) override{
    return cvt_copy(src, count, Base::_it);
  }
};

Header::CompressionType Header::compression = Header::CT_LOSSLESS;
float Header::image_data_ratio = 1.5;
float Header::meta_data_ratio = 1.5;
Header::Unzip  Header::unzip_front = cvt_12to16<Header::DstIterator>;
Header::UnzipR Header::unzip_back = cvt_12to16<Header::DstIteratorR>;

void Header::SetCompression(Header::CompressionType type){
  switch(type){
    case CT_NONE:{
      Header::compression = Header::CT_NONE;
      Header::image_data_ratio = 2;
      Header::unzip_front = cvt_copy<Header::DstIterator>;
      Header::unzip_back = cvt_copy<Header::DstIteratorR>;
    }break;
    default:{
      Header::compression = Header::CT_LOSSLESS;
      Header::image_data_ratio = 1.5;
      Header::unzip_front = cvt_12to16<Header::DstIterator>;
      Header::unzip_back = cvt_12to16<Header::DstIteratorR>;
    }break;
  }
}

inline const void* extract(const void* data, Header& hdr){
  const unsigned char* p = (const unsigned char*)data;
  hdr.data_count = ((p[0]<<8)+p[1]);
  p+=2;
  if((*p&0xF0)==0x10) hdr.image_type = Header::IT_SWIR;
  else hdr.image_type = Header::IT_VNIR;
  if((*p&0x0F)==7) hdr.data_type = Header::DT_IMAGE;
  else hdr.data_type = Header::DT_META;
  ++p;
  // hdr.compression = (Header::CompressionType)(p[0]>>6);
  if(p[0]>0x02)
    hdr.band_id = (int)p[1];
  else
    hdr.band_id = cvt_12to16_odd( p[0], p[1] );
  p+=2;
  hdr.frame_id = (p[0]<<16)+(p[1]<<8)+p[2];
  p+=3;
  if( hdr.data_type == Header::DT_IMAGE){
    hdr.data_count = int(hdr.data_count*hdr.image_data_ratio);
  }else{
    hdr.data_count = int(hdr.data_count*hdr.meta_data_ratio);
  }
  return p;
}

class SecFile{
 private:
  size_t _buffer_capacity;
  HANDLE _fp;
  char* _data;
 public:
  SecFile(size_t capacity) : _buffer_capacity(capacity), _fp(INVALID_HANDLE_VALUE), _data(new char[capacity]) {}
  ~SecFile(){
    if(_data) delete[] _data;
    if(_fp) CloseHandle(_fp);
  }
  bool open(const char* filepath){
    if(_fp) CloseHandle(_fp);
    _fp = CreateFileE(filepath, GENERIC_READ);
    return _fp&&_fp!=INVALID_HANDLE_VALUE;
  }
  template<class Op>
  void run(Op& op){
    DWORD size = _buffer_capacity;
    char* pdata = _data;

    op.begin();

    while( ReadFile( _fp, pdata, size, &size, nullptr) && size>0 ){
      size = pdata+size-_data;
      pdata = _data;

      if(op.apply(pdata, size)){
        size = _data+size-pdata;
        if(size!=0) memcpy( _data, pdata, size);
        pdata = _data+size;
        size = _buffer_capacity-size;
      }else{
        pdata = _data;
        size = _buffer_capacity;
      }

    }

    op.end();
  }

  void reset(){
    SetFilePointer(_fp, 0, NULL, FILE_BEGIN);
  }
};

class SecOp{
 private:
  tokenizer::Token* _it_token;
 public:
  SecOp(tokenizer::Token* it) : _it_token(it) {}
  virtual void begin() =0;
  virtual void end() = 0;
  virtual void process( Header& hdr, const char* p, bool valid) = 0;
  virtual bool apply(char* &data, size_t size){
    std::vector<tokenizer::Section> sections;
    tokenizer::apply(data, data+size, _it_token, _it_token+1, std::back_inserter(sections));
    if(sections.size()==0){
      return false;
    }

    auto it1 = sections.begin();
    auto it0 = it1++;
    while(it1!=sections.end()){
      Header hdr;
      char* p = (char*)extract(it0->data, hdr);
      bool valid = (p+hdr.data_count+_it_token->size()<=it1->data);
      process(hdr, p, valid);
      it0 = it1++;
    }
    {
      Header hdr;
      char* p = (char*)extract(it0->data, hdr);
      bool valid = (p+hdr.data_count<=data+size);
      if(valid){
        process(hdr, p, valid);
        data = p+hdr.data_count;
      }else{
        data = (char*)(it0->data)-_it_token->size();
      }
    }
    return true;
  }
};

class CountOp : public SecOp{
 public:
  struct Info{
    int band_count;
    Info(int cnt) : band_count(cnt) {}
  };
  typedef std::vector<Info> InfoContainer;
 private:
  InfoContainer _infos;
  int _band_max_count;
  int _band_count;
  int _width;
  bool _meta_valid;
  Header::ImageType _image_type;
  // std::vector<int> _band_list;
 public:
  CountOp(tokenizer::Token* it) : SecOp(it) {}
  InfoContainer& infos() { return _infos; }
  int band_count() const { return _band_max_count+1; }
  int width() const { return _width; }
  Header::ImageType image_type() const { return _image_type; };
  void begin() override{
    _band_max_count = 0;
    _width = -1;
    _band_count = -1;
    _meta_valid = false;
    _infos.clear();
  }
  void process(Header& hdr, const char* p, bool valid) override{
    _image_type = hdr.image_type;
    if(hdr.data_type==Header::DT_META){
      if(_band_count>=0){
        _infos.emplace_back(_band_count);
      }
      _band_count = 0;
      // _band_list.clear();
      _meta_valid = valid;
    }else if(_meta_valid&&valid){
      if(_width<=0) _width = int(hdr.data_count/hdr.image_data_ratio);
      ++_band_count;
      if(hdr.band_id>_band_max_count) _band_max_count = hdr.band_id;
    }
    // if(hdr.data_type==Header::DT_IMAGE)
    //   _band_list.push_back(hdr.band_id);
  }
  void end() override{
    if(_band_count>0){
      _infos.emplace_back(_band_count);
    }
  }
};

class SaveOp : public SecOp{
 private:
  int _width;
  int _height;
  int _frame;
  GDALDataset* _image;
  HANDLE _meta;
  CountOp::InfoContainer& _infos;
  CountOp::InfoContainer::iterator _it_info;
  int _frame_id;
  bool _meta_valid;

  bool _flipping;
 public:
  SaveOp(
      int width, int height, int frame,
      const char* imagepath, const char* metapath,
      CountOp::InfoContainer& infos, tokenizer::Token* it
         )
      : _width(width), _height(height), _frame(frame), _infos(infos), _flipping(false), SecOp(it) {
    _image = GDALCreate(imagepath, _width, _frame, _height, GDT_UInt16);
    _meta = CreateFileE(metapath, GENERIC_WRITE);
  }
  ~SaveOp(){
    if(_image) GDALClose(_image);
    if(_meta&&_meta!=INVALID_HANDLE_VALUE) CloseHandle(_meta);
  }

  bool valid() const { return _image&&_meta&&_meta!=INVALID_HANDLE_VALUE; }
  void set_flipping(bool f) { _flipping = f; }

  void begin() override {
    _it_info = _infos.begin();
    _meta_valid = false;
    _frame_id = -1;
  }
  void end() override {}

  void process( Header& hdr, const char* p, bool ) override{
    if(hdr.data_type==Header::DT_META){
      _meta_valid = (_it_info->band_count==_height);
      if(_meta_valid){
        ++_frame_id;
        if(_meta){
          // fprintf(_meta, "%d\t", hdr.frame_id);
          // std::vector<char> data(hdr.data_count*2.0/3);
          // cvt_12to16( p, hdr.data_count, data.begin());
          // fprintf(_meta, "%s\n", (char*)data.data()+1);
          DWORD rw;
          WriteFile(_meta, &hdr.frame_id, sizeof(int), &rw, nullptr);
          WriteFile(_meta, &hdr.data_count, sizeof(int), &rw, nullptr);
          WriteFile(_meta, p, hdr.data_count, &rw, nullptr);
        }
      }
      ++_it_info;
    }else if(_meta_valid&&_image){
      std::vector<unsigned short> data(_width);

      if(_flipping){
        if(hdr.image_type==Header::IT_VNIR)
          hdr.unzip_back( p, std::ceil(_width*hdr.image_data_ratio), data.rbegin());
        else{
          int num_module = 4;
          int pixels_module = _width/num_module;
          int bytes_module = pixels_module*hdr.image_data_ratio;
          for(int i=0; i<num_module; ++i)
            hdr.unzip_back(p+(num_module-1-i)*bytes_module, bytes_module, data.rbegin()+pixels_module*i );
        }
      } else{
        hdr.unzip_front( p, std::ceil(_width*hdr.image_data_ratio), data.begin());
      }
      if(_image->GetRasterBand(hdr.band_id+1)->RasterIO(GF_Write, 0, _frame_id, _width, 1, &data[0], _width, 1, GDT_UInt16, 0, 0) != CE_None) {}
    }
  }
};

bool Decode::apply(const char* filepath, size_t buffer_capacity, const char* imagepath, const char* metapath){

    if(_verbose){
      std::cout << "[DECODE] Processing " << filepath << " ..." << std::endl;
      printf("[DECODE] Token= %x\n", *(int*)(&_token.key[0]));
    }

    SecFile file(buffer_capacity);
    if(!file.open(filepath)) {
      std::cout << "[DECODE] File open FAILED." << std::endl;
      return false;
    }

    CountOp co(&_token);
    file.run(co);

    int band_count = _height>0?_height:co.band_count();
    int width = co.width();
    if(_width>0&&_width<width) width = _width;
    size_t frame_count = 0;
    {
      auto& infos = co.infos();
      for(auto it=infos.begin(); it!=infos.end(); ++it){
        if (it->band_count==band_count) ++frame_count;
      }
    }

    if(_verbose){
      const char* type = (co.image_type()==Header::IT_SWIR?"SWIR":"VNIR");
      std::cout << "[DECODE] type = " << type << "\twidth = " << width << "\tband = " << band_count << "\tframe = " << frame_count << "/" << co.infos().size() << std::endl;
    }

    if(frame_count<1) return false;

    SaveOp so( width, band_count, frame_count, imagepath, metapath, co.infos(), &_token);

    {
      int flipping = 0;
      if(flipping == 0){
        if (co.image_type()==Header::IT_SWIR) flipping = 1;
        else flipping = -1;
      }
      so.set_flipping(flipping>0);
    }

    if(!so.valid()){
      if(_verbose){
        std::cout << "[DECODE] File created FAILED." << std::endl;
      }
    }

    if(_verbose){
      std::cout << "[DECODE] Output image to " << imagepath << std::endl;
      std::cout << "[DECODE] Output auxdata to " << metapath << std::endl;
    }

    file.reset();
    file.run(so);

    return true;
  }
};

/*
template<typename T>
bool check_segment_length(T a, T b){
  return a == b;
}

class SegOp{
 protected:
  int _length;
 public:
  SegOp(int segment_length) : _length(segment_length) {}
  virtual size_t apply(SegInfoContainer& segs){
    size_t count = 0;
    int len = _length;
    for(auto it=segs.begin(); it!=segs.end(); ++it){
      if(check_segment_length(it->second,len)){
        count += apply(count, it->first);
      }
    }
    return count;
  }
  virtual size_t apply(size_t , Segment::DataType* ){
    return 1;
  }
};

class SaveOp : public SegOp{
 protected:
  GDALDataset* _dataset;
  int _band_header_length;
 public:
  SaveOp(GDALDataset* dataset, int segment_length, int band_header_length = -1) : _dataset(dataset), SegOp(segment_length) {
    if(band_header_length>0) _band_header_length = band_header_length;
    else{
      Segment::Size cols = dataset->GetRasterXSize();
      Segment::Size bands = dataset->GetRasterCount();
      _band_header_length = segment_length/bands-cols*1.5;
    }
  }
  size_t apply(size_t id, Segment::DataType* seg_data){
    Segment::Size cols = dataset->GetRasterXSize();
    Segment::Size bands = dataset->GetRasterCount();

    std::vector<unsigned short> data(cols*bands);
    Segment::Size count = std::ceil(cols*1.5);
    Segment::Size linespace = count+_band_header_length;
    Segment::DataType* org = seg_data+_band_header_length;
    for(int i=0; i<bands; ++i){
      convert( org+i*linespace, count, &data[i*cols]);
    }

    if(_dataset->RasterIO(GF_Read, 0, id, cols, 1, data.data(), cols, 1, GDT_UInt16, bands, nullptr, sizeof(unsigned short), 0, sizeof(unsigned short)*cols) != CE_None) {}
    return 1;
  }
  virtual size_t convert(unsigned char* src, size_t count, unsigned short* dst){
    return cvt_12to16( src, count, dst);
  }
};

class Decode{
 public:
  // typedef std::vector<std::string> PathContainer;
  typedef SegFile::Keyword::DataType DataType;
  typedef SegFile::Keyword::IndexT   SegSize;
  enum KEYPOS{
    KP_DATA = 0,
    KP_META
  };
 private:
  // PathContainer _files;
  int _width;
  int _height;

  int _verbose;
  size_t _buffer_capacity;
  Keywords _keywords;
  std::vector<int> _keyword_segment_length;
 public:
  Decode(size_t buffer_size) : _buffer_capacity(buffer_size){
  }
  Decode(int width, int height, size_t buffer_size)
      : _width(width), _height(height), _buffer_capacity(buffer_size) {
  }
  void set_width(int w){
    _width = w;
  }
  void set_height(int h){
    _height = h;
  }
  int width() const { return _width; }
  int height() const { return _height; }

  template<typename KeyIterator>
  void InitSegKeys(
      KeyIterator data_first, KeyIterator data_last, int data_len,
      KeyIterator meta_first, KeyIterator meta_last, int meta_len
                   ){
    _keywords.clear();
    _keywords.emplace_back(data_first, data_last);
    _keywords.emplace_back(meta_first, meta_last);
    _keyword_segment_length.clear();
    _keyword_segment_length.push_back(data_len);
    _keyword_segment_length.push_back(meta_len);
  }
  template<typename KeyIterator>
  void AddSegKey(KeyIterator first, KeyIterator last, int len){
    _keywords.emplace_back(first, last);
    _keyword_segment_length.push_back(len);
  }
  bool apply(const char* filepath, const char* dstpath[2]){
    FILE* fp = fopen(filepath, "r");
    if(fp==nullptr) return false;

    char* data = new char[_capacity];
    size_t size = _capacity;

    std::vector<int> seg_lengths=_keyword_segment_length;

    size_t count[2] = {0, 0};
    char* pdata = data;
    while( size = fread( pdata, 1, size, fp) ){
      size = pdata+size-data;
      Segments segs_all = parse(data, size, _keywords);
      SegInfoContainer segs[2];
      segs[KP_DATA] = find(segs_all, KP_DATA);
      segs[KP_META] = find(segs_all, KP_META);
      if(seg_lengths[KP_DATA]<=0) seg_lengths[KP_DATA] = segs[KP_DATA].front().second;
      if(seg_lengths[KP_META]<=0) seg_lengths[KP_META] = segs[KP_META].front().second;
      SegOp op_data(seg_lengths[KP_DATA]), op_meta(seg_lengths[KP_META]);
      count[KP_DATA] += op_data.apply(segs[KP_DATA]);
      count[KP_META] += op_meta.apply(segs[KP_META]);
      auto id = segs_all.back().id;

      if(check_segment_length(segs_all.back().count, seg_lengths[id])){
        pdata = segs_all.back().data+segs_all.back().count;
      }else{
        pdata = segs_all.back().data-keywords[id].key.size();
      }
      size = data+size-pdata;
      if(size!=0) memcpy( data, pdata, size);
      pdata = data+size;
      size = _capacity-size;
    }

    GDALDataset* dataset[2];
    dataset[KP_DATA] = GDALCreate(dstpath[KP_DATA], _width, count[KP_DATA], _height, GDT_UInt16);
    dataset[KP_META] = GDALCreate(dstpath[KP_META], _width, count[KP_META], _height, GDT_UInt16);

    if(dataset[KP_DATA]&&dataset[KP_META]){
      pdata = data;
      size = _capacity;
      rewind(fp);
      while( size = fread( pdata, 1, size, fp) ){
        size = pdata+size-data;
        Segments segs_all = parse(data, size, _keywords);
        SegInfoContainer segs[2];
        segs[KP_DATA] = find(segs_all, KP_DATA);
        segs[KP_META] = find(segs_all, KP_META);
        SaveOp op_data(seg_lengths[KP_DATA]), op_meta(seg_lengths[KP_META]);
        op_data.apply(segs[KP_DATA]);
        op_meta.apply(segs[KP_META]);
        auto id = segs_all.back().id;
        if(check_segment_length(segs_all.back().count, seg_lengths[id])){
          pdata = segs_all.back().data+segs_all.back().count;
        }else{
          pdata = segs_all.back().data-keywords[id].key.size();
        }
        size = data+size-pdata;
        if(size!=0) memcpy( data, pdata, size);
        pdata = data+size;
        size = _capacity-size;
      }
    }

    if(dataset[KP_DATA]) GDALClose(dataset[KP_DATA]);
    if(dataset[KP_META]) GDALClose(dataset[KP_META]);

    delete[] data;
    fclose(fp);
    return true;
  }
  };
*/
