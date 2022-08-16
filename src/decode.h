#ifndef HSP_DECODE_H
#define HSP_DECODE_H

#include <vector>

#include "tokenizer.hpp"

namespace HSP{

inline unsigned short cvt_12to16_even(unsigned char a1, unsigned char a2){
  return (a1<<4)+((a2&0xF0)>>4);
}
inline unsigned short cvt_12to16_odd(unsigned char a1, unsigned char a2){
  return ((a1&0x0F)<<8)+a2;
}

template<typename DstIterator>
size_t cvt_12to16(const void* src, size_t count, DstIterator dst){
  const unsigned char* data = (const unsigned char*)src;
  size_t i = 1;
  auto it = dst;
  while(i<count){
    *it = cvt_12to16_even(data[i-1], data[i]);
    ++i; ++it;
    if(i>=count) break;
    *it = cvt_12to16_odd(data[i-1], data[i]);
    i+=2; ++it;
  }
  return std::distance(dst, it);
}

struct Header{
  typedef std::vector<unsigned short>::iterator DstIterator;
  typedef std::vector<unsigned short>::reverse_iterator DstIteratorR;
  typedef size_t (*Unzip)(const void*, size_t, DstIterator);
  typedef size_t (*UnzipR)(const void*, size_t, DstIteratorR);
  enum DataType{
    DT_IMAGE = 0,
    DT_META
  };
  enum ImageType{
    IT_SWIR = 0,
    IT_VNIR
  };
  enum CompressionType{
    CT_LOSSLESS = 0,
    CT_LOSS8 = 1,
    CT_LOSS4 = 2,
    CT_NONE = 3
  };
  int data_count;
  DataType data_type;
  ImageType image_type;
  static CompressionType compression;
  static float image_data_ratio;
  static float meta_data_ratio;
  static Unzip  unzip_front;
  static UnzipR unzip_back;

  int band_id;
  int frame_id;

  static void SetCompression(CompressionType);
};

class Decode{
 private:
  int _width;
  int _height;

  int _verbose;
  tokenizer::Token _token;
 public:
  Decode(int width, int height, int verbose) : _width(width), _height(height), _verbose(verbose) {}
  ~Decode() {}
  template<typename Iterator>
  void add_token(Iterator first, Iterator last){
    _token = tokenizer::Token(first, last);
  }
  bool apply(const char* filepath, size_t buffer_capacity, const char* imagepath, const char* metapath);
};

bool Aux2Pos(const char* auxpath, const char* pospath);

};

#endif
