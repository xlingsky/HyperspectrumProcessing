#ifndef GDALEX_HPP
#define GDALEX_HPP

#include <gdal_priv.h>

inline char* strlwr( char *str ){
  char *orig = str; char d = 'a'-'A';
  for( ;*str!='\0';str++ ){ if ( *str>='A' && *str<='Z' ) *str = *str+d; }
  return orig;
}
inline const char* GetGDALDescription(const char* lpstrExt)
{
  static const char* strExt[] = { "tif", "tiff", "bmp", "jpg", "png", "img", "bt", "ecw", "fits", "gif", "hdf", "hdr","pix", "dat"};
  static const char* strDescrip[] = { "GTiff","GTiff", "BMP", "JPEG","PNG", "HFA", "BT", "ECW", "FITS", "GIF", "HDF4", "EHdr","PCIDSK", "ENVI" };

  char suffix[10];	if (*lpstrExt == '.') strcpy(suffix, lpstrExt + 1); else strcpy(suffix, lpstrExt);
  strlwr(suffix);

  for (unsigned long i = 0; i < sizeof(strExt) / sizeof(const char*);i++){
    if (!strcmp(strExt[i], suffix)) {
      return strDescrip[i];
    }
  }
  return nullptr;
}

inline GDALDataset* GDALCreate(const char* filepath, int cols, int rows, int bands, GDALDataType type){
  GDALDataset* dataset;
  const char* ext = strrchr(filepath,'.');
  if(ext==nullptr) ext = "tif";
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(GetGDALDescription(ext));
  dataset = poDriver->Create(filepath, cols, rows, bands, type, nullptr);
  return dataset;
}


#endif
