#ifndef GDALEX_HPP
#define GDALEX_HPP

#include <gdal_priv.h>

inline char* strlwr( char *str ){
  char *orig = str; char d = 'a'-'A';
  for( ;*str!='\0';str++ ){ if ( *str>='A' && *str<='Z' ) *str = *str+d; }
  return orig;
}

inline const char* GetGDALDescription(const char* lpstrExt, const char* default_ret)
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
  return default_ret;
}

inline bool IsRasterDataset(const char* filepath){
  FILE* fp = fopen(filepath,"r");
  if (fp==nullptr) return false;
  fclose(fp);
  char path[512];
  strcpy(path,filepath);
  char* p = strrchr(path, '.');
  if( p && GetGDALDescription(p+1, nullptr) ) return true;
  if(p==nullptr) p = path+strlen(path);
  strcpy(p,".hdr");
  fp = fopen(path, "r");
  if(fp==nullptr) return false;
  fclose(fp);
  return true;
}

inline std::pair<double, double> GetDataTypeMinMax(GDALDataType type){
  switch(type){
    case GDT_Byte:
      return std::make_pair(0, 255);
    case GDT_Int16:
      return std::make_pair(std::numeric_limits<short>::min(), std::numeric_limits<short>::max());
    case GDT_UInt16:
      return std::make_pair(std::numeric_limits<unsigned short>::min(), std::numeric_limits<unsigned short>::max());
    case GDT_Float64:
    default:
      return std::make_pair(std::numeric_limits<float>::min(), std::numeric_limits<float>::max());
  }
}

inline GDALDataset* GDALCreate(const char* filepath, int cols, int rows, int bands, GDALDataType type){
  GDALDataset* dataset;
  const char* ext = strrchr(filepath,'.');
  if(ext==nullptr) ext = "tif";
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(GetGDALDescription(ext,"ENVI"));
  dataset = poDriver->Create(filepath, cols, rows, bands, type, nullptr);
  return dataset;
}


#endif
