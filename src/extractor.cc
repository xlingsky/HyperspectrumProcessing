#include <gdal_priv.h>
#include <iostream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <boost/algorithm/string.hpp>
#include <set>

#include "util/pathadaptor.hpp"
#include "util/bigfileio.h"
#include "decode/decode.h"
#include "raster/gdal_traits.hpp"
#include "raster/gdalex.hpp"
#include "geometry/pos.h"

DEFINE_string(o,"", "output directory");
DEFINE_int32(t,0x00C01509, "header tag");//0x0915C000
DEFINE_int32(w, -1, "input image width");
DEFINE_int32(b, -1, "input image bands");
DEFINE_int32(splice, 1, "unite images into one");
DEFINE_int32(c, 0, "compression type: 0=LOSSLESS, 1=LOSS8, 2=LOSS4, 3=NONE");
DEFINE_string(ext, ".tif", "default extension for output decoded image");
DEFINE_int32(buffer, 2048, "buffer size/MB");//(std::numeric_limits<int>::max)()

bool IsRaw(boost::filesystem::path& file){
  std::string ext;
  ext = boost::algorithm::to_lower_copy( file.extension().string());
  if(ext==".dat"){
    boost::filesystem::path t(file);
    t.replace_extension(".hdr");
    if(boost::filesystem::exists(t))
      return false;
    else
      return true;
  }
  return false;
}

void RemoveNewline(char* str){
  const char r = 0;
  if(*str=='\n') *str=r;
  while(*str){
    if(*str=='\n'){
      if(*(str-1)=='\r') *(str-1) = r;
      *str = r;
    }
    ++str;
  }
}

struct SegFrm{
  int id;
  int st;
  int ed;
  friend bool operator < (const SegFrm& s1, const SegFrm& s2){
    return s1.st<s2.st;
  }
};

std::vector< std::vector<std::string> > Group(std::vector<std::string>& list){
  std::vector< std::vector<std::string> > grps(1);
  std::set<SegFrm> segs;
  for(int i=0; i<list.size(); ++i ){
    boost::filesystem::path t(list[i]);
    t.replace_extension(".aux");
    if(boost::filesystem::exists(t)){
      HSP::Pos pos;
      pos.load(t.string().c_str());
      if(pos.size()==0)
        grps[0].push_back(list[i]);
      else{
        SegFrm seg;
        seg.id = i;
        seg.st = pos.data().begin()->second.frame;
        seg.ed = pos.data().rbegin()->second.frame;
        segs.insert(seg);
      }
    }
  }

  if(grps[0].size() != list.size()) {
    std::vector<const SegFrm*> segrp;
    for (auto it=segs.begin(); it!=segs.end(); ++it) {
      auto& seg = *it;
      auto it_grp = segrp.rbegin();
      while(it_grp!=segrp.rend()){
        if((*it_grp)->ed<seg.st) break;
        ++it_grp;
      }
      if(it_grp!=segrp.rend()){
        *it_grp = &seg;
        (grps.rbegin()+std::distance(segrp.rbegin(), it_grp))->push_back(list[seg.id]);
      }else{
        segrp.push_back(&seg);
        std::vector<std::string> t;
        t.push_back(list[seg.id]);
        grps.push_back(t);
      }
    }
  }

  return grps;
}

std::vector<char> ReadFileToVector(const char* filepath){
  std::vector<char> ret;
  // FILE* fp = fopen(filepath, "r");
  // if(fp) {
  //   fseek(fp, 0, SEEK_END);
  //   ret.resize(ftell(fp));
  //   if(ret.size())
  //   {
  //     rewind(fp);
  //     fread(ret.data(), sizeof(char), ret.size(), fp);
  //   }
  //   fclose(fp);
  // }
  HANDLE fp = CreateFileE(filepath, GENERIC_READ);
  if(fp&&fp!=INVALID_HANDLE_VALUE) {
    ret.resize(SetFilePointer(fp, 0, nullptr, FILE_END));
    SetFilePointer(fp, 0, nullptr, FILE_BEGIN);
    DWORD rw;
    ReadFile(fp, ret.data(), ret.size(), &rw, nullptr);
    CloseHandle(fp);
  }
  return ret;
}

template<typename PathIterator, typename T = unsigned short>
bool Splice(PathIterator first, PathIterator last, const char* dstpath){
  bool aux = true;
  int width = 0, height = 0, band = 0;
  for(auto it=first; it!=last; ++it){
    boost::filesystem::path file(*it);
    GDALDataset* src = (GDALDataset*)GDALOpen(file.string().c_str(), GA_ReadOnly);
    if(src==nullptr) continue;
    if(width<src->GetRasterXSize()) width = src->GetRasterXSize();
    if(band<src->GetRasterCount()) band = src->GetRasterCount();
    height += src->GetRasterYSize();
    GDALClose(src);
    file.replace_extension(".aux");
    if(!boost::filesystem::exists(file)) aux = false;
  }

  if(width==0||height==0||band==0) return false;

  GDALDataType type = gdal::DataType<T>::type();
  GDALDataset* dst = GDALCreate( dstpath, width, height, band, type);
  if(dst==nullptr) return false;

  HANDLE fp = 0;
  boost::filesystem::path auxfile;
  if(aux)
  {
    auxfile = dstpath;
    auxfile.replace_extension(".aux");
    fp = CreateFileE(auxfile.string().c_str() , GENERIC_WRITE);
    if(fp==INVALID_HANDLE_VALUE) fp = 0;
  }

  std::vector<T> data;
  int rowid = 0;
  for(auto it=first; it!=last; ++it){
    boost::filesystem::path file(*it);
    GDALDataset* src = (GDALDataset*)GDALOpen(file.string().c_str(), GA_ReadOnly);
    if(src==nullptr) continue;
    data.resize((size_t)src->GetRasterXSize()*src->GetRasterYSize()*src->GetRasterCount());
    if(src->RasterIO(GF_Read, 0, 0, src->GetRasterXSize(), src->GetRasterYSize(), data.data(), src->GetRasterXSize(), src->GetRasterYSize(), type, src->GetRasterCount(), nullptr, 0, 0, 0)==CE_None && dst->RasterIO(GF_Write, 0, rowid, src->GetRasterXSize(), src->GetRasterYSize(), data.data(), src->GetRasterXSize(), src->GetRasterYSize(), type, src->GetRasterCount(), nullptr, 0, 0, 0) == CE_None){
    }
    rowid += src->GetRasterYSize();
    GDALClose(src);

    if(fp){
      file.replace_extension(".aux");
      auto str = ReadFileToVector(file.string().c_str());
      if(str.size()){
        DWORD rw;
        WriteFile(fp, str.data(), str.size(), &rw, nullptr);
      }
    }
  }

  GDALClose(dst);
  if(fp) {
    CloseHandle(fp);
    boost::filesystem::path pos(auxfile);
    pos.replace_extension();
    pos += ".pos";
    HSP::Aux2Pos(auxfile.string().c_str(), pos.string().c_str());
  }
  return true;
}

int main(int argc, char* argv[]){

  FLAGS_logtostderr = 1;

  GDALAllRegister();
  CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

  std::string usage("This program extracts hyper-spectrum raw data. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + name + " <raw data file>*\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("1.0");

  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if(argc<2){
    gflags::ShowUsageWithFlagsRestrict(argv[0], "main");
    return 1;
  }

  int verbose = FLAGS_v;

  if(argc==2) {
    boost::filesystem::path path(argv[1]);
    if(path.extension() == ".aux"){
      auto pospath = path;
      pospath.replace_extension(".pos");
      HSP::Aux2Pos(path.string().c_str(), pospath.string().c_str());
      return 0;
    }
  }

  boost::filesystem::path tempdir;

  std::vector<std::string> rawlist, imagelist;

  for(int i=1; i<argc; ++i) {
    boost::filesystem::path file(argv[i]);
    if(!boost::filesystem::exists(file)) continue;
    if(file.extension()==".txt"){
      FILE* fp = fopen(file.string().c_str(), "r");
      if(fp){
        boost::filesystemEx::pathadaptor adaptor(file);
        char strline[512];
        while(fgets(strline, 512, fp)){
          RemoveNewline(strline);
          if(strlen(strline)>1){
            file = adaptor.absolutepath(strline);
            if(!file.empty()){
              if(IsRaw(file))
                rawlist.push_back(file.string());
              else
                imagelist.push_back(file.string());
            }
          }
        }
        fclose(fp);
      }
    }else if(IsRaw(file))
      rawlist.push_back(file.string());
    else
      imagelist.push_back(file.string());
  }

  if(rawlist.size()>0){
    char* token = (char*)&FLAGS_t;
    int token_length = 4;

    HSP::Header::SetCompression((HSP::Header::CompressionType)FLAGS_c);

    HSP::Decode decode( FLAGS_w, FLAGS_b, FLAGS_v);
    decode.add_token(token, token+token_length);

    if(verbose){
      std::cout << "[DECODE] Total to-decoded raw data number= " << rawlist.size() << std::endl;
    }

    boost::filesystem::path dirpath;
    if(FLAGS_o.empty()){
    } else if(boost::filesystem::is_directory(FLAGS_o)) {
      dirpath = FLAGS_o;
    }else if(boost::filesystem::path(FLAGS_o).extension().empty()){
      boost::filesystem::create_directories(FLAGS_o);
      dirpath = FLAGS_o;
    }else if(rawlist.size()>1){
      if(FLAGS_splice){
        tempdir = FLAGS_o;
        tempdir.replace_extension();
        tempdir += "_temp";
        boost::filesystem::create_directories(tempdir);
        dirpath = tempdir;
      }else{
        dirpath = boost::filesystem::path(FLAGS_o).parent_path();
      }
    }else
      dirpath = FLAGS_o;

    for(auto it=rawlist.begin(); it!=rawlist.end(); ++it){
      boost::filesystem::path srcpath(*it), dstpath;
      if(dirpath.empty()){
        dstpath = srcpath;
        dstpath.replace_extension();
        dstpath += "_decode";
      }else if(boost::filesystem::is_directory(dirpath)){
        dstpath = dirpath;
        dstpath /= srcpath.filename();
        dstpath.replace_extension();
      }else {
        dstpath = dirpath;
        FLAGS_ext = dstpath.extension().string();
        dstpath.replace_extension();
      }

      boost::filesystem::path image(dstpath);
      image += FLAGS_ext;
      boost::filesystem::path meta(dstpath);
      meta  += ".aux";

      decode.apply(srcpath.string().c_str(), FLAGS_buffer, image.string().c_str(), meta.string().c_str());
      {
        boost::filesystem::path pos(dstpath);
        pos += ".pos";
        HSP::Aux2Pos(meta.string().c_str(), pos.string().c_str());
      }
      imagelist.push_back(image.string());
    }
  }

  if(FLAGS_splice > 0 && imagelist.size() >1){
    if(verbose){
      std::cout << "[SPLICE] Total image number= " << imagelist.size() << std::endl;
    }
    std::vector< std::vector<std::string> > grps;
    if(FLAGS_splice == 1)
      grps.push_back(imagelist);
    else
      grps = Group(imagelist);
    if(verbose){
      std::cout << "[SPLICE] Group number= " << grps.size() << std::endl;
    }

    boost::filesystem::path dirpath(FLAGS_o);

    for(int i=0; i<grps.size(); ++i){
      auto& g = grps[i];
      if(verbose){
        std::cout << " [SPLICE] group["<< i << "]  number= " << g.size() << std::endl;
      }
      if(g.size()==0) continue;
      boost::filesystem::path srcpath(g[0]), dstpath;
      if(dirpath.empty()){
        dstpath = srcpath;
        dstpath.replace_extension();
        dstpath += "_splice";
        dstpath += srcpath.extension();
      }else if(boost::filesystem::is_directory(dirpath)){
        dstpath = dirpath;
        char name[512];
        snprintf(name,512,"%s_splice_%d%s", srcpath.stem().string().c_str(), i, srcpath.extension().string().c_str() );
        dstpath /= name;
      }else {
        dstpath = dirpath;
        if(grps.size()>1){
          char post[512];
          snprintf(post,512,"_%d%s", i, dstpath.extension().string().c_str());
          dstpath.replace_extension();
          dstpath /= post;
        }
      }
      Splice(g.begin(), g.end(), dstpath.string().c_str());
    }
  }

  if(boost::filesystem::exists(tempdir))
    boost::filesystem::remove_all(tempdir);

  return 0;
}
