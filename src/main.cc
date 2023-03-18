#include <gdal_priv.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <set>
#include <iostream>
#include <limits>
#include <cctype>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "pathadaptor.hpp"
#include "decode/decode.h"
#include "raster.hpp"
#include "ipf.hpp"
#include "bigfileio.h"

#include "pos.h"

DEFINE_string(file, "", "image list");
DEFINE_string(o,"", "output directory");
DEFINE_string(task,"", "file for tasklist");
DEFINE_int32(t,0x00C01509, "header tag");//0x0915C000
DEFINE_int32(w, -1, "input image width");
DEFINE_int32(b, -1, "input image bands");
DEFINE_int32(buffer, 2048, "buffer size/MB");//(std::numeric_limits<int>::max)()
DEFINE_int32(splice, 1, "unite images into one");
DEFINE_string(ext, ".tif", "default extension for output decoded image");
DEFINE_int32(c, 0, "compression type: 0=LOSSLESS, 1=LOSS8, 2=LOSS4, 3=NONE");
DEFINE_string(gcp, "", "pos to add gcp to tif");
DEFINE_string(nodata, "", "specify the nodata value.");
DEFINE_int32(fc, 1, "flush cache mode: 0, 1, 2");
DEFINE_string(ot, "", "output data type {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}");
DEFINE_bool(ow, false, "whether to overwrite the destination image");
DEFINE_bool(oc, false, "Only save the window of source image.");
DEFINE_string(srcwin, "", "the window of source image");//extend nuc

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  os << "[";
  for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
    if (it != v.begin()) {
      os << ", ";
    }
    os << *it;
  }
  os << "]";
  return os;
}

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

  std::string usage("This program decodes hyper-spectrum raw data. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + "[Decode]: " + name + " <raw data file>*\n" +
            "[task]: " + name + " -task <task xml file> <image file>\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("2.3");

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
    }else if(!FLAGS_task.empty()){
      boost::property_tree::ptree tree;
      try {
        boost::property_tree::read_xml(FLAGS_task, tree, boost::property_tree::xml_parser::trim_whitespace);
      }
      catch (const boost::property_tree::ptree_error& e) {
        LOG(ERROR) << e.what();
        return 2;
      }
      boost::filesystemEx::pathadaptor adaptor(FLAGS_task);
      GDALDataset* src = (GDALDataset*)GDALOpen(path.string().c_str(), GA_ReadOnly);
      if (src == nullptr) return 1;
      std::vector<double> nodata;
      if(!FLAGS_nodata.empty()){
        std::vector<std::string> vs;
        boost::split( vs, FLAGS_nodata, boost::is_any_of(" ,;|"), boost::token_compress_on);
        for(auto& v : vs){
          nodata.push_back(std::stod(v));
        }
      }else{
        int nodata_success;
        double v = src->GetRasterBand(1)->GetNoDataValue(&nodata_success);
        if(nodata_success)
          nodata.push_back(v);
      }

      if(nodata.size()>0){
        VLOG(2) << "NoData# = " << nodata.size();
        VLOG(3) << "NoData Value = " << nodata;
      }

      int src_size[3] = {src->GetRasterXSize(), src->GetRasterYSize(), src->GetRasterCount()};
      int src_win[6] = {0, src->GetRasterXSize(), 0, src->GetRasterYSize(), 0, src->GetRasterCount()};
      if(!FLAGS_srcwin.empty()){
        std::vector<std::string> vs;
        boost::split( vs, FLAGS_srcwin, boost::is_any_of(" ,;|"), boost::token_compress_on);
        for(int i=0; i<vs.size()-1; i+=2){
          int vt = std::stoi(vs[i]);
          src_win[i] = (std::max)(vt, src_win[i]);
          if(src_win[i]>=src_win[i+1]) {
            LOG(ERROR) << "Source window size is zero";
            GDALClose(src);
            return 1;
          }
          vt = std::stoi(vs[i+1]);
          if(vt>0) src_win[i+1] = (std::min)(vt, src_win[i+1]-src_win[i]);
          else if(vt<0) src_win[i+1] += vt-src_win[i];
        }
      }
      int dst_size[3];
      int dst_win[6];
      if(FLAGS_oc){
        dst_win[0] = dst_win[2] = dst_win[4] = 0;
        dst_size[0] = dst_win[1] = src_win[1];
        dst_size[1] = dst_win[3] = src_win[3];
        dst_size[2] = dst_win[5] = src_win[5];
      }else{
        memcpy(dst_size, src_size, sizeof(src_size));
        memcpy(dst_win, src_win, sizeof(src_win));
      }
      VLOG(1) << "Source size = (" << src_size[0] << ", " << src_size[1] << ", " << src_size[2] << ")";
      VLOG_IF(1, (src_win[1]!=src_size[0])||(src_win[3]!=src_size[1])||(src_win[5]!=src_size[2]))
          << "Source win = (<" << src_win[0] << "," << src_win[0]+src_win[1]-1
          << ">, <" << src_win[2] << "," << src_win[2]+src_win[3]-1
          << ">, <" << src_win[4] << "," << src_win[4]+src_win[5]-1 <<">)";

      GDALDataset* dst = nullptr;

      xlingsky::raster::ComboOperator* ops = new xlingsky::raster::ComboOperator;
      xlingsky::raster::Processor frame( ops,GDT_Float32, FLAGS_fc);
      int store_prior[3] = {0,2,1};
      boost::filesystem::path outpath;
      int boutput = 0;
      {
        auto prior = tree.get<std::string>("HSP.dim_prior", "0,1,2");
        int boutput = tree.get<int>("HSP.output", 0);
        std::vector<std::string> result;
        boost::split(result, prior, boost::is_any_of(","));
        if (result.size() == 3) {
          store_prior[0] = std::stoi(result[0]);
          store_prior[1] = std::stoi(result[1]);
          store_prior[2] = std::stoi(result[2]);
        }
      }
      std::vector<std::string> post_tasks;

      BOOST_FOREACH(boost::property_tree::ptree::value_type & v, tree.get_child("HSP")) {
        if (v.first != "task") continue;
        auto name = v.second.get_child("<xmlattr>.name").get_value<std::string>();
        if (name == "uniform") {
          std::string a, b;
          try{
            a = v.second.get<std::string>("a");
            b = v.second.get<std::string>("b");
          }catch(const boost::property_tree::ptree_error& e){
            std::cerr << e.what() << std::endl;
            return 2;
          }
          float dst_min = v.second.get<float>("dst_min", (std::numeric_limits<float>::min)());
          float dst_max = v.second.get<float>("dst_max", (std::numeric_limits<float>::max)());
          xlingsky::raster::radiometric::NonUniformCorrection* op = new xlingsky::raster::radiometric::NonUniformCorrection(src_size[store_prior[0]],src_size[store_prior[1]], dst_min, dst_max);
          if (!op->load(adaptor.absolutepath(a).string().c_str(), adaptor.absolutepath(b).string().c_str())) {
            std::cout << "ERROR:uniform file not loaded a and b" << std::endl;
            return 1;
          }
          if(nodata.size()>0) op->AddNoDataValue(nodata.begin(), nodata.end());
          ops->Add(op);
          boutput = 1;
        }
        else if (name == "dpc") {
          std::string b = v.second.get<std::string>("file");
          xlingsky::raster::radiometric::DefectivePixelCorrection* op = nullptr;
          xlingsky::raster::radiometric::DefectivePixelCorrectionV2* p1 = new xlingsky::raster::radiometric::DefectivePixelCorrectionV2;
          if( p1->load(adaptor.absolutepath(b).string().c_str(), src_size[0], src_size[1], src_size[2])){
            op = p1;
          }else {
            delete p1;
            xlingsky::raster::radiometric::DefectivePixelCorrectionV1* p2 = new xlingsky::raster::radiometric::DefectivePixelCorrectionV1;
            int pt[3];
            {
              std::string prior = v.second.get<std::string>("dim_prior", "");
              std::vector<std::string> result;
              boost::split(result, prior, boost::is_any_of(","));
              if (result.size() == 3) {
                pt[0] = std::stoi(result[0]);
                pt[1] = std::stoi(result[1]);
                pt[2] = std::stoi(result[2]);
              }else{
                memcpy(pt, store_prior, sizeof(int)*3);
              }
            }
            if(!p2->load(adaptor.absolutepath(b).string().c_str(), src_size[pt[0]], src_size[pt[1]], pt[2])){
              delete p2;
              std::cout << "ERROR:dpc file not loaded!" << std::endl;
              return 1;
            }
            op = p2;
          }
          ops->Add(op);
          boutput = 1;
        }
        else if (name == "extend"){
          xlingsky::raster::common::Extend* op = new xlingsky::raster::common::Extend();
          std::string b = v.second.get<std::string>("file");
          if (!op->load(adaptor.absolutepath(b).string().c_str(), src_size[store_prior[0]], src_size[store_prior[1]])) {
            std::cout << "ERROR:extend file not loaded " << std::endl;
            return 1;
          }
          ops->Add(op);
          boutput = 1;
        }
        else if (name == "gauss") {
          int band = v.second.get<int>("band");
          int ksize = v.second.get<int>("ksize");
          int dim = v.second.get<int>("dim", 1);
          xlingsky::raster::filter::LinearFilter* op = new xlingsky::raster::filter::LinearFilter(ksize, band, dim);
          ops->Add(op);
          boutput = 1;
        }else if (name == "median") {
          int ksize = v.second.get<int>("ksize", 3);
          std::string method = v.second.get<std::string>("method", "adapt");
          xlingsky::raster::Operator* op = nullptr;
          if(method == "adapt"){
            float start = v.second.get<float>("hist_min", 7);
            float end = v.second.get<float>("hist_max", 4000);
            float step = v.second.get<float>("hist_step", 1);
            xlingsky::raster::enhancement::Despeckle* p = new xlingsky::raster::enhancement::Despeckle(ksize, start, end, step);
            op = p;
          }else{
            xlingsky::raster::filter::MedianBlur* p = new xlingsky::raster::filter::MedianBlur(ksize);
            op = p;
          }
          ops->Add(op);
          boutput = 1;
        }else if(name == "dark"){
          std::string b = v.second.get<std::string>("b");
          std::string a = v.second.get<std::string>("a","");
          std::string index = v.second.get<std::string>("index","");
          float dst_min = v.second.get<float>("dst_min", (std::numeric_limits<float>::min)());
          float dst_max = v.second.get<float>("dst_max", (std::numeric_limits<float>::max)());
          xlingsky::raster::radiometric::PixelCorrection* op = nullptr;
          if(a.empty() || index.empty()){
            xlingsky::raster::radiometric::DarkBackgroundCorrection* p = new xlingsky::raster::radiometric::DarkBackgroundCorrection(src_size[store_prior[0]],src_size[store_prior[1]], dst_min, dst_max);
            if (!p->load(adaptor.absolutepath(b).string().c_str())) {
              std::cout << "ERROR:dark file not loaded b" << std::endl;
              return 1;
            }
            op = p;
          }else{
            xlingsky::raster::radiometric::DarkBackgroundLinear* p = new xlingsky::raster::radiometric::DarkBackgroundLinear(src_size[store_prior[0]],src_size[store_prior[1]],src_size[store_prior[2]], dst_min, dst_max);
            if (!p->load(adaptor.absolutepath(a).string().c_str(),adaptor.absolutepath(b).string().c_str(),adaptor.absolutepath(index).string().c_str())) {
              std::cout << "ERROR:linear dark file not loaded a,b,index" << std::endl;
              return 1;
            }
            op = p;
          }
          if(nodata.size()>0) op->AddNoDataValue(nodata.begin(), nodata.end());
          ops->Add(op);
          boutput = 1;
        }else if(name == "statistic"){
          std::string method = v.second.get<std::string>("method", "mean");
          boost::filesystem::path dstpath;
          if(FLAGS_o.empty()){
            dstpath = path;
          }else{
            dstpath = FLAGS_o;
            if(boost::filesystem::is_directory(dstpath)){
              dstpath /= path.filename();
            }
          }
          dstpath.replace_extension();
          dstpath += "_"+method+".txt";
          xlingsky::raster::FrameIterator* op = nullptr;
          if(method=="mean"){
            float cut_ratio_lower = v.second.get<float>("cut_lower", 0);
            float cut_ratio_upper = v.second.get<float>("cut_upper", 0);
            xlingsky::raster::radiometric::MeanStdCalculator* p = new xlingsky::raster::radiometric::MeanStdCalculator(src_size[store_prior[1]], cut_ratio_lower, cut_ratio_upper);
            bool nuc = v.second.get<bool>("nuc", false);
            if(nuc){
              boost::filesystem::path xmlpath = dstpath;
              xmlpath.replace_extension(".xml");
              p->SetXmlPath(xmlpath.string().c_str());
              p->SetDimOrder(store_prior);
              bool apply = v.second.get<bool>("nuc.<xmlattr>.apply", true);
              if(apply)
              {
                std::string cmd(argv[0]);
                cmd = cmd + " -task " + xmlpath.string() + " " + path.string();
                if(!FLAGS_o.empty())
                  cmd += " -o " + FLAGS_o;
                post_tasks.push_back(cmd);
              }
            }
            p->SetFilePath(dstpath.string().c_str());
            op = p;

          } else{
            xlingsky::raster::radiometric::MedianCalculator* p = new xlingsky::raster::radiometric::MedianCalculator(src_size[store_prior[1]]);
            p->SetFilePath(dstpath.string().c_str());
            op = p;
          }
          if(nodata.size()>0) op->AddNoDataValue(nodata.begin(), nodata.end());
          ops->Add(op);
        }else if(name == "sort"){
          xlingsky::raster::common::Sort* op = new xlingsky::raster::common::Sort();
          ops->Add(op);
          boutput = 1;
        }
        else if (name == "nuc"){
          float cut_ratio_dark = v.second.get<float>("cut_dark", 0.03);
          float cut_ratio_bright = v.second.get<float>("cut_bright", 0.1);
          float ratio_threshold_dark = v.second.get<float>("threshold_dark", 0.03);
          float ratio_threshold_bright = v.second.get<float>("threshold_bright", 0.03);
          int sample_num_dark = v.second.get<int>("sample_dark", 40);
          int sample_num_bright = v.second.get<int>("sample_bright", 40);
          int tile_size = v.second.get<int>("tile_size", (std::numeric_limits<int>::max)());
          int tile_overlap = v.second.get<int>("tile_overlap", -1);
          bool preferred_a = v.second.get<bool>("a", true);
          bool debug = v.second.get<bool>("debug", true);
          std::string t = v.second.get<std::string>("interp", "pchip");
          xlingsky::raster::radiometric::NucCalculator::InterpType type;
          if(t=="makima")
            type = xlingsky::raster::radiometric::NucCalculator::MAKIMA;
          else if(t=="barycentric")
            type = xlingsky::raster::radiometric::NucCalculator::BARYCENTRIC;
          else
            type = xlingsky::raster::radiometric::NucCalculator::PCHIP;
          if(tile_overlap<0) tile_overlap = tile_size/2;
          xlingsky::raster::radiometric::NucCalculator* op = new xlingsky::raster::radiometric::NucCalculator(src_size[store_prior[1]], cut_ratio_dark, cut_ratio_bright, ratio_threshold_dark, ratio_threshold_bright, sample_num_dark, sample_num_bright, tile_size, tile_overlap, type, preferred_a?xlingsky::raster::radiometric::NucCalculator::SCALE:xlingsky::raster::radiometric::NucCalculator::OFFSET);
          if(nodata.size()>0) op->AddNoDataValue(nodata.begin(), nodata.end());
          char apath[512], bpath[512], bppath[512], xmlpath[512], hipath[512], lopath[512];
          boost::filesystem::path dstpath;
          if(FLAGS_o.empty()){
            dstpath = path;
          }else{
            dstpath = FLAGS_o;
            if(boost::filesystem::is_directory(dstpath)){
              dstpath /= path.filename();
            }
          }
          dstpath.replace_extension();
          dstpath += "_";
          {
            boost::filesystem::path t;
            t = dstpath;
            t += "a.txt";
            strcpy(apath, t.string().c_str());
            t = dstpath;
            t += "b.txt";
            strcpy(bpath, t.string().c_str());
            t = dstpath;
            t += "badpixels.txt";
            strcpy(bppath, t.string().c_str());
            t = dstpath;
            t += "uniform.xml";
            strcpy(xmlpath, t.string().c_str());
            if (debug){
              t = dstpath;
              t += "hi.txt";
              strcpy(hipath, t.string().c_str());
              t = dstpath;
              t += "lo.txt";
              strcpy(lopath, t.string().c_str());
            }else{
              hipath[0] = lopath[0] = 0;
            }
          }
          op->SetFilePath(apath, bpath, bppath, xmlpath, hipath, lopath);
          op->SetDimOrder(store_prior);
          ops->Add(op);

          bool apply = v.second.get<bool>("apply", true);
          if(apply && xmlpath[0])
          {
            std::string cmd(argv[0]);
            cmd = cmd + " -task " + xmlpath + " " + path.string();
            if(!FLAGS_o.empty())
              cmd += " -o " + FLAGS_o;
            post_tasks.push_back(cmd);
          }
        }
        else if (name == "interp") {
          std::vector<double> wl_old, wl_new;
          xlingsky::raster::spectrum::Interpolator::InterpType type;
          int num_lines = src->GetRasterXSize(), num_samples_old = src->GetRasterCount(), num_samples_new;
          {
            std::string w0 = v.second.get<std::string>("wl0");
            std::string w1 = v.second.get<std::string>("wl1");

            std::vector<std::string> result;
            boost::split(result, w0, boost::is_any_of(", ;"));
            for (auto& r : result)
              wl_old.push_back(std::stod(r));
            boost::split(result, w1, boost::is_any_of(", ;"));
            for (auto& r : result)
              wl_new.push_back(std::stod(r));

            bool flag = false;
            if(wl_old.size()==num_samples_old){
              num_samples_new = wl_new.size();
              flag = true;
            }else if(wl_old.size()==(size_t)num_samples_old*num_lines){
              num_samples_new = wl_new.size()/num_lines;
              if(num_samples_new>0) flag = true;
            }
            if (!flag) {
              std::cout << "ERROR: interp wl length MISMATCH" << std::endl;
              GDALClose(src);
              return 1;
            }

            std::string t = v.second.get<std::string>("type", "pchip");
            if (t == "spline_cubic")
              type = xlingsky::raster::spectrum::Interpolator::BSPLINE_CUBIC;
            else if(t=="spline_quintic")
              type = xlingsky::raster::spectrum::Interpolator::BSPLINE_QUINTIC;
            else if(t=="spline_quadratic")
              type = xlingsky::raster::spectrum::Interpolator::BSPLINE_QUADRATIC;
            else if(t=="makima")
              type = xlingsky::raster::spectrum::Interpolator::MAKIMA;
            else if(t=="barycentric")
              type = xlingsky::raster::spectrum::Interpolator::BARYCENTRIC;
            else
              type = xlingsky::raster::spectrum::Interpolator::PCHIP;
          }
          dst_size[2] = dst_win[5] = num_samples_new;
          xlingsky::raster::spectrum::Interpolator* op = new xlingsky::raster::spectrum::Interpolator(num_lines, wl_old, num_samples_old, wl_new, num_samples_new);
          op->SetInterpType(type);
          ops->Add(op);
          boutput = 1;
        } else if (name == "destripe") {
          int tile_size = v.second.get<int>("tile_size", 31);
          xlingsky::raster::enhancement::Destripe* op =
              new xlingsky::raster::enhancement::Destripe(tile_size);
          ops->Add(op);
          boutput = 1;
        }else if(name == "render"){
          float src_min = v.second.get<float>("src_min", 7);
          float src_umax = v.second.get<float>("src_umax", 4000);
          float dst_min = v.second.get<float>("dst_min", 0);
          float dst_max = v.second.get<float>("dst_max", 255);
          if(FLAGS_ot.empty()){
            if(dst_max<(int)(std::numeric_limits<unsigned char>::max)()+1)
              FLAGS_ot = "byte";
          }
          std::string m = v.second.get<std::string>("mode", "WALLIS|MINMAX|CUT|TILE");
          std::transform(m.begin(), m.end(), m.begin(),
                         [](unsigned char c){ return std::tolower(c); });
          int mode = 0;
          if(m.find("minmax")!=std::string::npos)
            mode |= 0x01;
          if(m.find("cut")!=std::string::npos)
            mode |= 0x02;
          if(m.find("hist")!=std::string::npos)
            mode |= 0x04;
          if(m.find("wallis")!=std::string::npos)
            mode |= 0x08;
          if(m.find("tile")!=std::string::npos)
            mode |= 0x40;
          if(m.find("global")!=std::string::npos)
            mode |= 0x80;

          xlingsky::raster::LutCreator* lut = nullptr;
          if ((mode & 0x02) || (mode & 0x04) || (mode & 0x08)) {
            float src_step = v.second.get<float>("src_step", 1);
            if(mode&0x08){
              float dst_mean = v.second.get<float>("dst_mean", 127);
              float dst_std = v.second.get<float>("dst_std", 70);//40-70
              float c = v.second.get<float>("c", 0.8);//0-1
              float b = v.second.get<float>("b", 0.9);//0-1
              xlingsky::raster::LutWallis* l = new xlingsky::raster::LutWallis(src_step, src_min, src_umax, dst_min, dst_max);
              l->Setup(dst_mean, dst_std, c, b);
              lut = l;
            }else{
              float cut_ratio_lower = v.second.get<float>("cut_lower", 0.002);
              float cut_ratio_upper = v.second.get<float>("cut_upper", 0.002);
              if(mode&0x02){
                xlingsky::raster::LutLinear* l = new xlingsky::raster::LutLinear(src_step, src_min, src_umax, dst_min, dst_max);
                l->set_cut_ratio(cut_ratio_lower, cut_ratio_upper);
                lut = l;
              }else{
                float hist_clip = v.second.get<float>("hist_clip", 10);
                xlingsky::raster::LutClahe* l = new xlingsky::raster::LutClahe(hist_clip,src_step, src_min, src_umax, dst_min, dst_max);
                l->set_cut_ratio(cut_ratio_lower, cut_ratio_upper);
                lut = l;
              }
            }
          }else{
            xlingsky::raster::LutCreator* l = new xlingsky::raster::LutCreator(src_min, src_umax, dst_min, dst_max);
            lut = l;
          }

          int resample_col_step = v.second.get<float>("resample_col_step", 3);
          int resample_row_step = v.second.get<float>("resample_row_step", 3);
          if(mode&0x40){
            unsigned int tile_cols = v.second.get<unsigned int>("tile_cols", 0);
            unsigned int tile_rows = v.second.get<unsigned int>("tile_rows", 0);
            unsigned int grid_x = v.second.get<unsigned int>("grid_x", 8);
            unsigned int grid_y = v.second.get<unsigned int>("grid_y", 8);
            if (tile_cols == 0) {
              if (grid_x == 0)
                tile_cols = 512;
              else {
                tile_cols = src_size[store_prior[0]] / grid_x;
                if (tile_cols < 1) tile_cols = 1;
              }
            }
            if (tile_rows == 0) {
              if (grid_y == 0)
                tile_rows = 512;
              else {
                tile_rows = src_size[store_prior[1]] / grid_y;
                if (tile_rows < 1) tile_rows = 1;
              }
            }
            xlingsky::raster::BiTileLut* op =
                new xlingsky::raster::BiTileLut( tile_cols, tile_rows, (mode&0x80)?xlingsky::raster::BiTileLut::GLOBAL:xlingsky::raster::BiTileLut::NONE);
            op->set_lut_creator(lut);
            op->set_resample_interval(resample_col_step, resample_row_step);
            ops->Add(op);
          }else{
            xlingsky::raster::enhancement::Render* op =
                new xlingsky::raster::enhancement::Render((mode&0x80)?xlingsky::raster::enhancement::Render::GLOBAL:xlingsky::raster::enhancement::Render::NONE);
            op->set_lut_creator(lut);
            op->set_resample_interval(resample_col_step, resample_row_step);
            ops->Add(op);
          }
          boutput = 1;
        }
      }

      VLOG(1) << "Total number of tasks is " << ops->size();
      if(ops->size()<1){
        delete ops;
        GDALClose(src);
        return 1;
      }

      if (boutput) {
        if (!FLAGS_o.empty()) outpath = FLAGS_o;
        else {
          outpath = path;
          outpath.replace_extension();
          outpath += "_mod";
          outpath += path.extension();
        }
      }
      if (!outpath.empty()) {
        GDALDataType datatype =src->GetRasterBand(1)->GetRasterDataType();
        {
          if(!FLAGS_ot.empty()){
            const char* tag_name[] = {"byte", "uint16", "int16", "uint32", "int32", "float32", "float64", "cint16", "cint32", "cfloat32", "cfloat64"};
            GDALDataType tag_type[] = {GDT_Byte, GDT_UInt16, GDT_Int16, GDT_UInt32, GDT_Int32, GDT_Float32, GDT_Float64, GDT_CInt16, GDT_CInt32, GDT_CFloat32, GDT_CFloat64};
            std::transform(FLAGS_ot.begin(), FLAGS_ot.end(), FLAGS_ot.begin(),
                           [](unsigned char c){ return std::tolower(c); });
            for(int i=0; i<sizeof(tag_name)/sizeof(tag_name[0]); ++i)
              if(FLAGS_ot==tag_name[i]) { datatype = tag_type[i]; break; }
          }
        }
        VLOG(1) << "Output path is " << outpath.string();
        VLOG(1) << "Target size = (" << dst_size[0] << ", " << dst_size[1] << ", " << dst_size[2] << ")";
        if(boost::filesystem::exists(outpath) && !FLAGS_ow){
          dst = (GDALDataset*)GDALOpen(outpath.string().c_str(), GA_Update);
          if(dst==nullptr || dst->GetRasterXSize()!=dst_size[0] || dst->GetRasterYSize()!=dst_size[1] || dst->GetRasterCount()!=dst_size[2]){
            LOG(ERROR) << "The existing target size doesn't match the estimation";
            delete ops;
            GDALClose(src);
            GDALClose(dst);
            return 3;
          }
        }else
          dst = GDALCreate(outpath.string().c_str(), dst_size[0], dst_size[1], dst_size[2], datatype);
        if (dst == nullptr) {
          delete ops;
          GDALClose(src);
          return 1;
        };
      }
      frame.SetStoreOrder(store_prior);
      frame.SetSource(src, src_win);
      if (dst) frame.SetDestination(dst, dst_win);

      int buffer_size[3] = {dst_win[1], dst_win[3], dst_win[5]};
      {
        size_t max_buffer_size = (size_t)FLAGS_buffer<<20;
        size_t l = max_buffer_size / ((size_t)frame.GetDataTypeSize() * buffer_size[store_prior[0]] * buffer_size[store_prior[1]]);
        if (l == 0) l = 1;
        if (buffer_size[store_prior[2]]>l)
          buffer_size[store_prior[2]] = l;
      }
      frame.ReserveBufferSize((size_t)buffer_size[0] * buffer_size[1] * buffer_size[2]);
      if (src_win[5] < buffer_size[2])  buffer_size[2] = src_win[5];

      xlingsky::ipf::TileProcessing( frame.source().win+3, buffer_size, &frame);

      delete ops;
      if (dst) {
        double aop6[6];
        if(src->GetGeoTransform(aop6)==CE_None)
          dst->SetGeoTransform(aop6);
        dst->SetProjection(src->GetProjectionRef());
        GDALClose(dst);
      }
      GDALClose(src);

      frame.ReserveBufferSize(0);
      if(post_tasks.size() > 0){
        for(auto& cmd : post_tasks){
          VLOG(1) << "CALL: " << cmd;
          system(cmd.c_str());
        }
      }
      return 0;
    } else if (!FLAGS_gcp.empty()) {
      HSP::Pos pos;
      pos.load(FLAGS_gcp.c_str());
      HSP::PinholeCamera cam;
      cam.SetFocalLength(59750/30.0);
      cam.SetPrincipalPoint(-724,0);
      double a[3] = {0,0,0};
      cam.SetAngles(a);
      double t[3] = {0.006, 0.25768, -0.6773 };
      cam.SetTranslation(t);
      HSP::LinescanModel model;
      model.SetCamera(&cam);
      model.SetPos(&pos);
      boost::filesystem::path rpcpath(path);
      rpcpath.replace_extension(".rpc");
      model.test();
      double range_samp[] = {10, 1300};
      double range_line[] = {8, (double)pos.data().rbegin()->first-8};
      double range_height[] = {1500, 2500};
      model.GenerateRPC(range_samp, range_line, range_height, path.string().c_str());
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
