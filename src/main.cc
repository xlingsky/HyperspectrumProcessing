#include <gdal_priv.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <set>
#include <iostream>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "pathadaptor.hpp"
#include "gdalex.hpp"
#include "decode.h"
#include "radiometric_correction.hpp"
#include "pos.h"
#include "bigfileio.h"

DEFINE_string(file, "", "image list");
DEFINE_string(o,"", "output directory");
DEFINE_string(task,"", "file for tasklist");
DEFINE_int32(t,0x00C01509, "header tag");//0x0915C000
DEFINE_int32(w, -1, "input image width");
DEFINE_int32(b, -1, "input image bands");
DEFINE_int32(buffer, (100*1024*1024), "buffer size");
DEFINE_int32(splice, 1, "unite images into one");
DEFINE_string(m, "", "method for processing: mean or median");
DEFINE_bool(f, false, "operator with force");
DEFINE_string(ext, ".tif", "default extension for output decoded image");
DEFINE_int32(c, 0, "compression type: 0=LOSSLESS, 1=LOSS8, 2=LOSS4, 3=NONE");
DEFINE_string(gcp, "", "pos to add gcp to tif");

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
    GDALClose(src);
    rowid += src->GetRasterYSize();

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

#ifndef NDEBUG
  FLAGS_logtostderr = 1;
#endif

  GDALAllRegister();
  CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

  std::string usage("This program decodes hyper-spectrum raw data. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + "[Decode]: " + name + " <raw data file>*\n"
            +"[DarkCorrect]: " + name + " -r0 <dark level file> <image file>\n"
            +"[BadpixelCorrect]: "+name+" -r1 <bad pixel file> <image file>\n"
            +"[NonuniformCorrect]: "+name+" -r0 <file b> -r1 <file a> <image file>\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("1.0");

  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if(argc<2){
    gflags::ShowUsageWithFlagsRestrict(argv[0], "decode");
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
            std::cerr << e.what() << std::endl;
            return 2;
        }
        GDALDataset* src = (GDALDataset*)GDALOpen(path.string().c_str(), GA_ReadOnly);
        if (src == nullptr) return 1;
        GDALDataset* dst = nullptr;
        int dst_cols = src->GetRasterXSize();
        int dst_rows = src->GetRasterYSize();
        int dst_bands = src->GetRasterCount();
        
        radiometric::ComboOperator* ops = new radiometric::ComboOperator;
        ipf::RasterOperator frame( ops,GDT_Float32);
        int store_prior[3] = {0,2,1};
        int buffer_size[3] = { src->GetRasterXSize(), src->GetRasterYSize(), src->GetRasterCount()};
		boost::filesystem::path outpath;
        int boutput = 0;
        {
            auto prior = tree.get<std::string>("decode.band_prior", "0,1,2");
            int boutput = tree.get<int>("decode.output", 0);
			std::vector<std::string> result;
			boost::split(result, prior, boost::is_any_of(","));
			if (result.size() == 3) {
				store_prior[0] = std::stoi(result[0]);
				store_prior[1] = std::stoi(result[1]);
				store_prior[2] = std::stoi(result[2]);
			}
        }
    
        BOOST_FOREACH(boost::property_tree::ptree::value_type & v, tree.get_child("decode")) {
            if (v.first != "task") continue;
            auto name = v.second.get_child("<xmlattr>.name").get_value<std::string>();
            if (name == "uniform") {
                std::string a = v.second.get<std::string>("a");
                std::string b = v.second.get<std::string>("b");
                radiometric::NonUniformCorrection* op = new radiometric::NonUniformCorrection(src->GetRasterXSize(), src->GetRasterYSize(), src->GetRasterCount());
                if (!op->load(a.c_str(), b.c_str())) {
                    std::cout << "ERROR:uniform file not loaded a and b" << std::endl;
                    return 1;
                }
                ops->Add(op);
                boutput = 1;
            }
            else if (name == "badpixel") {
                std::string b = v.second.get<std::string>("file");
                radiometric::BadPixelCorrection* op = new radiometric::BadPixelCorrection(src->GetRasterXSize(), src->GetRasterYSize(), src->GetRasterCount());
                if (!op->load( b.c_str())) {
                    std::cout << "ERROR:badpixel file not loaded a and b" << std::endl;
                    return 1;
                }
                ops->Add(op);
                boutput = 1;
            }
            else if (name == "gauss") {
                int band = v.second.get<int>("band");
                int ksize = v.second.get<int>("ksize");
                radiometric::GaussianBlur* op = new radiometric::GaussianBlur(ksize, band);
                ops->Add(op);
                boutput = 1;
            }else if(name == "dark"){
              std::string b = v.second.get<std::string>("b");
              radiometric::DarkLevelCorrection* op = new radiometric::DarkLevelCorrection(src->GetRasterXSize(), src->GetRasterYSize(), src->GetRasterCount());
              if (!op->load(b.c_str())) {
                std::cout << "ERROR:dark file not loaded a and b" << std::endl;
                return 1;
              }
              ops->Add(op);
			  boutput = 1;
            }else if(name == "statistic"){
              std::string method = v.second.get<std::string>("method", "mean");
              boost::filesystem::path dstpath(path);
              dstpath.replace_extension();
              dstpath += "_"+method+".txt";
              radiometric::FrameIterator* op = nullptr;
              if(method=="mean"){
                radiometric::MeanStdCalculator* p = new radiometric::MeanStdCalculator(buffer_size[store_prior[1]]);
                p->SetFilePath(dstpath.string().c_str());
                op = p;
              } else{
                radiometric::MedianCalculator* p = new radiometric::MedianCalculator(buffer_size[store_prior[1]]);
                p->SetFilePath(dstpath.string().c_str());
                op = p;
              }
              ops->Add(op);
            }
            else if (name == "interp" && !outpath.empty()) {
                std::vector<double> wl_old, wl_new;
                radiometric::BandInterpolator::InterpType type;
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
                    if (wl_old.size() != src->GetRasterCount()) {
                        std::cout << "ERROR: interp wl length MISMATCH" << std::endl;
						GDALClose(src);
						return 1;
                    }

					std::string t = v.second.get<std::string>("type", "pchip");
                    if (t == "spline_cubic")
                        type = radiometric::BandInterpolator::BSPLINE_CUBIC;
                    else if(t=="spline_quintic")
                        type = radiometric::BandInterpolator::BSPLINE_QUINTIC;
                    else if(t=="spline_quadratic")
                        type = radiometric::BandInterpolator::BSPLINE_QUADRATIC;
                    else if(t=="makima")
                        type = radiometric::BandInterpolator::MAKIMA;
                    else if(t=="barycentric")
                        type = radiometric::BandInterpolator::BARYCENTRIC;
                    else
                        type = radiometric::BandInterpolator::PCHIP;
                }
                dst_bands = wl_new.size();
                radiometric::BandInterpolator* op = new radiometric::BandInterpolator(std::move(wl_old), std::move(wl_new));
                op->SetInterpType(type);
                ops->Add(op);
                boutput = 1;
            }
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
			dst = GDALCreate(outpath.string().c_str(), dst_cols, dst_rows, dst_bands, src->GetRasterBand(1)->GetRasterDataType());
			if (dst == nullptr) {
				GDALClose(src);
				return 1;
			}
        }
        frame.SetStoreOrder(store_prior);
        frame.SetSource(src);
        if (dst) frame.SetDestination(dst);
        if (buffer_size[2] < dst_bands) buffer_size[2] = dst_bands;
        {
            size_t max_buffer_size = (size_t)1024 * 1024 * 1024;
            buffer_size[store_prior[2]] = max_buffer_size / (frame.GetDataTypeSize() * buffer_size[store_prior[0]] * buffer_size[store_prior[1]]);
            if (buffer_size[store_prior[2]] == 0) buffer_size[store_prior[2]] = 1;
        }
		frame.ReserveBufferSize((size_t)buffer_size[0] * buffer_size[1] * buffer_size[2]);
        if (src->GetRasterCount() < buffer_size[2])  buffer_size[2] = src->GetRasterCount();
        ipf::TileProcessing( frame.source().win+3, buffer_size, &frame);

        delete ops;
        if (dst) {
            double aop6[6];
            if(src->GetGeoTransform(aop6)==CE_None)
                dst->SetGeoTransform(aop6);
            dst->SetProjection(src->GetProjectionRef());
            GDALClose(dst);
        }
        GDALClose(src);
        return 0;
    } else if (!FLAGS_gcp.empty()) {
        HSP::Pos pos;
        pos.load(FLAGS_gcp.c_str());
        HSP::PinholeCamera cam;
        cam.Load("G:\\jincang\\camera_vnir.txt");
        HSP::LinescanModel model;
        model.SetCamera(&cam);
        model.SetPos(&pos);
//        boost::filesystem::path rpcpath(path);
//        rpcpath.replace_extension(".rpc");
        double range_samp[] = {10, 1520};
        double range_line[] = {8, pos.data().rbegin()->first-8};
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
    decode.add_token(token, token+4);

    if(verbose){
      std::cout << "[DECODE] Total to-decoded raw data number= " << rawlist.size() << std::endl;
    }

    boost::filesystem::path dirpath;
    if(FLAGS_o.empty()){
    } else if(boost::filesystem::is_directory(FLAGS_o)) {
      dirpath = FLAGS_o;
    }else if(boost::filesystem::extension(FLAGS_o).empty()){
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
        sprintf(name,"%s_splice_%d%s", srcpath.stem().string().c_str(), i, srcpath.extension().string().c_str() );
        dstpath /= name;
      }else {
        dstpath = dirpath;
        if(grps.size()>1){
          char post[512];
          sprintf(post,"_%d%s", i, dstpath.extension().string().c_str());
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
