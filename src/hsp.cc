#include <gdal_priv.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include <cctype>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "raster.hpp"
#include "ipf.hpp"
#include "util/config.hpp"

DEFINE_string(o,"", "output directory");
DEFINE_string(task,"", "file for tasklist");
DEFINE_int32(buffer, 2048, "buffer size/MB");//(std::numeric_limits<int>::max)()
DEFINE_string(nodata, "", "specify the nodata value.");
DEFINE_int32(fc, 1, "flush cache mode: 0, 1, 2");
DEFINE_string(ot, "", "output data type {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}");
DEFINE_bool(ow, false, "whether to overwrite the destination image");
DEFINE_bool(oc, false, "Only save the window of source image.");
DEFINE_string(srcwin, "", "the window of source image");//extend nuc

static bool ValidateTask(const char* flagname, const std::string& value) {
  if(!boost::filesystem::exists(value)){
    return false;
  }
  return true;
}
DEFINE_validator(task, &ValidateTask);

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

int main(int argc, char* argv[]){

  FLAGS_logtostderr = 1;

  GDALAllRegister();
  CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

  std::string usage("This program processes temporal raw data. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + name + " -task <task file> <image file>\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("1.4");

  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if(argc<2){
    gflags::ShowUsageWithFlagsRestrict(argv[0], "hsp");
    return 1;
  }

  int verbose = FLAGS_v;

  boost::filesystem::path path(argv[1]);
  boost::property_tree::ptree tree;
  bool flag_xml = true;
  try {
    if ( FLAGS_task.rfind(".xml") != std::string::npos || FLAGS_task.rfind(".XML") != std::string::npos ){
      boost::property_tree::read_xml(FLAGS_task, tree, boost::property_tree::xml_parser::trim_whitespace);
    }else{
      boost::property_tree::read_json(FLAGS_task, tree);
      flag_xml = false;
    }
  }
  catch (const boost::property_tree::ptree_error& e) {
    LOG(ERROR) << e.what();
    return 2;
  }
  boost::filesystem::path taskpath(FLAGS_task.c_str());

  GDALDataset* src = (GDALDataset*)GDALOpen(path.string().c_str(), GA_ReadOnly);
  if (src == nullptr){
    LOG(ERROR) << "CANNOT OPEN input file " << path.string() ;
    return 1;
  }

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
    while (nodata.size()>0 && std::isnan(nodata.back())) nodata.pop_back();
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
      else src_win[i+1] += vt-src_win[i];
    }
  }

  int dst_size[3], dst_win[6], buffer_size[3];
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
  xlingsky::raster::ComboOperator* ios = new xlingsky::raster::ComboOperator;
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
  int retcode = 0;

  boost::property_tree::ptree tasks = tree.get_child("HSP.task");
  for (const auto& task : tasks) {
    auto& v = task.second;
    auto name = v.get<std::string>( flag_xml?"<xmlattr>.name":"name");
    xlingsky::util::config config( v, taskpath);
    if (name == "dpc") {
      xlingsky::raster::radiometric::DefectivePixelCorrection* op = xlingsky::raster::radiometric::DefectivePixelCorrection::Create(config, src_size, store_prior);
      if (op==nullptr) {
        std::cout << "DefectivePixelCorrection: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      ops->Add(op);
      boutput = 1;
    }
    else if(name == "statistic"){
      std::string method = v.get<std::string>("method", "mean");
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
        float cut_ratio_lower = v.get<float>("cut_lower", 0);
        float cut_ratio_upper = v.get<float>("cut_upper", 0);
        xlingsky::raster::radiometric::MeanStdCalculator* p = new xlingsky::raster::radiometric::MeanStdCalculator(src_size[store_prior[1]], cut_ratio_lower, cut_ratio_upper);
        bool nuc = v.get<bool>("nuc", false);
        if(nuc){
          boost::filesystem::path xmlpath = dstpath;
          xmlpath.replace_extension(".xml");
          p->SetXmlPath(xmlpath.string().c_str());
          p->SetDimOrder(store_prior);
          bool apply = v.get<bool>("nuc.<xmlattr>.apply", true);
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
    } else if(name == "extractor"){
      std::string type = v.get<std::string>("type", "peak");
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
      dstpath += "_"+type;
      xlingsky::raster::Operator* op = nullptr;
      if(type=="peak"){
        xlingsky::raster::extractor::PeakDetection* p = 
        new xlingsky::raster::extractor::PeakDetection(src_size[store_prior[2]], v);
        std::string pointpath = dstpath.string() + "_points.txt";
        std::string linepath = dstpath.string() + "_lines.txt";
        p->SetFilePath(pointpath.c_str(), linepath.c_str());
        op = p;
        boutput = 1;
      } else{
        dstpath += ".txt";
        auto *p = new xlingsky::raster::extractor::LSDExtractor(
            src_size[store_prior[2]], v);
        p->SetFilePath(dstpath.string().c_str());
        op = p;
      }
      ops->Add(op);
      if (FLAGS_v<2) {
        boutput = 0;
      }
    }else if (name == "background") {
      auto * op = new xlingsky::raster::BackgroundExtraction<xlingsky::raster::ImageWriter>();
      if (!op->load_config(config)) {
        std::cout << "BACKGROUND: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }

      std::string dstpath = "in_memory";
      bool removal = v.get<bool>("removal", true);
      if (!removal) {
        boost::filesystem::path temp;
        if (FLAGS_o.empty()) {
          temp = path;
        } else {
          temp = FLAGS_o;
          if (boost::filesystem::is_directory(dstpath)) {
            temp /= path.filename();
          }
        }
        temp.replace_extension();
        temp += "_background.tif";
        dstpath = temp.string();
      }

      auto* imwriter = new xlingsky::raster::ImageWriter(
          dstpath.c_str(), src_size[store_prior[0]], src_size[store_prior[1]], op->get_blknum(src_size[store_prior[2]]), GDT_Float32);
      op->set_exporter(imwriter);
      ios->Add(imwriter);
      ops->Add(op);
      if (removal) {
        auto* t = new xlingsky::raster::common::Minus(imwriter->ptr(), op->get_ksize());
        ops->Add(t);
      }
    }else if (name == "removal") {
      auto* op = new xlingsky::raster::common::Minus();
      ops->Add(op);
      if (!op->load_config(config)) {
        std::cout << "REMOVAL: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      if (!op->good()) {
        if (argc < 3 || !op->set_dataset(argv[2]) ) {
          std::cout << "REMOVAL: need a background dataset !" << std::endl;
          retcode = 1;
          goto finished;
        }
      }
      boutput = 1;
      if (FLAGS_ot.empty()){
        FLAGS_ot = "float32";
      }
    }else if (name == "anomaly") {
      xlingsky::raster::RXDetector* op = new xlingsky::raster::RXDetector();
      if (!op->load_config(config)) {
        std::cout << "ANOMALYDETECTION: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      dst_size[store_prior[2]] = dst_win[store_prior[2]*2+1] = op->get_blknum(src_size[store_prior[2]]);
      ops->Add(op);
      boutput = 1;
      if (FLAGS_ot.empty()){
        FLAGS_ot = "float32";
      }
    }
  }

  VLOG(1) << "Total number of tasks is " << ops->size();
  if(ops->size()<1){
    retcode = 1;
    goto finished;
  }

  if (boutput) {
    if (!FLAGS_o.empty()) outpath = FLAGS_o;
    else {
      outpath = path;
      outpath.replace_extension();
      outpath += "_mod";
      if (path.extension() == ".png" || path.extension() == ".jpg") {
        outpath += ".tif";
      }else {
        outpath += path.extension();
      }
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
        retcode = 3;
        goto finished;
      }
    }else
      dst = GDALCreate(outpath.string().c_str(), dst_size[0], dst_size[1], dst_size[2], datatype);
    if (dst == nullptr) {
      retcode = 3;
      goto finished;
    };
  }
  frame.SetStoreOrder(store_prior);
  frame.SetSource(src, src_win);
  if (dst) frame.SetDestination(dst, dst_win);

  {
    buffer_size[0] = std::max(src_win[1], dst_win[1]);
    buffer_size[1] = std::max(src_win[3], dst_win[3]);
    buffer_size[2] = std::max(src_win[5], dst_win[5]);
    size_t max_buffer_size = (size_t)FLAGS_buffer<<20;
    size_t l = max_buffer_size / ((size_t)frame.GetDataTypeSize() * buffer_size[store_prior[0]] * buffer_size[store_prior[1]]);
    if (l == 0) l = 1;
    if (buffer_size[store_prior[2]]>l)
      buffer_size[store_prior[2]] = l;
  }
  frame.ReserveBufferSize((size_t)buffer_size[0] * buffer_size[1] * buffer_size[2]);
  if (src_win[5] < buffer_size[2])  buffer_size[2] = src_win[5];

  xlingsky::ipf::TileProcessing( frame.source().win+3, buffer_size, &frame);

  finished:

  delete ops;
  delete ios;
  if (dst) {
    if (FLAGS_oc) {
      double aop6[6];
      if (src->GetGeoTransform(aop6) == CE_None)
        dst->SetGeoTransform(aop6);
      dst->SetProjection(src->GetProjectionRef());
    }
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
  return retcode;
}
