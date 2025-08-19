#include <gdal_priv.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <iostream>
#include <cctype>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "raster.hpp"
#include "ipf.hpp"
#include "raster/anomalydetection.hpp"
#include "raster/common.hpp"
#include "raster/extractor.hpp"
#include "raster/filter.hpp"
#include "raster/io.hpp"
#include "raster/operator.h"
#include "util/config.hpp"

//////////////////////////////////////////////////////////////////////////
/*
new framework:
1. flowing arguments: tile info like rasterio
2. change tilemanager|buffer & io to operator
*/
//////////////////////////////////////////////////////////////////////////

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

  std::string usage("This program processes hyper-spectrum raw data. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + name + " -task <task xml file> <image file>\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("2.4");

  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if(argc<2){
    gflags::ShowUsageWithFlagsRestrict(argv[0], "hsp");
    return 1;
  }

  int verbose = FLAGS_v;

  boost::filesystem::path path(argv[1]);
  boost::property_tree::ptree tree;
  try {
    boost::property_tree::read_xml(FLAGS_task, tree, boost::property_tree::xml_parser::trim_whitespace);
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

  BOOST_FOREACH(boost::property_tree::ptree::value_type & v, tree.get_child("HSP")) {
    if (v.first != "task") continue;
    auto name = v.second.get_child("<xmlattr>.name").get_value<std::string>();
    xlingsky::util::config config( v.second, taskpath);
    if (name == "uniform") {
      xlingsky::raster::radiometric::NonUniformCorrection* op = new xlingsky::raster::radiometric::NonUniformCorrection(src_size[store_prior[0]],src_size[store_prior[1]]);
      if (!op->load_config(config)) {
        std::cout << "NonUniformCorrection: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      if(nodata.size()>0) op->AddNoDataValue(nodata.begin(), nodata.end());
      ops->Add(op);
      boutput = 1;
    }
    else if (name == "dpc") {
      xlingsky::raster::radiometric::DefectivePixelCorrection* op = xlingsky::raster::radiometric::DefectivePixelCorrection::Create(config, src_size, store_prior);
      if (op==nullptr) {
        std::cout << "DefectivePixelCorrection: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      ops->Add(op);
      boutput = 1;
    }
    else if (name == "extend"){
      xlingsky::raster::common::Extend* op = new xlingsky::raster::common::Extend();
      if (!op->load( config.filepath("file").string().c_str(), src_size[store_prior[0]], src_size[store_prior[1]])) {
        std::cout << "EXTEND: FALIED to load file!" << std::endl;
        retcode = 1;
        goto finished;
      }
      ops->Add(op);
      boutput = 1;
    }
    else if (name == "gauss") {
      xlingsky::raster::filter::LinearFilter* op = new xlingsky::raster::filter::LinearFilter();
      if (!op->load_config(config)) {
        std::cout << "GAUSSIANFILTER: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      ops->Add(op);
      boutput = 1;
    }else if (name == "median") {
      std::string method = v.second.get<std::string>("method", "else");
      xlingsky::raster::Operator* op = nullptr;
      if(method == "adapt"){
        xlingsky::raster::filter::Despeckle* p = new xlingsky::raster::filter::Despeckle();
        if (!p->load_config(config)) {
          std::cout << "DESPECKLE: config not loaded correctly!"
                    << std::endl;
          retcode = 1;
          goto finished;
        }
        op = p;
      }else{
        xlingsky::raster::filter::MedianBlur* p = new xlingsky::raster::filter::MedianBlur();
        if (!p->load_config(config)) {
          std::cout << "MEDIANBLUR: config not loaded correctly!"
                    << std::endl;
          retcode = 1;
          goto finished;
        }
        op = p;
      }
      ops->Add(op);
      boutput = 1;
    }else if(name == "dark"){
      std::string a = v.second.get<std::string>("a","");
      float dst_min = v.second.get<float>("dst_min", (std::numeric_limits<float>::min)());
      float dst_max = v.second.get<float>("dst_max", (std::numeric_limits<float>::max)());
      xlingsky::raster::radiometric::PixelCorrection* op = nullptr;
      if(a.empty()){
        xlingsky::raster::radiometric::DarkBackgroundCorrection* p = new xlingsky::raster::radiometric::DarkBackgroundCorrection(src_size[store_prior[0]],src_size[store_prior[1]], dst_min, dst_max);
        if (!p->load(config.filepath("b").string().c_str())) {
          std::cout << "ERROR:dark file not loaded b" << std::endl;
          retcode = 1;
          goto finished;
        }
        op = p;
      }else{
        xlingsky::raster::radiometric::DarkBackgroundLinear* p = new xlingsky::raster::radiometric::DarkBackgroundLinear(src_size[store_prior[0]],src_size[store_prior[1]],src_size[store_prior[2]], dst_min, dst_max);
        if (!p->load(config.filepath("a").string().c_str(),
                     config.filepath("b").string().c_str(),
                     config.filepath("index", "").string().c_str())) {
          std::cout << "ERROR:linear dark file not loaded a,b,index" << std::endl;
          retcode = 1;
          goto finished;
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
      xlingsky::raster::radiometric::NucCalculator* op = new xlingsky::raster::radiometric::NucCalculator(src_size[store_prior[1]]);
      if (!op->load_config(config)) {
        std::cout << "NUCCALCULATOR: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
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
        if (v.second.get<bool>("debug", false)) {
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
      xlingsky::raster::spectrum::Interpolator* op = new xlingsky::raster::spectrum::Interpolator(src->GetRasterXSize(), src->GetRasterCount());
      if (!op->load_config(config)) {
        std::cout << "SPECTRUMINTERPOLATOR: config not loaded correctly!" << std::endl;
        retcode = 1;
        goto finished;
      }
      dst_size[2] = dst_win[5] = op->num_samples_new();
      ops->Add(op);
      boutput = 1;
    } else if (name == "destripe") {
      int tile_size = v.second.get<int>("tile_size", 31);
      xlingsky::raster::enhancement::Destripe* op =
          new xlingsky::raster::enhancement::Destripe(tile_size);
      ops->Add(op);
      boutput = 1;
    }else if(name == "render"){
      std::string m = v.second.get<std::string>("mode", "WALLIS|MINMAX|CUT|TILE");
      std::transform(m.begin(), m.end(), m.begin(),
                     [](unsigned char c){ return std::tolower(c); });
      bool bglobal = (m.find("global")!=std::string::npos);

      xlingsky::raster::LutCreator* lut = xlingsky::raster::LutCreate(config);
      if(FLAGS_ot.empty()){
        if(lut->dst_maximum()<(int)(std::numeric_limits<unsigned char>::max)()+1)
          FLAGS_ot = "byte";
      }

      int resample_col_step = v.second.get<float>("resample_col_step", 3);
      int resample_row_step = v.second.get<float>("resample_row_step", 3);
      if(m.find("tile")!=std::string::npos){
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
            new xlingsky::raster::BiTileLut( tile_cols, tile_rows, bglobal?xlingsky::raster::BiTileLut::GLOBAL:xlingsky::raster::BiTileLut::NONE);
        op->set_lut_creator(lut);
        op->set_resample_interval(resample_col_step, resample_row_step);
        ops->Add(op);
      }else{
        xlingsky::raster::enhancement::Render* op =
            new xlingsky::raster::enhancement::Render(bglobal?xlingsky::raster::enhancement::Render::GLOBAL:xlingsky::raster::enhancement::Render::NONE);
        op->set_lut_creator(lut);
        op->set_resample_interval(resample_col_step, resample_row_step);
        ops->Add(op);
      }
      boutput = 1;
    } else if(name == "extractor"){
      std::string type = v.second.get<std::string>("type", "peak");
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
        new xlingsky::raster::extractor::PeakDetection(src_size[store_prior[2]], v.second);
        std::string pointpath = dstpath.string() + "_points.txt";
        std::string linepath = dstpath.string() + "_lines.txt";
        p->SetFilePath(pointpath.c_str(), linepath.c_str());
        op = p;
        boutput = 1;
      } else{
        dstpath += ".txt";
        auto *p = new xlingsky::raster::extractor::LSDExtractor(
            src_size[store_prior[2]], v.second);
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
      bool removal = v.second.get<bool>("removal", true);
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
