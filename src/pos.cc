#include "pos.h"
#include "decode.h"
#include "lx_geometry_rpc.h"
#include "InterpolatorAdaptor.hpp"

#include <gdal_priv.h>
#include <gdal_alg.h>
#include <gdal_alg_priv.h>
#include <cpl_error.h>
#include <ogr_spatialref.h>
#include <ogr_srs_api.h>
#include <vector>

#include <eigen3/Eigen/Dense>

#ifndef SRS_WKT_WGS84_LAT_LONG
#define SRS_WKT_WGS84_LAT_LONG SRS_WKT_WGS84
#endif

#ifdef _DEBUG
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>
template<class DataIt>
bool WriteVector(const char* filepath, DataIt first, DataIt last) {
    std::ofstream file(filepath, std::ios::out);
    file << std::setiosflags(std::ios::fixed);
    while (first != last) {
        file << std::setprecision(6) << *first << "\t";
        ++first;
    }
    file.close();
    return true;
}
#endif

namespace PosParse{

bool MARK1PVAA(const char* msg, HSP::Pos::record& rec){
  const char* p = strrchr(msg, ';');
  if (p==nullptr) return false;

  if( sscanf(p+1, "%u%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf"
             , &rec.week, &rec.time
             , rec.blh+1, rec.blh, rec.blh+2
             , rec.vx+1, rec.vx, rec.vx+2
             , rec.a, rec.a+1, rec.a+2
             ) != 11){
    return false;
  }

  if (true) {
      rec.a[0] *= (M_PI / 180);
      rec.a[1] *= (M_PI / 180);
      rec.a[2] *= (M_PI / 180);
      rec.a[0] = rec.a[1] = 0;
  }

  return true;
}

bool MARKPOSA(const char* msg, HSP::Pos::record& rec){
  const char* p = strstr(msg, "SINGLE");
  if( sscanf(p+strlen("SINGLE"), "%*c%lf%*c%lf%*c%lf"
             , rec.blh+1, rec.blh, rec.blh+2
             ) != 3){
    return false;
  }
  return true;
}

bool Parse(const char* msg, HSP::Pos::record& rec){
  if(strstr(msg,"MARK1PVAA")){
    return MARK1PVAA(msg, rec);
  }
  if(strstr(msg, "MARKPOSA"))
    return MARKPOSA(msg, rec);
  return false;
}

};

namespace HSP{

    typedef pchip D1Interp;
    typedef cubic_hermite D2Interp;

class Interpolator{
 public:
  InterpolatorAdaptor* _x[3];
  InterpolatorAdaptor* _a[3];
 public:
  Interpolator() {
    _x[0] = _x[1] = _x[2] = nullptr;
    _a[0] = _a[1] = _a[2] = nullptr;
  }
  ~Interpolator() {
    for(int i=0; i<3; ++i)
      if(_x[i]) delete _x[i];
    for(int i=0; i<3; ++i)
      if(_a[i]) delete _a[i];
  }
  int GetPosition(double lineid, double x[3]){
    x[0] = (*_x[0])(lineid);
    x[1] = (*_x[1])(lineid);
    x[2] = (*_x[2])(lineid);
    return 1;
  }
  int GetAngle(double lineid, double a[3]){
    a[0] = (*_a[0])(lineid);
    a[1] = (*_a[1])(lineid);
    a[2] = (*_a[2])(lineid);
    return 1;
  }
};

Interpolator* ToInterpolator(void* t){
  return (Interpolator*)t;
}

void DestroyInterpolator(void* t){
  delete (Interpolator*)t;
}

Interpolator* CreateInterpolator(std::map<int, Pos::record>& data){
  if(data.size()<3) return nullptr;
  ValueContainer t, x[3], vx[3], a[3];
  bool vx_flag = true;
  t.reserve(data.size());
  for(int i=0; i<3; ++i){
    x[i].reserve(data.size());
    vx[i].reserve(data.size());
    a[i].reserve(data.size());
  }
  auto data_it1 = data.begin();
  ++data_it1;
  auto data_it0 = data_it1++;
  double ratio = (data_it1->second.time-data_it0->second.time)/(data_it1->first-data_it0->first);
  for(auto it=data.begin(); it!=data.end(); ++it){
    t.push_back(it->first);
    if(vx_flag && fabs(it->second.vx[0])<1e-4&&fabs(it->second.vx[1])<1e-4&&fabs(it->second.vx[2])<1e-4) vx_flag = false;
    for(int i=0; i<3; ++i){
      x[i].push_back(it->second.x[i]);
      if(vx_flag) vx[i].push_back(it->second.vx[i]*ratio);
      a[i].push_back(it->second.a[i]);
    }
  }

  Interpolator* interp = new Interpolator;
  for(int i=0; i<3; ++i){
    ValueContainer t1 = t;
    if(vx_flag) interp->_x[i] = new D2Interp(std::move(t1), std::move(x[i]), std::move(vx[i]));
    else interp->_x[i] = new D1Interp(std::move(t1), std::move(x[i]));
    t1 = t;
    interp->_a[i] = new D1Interp(std::move(t1), std::move(a[i]));
  }

  return interp;
}

Eigen::Quaterniond Angle2Quat(double ax, double ay, double az, const int order[3]) {
	Eigen::AngleAxisd aa[3];
	aa[0] = Eigen::AngleAxisd(ax, Eigen::Vector3d::UnitX());
	aa[1] = Eigen::AngleAxisd(ay, Eigen::Vector3d::UnitY());
	aa[2] = Eigen::AngleAxisd(az, Eigen::Vector3d::UnitZ());
    return aa[order[2]] * aa[order[1]] * aa[order[0]];
}
Eigen::Quaterniond rpy2Quat(double ax, double ay, double az) {
    const int order[3] = {1,0,2};
    return Angle2Quat(ax, ay, az, order);
}


Pos::~Pos(){
  if(_reprojector)
    GDALDestroyTransformer(_reprojector);
  if(_interpolator)
    DestroyInterpolator(_interpolator);
  if (_utm_wkt)
      CPLFree(_utm_wkt);
}

bool Pos::load(const char *filepath){
  const char* ext = strrchr(filepath, '.');

  if(ext==nullptr)
    return false;

  bool ret = false;

  if(strcmp(ext+1, "aux")==0)
    ret = LoadAux(filepath);

  if (ret) {
      Reprojection();
      BuildInterpolator();
  }

  return ret;
}

bool Pos::LoadAux(const char *filepath) {
  FILE* fp = fopen(filepath, "rb");

  std::map<int, record> pos_data;

  int line = 0;
  double last_time; last_time = -1;
  while(1){
    int frame, count;
    if(fread(&frame, sizeof(int), 1, fp)!=1) break;
    if(fread(&count, sizeof(int), 1, fp)!=1) break;
    if(count==0) continue;
    std::vector<char> buffer(count);
    if(fread(buffer.data(), sizeof(char), count, fp)!=count) break;

    std::vector<char> pos(count*2.0/3);
    HSP::cvt_12to16( buffer.data(), count, pos.begin());

    record rec;
    if (PosParse::Parse(pos.data(), rec)
        && fabs(rec.time-last_time)>1e-5) {
      last_time = rec.time;
      rec.frame = frame;
      pos_data.insert({line,rec});
    }
    ++line;
  }

  if(pos_data.size()>0) std::swap(_data, pos_data);

  fclose(fp);
  return true;
}

void Pos::test() {
  Reprojection();
  BuildInterpolator();
  Check();
}

void Pos::Reprojection(){
  if(_data.size()<1) return;
  if(!_reprojector){
    double lon,lat;
    auto& item = _data.begin()->second;
    lon = item.blh[0];
    lat = item.blh[1];
    OGRSpatialReference srcSRS, dstSRS;
    srcSRS.SetFromUserInput(SRS_WKT_WGS84_LAT_LONG);
    dstSRS.SetUTM(int((lon+180)/6)+1,lat>0);

    if(srcSRS.IsCompound()){
      srcSRS.StripVertical();
    }
    if(dstSRS.IsCompound()){
      dstSRS.StripVertical();
    }

#if GDAL_VERSION_MAJOR >= 3
    srcSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dstSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    _reprojector = GDALCreateReprojectionTransformerEx( reinterpret_cast<OGRSpatialReferenceH>(&srcSRS), reinterpret_cast<OGRSpatialReferenceH>(&dstSRS), nullptr );
#else
    char* wkt_src = nullptr;
    char* wkt_dst = nullptr;
    srcSRS.exportToWkt(&(wkt_src));
    dstSRS.exportToWkt(&(wkt_dst));
    _reprojection = GDALCreateReprojectionTransformer(wkt_src,wkt_dst);
    if(wkt_src!=nullptr) CPLFree(wkt_src);
    if(wkt_dst!=nullptr) CPLFree(wkt_dst);
#endif
    dstSRS.exportToWkt(&_utm_wkt);
  }
  for(auto it=_data.begin(); it!=_data.end(); ++it){
    auto& item = it->second;
    memcpy(item.x, item.blh, sizeof(double)*3);
    int valid;
    GDALUseTransformer(_reprojector, 0, 1, item.x, item.x+1, item.x+2, &valid);
    if (_offset[0] == 0) {
        _offset[0] = item.x[0];
        _offset[1] = item.x[1];
        _offset[2] = item.x[2];
        item.x[0] = item.x[1] = item.x[2] = 0;
    }
    else {
        item.x[0] -= _offset[0];
        item.x[1] -= _offset[1];
        item.x[2] -= _offset[2];
    }
  }
}

void Pos::BuildInterpolator(){
  if(_interpolator) DestroyInterpolator(_interpolator);
  _interpolator = CreateInterpolator(_data);
}

int Pos::GetPosition(double lineid, double x[3]){
  return ToInterpolator(_interpolator)->GetPosition(lineid, x);
}
int Pos::GetAngle(double lineid, double a[3]){
  return ToInterpolator(_interpolator)->GetAngle(lineid, a);
}

int Pos::GetPos(double lineid, double x[3], double a[3]){
  int ret = 0;
  ret += ToInterpolator(_interpolator)->GetPosition(lineid, x);
  ret += ToInterpolator(_interpolator)->GetAngle(lineid, a);
  return ret;
}

Eigen::Quaterniond Pos::GetQuaternion(double lineid){
  double a[3];
  GetAngle(lineid, a);
  return rpy2Quat(a[1], a[0], a[2]);
//  auto pitch = Eigen::AngleAxisd(a[1], Eigen::Vector3d::UnitX());
//  auto roll = Eigen::AngleAxisd(a[0], Eigen::Vector3d::UnitY());
//  auto yaw = Eigen::AngleAxisd(a[2], Eigen::Vector3d::UnitZ());
//  return yaw*pitch*roll;
}

int Pos::Check(){
  if(_data.size()<2) return 0;
  std::vector<std::pair<double, double> > err;
  err.reserve(_data.size()-1);
  int count = 0;
  auto it1 = _data.begin();
  auto it0 = it1++;
  while (it1 != _data.end()) {
    auto& t0 = it0->second;
    auto& t1 = it1->second;
    double dt = (t1.time-t0.time)/2;
    double dx[3], dc[3];
    dc[0] = (t1.vx[0]+t0.vx[0])*dt;
    dc[1] = (t1.vx[1]+t0.vx[1])*dt;
    dc[2] = (t1.vx[2]+t0.vx[2])*dt;
    dx[0] = t1.x[0]-t0.x[0];
    dx[1] = t1.x[1]-t0.x[1];
    dx[2] = t1.x[2]-t0.x[2];
    double norm_dx = std::sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
    double norm_dc = std::sqrt(dc[0]*dc[0]+dc[1]*dc[1]+dc[2]*dc[2]);
    double dp = dc[0]*dx[0]+dc[1]*dx[1]+dc[2]*dx[2];

    err.push_back(std::make_pair(dp/norm_dc/norm_dx, norm_dc/norm_dx));

    it0 = it1++;
  }

  it1 = _data.begin();
  it0 = it1++;
  double x[3], a[3];
  int rx = GetPosition(it1->first, x);
  int ra = GetAngle(it1->first, a);

  rx = GetPosition((it0->first+it1->first)/2, x);
  ra = GetAngle((it0->first+it1->first)/2, a);

  return 0;
}

int Pos::Cvt_BLH2Local(double x[3]) {
    int valid;
    GDALUseTransformer(_reprojector, 0, 1, x, x + 1, x + 2, &valid);
    x[0] -= _offset[0];
    x[1] -= _offset[1];
    x[2] -= _offset[2];
    return 1;
}

int Pos::Cvt_Local2BLH(double x[3]) {
    x[0] += _offset[0];
    x[1] += _offset[1];
    x[2] += _offset[2];
    int valid;
    GDALUseTransformer(_reprojector, 1, 1, x, x + 1, x + 2, &valid);
    return 1;
}

Eigen::Matrix<double, 4, Eigen::Dynamic> BackProjectionToPlane(const CameraMatrixType& camera, const Eigen::Matrix<double, 3, Eigen::Dynamic>& im, const Eigen::Matrix<double, 4, 3>& plane_transformation) {
    Eigen::Matrix<double, 3, 3> m = camera * plane_transformation;
    Eigen::Matrix<double, 3, Eigen::Dynamic> planex = m.inverse() * im;
    return plane_transformation * planex;
}

#define EPSILON std::numeric_limits<double>::epsilon()
bool plane_transformation_2pto3p(const Eigen::Matrix<double,4,1>& plane_equation, Eigen::Matrix<double,4,3>& m){
  if(fabs(plane_equation[0])>EPSILON){
    m << plane_equation[1],plane_equation[2],plane_equation[3],
        -plane_equation[0],0,0,
        0,-plane_equation[0],0,
        0,0,-plane_equation[0];
  }else if(fabs(plane_equation[1])>EPSILON){
    m << -plane_equation[1],0,0,
        plane_equation[0],plane_equation[2],plane_equation[3],
        0,-plane_equation[1],0,
        0,0,-plane_equation[1];
  }else if(fabs(plane_equation[2])>EPSILON){
    m << 1,0,0,
        0,1,0,
        0,0,-plane_equation[3]/plane_equation[2],
        0,0,1;
  }else if(fabs(plane_equation[3])>EPSILON){
    m<<1,0,0,
        0,1,0,
        0,0,1,
        0,0,0;
  }else return false;
  return true;
}

Eigen::Matrix<double,4,Eigen::Dynamic> BackProjectionToPlane(const CameraMatrixType& camera, const Eigen::Matrix<double,3,Eigen::Dynamic>& im, const Eigen::Matrix<double,4,1>& plane_equation){
  Eigen::Matrix<double,4,3> plane_transformation;
  if(plane_transformation_2pto3p(plane_equation, plane_transformation)){
    return BackProjectionToPlane(camera, im, plane_transformation);
  }
  return Eigen::Matrix<double,4,Eigen::Dynamic>::Zero(4,1);
}

bool PinholeCamera::Load(const char* filepath) {
    FILE* fp = fopen(filepath, "r");
    if (fp == nullptr) return false;
    bool ret = false;
    char line[512];
    double f, pixel, x0, y0;
    if (fgets(line, 512, fp) && sscanf(line, "%lf%lf%lf%lf", &f, &pixel, &x0, &y0) == 4) {
        SetFocalLength(f / pixel);
        SetPrincipalPoint(x0, y0);
        double a[3], offset[3];
        if (fgets(line, 512, fp) && 
            sscanf(line, "%lf%lf%lf%lf%lf%lf", a, a+1, a+2, offset, offset+1, offset+2) == 6) {
            SetAngles(a);
            SetTranslation(offset);
			ret = true;
        }
    }
    fclose(fp);
    return ret;
}

void PinholeCamera::SetAngles(double a[3]) {
    _pose = rpy2Quat(a[0], a[1], a[2]);
}
void PinholeCamera::SetTranslation(double offset[3]) {
    _translation << offset[0], offset[1], offset[2];
}

CameraMatrixType LinescanModel::CameraMatrix(double linenumber) const
{
  CameraMatrixType m = _camera->CameraMatrix();
  Eigen::Vector3d x;
  _pos->GetPosition(linenumber, x.data());
  auto pose = _pos->GetQuaternion(linenumber);
  Eigen::Matrix4d h = Eigen::Matrix4d::Identity();
  h.block<3, 3>(0,0) = pose.matrix();
  h.block<3, 1>(0,3) << -h.block<3, 3>(0, 0)*x;
  Eigen::Matrix4d hcvt;
  hcvt << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, -1, 0,
      0, 0, 0, 1;
  return m*hcvt*h;
}

std::vector<double> GenerateSample(double range[2], int num) {
    std::vector<double> samps;
    double interval = (range[1] - range[0]) / (num-1);
    for (int i = 0; i < num; ++i)
        samps.push_back(range[0] + i * interval);
    return samps;
}

bool LinescanModel::GenerateRPC(double range_samp[2], double range_line[2], double range_height[2], const char* rpcpath) {
    std::vector<double> samp1d =  GenerateSample(range_samp, 3);
    std::vector<double> height1d = GenerateSample(range_height, 3);
    std::vector<double> line1d ;//= GenerateSample(range_line, 100);
   auto& posdata = _pos->_data;
   for (auto& pos : posdata)
       if (pos.first >= range_line[0] && pos.first <= range_line[1]) line1d.push_back(pos.first);
    for (auto& h : height1d)
        h -= _pos->_offset[2];
    size_t count = samp1d.size() * line1d.size() * height1d.size();
    std::vector<double> grid_samp, grid_line, grid_height, grid_lon, grid_lat;
	grid_samp.reserve(count), grid_line.reserve(count), grid_height.reserve(count), grid_lon.reserve(count), grid_lat.reserve(count);
    for (auto& l : line1d) {
        auto cam = CameraMatrix(l);
		for (auto& h : height1d) {
            Eigen::Matrix<double, 4, 1> height_plane;
            height_plane << 0, 0, 1, -h;
			for (auto& s : samp1d) {
                grid_samp.push_back(s);
                grid_line.push_back(l);
                Eigen::Matrix<double, 3, 1> im;
                im << s, 0, 1;
                Eigen::Matrix<double,4,1> t = BackProjectionToPlane(cam, im, height_plane); 
                Eigen::Matrix<double, 3, 1> x = t.hnormalized();
                _pos->Cvt_Local2BLH(x.data());
                grid_lon.push_back(x[0]);
                grid_lat.push_back(x[1]);
                grid_height.push_back(x[2]);
            }
        }
    }
#ifdef _DEBUG
    char debug_path[512] = "H:\\jinchang\\ws\\vnir\\debug\\" ;
    char* ps = debug_path + strlen(debug_path);
    strcpy(ps, "samp.txt");
    WriteVector(debug_path, grid_samp.begin(), grid_samp.end());
    strcpy(ps, "line.txt");
    WriteVector(debug_path, grid_line.begin(), grid_line.end());
    strcpy(ps, "lon.txt");
    WriteVector(debug_path, grid_lon.begin(), grid_lon.end());
    strcpy(ps, "lat.txt");
    WriteVector(debug_path, grid_lat.begin(), grid_lat.end());
    strcpy(ps, "height.txt");
    WriteVector(debug_path, grid_height.begin(), grid_height.end());
#endif
  //   GDALDataset* src = (GDALDataset*)GDALOpen(rpcpath, GA_Update);
  //   GDAL_GCP* pasGCPs = new GDAL_GCP[count];
  //   GDALInitGCPs(count, pasGCPs);
  //   for (int i = 0; i < count; ++i) {
  //       pasGCPs[i].dfGCPPixel= grid_samp[i];
  //       pasGCPs[i].dfGCPLine = grid_line[i];
  //       pasGCPs[i].dfGCPX = grid_lon[i];
  //       pasGCPs[i].dfGCPY = grid_lat[i];
  //       pasGCPs[i].dfGCPZ = grid_height[i];
  //   }
	// src->SetGCPs(count, pasGCPs, SRS_WKT_WGS84_LAT_LONG);
  //   GDALDeinitGCPs(count, pasGCPs);
  //   delete[] pasGCPs;
    
  //   return true;
    xlingeo::Rpc rpc;
    rpc.Solve(grid_samp.data(), grid_line.data(), grid_lat.data(), grid_lon.data(), grid_height.data(), grid_samp.size());
    rpc.Save(rpcpath);
    return true;
}

void LinescanModel::test(){
    auto cam = CameraMatrix(85);
    Eigen::Vector4d X;

    X << 101.8514647,	38.50230176,	2046.9854	, 1;

    _pos->Cvt_BLH2Local(X.data());
    Eigen::Vector2d x = (cam * X).hnormalized();

    std::cout << x << std::endl;
	return;
}


};
