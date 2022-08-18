#include "pos.h"
#include "decode.h"

#include <gdal_alg.h>
#include <gdal_alg_priv.h>
#include <cpl_error.h>
#include <ogr_spatialref.h>
#include <ogr_srs_api.h>

#include <eigen3/Eigen/Dense>

#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/pchip.hpp>

#ifndef SRS_WKT_WGS84_LAT_LONG
#define SRS_WKT_WGS84_LAT_LONG SRS_WKT_WGS84
#endif

namespace PosParse{

bool MARK1PVAA(const char* msg, HSP::Pos::record& rec){
  const char* p = strrchr(msg, ';');
  if (p==nullptr) return false;

  if( sscanf(p+1, "%u%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf"
             , &rec.week, &rec.time
             , rec.blh, rec.blh+1, rec.blh+2
             , rec.vx, rec.vx+1, rec.vx+2
             , rec.a, rec.a+1, rec.a+2
             ) != 11){
    return false;
  }

  return true;
}

bool MARKPOSA(const char* msg, HSP::Pos::record& rec){
  const char* p = strstr(msg, "SINGLE");
  if( sscanf(p+strlen("SINGLE"), "%*c%lf%*c%lf%*c%lf"
             , rec.blh, rec.blh+1, rec.blh+2
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

typedef std::vector<double> ValueContainer;

class InterpolatorBase {
 public:
  typedef ValueContainer::value_type value_type;
  virtual ~InterpolatorBase() {}
  virtual value_type operator()(value_type x) const = 0;
};

class D1Interp : public InterpolatorBase, public boost::math::interpolators::pchip<ValueContainer>{
 public:
  using Real = InterpolatorBase::value_type;
  using Base = boost::math::interpolators::pchip<ValueContainer>;
  using Container = ValueContainer;
  D1Interp(Container&& abscissas, Container&& ordinates) : Base(std::move(abscissas),std::move(ordinates)) {}
  Real operator()(Real x) const override{
    return Base::operator()(x);
  }
};

class D2Interp : public InterpolatorBase, public boost::math::interpolators::cubic_hermite<ValueContainer>{
 public:
  using Real = InterpolatorBase::value_type;
  using Base = boost::math::interpolators::cubic_hermite<ValueContainer>;
  using Container = ValueContainer;
  D2Interp(Container&& abscissas, Container&& ordinates, Container&& derivatives) : Base(std::move(abscissas),std::move(ordinates),std::move(derivatives)) {}
  Real operator()(Real x) const override{
    return Base::operator()(x);
  }
};

class Interpolator{
 public:
  InterpolatorBase* _x[3];
  InterpolatorBase* _a[3];
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
    lon = item.blh[1];
    lat = item.blh[0];
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
    GDALUseTransformer(_reprojector, 0, 1, item.x+1, item.x, item.x+2, &valid);
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
  auto pitch = Eigen::AngleAxisd(a[1]*M_PI/180, Eigen::Vector3d::UnitY());
  auto roll = Eigen::AngleAxisd(a[0]*M_PI/180, Eigen::Vector3d::UnitX());
  auto yaw = Eigen::AngleAxisd(a[2]*M_PI/180, Eigen::Vector3d::UnitZ());
  return yaw*pitch*roll;
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
    double f, pixel, x0, y0;
    if (fscanf(fp, "%lf%lf%lf%lf", &f, &pixel, &x0, &y0) == 4) {
        SetFocalLength(f / pixel);
        SetPrincipalPoint(x0, y0);
        ret = true;
    }
    fclose(fp);
    return ret;
}

CameraMatrixType LinescanModel::CameraMatrix(double linenumber) const
{
  CameraMatrixType m = _camera->CameraMatrix();
  Eigen::Vector3d x;
  _pos->GetPosition(linenumber, x.data());
  auto pose = _pos->GetQuaternion(linenumber);
  Eigen::Matrix4d h = Eigen::Matrix4d::Identity();
  h.block<3, 3>(0,0) = pose.matrix().transpose();
  h.block<3, 1>(0,3) << -h.block<3, 3>(0, 0)*x;
  return m*h;
}

void LinescanModel::test(const char* gcppath, const char* prjpath){
  auto& posdata = _pos->_data;
  auto it1 = posdata.begin();
  ++it1;
  auto it2 = it1;//posdata.rbegin();
  ++it2;

  CameraMatrixType cam1 = _camera->CameraMatrix();
  Eigen::Vector3d c;
  _pos->GetPosition(it1->first, c.data());
  auto pose = _pos->GetQuaternion(it1->first);
  Eigen::Matrix4d h = Eigen::Matrix4d::Identity();
  h.block<3, 3>(0, 0) = pose.matrix().transpose();
  h.block<3, 1>(0, 3) << -h.block<3, 3>(0, 0) * c;
  auto y = h * c.homogeneous();
  auto z = cam1 * y;
  std::cout << z;
  

  Eigen::Matrix<double, 4, 1> x;
  x << it2->second.x[0], it2->second.x[1], it2->second.x[2] - 1500, 1;
  auto k = h*x;
  auto im = (cam1 * x).hnormalized();
  printf("(%lf,%lf)<->%d", im[0], im[1], it2->first);

  return;
  double gsd = std::fabs(it1->second.x[0]-it2->second.x[0])/((double)it2->first-it1->first);
  int x_range = 500;
  int x_center = 1024-300;

  int flag = std::abs(it1->second.a[2] - 90) > std::abs(it1->second.a[2] - 270) ? 1 : -1;

  FILE* fp = fopen(prjpath, "w");
  if (fp) {
      fprintf(fp, "%s", _pos->_utm_wkt);
      fclose(fp);
  }
  fp = fopen(gcppath, "w");
  if (fp) {
      fprintf(fp, "%d %d %lf %lf\n", x_center - x_range, it1->first, it1->second.x[0], it1->second.x[1] + flag * gsd * (-x_range));
      fprintf(fp, "%d %d %lf %lf\n", x_center + x_range, it1->first, it1->second.x[0], it1->second.x[1] + flag * gsd * x_range);
      fprintf(fp, "%d %d %lf %lf\n", x_center - x_range, it2->first, it2->second.x[0], it2->second.x[1] + flag * gsd * (-x_range));
      fprintf(fp, "%d %d %lf %lf\n", x_center + x_range, it2->first, it2->second.x[0], it2->second.x[1] + flag * gsd * x_range);
      fclose(fp);
  }
  
  return;
}


};
