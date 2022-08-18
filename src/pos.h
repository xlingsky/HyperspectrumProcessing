#ifndef POS_H
#define POS_H

#include <vector>
#include <map>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

namespace HSP{

class LinescanModel;

class Pos{
 public:
  friend class LinescanModel;
  struct record{
    unsigned int frame;
    unsigned int week;
    double time;
    double blh[3];
    double x[3];
    double vx[3];
    double a[3];
  };

 protected:
  std::map<int, record> _data;
  void* _reprojector;
  void* _interpolator;
  char* _utm_wkt;
  double _offset[3];
 public:
     Pos() : _reprojector(nullptr), _interpolator(nullptr), _utm_wkt(nullptr) {
         _offset[0] = _offset[1] = _offset[2] = 0;
  }
  ~Pos();
  bool load(const char* filepath);
  size_t size() const { return _data.size(); }
  record& front() { return _data.begin()->second; }
  record& back() { return  _data.rbegin()->second; }
  void test();
  int GetPosition(double lineid, double x[3]);
  int GetAngle(double lineid, double a[3]);
  int GetPos(double lineid, double x[3], double a[3]);
  Eigen::Quaterniond GetQuaternion(double lineid);
 private:
  void Reprojection();
  void BuildInterpolator();
  int Check();
  bool LoadAux(const char* filepath);
};

typedef Eigen::Matrix<double, 3, 4> CameraMatrixType;

Eigen::Matrix<double, 4, Eigen::Dynamic> BackProjectionToPlane(const CameraMatrixType& camera, const Eigen::Matrix<double,3,Eigen::Dynamic>& im, const Eigen::Matrix<double,4,3>& plane_transformation);
Eigen::Matrix<double,4,Eigen::Dynamic> BackProjectionToPlane(const CameraMatrixType& camera,const Eigen::Matrix<double,3,Eigen::Dynamic>& im, const Eigen::Matrix<double,4,1>& plane_equation);

class PinholeCamera{
 protected:
  Eigen::Matrix3d _intrinsic;
  Eigen::Quaterniond _pose;
  Eigen::Vector3d _translation;
 public:
  PinholeCamera() : _intrinsic(Eigen::Matrix3d::Identity()), _pose(Eigen::Quaterniond::Identity()), _translation(Eigen::Vector3d::Zero()) {};
  bool Load(const char* filepath);
  CameraMatrixType CameraMatrix() const{
    CameraMatrixType m;
    m.block<3,3>(0,0).noalias() = _intrinsic*_pose.matrix();
    m.block<3,1>(0,3).noalias() = -m.block<3,3>(0,0)*_translation;
    return m;
  }
  void SetFocalLength(double f){
    _intrinsic(0,0) = _intrinsic(1,1) = f;
  }
  void SetPrincipalPoint(double x, double y){
    _intrinsic(0, 2) = x;
    _intrinsic(1, 2) = y;
  }
};

class LinescanModel{
 protected:
  PinholeCamera* _camera;
  Pos* _pos;
 public:
  LinescanModel() : _camera(nullptr), _pos(nullptr) {
  }
  ~LinescanModel(){
  }
  void SetPos(Pos* pos){
    _pos = pos;
  }
  void SetCamera(PinholeCamera* cam){
    _camera = cam;
  }
  CameraMatrixType CameraMatrix(double linenumber) const;
  void test( const char* gcppath, const char* prjpath);
};

};

#endif
