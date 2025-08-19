#include "geometry/xsofa.h"
#include "geometry/ray2ellipsoid.hpp"

#include <sstream>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <boost/filesystem.hpp>
#include <gflags/gflags.h>
#include <glog/logging.h>

#include <ogr_spatialref.h>

DEFINE_string(eopfile, "EOP.txt", "EOP data file");
DEFINE_string(time, "", "time");  
DEFINE_bool(jd, false, "J2000 or UTC time");  
DEFINE_bool(geo, true, "convert to geodetic coordinates");
DEFINE_string(ee, "", "exterior elements: x, y, z, qx, qy, qz, qw"); 
DEFINE_string(ie, "100,0,0", "interior elements: fx, fy, cx, cy");
DEFINE_int32(rs, 0, "exterior rotation status: 0: no convert, 1: conjugate rotation, 2: opengl rotation");
DEFINE_double(eps, 0.01, "epsilon in latitude/longtitude");

void test(){
  OGRSpatialReference geocentricCS;
  // Set CGCS2000 geocentric coordinate system
  geocentricCS.importFromEPSG(4479); // EPSG:4479 - CGCS2000 geocentric
  OGRSpatialReference geographicCS;

  // Set CGCS2000 geographic coordinate system
  geographicCS.importFromEPSG(4490); // EPSG:4490 - CGCS2000 geographic

  // Create coordinate transformation
  OGRCoordinateTransformation *c2g =  OGRCreateCoordinateTransformation(&geocentricCS, &geographicCS);
  OGRCoordinateTransformation *g2c =  OGRCreateCoordinateTransformation(&geographicCS, &geocentricCS);

  Eigen::Vector3d loc( 34, 34, 0.0);
  Eigen::Vector3d sat( 0, 30, 36000000.0);

  g2c->Transform(1, &loc[0], &loc[1], &loc[2]);
  g2c->Transform(1, &sat[0], &sat[1], &sat[2]);

  Eigen::Vector3d origin = sat;
  Eigen::Vector3d direction = loc - sat;
  direction.normalize();

  double t;
  Eigen::Vector3d intersection;
  if (!s2r::utils::RayToEllipsoid(origin, direction, geocentricCS.GetSemiMajor(), geocentricCS.GetSemiMajor(), geocentricCS.GetSemiMinor(), &t)){
    std::cout << "ERROR!" << std::endl;
  } else{
    intersection = origin + direction * t;
    Eigen::Vector3d diff = intersection-loc;
    std::cout << "diff: " << diff << std::endl;
  }

  OCTDestroyCoordinateTransformation(c2g);
  OCTDestroyCoordinateTransformation(g2c);
}

inline int ValidateDoubleArray( const std::string& value, double* m, int sz){
  const char* p = value.c_str();
  if (*p != '-' && (*p < '0' || *p > '9') ) ++p;
  int n; int i = 0;
  char token;
  while ( i<sz ) {
    int t = sscanf(p, "%lf%c%n", m++, &token, &n);
    if(t>0) i++;
    if(t<2) break;
    p += n;
  }
  return i;
}

// Function to display tm struct
void printTimeStruct(const std::tm& tm) {
    std::cout << "Year: " << tm.tm_year + 1900 << "\n"
              << "Month: " << tm.tm_mon + 1 << "\n"
              << "Day: " << tm.tm_mday << "\n"
              << "Hour: " << tm.tm_hour << "\n"
              << "Minute: " << tm.tm_min << "\n"
              << "Second: " << tm.tm_sec << std::endl;
}

// Function to parse time from string in various formats
bool ParseTimeFromString(const std::string& timeStr, std::tm& timeStruct) {
    std::istringstream ss(timeStr);
    
    // Try different formats one by one
    if (ss >> std::get_time(&timeStruct, "%Y-%m-%d %H:%M:%S")) {
        return true;
    }
    ss.clear();
    ss.str(timeStr);
    
    if (ss >> std::get_time(&timeStruct, "%Y/%m/%d %H:%M:%S")) {
        return true;
    }
    ss.clear();
    ss.str(timeStr);
    
    if (ss >> std::get_time(&timeStruct, "%d-%m-%Y %H:%M:%S")) {
        return true;
    }
    ss.clear();
    ss.str(timeStr);
    
    if (ss >> std::get_time(&timeStruct, "%m/%d/%Y %H:%M:%S")) {
        return true;
    }
    ss.clear();
    ss.str(timeStr);
    
    if (ss >> std::get_time(&timeStruct, "%H:%M:%S")) {
        // If only time is provided, use current date
        std::time_t now = std::time(nullptr);
        std::tm* now_tm = std::localtime(&now);
        timeStruct.tm_year = now_tm->tm_year;
        timeStruct.tm_mon = now_tm->tm_mon;
        timeStruct.tm_mday = now_tm->tm_mday;
        return true;
    }
    
    ss.clear();
    ss.str(timeStr);
    double seconds;
    if (ss >> seconds) {
      int year, mon, day, hour, min;
      double sec;
      s2r::xsofa::jdTime2utc(seconds/24.0/3600.0, 
        &year, &mon, &day, 
        &hour, &min, &sec
       );
      timeStruct.tm_year = year - 1900;
      timeStruct.tm_mon = mon - 1;
      timeStruct.tm_mday = day;
      timeStruct.tm_hour = hour;
      timeStruct.tm_min = min;
      timeStruct.tm_sec = static_cast<int>(sec);
    }
    
    return false;
}

int main(int argc, char* argv[]) {
  FLAGS_logtostderr = 1;

  std::string usage("This program converts J2000 to CGCS2000. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + name + " <flags> <J2000 points>\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("0.1");

  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (argc < 2) {
    gflags::ShowUsageWithFlagsRestrict(argv[0], "j2cgcs");
    return 1;
  }

  if (!boost::filesystem::exists(FLAGS_eopfile)) {
    FLAGS_eopfile = boost::filesystem::path(argv[0]).parent_path().string() + "/" + FLAGS_eopfile;
  } 
  if (boost::filesystem::exists(FLAGS_eopfile)) {
    FLAGS_eopfile = boost::filesystem::absolute(FLAGS_eopfile).string();
  } else {
    LOG(ERROR) << "EOP file not found: " << FLAGS_eopfile;
    return 1;
  }

  bool flag_imgcoords = false;
  Eigen::Matrix3d k;
  Eigen::Matrix4d cam2j2000 = Eigen::Matrix4d::Identity();
  if (!FLAGS_ee.empty() && !FLAGS_ie.empty()) {
    double ie[4];
    int n = ValidateDoubleArray(FLAGS_ie, ie,4);
    switch (n) {
      case 3:
        k << ie[0], 0, ie[1],
             0, ie[0], ie[2],
             0, 0, 1;
        break;
      case 4:
        k << ie[0], 0, ie[2],
             0, ie[1], ie[3],
             0, 0, 1;
        break;
      default:
        LOG(ERROR) << "Invalid interior elements: " << FLAGS_ie;
        return 1;
    }
    double ee[7];
    n = ValidateDoubleArray(FLAGS_ee, ee, 7);
    if (n != 7) {
      LOG(ERROR) << "Invalid exterior elements: " << FLAGS_ee;
      return 1;
    }
    Eigen::Quaterniond q( ee[6], ee[3], ee[4], ee[5]);
    if (FLAGS_rs & 0x01) {
      q = q.conjugate();
    }
    cam2j2000.topLeftCorner<3,3>() = q.toRotationMatrix();
    cam2j2000.rightCols<1>() = Eigen::Vector4d(ee[0], ee[1], ee[2], 1);
    if (FLAGS_rs & 0x02) {
      cam2j2000.block<3,2>(0,1) = -cam2j2000.block<3,2>(0,1);
    }
    flag_imgcoords = true;

    VLOG(1) << "K: \n" << k;
    VLOG(2) << "Euler angles in roll, pitch, yaw (degree): \n"
            << Eigen::Quaterniond(cam2j2000.topLeftCorner<3, 3>())
                   .toRotationMatrix()
                   .eulerAngles(0, 1, 2).transpose() * 180.0 / M_PI;
    VLOG(1) << "cam2j2000: \n" << cam2j2000;
  }

  s2r::xsofa::Transformation sofa;
  if ( !sofa.Load(FLAGS_eopfile.c_str()) ) {
    LOG(ERROR) << "Failed to load EOP data from: " << FLAGS_eopfile;
    return 1;
  }
  VLOG(1) << "EOP file: " << FLAGS_eopfile;

  std::tm timeStruct = {};
  if( !ParseTimeFromString(FLAGS_time, timeStruct) ) {
    LOG(ERROR) << "Invalid time: " << FLAGS_time;
    return 1;
  }
  VLOG(1) << "UTC Time: " << std::put_time(&timeStruct, "%Y-%m-%d %H:%M:%S");

  OGRSpatialReference geocentricCS;
  // Set CGCS2000 geocentric coordinate system
  geocentricCS.importFromEPSG(4479); // EPSG:4479 - CGCS2000 geocentric

  OGRCoordinateTransformation *centric2graphic = nullptr;
  OGRCoordinateTransformation *graphic2centric = nullptr;
  if (FLAGS_geo) {
    OGRSpatialReference geographicCS;
    // Set CGCS2000 geographic coordinate system
    geographicCS.importFromEPSG(4490); // EPSG:4490 - CGCS2000 geographic
    // Create coordinate transformation
    centric2graphic =
        OGRCreateCoordinateTransformation(&geocentricCS, &geographicCS);
    graphic2centric = 
        OGRCreateCoordinateTransformation(&geographicCS, &geocentricCS);
        // centric2graphic->GetInverse();
  }

  Eigen::Matrix3d rotation;
  {
    double r[3][3];
    sofa.ComputeRotation(timeStruct.tm_year + 1900, timeStruct.tm_mon + 1,
                         timeStruct.tm_mday, timeStruct.tm_hour,
                         timeStruct.tm_min, timeStruct.tm_sec, r);

    rotation << r[0][0], r[0][1], r[0][2],
                r[1][0], r[1][1], r[1][2],
                r[2][0], r[2][1], r[2][2];
  }
  VLOG(1) << "J2000 to CGCS2000 matrix: \n" << rotation;

  if (flag_imgcoords) {
    Eigen::Matrix4d cam2cgcs = Eigen::Matrix4d::Identity();
    cam2cgcs.topLeftCorner<3,3>() = rotation*cam2j2000.topLeftCorner<3,3>();
    cam2cgcs.topRightCorner<3,1>() = rotation*cam2j2000.topRightCorner<3,1>();
    VLOG(2) << "cam2cgcs: \n" << cam2cgcs;
    Eigen::Vector3d rayorigin = cam2cgcs.topRightCorner<3,1>();
    VLOG(1) << "rayorigin: " << rayorigin.transpose();

    for (int i = 1; i < argc; ++i) {
      double xyz[3] = {0, 0, 0};
      int n = ValidateDoubleArray(argv[i], xyz, 3);
      if ( n < 2) {
        LOG(ERROR) << "Invalid image coordinates: " << argv[i];
        continue;
      }
      Eigen::Vector3d x(xyz[0], xyz[1], 1);
      Eigen::Vector3d raydirection = cam2cgcs.topLeftCorner<3, 3>() * k.inverse() * x;

      VLOG(1) << "raydirection: " << raydirection.transpose();

      std::cout << std::fixed << std::setprecision(2);
      std::cout << "[" << xyz[0] << ", " << xyz[1] << "] \t";

      double t;
      if (!s2r::utils::RayToEllipsoid(rayorigin, raydirection, 
        geocentricCS.GetSemiMajor()+xyz[2], geocentricCS.GetSemiMajor()+xyz[2], geocentricCS.GetSemiMinor()+xyz[2], &t)){
        std::cout << "No intersection" << std::endl;
        continue;
      }
      Eigen::Vector3d intersection;
      intersection = rayorigin + raydirection * t;

      // if (n == 3 && centric2graphic) {
      //   centric2graphic->Transform(1, &intersection[0], &intersection[1], &intersection[2]);
      //   Eigen::Vector3d p1, p2, p3;
      //   p1 << intersection[0], intersection[1], xyz[2];
      //   p2 << intersection[0] + FLAGS_eps, intersection[1], xyz[2];
      //   p3 << intersection[0], intersection[1] + FLAGS_eps, xyz[2];
      //   graphic2centric->Transform(1, &p1[0], &p1[1], &p1[2]);
      //   graphic2centric->Transform(1, &p2[0], &p2[1], &p2[2]);
      //   graphic2centric->Transform(1, &p3[0], &p3[1], &p3[2]);
      //   Eigen::Vector3d v1 = p2 - p1;
      //   Eigen::Vector3d v2 = p3 - p1;
      //   Eigen::Vector3d normal = v1.cross(v2);
      //   normal.normalize();
      //   Eigen::Vector4d plane;
      //   plane << normal[0], normal[1], normal[2], -normal.dot(p1);

      //   if (!s2r::utils::RayToPlane(rayorigin, raydirection, plane, &t)) {
      //     std::cout << "No intersection" << std::endl;
      //     continue;
      //   }
      //   intersection = rayorigin + raydirection * t;
      // }

      std::cout << intersection[0] << ", " << intersection[1] << ", " << intersection[2];

      if (centric2graphic) {
        centric2graphic->Transform(1, &intersection[0], &intersection[1], &intersection[2]);
        // lat lon height
        std::cout << std::fixed << std::setprecision(6) << "\t" << intersection[0] << ", " << intersection[1] << ", "
                  << intersection[2];
      }
      std::cout << std::endl;
    }
  }else {
    for (int i = 1; i < argc; ++i) {
      double xyz[3];
      if (ValidateDoubleArray(argv[i], xyz, 3) < 3) {
        LOG(ERROR) << "Invalid J2000 xyz points: " << argv[i];
        continue;
      }
      Eigen::Vector3d j2000(xyz[0], xyz[1], xyz[2]);
      Eigen::Vector3d result = rotation * j2000;

      std::cout << std::fixed << std::setprecision(2);
      std::cout << result[0] << ", " << result[1] << ", " << result[2];

      if (centric2graphic) {
        centric2graphic->Transform(1, &result[0], &result[1], &result[2]);
        // lat lon height
        std::cout << std::fixed << std::setprecision(6) << "\t" << result[0] << ", " << result[1] << ", "
                  << result[2];
      }
      std::cout << std::endl;
    }
  }

  if (centric2graphic) {
    OCTDestroyCoordinateTransformation(centric2graphic);
  }
  if (graphic2centric) {
    OCTDestroyCoordinateTransformation(graphic2centric);
  }

  return 0;
}