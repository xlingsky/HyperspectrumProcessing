#include <cassert>
#include <cstddef>
#include <cstdio>
#include <ostream>
#include <string>
#include <fstream>

#include <opencv2/core.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <nanoflann.hpp>

#include "tracker.h"

DEFINE_string(sp, "", "point seeds for tracking");
DEFINE_string(sl, "NA", "line seeds for tracking");
DEFINE_string(o, "", "output prefix");
DEFINE_string(config, "config.xml", "configuration file");
DEFINE_bool(ow, false, "overwriting mode");

typedef float DataType;
typedef cv::Point2f Point2;
typedef cv::Point3f Point;
typedef std::vector<Point> PointSeeds;
typedef std::pair<Point, Point> Line;
typedef std::vector<Line> LineSeeds;

std::ostream& operator<<(std::ostream& os, const Point& pt) {
  os << pt.x << " " << pt.y << " " << pt.z;
  return os;
}
std::ostream& operator<<(std::ostream& os, const Line& obj) {
  os << obj.first << " " << obj.second;
  return os;
}

std::istream& operator>>(std::istream& is, Point& pt) {
  is >> pt.x >> pt.y >> pt.z;
  return is;
}

std::istream& operator>>(std::istream& is, Line& line) {
  is >> line.first >> line.second;
  return is;
}

template<class Node>
std::istream& operator>>(std::istream& is, std::vector<Node>& seeds) {
  int n;
  is >> n;
  seeds.resize(n);
  for (int i = 0; i < n; ++i) {
    is >> seeds[i];
  }
  return is;
}

template<class _PointSet>
struct PointKNN{
  _PointSet& seeds_;
  PointKNN(_PointSet& seeds) : seeds_(seeds) {}
  inline size_t kdtree_get_point_count() const {
    return seeds_.size();
  }
  inline int kdtree_get_pt(const size_t idx, int dim) const {
    if (dim == 0) {
      return seeds_[idx].x;
    } else {
      return seeds_[idx].y;
    }
  }
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const {
    return false; // no bounding box search
  }
};

template<typename ObjectList>
bool read_seeds_file(const std::string& filename, std::vector<ObjectList>& seeds) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        return false;
    }
    int n;
    fin >> n;
    seeds.resize(n);
    for (int i = 0; i < n; ++i) {
      fin >> seeds[i];
    }

    return true;
}

template<class Point, class Tracker, class Params, int NUM=1>
void tracking(std::vector<Point>& seeds, std::vector<Tracker>& trackers, int frameid, Params& params) {
  struct _Point{
    int id;
    DataType x, y;
    _Point(int id, DataType x, DataType y) : id(id), x(x), y(y) {}
  };
  std::vector<_Point> predicted_points;
  predicted_points.reserve(trackers.size());
  for (size_t i=0; i<trackers.size(); ++i) {
    auto& tracker = trackers[i];
    if (!tracker.good()){
      continue;
    }
    auto pt = tracker.estimate();
    predicted_points.push_back(_Point(i, pt.x, pt.y));
  }

  if (seeds.size() == 0) {
    for (auto& predicted : predicted_points) {
      auto &tracker = trackers[predicted.id];
      tracker.missing({predicted.x, predicted.y});
    }
    return;
  }

   // Construct the KD-tree
   typedef std::vector<Point> PointSet;
    using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<DataType, PointKNN<PointSet>>,
        PointKNN<PointSet>,
        2 /* 2D space */
    >;

    KDTree index(2 /*dim*/, seeds, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    index.buildIndex();
    int n = seeds.size()-predicted_points.size();
    if (n>0) 
      trackers.reserve(trackers.size()+n);

    std::vector<int> status(seeds.size(), 0);
    for (auto& predicted : predicted_points) {

      size_t ret_index[NUM];
      DataType out_dist_sqr[NUM];
      nanoflann::KNNResultSet<DataType> result_set(NUM);
      result_set.init(ret_index, out_dist_sqr);

      index.findNeighbors(result_set,
                          &predicted.x,                // Query point
                          nanoflann::SearchParameters() // Search params
      );
      auto& tracker = trackers[predicted.id];
      if ( status[ret_index[0]] == 0 && out_dist_sqr[0] < tracker.search_radius() * tracker.search_radius()) {
        tracker.update( ret_index[0], {seeds[ret_index[0]].x, seeds[ret_index[0]].y});
        status[ret_index[0]] = 1; // Mark as used
      } else {
        tracker.missing({predicted.x, predicted.y});
      }
    }
    for (size_t i=0; i<seeds.size(); ++i) {
      if (status[i] == 0) {
        trackers.push_back(Tracker(frameid, i, {seeds[i].x, seeds[i].y}, params));
      }
    }
}

template<class Node>
struct Trajectory2D{
    int frame_st;
    std::vector<Node> nodes;
    Trajectory2D() : frame_st(-1) {}
};

struct PointEx : public Point {
  int status;
  PointEx(float x, float y, float z, int status = 0) : Point(x, y, z), status(status) {}
  PointEx(const Point& pt, int status = 0) : Point(pt), status(status) {}
};
struct LineEx : public Line {
  int status;
  LineEx(const Line& line, int status = 0) : Line(line), status(status) {}
};

typedef Trajectory2D<PointEx> PointTrajectory;
typedef std::vector<PointTrajectory> PointTrajectories;
typedef Trajectory2D<LineEx> LineTrajectory;
typedef std::vector<LineTrajectory> LineTrajectories;

std::ostream& operator<<(std::ostream& os, const PointEx& obj) {
  os << *((Point*)&obj) << " " << obj.status;
  return os;
}

std::ostream& operator<<(std::ostream& os, const LineEx& obj) {
  os << *((Line*)&obj) << " " << obj.status;
  return os;
}

template<class Node>
std::ostream& operator<<(std::ostream& os, const std::vector<Node>& obj) {
  os << obj.size();
  for (const auto& node : obj) {
    os << "\n" << node ;
  }
  if (obj.size()==0) os << "\n";
  return os;
}

template<class Node>
std::ostream& operator<<(std::ostream& os, const Trajectory2D<Node>& obj) {
  os << obj.frame_st << "\n" << obj.nodes;
  return os;
}

int main(int argc, char* argv[]) {
  FLAGS_logtostderr = 1;
  std::string usage("This program tracks objects from frames. Usage:\n");
  {
    std::string name = boost::filesystem::path(argv[0]).filename().string();
    usage = usage + name + " <flags> <frames>\n";
  }
  gflags::SetUsageMessage(usage);
  gflags::SetVersionString("0.3");

  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (argc < 1) {
    gflags::ShowUsageWithFlagsRestrict(argv[0], "tracker");
    return 1;
  }

  std::vector<PointSeeds> pointseeds;
  std::vector<LineSeeds> lineseeds;

  if (FLAGS_sp != "NA") {
    if ( boost::filesystem::exists(boost::filesystem::path(FLAGS_sp)) && !read_seeds_file(FLAGS_sp, pointseeds)) {
      LOG(ERROR) << "Failed to read point seeds from file: " << FLAGS_sp;
      return 3;
    }
  } 
  if (FLAGS_sl != "NA") {
    if ( boost::filesystem::exists(boost::filesystem::path(FLAGS_sl)) && !read_seeds_file(FLAGS_sl, lineseeds)) {
      LOG(ERROR) << "Failed to read line seeds from file: " << FLAGS_sl;
      return 3;
    }
  }

  if (FLAGS_o.empty()) {
    auto path = boost::filesystem::path(pointseeds.empty() ? FLAGS_sl : FLAGS_sp);
    path.replace_extension();
    FLAGS_o = path.string();
  }else{
    auto dir = boost::filesystem::path(FLAGS_o).parent_path();
    if ( !boost::filesystem::exists(dir) && !boost::filesystem::create_directories(dir) ) {
      LOG(ERROR) << "Failed to create directory: " << dir.string();
      return 3;
    }
  }

  if (!boost::filesystem::exists(FLAGS_config)) {
    FLAGS_config = boost::filesystem::path(argv[0]).parent_path().string() + "/" + FLAGS_config;
  } 
  if (boost::filesystem::exists(FLAGS_config)) {
    FLAGS_config = boost::filesystem::absolute(FLAGS_config).string();
  } else {
    LOG(ERROR) << "config not found: " << FLAGS_config;
    return 1;
  }
  boost::property_tree::ptree config;
  try {
    boost::property_tree::read_xml(FLAGS_config, config);
    config = config.get_child("HSP");
  } catch (const std::exception &e) {
    LOG(ERROR) << "Failed to read config file: " << FLAGS_config << ", error: " << e.what();
    return 1;
  }

  int frames_point = pointseeds.size();
  int frames_line = lineseeds.size();

  if (frames_point==0&&frames_line==0) {
    LOG(ERROR) << "NO data in point seed file and line seed file!";
    return 3;
  }
  if ((frames_point>0&&frames_line>0)&&frames_line!=frames_point){
    LOG(ERROR) << "MISMATCHING dimemsion between point seeds and line seeds!";
    return 3;
  }

  int frames = std::max(frames_point, frames_line);

  int total_points = 0;
  int total_lines = 0;
  {
    char message[512] = "Frame ";
    char* p = message + strlen(message);
    for (int i=0; i<frames; ++i) {
      snprintf(p, 100, "%d: ", i+1);
      if (pointseeds.size() > 0){
        total_points += pointseeds[i].size();
        snprintf(p+strlen(p), 100, "%ld points ", pointseeds[i].size());
      }
      if (lineseeds.size() > 0){
        total_lines += lineseeds[i].size();
        snprintf(p+strlen(p), 100, "%ld lines ", lineseeds[i].size());
      }
      VLOG(2) << message;
    }
    VLOG(1) << "Total points: " << total_points
             << ", Total lines: " << total_lines;
  }
  if (total_lines==0 && total_points==0) {
    LOG(ERROR) << "MISSING point seeds or line seeds!";
    return 3;
  }

  xlingsky::KalmanParams params(config);

  std::vector< xlingsky::KalmanTracker<Point2>> point_trackers, line_trackers;
  if (pointseeds.size()>0) {
    for (int i = 0; i < frames; ++i) {
      tracking(pointseeds[i], point_trackers, i, params);
    }
  }
  if (lineseeds.size()>0) {
    for (int i = 0; i < frames; ++i) {
      std::vector<Point2> seeds;
      for (auto& seed : lineseeds[i]) {
        auto pt = (seed.first + seed.second) / 2;
        seeds.push_back(Point2(pt.x, pt.y));
      }
      tracking(seeds, line_trackers, i,  params);
    }
  }
  for (auto& tracker  : point_trackers) {
    tracker.remove_all_missings();
  }
  for (auto& tracker   : line_trackers) {
    tracker.remove_all_missings();
  }

  int min_frame_number = config.get<int>("min_frame_number", 10);
  float min_speed = config.get<float>("min_speed", 0.1f);
  float max_acceleration = config.get<float>("max_acceleration", 1);

  PointTrajectories point_trajectories;
  point_trajectories.reserve(point_trackers.size());
  for (auto& tracker : point_trackers) {
    if (tracker.is_valid(min_frame_number, min_speed, max_acceleration)) {
      PointTrajectory trajectory;
      trajectory.frame_st = tracker._frame_start;
      for (size_t i=0; i<tracker._point_ids.size(); ++i) {
        auto& seeds = pointseeds[tracker._frame_start+i];
        if (tracker._point_ids[i] >= 0) {
          trajectory.nodes.emplace_back(seeds[tracker._point_ids[i]], 1);
        }else{
          trajectory.nodes.emplace_back(Point(tracker._points[i].x, tracker._points[i].y, 0), 0);
        }
      }
      point_trajectories.push_back(std::move(trajectory));
    }
  }

  LineTrajectories line_trajectories;
  line_trajectories.reserve(line_trackers.size());
  for (auto& tracker : line_trackers) {
    if (tracker.is_valid(min_frame_number, min_speed, max_acceleration)) {
      LineTrajectory trajectory;
      trajectory.frame_st = tracker._frame_start;
      for (size_t i=0; i<tracker._point_ids.size(); ++i) {
        auto& seeds = lineseeds[tracker._frame_start+i];
        if (tracker._point_ids[i] >= 0) {
          trajectory.nodes.emplace_back(seeds[tracker._point_ids[i]], 1);
        }else{
          Point pt(tracker._points[i].x, tracker._points[i].y, 0);
          trajectory.nodes.emplace_back(std::make_pair(pt, pt), 0);
        }
      }
      line_trajectories.push_back(std::move(trajectory));
    }
  }

  VLOG(1) << "Total point trajectories: " << point_trajectories.size()
             << ", Total line trajectories: " << line_trajectories.size();

  std::ofstream dstfile;
  if (pointseeds.size()>0) {
    std::string dstpath = FLAGS_o + "_point_trajectories.txt";
    dstfile.open(dstpath, std::ios::out|std::ios::trunc);
    if(dstfile.is_open()){
      dstfile << point_trajectories;
      dstfile.close();
    }
  }
  if (lineseeds.size()>0) {
    std::string dstpath = FLAGS_o + "_line_trajectories.txt";
    dstfile.open(dstpath, std::ios::out|std::ios::trunc);
    if(dstfile.is_open()){
      dstfile << line_trajectories;
      dstfile.close();
    }
  }

  return 0;
}

#ifdef _USE_OPENCV_TRACKER
  std::vector<Tracker> trackers;
  std::vector<float> data((size_t)width*height);
  for (int i=1; i<=frames; ++i) {
      VLOG(1) << "Tracking frame " << i << " ...";
      if (dataset->GetRasterBand(i)->RasterIO(GF_Read, 0, 0, width, height, data.data(), width, height, GDT_Float32, 0, 0) != CE_None) {
          LOG(ERROR) << "Failed to read frame " << i;
          GDALClose(dataset);
          return 4;
      }
      VLOG(1) << "\tUpdating existing #" << trackers.size() << " trackers";
      for (auto& tracker : trackers) {
          if (tracker.trajectory.ed>0 && tracker.trajectory.ed < i) continue;
          cv::Rect rc;
          try {
            if (tracker.tracker->update(cv::Mat(height, width, CV_32FC1, data.data()), rc)) {
                tracker.trajectory.bboxes.push_back(rc);
            }else  if (tracker.trajectory.ed <= -FLAGS_mmf) {
                  tracker.trajectory.ed += i-1;
            }else{
              tracker.trajectory.bboxes.push_back(cv::Rect(0,0,0,0));
              --tracker.trajectory.ed;
            }
          } catch (cv::Exception& e) {
            if (tracker.trajectory.ed <= -FLAGS_mmf) {
                  tracker.trajectory.ed += i-1;
            }else{
              tracker.trajectory.bboxes.push_back(cv::Rect(0,0,0,0));
              --tracker.trajectory.ed;
            }
          }
      }
      if (seeds.find(i) != seeds.end()) {
          const Seeds& frame_seeds = seeds[i];
          VLOG(1) << "\tAdding #" << frame_seeds.size() << " candidate trackers";
          for (const auto& seed : frame_seeds) {
              if (contain(trackers, i, seed)) {
                  continue;
              }
              Tracker tracker;
              if (tracker.init(FLAGS_algo)) {
                tracker.tracker->init(
                    cv::Mat(height, width, CV_32FC1, data.data()), seed);
                tracker.trajectory.st = i;
                tracker.trajectory.bboxes.push_back(seed);
                trackers.push_back(std::move(tracker));
              } else {
                  LOG(ERROR) << "Failed to create tracker: " << FLAGS_algo;
              }
          }
      }
  }
  #endif
