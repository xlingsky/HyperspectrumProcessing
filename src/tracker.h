#ifndef TRACKER_H
#define TRACKER_H

#include <cstddef>
#include <opencv2/core.hpp>
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <opencv2/video/tracking.hpp>
#include <opencv2/imgproc.hpp>

namespace xlingsky {

    struct KalmanParams {
        cv::Mat transition_matrix;
        float process_noise = 1e-2f;
        float measurement_noise = 1e-1f;
        float error_cov_post = 1.0f;
        int tolerance = 30;
        float search_radius = 10.0f;
        bool set_transition_matrix(std::string matrix_str) {
          std::istringstream iss(matrix_str);
          std::vector<float> values(16);
          for (int i = 0; i < 16; ++i) {
            iss >> values[i];
            if (iss.fail()) {
              return false;
            }
          }
          transition_matrix = cv::Mat(4, 4, CV_32F, values.data()).clone();
          return true;
        }
        bool load(const boost::property_tree::ptree& pt) {
          transition_matrix = (cv::Mat_<float>(4, 4) <<
              1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1);

          auto kf = pt.get_child("kalman");
          // Transition matrix (constant velocity model)
          std::string transition_matrix = kf.get<std::string>(
              "transition_matrix", "1 0 1 0 0 1 0 1 0 0 1 0 0 0 0 1");
          set_transition_matrix(transition_matrix);
          process_noise = kf.get<float>("process_noise", 1e-1f);
          measurement_noise = kf.get<float>("measurement_noise", 1e-1f);
          error_cov_post = kf.get<float>("error_cov_post", 1.0f);

          tolerance = kf.get<int>("max_missing_frames", 30);
          search_radius = kf.get<float>("search_radius", 2.0f);

          return true;
        }
        KalmanParams(const boost::property_tree::ptree& pt) {
            load(pt);
        }
    };

template<class Point>
class KalmanTracker : cv::KalmanFilter {
public:
  std::vector<int> _point_ids;
  std::vector<Point> _points;
  int _frame_start;

private:
  int _missings;
  int _tolerance;

  float _search_radius;

public:
  KalmanTracker(int frame_start, int ptid, const Point& pt, KalmanParams& params) : _frame_start(frame_start), _missings(0) {
    // State vector: [x, y, vx, vy] - position (x,y) and velocity (vx,vy)
    const int state_size = 4;
    // Measurement vector: [x, y] - we can only observe position
    const int measurement_size = 2;
    init(state_size, measurement_size, 0); 

    // Transition Matrix (A): Describes how the state changes from one time step to the next.
    // [1 0 dt 0]
    // [0 1 0 dt]
    // [0 0 1 0 ]
    // [0 0 0 1 ]
    // We assume constant velocity, so dt=1 for simplicity.
    transitionMatrix = params.transition_matrix;

    // Measurement Matrix (H): Maps the state to the measurement.
    // We can only measure position, not velocity.
    // [1 0 0 0]
    // [0 1 0 0]
    measurementMatrix = (cv::Mat_<float>(2, 4) << 1, 0, 0, 0, 0, 1, 0, 0);

    // Process Noise Covariance Matrix (Q): Represents uncertainty in our motion model.
    // A small value here means we trust our constant velocity model.
    cv::setIdentity(processNoiseCov, cv::Scalar::all(params.process_noise));

    // Measurement Noise Covariance Matrix (R): Represents uncertainty in our measurements.
    // A larger value here means we have a noisy sensor.
    cv::setIdentity(measurementNoiseCov, cv::Scalar::all(params.measurement_noise));

    // Error Covariance Matrix (P): Represents our initial uncertainty about the state.
    // Start with high uncertainty.
    cv::setIdentity(errorCovPost, cv::Scalar::all(params.error_cov_post));

    // Initial State Posterior (x'): Our initial guess of the state.
    // We'll set it with the first measurement.
    statePost = (cv::Mat_<float>(state_size, 1) << pt.x, pt.y, 0, 0);

    _tolerance = params.tolerance;
    _search_radius = params.search_radius;

    _point_ids.push_back(ptid);
    _points.push_back(pt);
  }
//   const Point& operator[](const size_t idx) const {
//     return _points[idx];
//   }

//   const Point& last() const {
//     return _points.back();
//   }

  bool good() const { 
    return _missings <= _tolerance;
  }

  float search_radius() const {
    return _missings==0?_search_radius:_search_radius*2;
  }

  void update(int id, const Point& pt) {
    _point_ids.push_back(id);
    _points.push_back(pt);
    cv::Mat m = (cv::Mat_<float>(2, 1) << pt.x, pt.y);
    correct(m);
    _missings = 0;
  }
  void missing(const Point& pt) {
    _point_ids.push_back(-1);
    _points.push_back(pt);
    cv::Mat m = (cv::Mat_<float>(2, 1) << pt.x, pt.y);
    correct(m);
    _missings++;
  }
  void remove_all_missings() {
    _point_ids.resize(_point_ids.size() - _missings);
    _points.resize(_points.size() - _missings);
    _missings = 0;
  }
  Point estimate(){
    // if (_point_ids.size()<=10)
    //   return _points.back();
    cv::Mat prediction = predict();
    return Point(prediction.at<float>(0), prediction.at<float>(1));
  }

  bool is_valid(int min_frame_number, float min_speed, float max_acceleration) const {
    if (_points.size() < min_frame_number) return false;
    float distance = 0.0f;
    std::vector<float> speeds(_points.size()-1);
    float displacement = 0;
    int min_displacement = 0;
    for (size_t i = 1; i < _points.size(); ++i) {
      auto& pt1 = _points[i - 1];
      auto& pt2 = _points[i];
      auto& v = speeds[i-1];
      auto dx = pt2.x - pt1.x;
      auto dy = pt2.y - pt1.y;
      v = std::abs(dx)>std::abs(dy)?dx:dy;
      distance += std::abs(v);
      if ( (distance-i*min_speed) < -1 ) return false;
      // displacement = std::max(std::abs(pt2.x-_points[0].x), std::abs(pt2.y-_points[0].y));
      // if (distance>_search_radius) {
      //   ++min_displacement; 
      //   if (min_displacement > max_acceleration) {
      //     return false;
      //   }
      // }
    }
    if (_points.size() > 5){
      cv::Mat data(1,speeds.size(), CV_32F, &speeds[0]);
      double minVal, maxVal;
      cv::minMaxIdx(data, &minVal, &maxVal);
      if (minVal < 0 && maxVal > 0){
        cv::Mat hist;
        int histSize = std::max(5, (int)std::ceil(maxVal-minVal));
        float range[] = {(float)minVal, (float)maxVal};
        const float* histRange = {range};
        cv::calcHist(&data, 1, 0, cv::Mat(), hist, 1, &histSize, &histRange, true, false);
        int id0 = int(-histSize*minVal/(maxVal-minVal));
        if (hist.at<float>(id0) > max_acceleration*_points.size()) return false;
        if (hist.at<float>(0)+hist.at<float>(histSize-1)>max_acceleration*_points.size()) return false;
      }
    }
    return true;
  }

};

};

#endif