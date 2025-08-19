#ifndef EXTRACTOR_HPP
#define EXTRACTOR_HPP

#include "raster/detail/lsd.h"
#include "raster/operator.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ios>
#include <opencv2/imgproc.hpp>
#include <boost/property_tree/ptree.hpp>

namespace xlingsky {
  namespace raster {
    namespace extractor {
      class PeakDetection : public FrameParallel {
        private:
          typedef float DataType;
          typedef cv::Point3_<DataType> Point;
          typedef std::pair<Point,Point> LineSegment;

          std::vector<DataType> _pdf;
          int _pdf_width;
          int _pdf_height;
          cv::Point2i _anchor;
          float _threshold;
          int _tolerance;

          int _length_min;
          int _length_max;

          std::vector<std::vector<LineSegment>> _lineseeds;
          std::vector<std::vector<Point> > _pointseeds;
          std::string _pointpath;
          std::string _linepath;

        public:
          PeakDetection( int bands, const boost::property_tree::ptree& config) {
            std::string line = config.get<std::string>("dimension", "5 5");
            std::istringstream iss(line);
            iss >> _pdf_width >> _pdf_height;
            if (iss.fail()) {
              _pdf_width = _pdf_height = 5;
            }
            iss.clear();

            _anchor.x = _pdf_width / 2;
            _anchor.y = _pdf_height / 2;

            _tolerance = std::max(_pdf_height, _pdf_width);

            _pdf.resize(_pdf_width * _pdf_height);
            line = config.get<std::string>(
                "matrix",
                "0.1 0.1 0.1 0.1 0.1 0.1   1   1   1 0.1 0.1   1 100   1 "
                "0.1 0.1   1   1   1 0.1 0.1 0.1 0.1 0.1 0.1");
            iss.str(line);
            for (int i = 0; i < _pdf_width * _pdf_height; ++i) {
              iss >> _pdf[i];
              if (iss.fail()) {
                _pdf[0] = _pdf[1] = _pdf[2] = _pdf[3] = _pdf[4] = _pdf[5] =
                    _pdf[9] = _pdf[10] = _pdf[14] = _pdf[15] = 0.1;
                _pdf[19] = _pdf[20] = _pdf[21] = _pdf[22] = _pdf[23] =
                    _pdf[24] = 0.1;
                _pdf[6] = _pdf[7] = _pdf[8] = _pdf[11] = _pdf[13] = _pdf[16] =
                    _pdf[17] = _pdf[18] = 1;
                _pdf[12] = 100;
                break;
              }
            }
            float &c = _pdf[_anchor.y * _pdf_width + _anchor.x];
            if (c != 0) {
              _threshold = c;
              c = 0;
            } else
              _threshold = std::numeric_limits<float>::lowest();

            _length_min = config.get<int>("length_min", 11);
            _length_max = config.get<int>("length_max", 20);

            _pointseeds.resize(bands);
            _lineseeds.resize(bands);
          }
          ~PeakDetection() {
            assert(!_pointpath.empty() && !_linepath.empty());
            std::ofstream ofs(_pointpath, std::ios::out | std::ios::trunc);
            if (ofs.good()) {
              ofs << _pointseeds.size() << "\n";
              ofs << std::fixed << std::setprecision(1);
              for (auto &seeds : _pointseeds) {
                ofs << seeds.size() << "\n";
                if (seeds.size() > 0) {
                  for (const auto &seed : seeds) {
                    ofs << seed.x << "\t" << seed.y << "\t" << seed.z << "\n";
                  }
                }
              }
              ofs.close();
            }
            ofs.open(_linepath, std::ios::out | std::ios::trunc);
            if (ofs.good()) {
              ofs << _lineseeds.size() << "\n";
              ofs << std::fixed << std::setprecision(1);
              for (const auto &bandseeds : _lineseeds) {
                ofs << bandseeds.size() << "\n";
                for (const auto &seed : bandseeds) {
                  ofs << seed.first.x << "\t" << seed.first.y << "\t" << seed.first.z << "\t"
                      << seed.second.x << "\t" << seed.second.y << "\t" << seed.second.z << "\n";
                }
              }
              ofs.close();
            }
          }
          void SetFilePath(const char *pointpath, const char *linepath) {
            _pointpath = pointpath;
            _linepath = linepath;
          }
          static int validate(const DataType *data, int rowstride,
                              const DataType *ratio, int xsize, int ysize,
                              const cv::Point2i &anchor, float threshold,
                              int tolerance) {
            DataType center = data[anchor.y * rowstride + anchor.x];
            if (center <= threshold)
              return -1;
            int count = 0;
            for (int y = 0; y < ysize; ++y) {
              int ro = y * rowstride;
              int rw = y * xsize;
              for (int x = 0; x < xsize; ++x) {
                if (ratio[rw + x] == 0)
                  continue; // skip empty pixels
                if (center * ratio[rw + x] <= data[ro + x]) {
                  ++count;
                  if (count > tolerance)
                    return -1;
                }
              }
            }
            return count;
          }
          bool operator()(int b, int xoff, int yoff, void *data, int cols,
                          int rows) override {
            std::vector<LineSegment> &bandlines = _lineseeds[b];
            std::vector<Point> &bandpoints = _pointseeds[b];

            DataType *pdata = (DataType *)data;
            cv::Mat mat(rows, cols, cv::DataType<DataType>::type, pdata);
            cv::Mat buffer(rows + _pdf_height, cols + _pdf_width,
                           cv::DataType<DataType>::type, cv::Scalar(0));
            mat.copyTo(buffer(cv::Rect(_anchor.x, _anchor.y, cols, rows)));
            for (int h = 0; h < rows; ++h) {
              for (int w = 0; w < cols; ++w) {
                int c = validate(buffer.ptr<DataType>(h, w), buffer.step1(),
                                 (DataType *)&_pdf[0], _pdf_width, _pdf_height,
                                 _anchor, _threshold, _tolerance);
                if (c < 0) {
                  pdata[h * cols + w] = 0;
                }else if(c==0){
                bandpoints.emplace_back(w + xoff, h + yoff,
                                      pdata[h * cols + w]);
                }
              }
            }

            cv::Mat binary = mat>0;
            std::vector<std::vector<cv::Point>> contours;
            std::vector<cv::Vec4i> hierarchy;
            cv::findContours(binary, contours, hierarchy, cv::RETR_EXTERNAL,
                             cv::CHAIN_APPROX_SIMPLE);
            for (auto &contour : contours) {
              if (contour.size() == 1){
                continue; 
              }
              double length = cv::arcLength(contour, true);
              std::vector<cv::Point> approx;
              cv::approxPolyDP(contour, approx, length * 0.01, true);
              if (approx.size() == 2) {
                double xx = std::abs(approx[0].x - approx[1].x) + 1;
                double yy = std::abs(approx[0].y - approx[1].y) + 1;
                double l = std::max(xx, yy);
                if ((_length_min >= 0 && l < _length_min) ||
                    (_length_max >= 0 && l > _length_max))
                  continue;
                bandlines.emplace_back(
                    Point(approx[0].x + xoff, approx[0].y + yoff,
                          mat.at<float>(approx[0])),
                    Point(approx[1].x + xoff, approx[1].y + yoff,
                          mat.at<float>(approx[1])));
              }
            }
            return true;
          }

      };
    
class LSDExtractor : public FrameParallel {
private:
  static const int dim = 7;
  double	_scale;
  double	_sigma_scale;
  double	_quant;
  double	_ang_th;
  double	_log_eps;
  double	_density_th;
  int _n_bins;
  int _length_min;
  int _length_max;
  int _thickness_min;
  int _thickness_max;

  typedef std::vector<std::array<double,dim>> Lines;
  std::vector<Lines> _seeds;
  std::string _filepath;

public:
  LSDExtractor(int bands, double scale, double sigma_scale, double quant,
               double ang_th, double log_eps, double density_th, int n_bins,
               int length_min, int length_max, int thickness_min,
               int thickness_max)
      : _scale(scale), _sigma_scale(sigma_scale), _quant(quant),
        _ang_th(ang_th), _log_eps(log_eps), _density_th(density_th),
        _n_bins(n_bins), _length_min(length_min), _length_max(length_max),
        _thickness_min(thickness_min), _thickness_max(thickness_max) {
    _seeds.resize(bands);
  }
  LSDExtractor(int bands, const boost::property_tree::ptree& line){
    _scale = line.get<double>("scale", 0.8);
    _sigma_scale = line.get<double>("sigma_scale", 0.6);
    _quant = line.get<double>("quant", 2.0);
    _ang_th = line.get<double>("ang_th", 22.5);
    _log_eps = line.get<double>("log_eps", 0);
    _density_th = line.get<double>("density_th", 0.7);
    _n_bins = line.get<int>("n_bins", 1024);
    _length_min = line.get<int>("length_min", 11);
    _length_max = line.get<int>("length_max", 20);
    _thickness_min = line.get<int>("thickness_min", 1);
    _thickness_max = line.get<int>("thickness_max", 1);
    _seeds.resize(bands);
  }
  ~LSDExtractor() {
    if (!_filepath.empty()) {
      std::ofstream ofs(_filepath, std::ios::out | std::ios::trunc);
      if (ofs.good()) {
        for (const auto &bandseeds : _seeds) {
          ofs << bandseeds.size() << "\n";
          for (const auto &seed : bandseeds) {
            ofs << seed[0] << " " << seed[1] << " "
                << seed[2] << " " << seed[3] << " "
                << seed[4] << " " << seed[5] << " "
                << seed[6] << "\n";
          }
        }
        ofs.close();
      }
    }
  }
  void SetFilePath(const char *filepath) { _filepath = filepath; }
  bool operator()(int b, int xoff, int yoff, void *data, int cols, int rows) override {
    double* segs;
    int n;
    segs = LineSegmentDetection(&n, (DataType*)data, cols, rows,
        _scale, _sigma_scale, _quant, _ang_th, _log_eps, _density_th, _n_bins, nullptr, nullptr, nullptr);

    Lines& bandseeds = _seeds[b];
    bandseeds.reserve(n);
    for (int i = 0; i < n; ++i) {
      double* seg = segs + i * dim;
      if ((_thickness_min>=0 && seg[4]<_thickness_min) || (_thickness_max>=0 && seg[4]>_thickness_max) ) continue;
      double xx = std::abs(seg[0]-seg[2])+1;
      double yy = std::abs(seg[1]-seg[3])+1;
      double l = std::max(xx, yy);
      if ( (_length_min>=0 && l<_length_min) || (_length_max>=0 && l>_length_max) ) continue;
      
      bandseeds.push_back({seg[0], seg[1], seg[2], seg[3], seg[4], seg[5], seg[6]});
    }

    free(segs);
    return true;
  }
};

    };
  };



};

#endif