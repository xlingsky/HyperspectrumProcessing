#include "xsofa.h"
#include "SOFA/sofa.h"
#include "SOFA/sofam.h"
#include <utility>
#include <cmath> 
#include <fstream>
#include <sstream>

namespace s2r {
  namespace xsofa {

void rotationECI2ECEF(double djmjd0, double tt, double date, double tut,
                      double xp, double yp, double r[3][3]){

  double j2000_to_GCRS[3][3], rc_to_it[3][3];

  double rb[3][3], rp[3][3], rbp[3][3];
  iauBp06(djmjd0, tt, rb, rp, rbp);
  iauTr(rb, j2000_to_GCRS);

  double x, y, s;
  iauXys06a(djmjd0, tt, &x, &y, &s);
  iauC2ixys(x, y, s, rc_to_it);

  double era = iauEra00(djmjd0 + date, tut);
  double rc2ti[3][3];
  iauCr(rc_to_it, rc2ti);
  iauRz(era, rc2ti);

  double rpom[3][3];
  iauPom00(xp * DAS2R, yp * DAS2R, iauSp00(djmjd0, tt), rpom);
  iauRxr(rpom, rc2ti, rc_to_it);

  iauRxr(rc_to_it, j2000_to_GCRS, r);
}

// Function to convert UTC time to Julian Date
double utc2jdTime(int year, int month, int day, 
                       int hour, int minute, double second) {
    // Handle January and February as months 13 and 14 of previous year
    if (month <= 2) {
        year -= 1;
        month += 12;
    }

    int A = year / 100;
    int B = 2 - A + (A / 4);

    // Calculate Julian Date
    double jd = floor(365.25 * (year + 4716)) + 
                floor(30.6001 * (month + 1)) + 
                day + B - 1524.5;

    // Add time of day
    double timeFraction = (hour + minute / 60.0 + second / 3600.0) / 24.0;
    jd += timeFraction;

    return jd;
}

void jdTime2utc(double jdTime, int* year, int* month, int* day,
    int* hour, int* minute, double* second){
      double fd;
      iauJd2cal(2451545.0, jdTime, year, month, day, &fd);
      *hour = static_cast<int>(fd * 24);
      *minute = static_cast<int>((fd * 24 - *hour) * 60);
      *second = ((fd * 24 - *hour) * 60 - *minute) * 60;
}

void time2jd(int year, int month, int day, int hour, int minute, double second, int eopDAT, double eopdut,
    double *djm0, double *tt, double *date, double *tut){
  iauCal2jd(year, month, day, djm0, date);
  *tt = *date + (hour * 3600 + minute * 60 + second) / 86400.0 +
              (eopDAT + TTMTAI) / 86400.0;
  *tut = (hour * 3600 + minute * 60 + second + eopdut) / 86400.0;
}

const std::vector<std::pair<int, int>> leapSeconds = {
    {44239, 19}, // 1980-01-01
    {44786, 20}, // 1981-07-01
    {45151, 21}, // 1982-07-01
    {45516, 22}, // 1983-07-01
    {46247, 23}, // 1985-07-01
    {47161, 24}, // 1988-01-01
    {47892, 25}, // 1990-01-01
    {48257, 26}, // 1991-01-01
    {48804, 27}, // 1992-07-01
    {49169, 28}, // 1993-07-01
    {49534, 29}, // 1994-07-01
    {50083, 30}, // 1996-01-01
    {50630, 31}, // 1997-07-01
    {51179, 32}, // 1999-01-01
    {53736, 33}, // 2006-01-01
    {54832, 34}, // 2009-01-01
    {56109, 35}, // 2012-07-01
    {57204, 36}, // 2015-07-01
    {57754, 37}, // 2017-01-01
};

double lagrangeInterp(double t, const double times[4], const double values[4]) {
    double result = 0.0;
    for (int i = 0; i < 4; ++i) {
        double term = values[i];
        for (int j = 0; j < 4; ++j) {
            if (i != j) term *= (t - times[j]) / (times[i] - times[j]);
        }
        result += term;
    }
    return result;
}

int getTAIUTCDiff(int mjd) {
    if (mjd < leapSeconds[0].first) {
        return 19; 
    }

    auto it = upper_bound(leapSeconds.begin(), leapSeconds.end(), mjd,
        [](int mjd, const std::pair<int, int>& p) { return mjd < p.first; });

    if (it == leapSeconds.begin()) {
        return leapSeconds.front().second;
    }
    return prev(it)->second;
}

Transformation::EopData Transformation::InterpolateEop(double targetMjd) {
    const std::vector<EopData>& eopList = eop_list_;
    size_t k = lower_bound(eopList.begin(), eopList.end(), targetMjd,
        [](const EopData& e, double mjd) { return e.mjd < mjd; }) - eopList.begin();
    k = std::max<size_t>(1, std::min(k, eopList.size() - 3));

    EopData interpEop;
    const double times[4]= {(double)eopList[k - 1].mjd, (double)eopList[k].mjd,(double)eopList[k + 1].mjd, (double) eopList[k + 2].mjd  };
    const double xpValues[] = { eopList[k - 1].xp, eopList[k].xp, eopList[k + 1].xp, eopList[k + 2].xp };
    interpEop.xp = lagrangeInterp(targetMjd, times, xpValues);
    const double ypValues[] = { eopList[k - 1].yp, eopList[k].yp, eopList[k + 1].yp, eopList[k + 2].yp };
    interpEop.yp = lagrangeInterp(targetMjd, times, ypValues);
    const double dutValues[] = { eopList[k - 1].dut, eopList[k].dut, eopList[k + 1].dut, eopList[k + 2].dut };
    interpEop.dut = lagrangeInterp(targetMjd, times, dutValues);
    interpEop.DAT = getTAIUTCDiff(targetMjd);  
    return interpEop;
}

Transformation::Transformation() {
}
Transformation::~Transformation() {
  eop_list_.clear();
}

bool Transformation::Load(const char* eopfile) {
    std::vector<EopData> eopList;
    std::ifstream file(eopfile);
    std::string line;

    if (file.is_open() == false) {
        return false;
    }

    while (getline(file, line)) {
        if (line.empty() || line.find("I") == std::string::npos) continue;

        EopData eop;
        line = line.substr(7);
        std::istringstream iss(line);
        std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{} };

        if (tokens.size() < 15) {
            break;
        }
        eop.mjd = stoi(tokens[0]);
        eop.xp = stod(tokens[2]);
        eop.yp = stod(tokens[4]);
        eop.dut = stod(tokens[7]);

        eopList.push_back(eop);
    }

    if (eopList.size() < 1) {
        return false;
    }

    std::swap(eop_list_, eopList);
    return true;
}

void Transformation::ComputeRotation(int year, int month, int day, int hour,
                                     int min, double sec,
                                     double r[3][3]) {

  double jdTime = utc2jdTime(year, month,
                             day, hour,
                             min, sec);
  double mjd = jdTime - 2400000.5;
  EopData eop = InterpolateEop(mjd);
  double djm0, tt, date, tut;
  time2jd(year, month, day,
          hour, min, sec, eop.DAT,
          eop.dut, &djm0, &tt, &date, &tut);

  rotationECI2ECEF(djm0, tt, date, tut, eop.xp, eop.yp, r);
}

  };
};
