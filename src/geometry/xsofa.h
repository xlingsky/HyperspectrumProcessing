#ifndef _XSOFA_H
#define _XSOFA_H

#include <vector>

namespace s2r {
namespace xsofa {

// J2000 to WGS84
void rotationECI2ECEF(double djmjd0, double tt, double date, double tut,
                      double xp, double yp, double r[3][3]);

double utc2jdTime(int year, int month, int day, 
                       int hour, int minute, double second) ;

void time2jd(int year, int month, int day, int hour, int minute, double second, int eopDAT, double eopdut,
    double *djm0, double *tt, double *date, double *tut);

void jdTime2utc(double jdTime, int* year, int* month, int* day,
    int* hour, int* minute, double* second);

class Transformation {
    private:
      struct EopData {
        int mjd, DAT;
        double xp, yp, dut;
      };
      std::vector<EopData> eop_list_;

    public:
    Transformation();
    virtual ~Transformation();
    bool Load(const char* eopfile);
    //The roation converting J2000 to CGCS2000
    void ComputeRotation(int year, int month, int day, int hour, int minute,
                         double second, double r[3][3]);

  private: 
    EopData InterpolateEop(double targetMjd);
};

};
};

#endif