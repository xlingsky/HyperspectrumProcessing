#if !defined LX_GEOMETRY_RPC_H
#define LX_GEOMETRY_RPC_H

namespace xlingeo{


//Make sure the longitude is between -180.00 .. 179.9
#define NORMALIZE_LONGTITUDE(lon) ((lon+180)-int((lon+180)/360)*360-180)

inline int      UTM_Zone(float lon){
  return int((lon + 180)/6) + 1;
}
inline int      UTM_Zone(float lon, float lat){
  int zonenumber = UTM_Zone(lon);
	if( lat >= 56.0 && lat < 64.0 && lon >= 3.0 && lon < 12.0 )
		zonenumber = 32;

	// Special zones for Svalbard
	if( lat >= 72.0 && lat < 84.0 ) 
	{
		if(      lon >= 0.0  && lon <  9.0 ) zonenumber = 31;
		else if( lon >= 9.0  && lon < 21.0 ) zonenumber = 33;
		else if( lon >= 21.0 && lon < 33.0 ) zonenumber = 35;
		else if( lon >= 33.0 && lon < 42.0 ) zonenumber = 37;
	}
  return zonenumber;
}
inline char UTM_LetterDesignator(double Lat)
{
	//This routine determines the correct UTM letter designator for the given latitude
	//returns 'Z' if latitude is outside the UTM limits of 84N to 80S
	//Written by Chuck Gantz- chuck.gantz@globalstar.com
	char LetterDesignator;

	if((84 >= Lat) && (Lat >= 72)) LetterDesignator = 'X';
	else if((72 > Lat) && (Lat >= 64)) LetterDesignator = 'W';
	else if((64 > Lat) && (Lat >= 56)) LetterDesignator = 'V';
	else if((56 > Lat) && (Lat >= 48)) LetterDesignator = 'U';
	else if((48 > Lat) && (Lat >= 40)) LetterDesignator = 'T';
	else if((40 > Lat) && (Lat >= 32)) LetterDesignator = 'S';
	else if((32 > Lat) && (Lat >= 24)) LetterDesignator = 'R';
	else if((24 > Lat) && (Lat >= 16)) LetterDesignator = 'Q';
	else if((16 > Lat) && (Lat >= 8)) LetterDesignator = 'P';
	else if(( 8 > Lat) && (Lat >= 0)) LetterDesignator = 'N';
	else if(( 0 > Lat) && (Lat >= -8)) LetterDesignator = 'M';
	else if((-8> Lat) && (Lat >= -16)) LetterDesignator = 'L';
	else if((-16 > Lat) && (Lat >= -24)) LetterDesignator = 'K';
	else if((-24 > Lat) && (Lat >= -32)) LetterDesignator = 'J';
	else if((-32 > Lat) && (Lat >= -40)) LetterDesignator = 'H';
	else if((-40 > Lat) && (Lat >= -48)) LetterDesignator = 'G';
	else if((-48 > Lat) && (Lat >= -56)) LetterDesignator = 'F';
	else if((-56 > Lat) && (Lat >= -64)) LetterDesignator = 'E';
	else if((-64 > Lat) && (Lat >= -72)) LetterDesignator = 'D';
	else if((-72 > Lat) && (Lat >= -80)) LetterDesignator = 'C';
	else LetterDesignator = 'Z'; //This is here as an error flag to show that the Latitude is outside the UTM limits

	return LetterDesignator;
}

typedef struct tagRpcPara
{
  double line_off;     // line offset in pixels
  double samp_off;     // sample offset in pixels
  double lat_off;      // latitude offset in degrees
  double long_off;     // longitude offset in degrees
  double height_off;   // height offset in meters
  double line_scale;   // line scale in pixels
  double samp_scale;   // sample scale in pixels
  double lat_scale;    // latitude scale in degrees
  double long_scale;   // longitude scale in degrees
  double height_scale; // height scale in meters
  double c[20];        // 20 line numerator coefficients
  double d[20];        // 20 line denominator coefficients
  double a[20];        // 20 sample numerator coefficients
  double b[20];        // 20 sample denominator coefficients

}RpcPara;

class Rpc
{
public:
  Rpc();
  virtual ~Rpc();
  void SetAop6(const double aop6[6]);
  void GetRpcPara(RpcPara* para) const;
  RpcPara* parameters();
  const RpcPara* parameters() const;
  //if dem is loaded, alt is relative height to dem 
  bool AttachDemFile(const char* lpstrPathName);

  bool Load(const char* lpstrPathName);
  bool Save(const char* lpstrPathName, double *sx = nullptr, double *sy = nullptr, double *mx = nullptr, double *my = nullptr) const;
  //use 2D and 3D points to solve RPC parameters
  bool Solve(const double* samp, const  double* line, const  double* lat, const double* lon, const  double* alt, int num, double *sx = nullptr, double *sy = nullptr, double *mx = nullptr, double *my = nullptr);

  void Ground2Photo(double lat, double lon, double alt, double *samp, double *line) const;
  void Ground2Photo(double lat, double lon, double alt, double *samp, double *line, bool bAop6) const
  {
    Ground2Photo(lat, lon, alt, samp, line);
    if (bAop6) {  pxy_to_ixy(*samp, *line, samp, line); }
  }
  void Ground2Photo(const double* lat, const double* lon, double* alt, int nPtNum, double *samp, double *line) const;
  void Ground2Photo(const double* lat, const double* lon, double* alt, int nPtNum, double *samp, double *line, bool bAop6) const
  {
    Ground2Photo(lat, lon, alt, nPtNum, samp, line);
    if (bAop6){
      while(nPtNum-->0){ pxy_to_ixy(*samp, *line, samp, line); ++samp; ++line; }
    }
  }

  void PhotoZ2Ground(double samp, double line, double gz, double *lat, double *lon, double *alt) const;
  void PhotoZ2Ground(double samp, double line, double gz, double *lat, double *lon, double *alt, bool bAop6) const
  {
    if (bAop6) { ixy_to_pxy(samp, line, &samp, &line); }
    PhotoZ2Ground(samp, line, gz, lat, lon, alt);
  }
  void PhotoZ2Ground(const double* samp, const  double* line, const  double* gz,int nPtNum, double *lat, double *lon, double *alt) const;
  void PhotoZ2Ground(const double* samp, const  double* line, const double* gz, int nPtNum, double *lat, double *lon, double *alt, bool bAop6) const;

  void pxy_to_ixy(double px, double py, double *ix, double *iy) const;
  void ixy_to_pxy(double ix, double iy, double *px, double *py) const;

  //x=samp y=line
  friend void RpcIntersection(
      const Rpc *pRpcL, double pxl, double pyl,
      const Rpc *pRpcR, double pxr, double pyr,
      double *lat, double *lon, double *alt,
      bool bAop6
                           );
protected:
    void Reset();

    int*    m_pDxyGrid; // It saved value is: dx= (realx-rpcx)*1000 , dy=(realy-rpcy)*1000;
    int     m_gridC, m_gridR, m_gridDx, m_gridDy;
    int*    NewGrid(int gridC, int gridR, int dx, int dy);
    void    GetDxy(double px, double py, double *dx, double *dy) const;
private:
  class RpcBase* rpc_;
  double      aop6_[6];
};

};


#endif

