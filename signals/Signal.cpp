#include "math.h"
#include "Signal.hpp"
#include "params.hpp"

parsSig parsSignal;
parsSig* getParsSrc(){ return &parsSignal; }; 

//Possible Regions of Boundary Source initials
inline bool RegionCircle(double x, double y) { if(x*x+y*y<parsSignal.Rh*parsSignal.Rh) return true; else return false; }
inline bool RegionLine(double x, double y) { if(fabs(y)<parsSignal.Rh) return true; else return false; }
inline bool RegionQuadrat(double x, double y) { if(fabs(x)<parsSignal.Rh && fabs(y)<parsSignal.Rh) return true; else return false; }
inline bool Region(double x, double y) { return RegionCircle(x,y); }

inline double EnvelopeR(double x, double y) {
  double r2=x*x+y*y; if(r2>=parsSignal.Rh*parsSignal.Rh) return 0.0;
/*Cosine,Cosine/2*/
    const double Kx=0.5*M_PI/parsSignal.Rh, Ky=0.5*M_PI/parsSignal.Rh; return
          0.5*(1.0+cos(2.0*Kx*sqrt(x*x+y*y)));
/*Gauss*/ //const double Kr=1.0/(parsSignal.Rh*parsSignal.Rh); double r2=Kr*(x*x+y*y); return exp(-r2)-sqrt(r2)*(exp(-9.)/3.0);
}
inline double Lambda2(float x) { return fabs(x)<1.5?(fabs(x)<0.5?3./4-x*x:pow((3./2-fabs(x)),2)/2):0 ; }
inline double Lambda3(float x) { return fabs(x)<2?(fabs(x)<1?2./3-x*x+pow(fabs(x),3)/2.:pow((2.-fabs(x)),3)/6):0; }
inline double Lambda3p1(float x) { return fabs(x)<2?(fabs(x)<1?-2*x+x*fabs(x)*3/2.:(x<0?(2.+x)*(2.+x):-(2.-x)*(2.-x))/2):0; }
inline double EnvelopeT(double t) {
  if(t<0 || t>parsSignal.Th) return 0.0;
  const double Kt=2.0*M_PI/parsSignal.Th; return
    //1.0;
    //sin(0.5*Kt*t);
    //cos(Kt*t);
    0.5*(1.0-cos(Kt*t));
    //1.0-cos(Kt*t);
    //const double Kt=3.0/(2.0*parsSignal.Th); return exp(-Kt*t)-(exp(-9.)/3.0)*Kt*t;
    //4./3.*Lambda2(8./3.*(t-1/parsSignal.F0)*parsSignal.F0);                            
    //3./2.*Lambda3p1(8./3.*(t-1/parsSignal.F0)*parsSignal.F0);                        
    //sig3(t,T)=3./2.*Lambda3(3*t/T)  
}

inline double Latency(double X, double Y) { return 0.; }

const double kx=0.*parsSignal.k0, ky=0.*parsSignal.k0, kz=sqrt(parsSignal.k0*parsSignal.k0-kx*kx-ky*ky);
inline double Phase(double X, double Y) {
  return 0.0;
}

inline double Kmode(double X, double Y) { return kx*X+ky*Y; }

float SrcDx(int ix, int iy, MainType deep, long int it) { 
  return parsSignal.A0*parsSignal.Px*EnvelopeR(ix*pars.dx-parsSignal.Xcnt,(iy+0.5)*pars.dy-parsSignal.Ycnt)*EnvelopeT(it*pars.dt)*cos(2*M_PI*parsSignal.F0*it*pars.dt);
}
float SrcDy(int ix, int iy, MainType deep, long int it) {
  return parsSignal.A0*parsSignal.Py*EnvelopeR((ix+0.5)*pars.dx-parsSignal.Xcnt,iy*pars.dy-parsSignal.Ycnt)*EnvelopeT(it*pars.dt)*cos(2*M_PI*parsSignal.F0*it*pars.dt);
}
float SrcDz(int ix, int iy, MainType deep, long int it) { return 0.; }

float SrcBx(int ix, int iy, MainType deep, long int it) { return 0.; }
float SrcBy(int ix, int iy, MainType deep, long int it) { return 0.; }
float SrcBz(int ix, int iy, MainType deep, long int it) { return 0.; }

