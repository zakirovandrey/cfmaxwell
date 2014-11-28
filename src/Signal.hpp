#ifndef SOURCE_HPP
#define SOURCE_HPP
#include "params.hpp"

//------------Signal params-------------------------------//
struct parsSig{
  double Xlength, Ylength, Zlength;
  double Xcnt, Ycnt, Zcnt;
  double Rx, Ry;
  double Rh, xQ, yQ;
  double Th;
  double F0;
  double A0;
  double omega, Omega;
  double Px, Py, Pz; 		//polarization
  double alpha;
  double kx, ky, kz, k0;
  double Lambda;

  void setPars( ) {
    Xlength = pars.dx*(1<<MaxRank)*pars.Nx;
    Ylength = pars.dy*(1<<MaxRank)*pars.Ny;
    Zlength = pars.dz*(1<<MaxRank)*pars.Nz; 
    Xcnt = 0.5*Xlength;
    Ycnt = 0.5*Ylength;
    Zcnt = 0.0;
    Rh = 1.0;
    xQ = Xlength - 10; 
    yQ = Ylength - 10; 
    A0 = 1.0; 
    Px = A0*cos(alpha); 	//polarization
    Py = 0.0; 			//polarization
    Pz = A0*sin(alpha); 	//polarization
    k0 = 1.0;
    kx = k0*sin(alpha);
    ky = 0.0;
    kz = sqrt(k0*k0 - kx*kx - ky*ky);
  }
  inline double Phase(double X, double Y) { return 0.0;}
  inline double Kmode(double X, double Y) { return kx*X + ky*Y;}
 };

extern parsSig parsSignal;
parsSig* getParsSrc(); 

namespace Plane{};
namespace Gauss{};
namespace Cosine{};

#ifndef SWIG
#define SRCS(Type) \ 
namespace Type{\ 
	float SrcDx_(int ix, int iy, MainType deep, long int it); \
	float SrcDy_(int ix, int iy, MainType deep, long int it); \
	float SrcDz_(int ix, int iy, MainType deep, long int it); \
	float SrcBx(int ix, int iy, MainType deep, long int it); \
	float SrcBy(int ix, int iy, MainType deep, long int it); \
	float SrcBz(int ix, int iy, MainType deep, long int it); \
	float SrcDx(int ix, int iy, MainType deep, long int it);  \
	float SrcDy(int ix, int iy, MainType deep, long int it); \
	float SrcDz(int ix, int iy, MainType deep, long int it); \
	float SrcInsideDx(MainType& val, int ix, int iy, int iz, int it);  \
	float SrcInsideDy(MainType& val, int ix, int iy, int iz, int it);  \
	float SrcInsideDz(MainType& val, int ix, int iy, int iz, int it);  \
	float SrcInsideBx(MainType& val, int ix, int iy, int iz, int it);  \
	float SrcInsideBy(MainType& val, int ix, int iy, int iz, int it);  \
	float SrcInsideBz(MainType& val, int ix, int iy, int iz, int it);  \
};

SRCS(Plane)
SRCS(Gauss)
SRCS(Cosine)
SRCS(continuousEmission)

using namespace Plane;

#endif //SWIG
#endif//SOURCE_HPP
