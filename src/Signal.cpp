#include "math.h"
#include "Signal.hpp"
#include "materials.hpp"
#include "params.hpp"

extern parsMaterial Material;
parsSig parsSignal;
parsSig* getParsSrc(){ return &parsSignal; }; 


//====================================Possible Regions of Boundary Source initials==================================
namespace Region{
	inline bool Circle(double x, double y){
	   	if(x*x+y*y < parsSignal.Rh*parsSignal.Rh) return true;
     	  	else return false;
	}
	inline bool Line(double x, double y) {
      	 	if(fabs(y)<parsSignal.Rh) return true;
       		else return false;
	}
 	inline bool Quadrate(double x, double y) {
    	   	if(fabs(x)<parsSignal.xQ && fabs(y)<parsSignal.yQ)
		     return true;
      	 	else return false;
	}
	inline bool Entire(double x, double y) {
		return true;
	}
}; 

//=================================================================================================================

//==========================================Possible Initial Signals===============================================
namespace Plane{
	inline double EnvelopeR(double x, double y, double t) {
		double Xcnt = parsSignal.Xcnt;
		double Ycnt = parsSignal.Ycnt;
		double xCurrent = !parsSignal.kx ? Material.Xlength : t*parsSignal.omega/parsSignal.kx;
		if( (x > xCurrent) && Region::Quadrate(x-Xcnt, y-Ycnt) ) return 0.0;
	//	if( x > xCurrent ) return 0.0;
		return 0.5*(1.0 + cos(-parsSignal.kx*x+parsSignal.omega*t));
	}; 
	inline double EnvelopeT(double t) {
  		if( (t < 0) || (t > 0.5*parsSignal.Th) ) return 0.0;
  		const double Omega = parsSignal.Omega;//1.0*M_PI/parsSignal.Th;
		const double omega = parsSignal.omega;//2.0*M_PI*parsSignal.F0;
  		return sin(Omega*t)*sin(Omega*t);//*sin(omega*t);
		//0.5 * (1.0 - cos(Kt*t)) * cos(2*M_PI*parsSignal.F0*t);
	};
	float SrcDx(int ix, int iy, MainType deep, long int it) { 
		return parsSignal.A0*parsSignal.Px*
		EnvelopeR(ix*pars.dx, iy*pars.dy - parsSignal.Ycnt, it*pars.dt)*
	 	EnvelopeT(it*pars.dt);
	};  
	float SrcDy(int ix, int iy, MainType deep, long int it) {
  		return parsSignal.A0*parsSignal.Py*
	 	EnvelopeR(ix*pars.dx - parsSignal.Xcnt, iy*pars.dy - parsSignal.Ycnt, it*pars.dt)*
	 	EnvelopeT(it*pars.dt);
	};
	float SrcDz(int ix, int iy, MainType deep, long int it) { 
  		return parsSignal.A0*parsSignal.Pz*
	 	EnvelopeR(ix*pars.dx, iy*pars.dy - parsSignal.Ycnt, it*pars.dt)*
	 	EnvelopeT(it*pars.dt);
	};

	float SrcBx(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBy(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBz(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcDx_(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcDy_(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcDz_(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcInsideDx(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcDx_(ix,iy,0,it); }
	float SrcInsideDy(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcDy_(ix,iy,0,it); }
	float SrcInsideDz(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcDz_(ix,iy,0,it); }
	float SrcInsideBx(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBx(ix,iy,0,it); }
	float SrcInsideBy(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBy(ix,iy,0,it); }
	float SrcInsideBz(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBz(ix,iy,0,it); }
};

namespace Gauss{
	inline double EnvelopeR(double x, double y) {
		if(x*x + y*y >= parsSignal.Rh*parsSignal.Rh) return 0.0;
//	       	const double Kx = 0.5*M_PI/parsSignal.Rh;
//		const double Ky = 0.5*M_PI/parsSignal.Rh;
//		return 0.5*(1.0 + cos(2.0*Kx*sqrt(x*x + y*y)));
	/*Gauss*/
		const double Kr=1.0/(parsSignal.Rh*parsSignal.Rh); 
		double r2 = (x-parsSignal.Xcnt) * (x-parsSignal.Xcnt) + (y-parsSignal.Ycnt) * (y-parsSignal.Ycnt); 
		return exp(-Kr*r2)-sqrt(Kr*r2)*(exp(-9.)/3.0);
	};
	inline double EnvelopeT(double t) {
		if( (t < 0) || (t > parsSignal.Th) ) return 0.0;
		const double Kt = 2.0*M_PI/parsSignal.Th;
		return cos(Kt*t); 
	};
	float SrcDx_(int ix, int iy, MainType deep, long int it) { 
		return parsSignal.A0*parsSignal.Px*
		EnvelopeR(ix*pars.dx - parsSignal.Xcnt, (iy+0.5)*pars.dy - parsSignal.Ycnt)*
		EnvelopeT(it*pars.dt);
	};
	float SrcDy_(int ix, int iy, MainType deep, long int it) {
  		return parsSignal.A0*parsSignal.Py*
	 	EnvelopeR(ix*pars.dx - parsSignal.Xcnt, iy*pars.dy - parsSignal.Ycnt)*
	 	EnvelopeT(it*pars.dt);
	};
	float SrcDz_(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBx(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBy(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBz(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcDx(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcDy(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcDz(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcInsideDx(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcDx_(ix,iy,0,it); }
	float SrcInsideDy(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcDy_(ix,iy,0,it); }
	float SrcInsideDz(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcDz_(ix,iy,0,it); }
	float SrcInsideBx(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBx(ix,iy,0,it); }
	float SrcInsideBy(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBy(ix,iy,0,it); }
	float SrcInsideBz(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBz(ix,iy,0,it); }
};

namespace Cosine{
	inline double EnvelopeR(double x, double y) {
		if(x*x + y*y >= parsSignal.Rh*parsSignal.Rh) return 0.0;
	       	const double Kx = 0.5*M_PI/parsSignal.Rh;
		const double Ky = 0.5*M_PI/parsSignal.Rh;
		return 0.5*(1.0 + cos(2.0*Kx*sqrt(x*x + y*y)));
	}
	inline double EnvelopeT(double t) {
		if( (t < 0) || (t > parsSignal.Th) ) return 0.0;
		const double Kt = 2.0*M_PI/parsSignal.Th;
		return 0.5 * (1.0 - cos(Kt*t)) * cos(2*M_PI*parsSignal.F0*t);
	}
	float SrcDx_(int ix, int iy, MainType deep, long int it) { 
		return  parsSignal.A0*parsSignal.Px*
			EnvelopeR(ix*pars.dx - parsSignal.Xcnt, (iy+0.5)*pars.dy - parsSignal.Ycnt)*
			EnvelopeT(it*pars.dt);
	}
	float SrcDy_(int ix, int iy, MainType deep, long int it) {
  		return  parsSignal.A0*parsSignal.Py*
			EnvelopeR(ix*pars.dx - parsSignal.Xcnt, iy*pars.dy - parsSignal.Ycnt)*
	 		EnvelopeT(it*pars.dt);
	}
	float SrcDz_(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBx(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBy(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBz(int ix, int iy, MainType deep, long int it) { return 0.; }
float SrcDx(int ix, int iy, MainType deep, long int it) { return SrcDx_(ix,iy,0,it); }
float SrcDy(int ix, int iy, MainType deep, long int it) { return SrcDy_(ix,iy,0,it); }
float SrcDz(int ix, int iy, MainType deep, long int it) { return SrcDz_(ix,iy,0,it); }
float SrcInsideDx(MainType& val, int ix, int iy, int iz, int it) { }
float SrcInsideDy(MainType& val, int ix, int iy, int iz, int it) { }
float SrcInsideDz(MainType& val, int ix, int iy, int iz, int it) { }
float SrcInsideBx(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBx(ix,iy,0,it); }
float SrcInsideBy(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBy(ix,iy,0,it); }
float SrcInsideBz(MainType& val, int ix, int iy, int iz, int it) { if(it*pars.dt<=parsSignal.Th) val= SrcBz(ix,iy,0,it); }
};
//=================================================================================================================
namespace continuousEmission{ 
	inline double EnvelopeR(double x, double y) {
		const double Rx2 = (parsSignal.Rx*parsSignal.Rx); 
		const double Ry2 = (parsSignal.Ry*parsSignal.Ry); 
		if(x*x/Rx2 + y*y/Ry2 >= 1) return 0.0;
		double x2 = x*x/Rx2;
		double y2 = y*y/Ry2;
		return cos(0.5*M_PI*sqrt(x2+y2))*cos(0.5*M_PI*sqrt(x2+y2));
	};
	inline double EnvelopeT(double t) {
		if( (t < 0) ) return 0.0;
		if(parsSignal.Omega*t > M_PI) return  cos(parsSignal.omega*t);
		return 0.5*cos(parsSignal.omega*t)*(1-cos(parsSignal.Omega*t));
	};
	float SrcDx(int ix, int iy, MainType deep, long int it) {
		double retval = parsSignal.A0*parsSignal.Px*
		EnvelopeR(ix*pars.dx - parsSignal.Xcnt, iy*pars.dy - parsSignal.Ycnt)*
		EnvelopeT(it*pars.dt);
//		if(fabs(ix-parsSignal.Xcnt/pars.dx)<3 && fabs(iy-parsSignal.Ycnt/pars.dy)<3) 
//			printf(" OK %d %d %g\n", ix, iy, EnvelopeR(ix*pars.dx - parsSignal.Xcnt, iy*pars.dy - parsSignal.Ycnt));
		return retval;
	}; 
	float SrcDy(int ix, int iy, MainType deep, long int it) {
  		return parsSignal.A0*parsSignal.Py*
	 	EnvelopeR(ix*pars.dx - parsSignal.Xcnt, iy*pars.dy - parsSignal.Ycnt)*
	 	EnvelopeT(it*pars.dt);
	}; 
	float SrcDz(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBx(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBy(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcBz(int ix, int iy, MainType deep, long int it) { return 0.; }
	float SrcInsideDx(MainType& val, int ix, int iy, int iz, int it) { }//if(it*pars.dt<=parsSignal.Th) val= SrcDx_(ix,iy,0,it); }
	float SrcInsideDy(MainType& val, int ix, int iy, int iz, int it) { }//if(it*pars.dt<=parsSignal.Th) val= SrcDy_(ix,iy,0,it); }
	float SrcInsideDz(MainType& val, int ix, int iy, int iz, int it) { }//if(it*pars.dt<=parsSignal.Th) val= SrcDz_(ix,iy,0,it); }
	float SrcInsideBx(MainType& val, int ix, int iy, int iz, int it) { }//if(it*pars.dt<=parsSignal.Th) val= SrcBx(ix,iy,0,it); }
	float SrcInsideBy(MainType& val, int ix, int iy, int iz, int it) { }//if(it*pars.dt<=parsSignal.Th) val= SrcBy(ix,iy,0,it); }
	float SrcInsideBz(MainType& val, int ix, int iy, int iz, int it) { }//if(it*pars.dt<=parsSignal.Th) val= SrcBz(ix,iy,0,it); }
};
