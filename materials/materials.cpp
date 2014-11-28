#include "params.hpp"
#include "materials.hpp"
#include <math.h>

parsMaterial parsMat;
parsMaterial* getParsMat(){ return &parsMat; }; 

const double sqrt3 = sqrt(3);
const double dsqrt3 = 1./sqrt(3);
const double sqrt2 = sqrt(2);

int Material (MainType x, MainType y, MainType z) { //return IndVac;
// periodic 3D ph. crystall
/*    if (x<parsMat.Xbeg || x>parsMat.Xend || y<parsMat.Ybeg || y>parsMat.Yend || z<parsMat.Zbeg || z>parsMat.Zend) return IndVac;
    MainType _x=x-parsMat.Xbeg, _y=y-parsMat.Ybeg, _z=z-parsMat.Zbeg;
    bool layer=false;
    for(int i=0;i<3;i++) {
      int mz = round((_y-i*parsMat.d)*parsMat.k_mz);
      int m1 = round((_x-i*parsMat.D1*dsqrt3)*parsMat.k_m1);
      int m2 = round((_z-(_x-i*parsMat.D1*dsqrt3)*dsqrt3)*parsMat.k_m2);
      MainType d_x = _x-i*parsMat.D1*dsqrt3-m1*parsMat.D1*sqrt3*0.5;
      MainType d_z = _z-m2*parsMat.D2-m1*parsMat.D1*0.5;
      layer = layer || (d_x*d_x + d_z*d_z <parsMat.r*parsMat.r && fabs(_y-i*parsMat.d-mz*3*parsMat.d)<0.5*parsMat.h);
    }
    layer = layer && (_x>=parsMat.a && x<=parsMat.Xend-parsMat.a && _y>=parsMat.a && y<=parsMat.Yend-parsMat.a && _z>=0.*parsMat.h && z<=parsMat.Zend-0.*parsMat.h);
    layer = layer || (_x>parsMat.xdel*parsMat.D1*sqrt3*0.5 && _x<(parsMat.xdel+1)*parsMat.D1*sqrt3*0.5 && _y>parsMat.ydel*(3*parsMat.d)-parsMat.h*0.5 && _y<parsMat.ydel*(3*parsMat.d)+parsMat.h*0.5-parsMat.d); // waveguide
    if(layer) return IndVac;
    else      return IndMat;*/

//particle
   if (x>parsMat.Xbeg && x<parsMat.Xend && y>parsMat.Ybeg && y<parsMat.Yend && z>parsMat.Zbeg && z<=parsMat.Zend) return IndMat;
   else if ( (x-parsMat.pCx)*(x-parsMat.pCx)+(y-parsMat.pCy)*(y-parsMat.pCy)+(z-parsMat.pCz)*(z-parsMat.pCz)<=parsMat.pR*parsMat.pR ) return IndMat;
   else return IndVac;
}

