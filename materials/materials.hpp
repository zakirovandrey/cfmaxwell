#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include "params.hpp"

enum {IndVac=0, IndMat=1};

int Material (MainType x, MainType y, MainType z);
inline int getMaterial (MainType coords[3]) {  
    MainType x=coords[0], y=coords[1], z=coords[2];
    return Material(x,y,z);
}

struct parsMaterial {
  MainType Xbeg,Xend,Ybeg,Yend,Zbeg,Zend;
  MainType a,h,d,r,D1,D2;
  MainType k_mz,k_m1,k_m2;
  int xdel;
  int ydel;

  MainType pCx,pCy,pCz,pR;
  void setPars(){
    double Xlength=pars.dx*(1<<MaxRank)*pars.Nx, Ylength=pars.dy*(1<<MaxRank)*pars.Ny, Zlength=pars.dz*(1<<MaxRank)*pars.Nz;
    Xbeg=10*pars.dx; Xend=Xlength-10*pars.dx; Ybeg=10*pars.dy; Yend=Ylength-10*pars.dy; Zbeg=10*pars.dz; Zend=Zlength*0.5;
    a=50*pars.dx; h=0.93*a; d=a/sqrt(3.);
    r=0.293*a;
    D1=a*0.5*sqrt(2.); D2=D1;
    k_mz=1./(3*d);
    k_m1=1./(D1*sqrt(3.)*0.5);
    k_m2=1./D2;
    xdel=7; ydel=3;

    pR=0.07;
    pCx=Xlength*0.5; pCy=Ylength*0.5; pCz=Zend+pR;
  }
};

parsMaterial* getParsMat();

#endif //MATERIAL_HPP

