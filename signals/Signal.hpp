#ifndef SOURCE_HPP
#define SOURCE_HPP
#include "params.hpp"

//------------Signal params--------------------------------------------------------//
struct parsSig{
  double Xlength, Ylength, Zlength;
  double Xcnt, Ycnt;

  double Rh;
  double Th;
  double F0;
  double A0, omega, k0;
  double Px, Py; //polarization
  void setPars() {
    Xlength=pars.dx*(1<<MaxRank)*pars.Nx; Ylength=pars.dy*(1<<MaxRank)*pars.Ny; Zlength=pars.dz*(1<<MaxRank)*pars.Nz; 
    Xcnt=0.5*Xlength; Ycnt=0.5*Ylength;
    F0=2.5;
    Rh=(Xlength<=Ylength?Xlength:Ylength)*0.35;
    Th=7.0/F0;
    A0=1.0; omega=2*M_PI*F0; k0=omega/1.;
    Px=1.0; Py=0.0; //polarization
  }
};
parsSig* getParsSrc();

float SrcDx(int ix, int iz, MainType deep, long int it);
float SrcDy(int ix, int iz, MainType deep, long int it);
float SrcDz(int ix, int iz, MainType deep, long int it);
float SrcBx(int ix, int iz, MainType deep, long int it);
float SrcBy(int ix, int iz, MainType deep, long int it);
float SrcBz(int ix, int iz, MainType deep, long int it);

#endif//SOURCE_HPP
