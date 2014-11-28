#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include "params.hpp"
#include <cmath>

enum {IndAir=0, IndGold=1, IndBIG=2, IndGGG=3};

int getMaterial (MainType coords[3]);

struct parsMaterial{
	double thinGGG, thinBIG, thinGold;
	double zGGG, zBIG, zGold;
	double X0, Y0;
	double Xlength, Ylength;
	double cell, period, shift;
  void setPars(){
   	Xlength = pars.dx*(1<<MaxRank)*pars.Nx;
   	Ylength = pars.dy*(1<<MaxRank)*pars.Ny;
	X0 = 0.5*Xlength;
	Y0 = 0.5*Ylength;
  };  
 inline bool inCell(double x){
   if (fmod(x+shift, period)<=cell) return true;
   return false;
 };
};
parsMaterial* getParsMat();
#endif //MATERIAL_HPP

