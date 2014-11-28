#include "params.hpp"
#include "materials.hpp"
#include "Signal.hpp"
#include <cmath>

parsMaterial Material;

parsMaterial* getParsMat(){ return &Material; }; 


int getMaterial (MainType coords[3]){
  MainType x = coords[0];
  MainType y = coords[1];
  MainType z = coords[2];

// Get Gold
  if (
    (fabs(x-Material.X0   ) <= Material.Xlength     ) &&
    (fabs(y-Material.Y0   ) <= Material.Ylength     ) &&
    (fabs(z-Material.zGold) <= 0.5*Material.thinGold) &&
    ( Material.inCell(x) )
  ) return IndGold;
// Get BIG
  else if(
    (fabs(x-Material.X0   ) <= Material.Xlength     ) &&
    (fabs(y-Material.Y0   ) <= Material.Ylength     ) &&
    (fabs(z-Material.zBIG) <= 0.5*Material.thinBIG)
  ) return IndBIG;
// Get GGG
  else if (
    (fabs(x-Material.X0   ) <= Material.Xlength     ) &&
    (fabs(y-Material.Y0   ) <= Material.Ylength     ) &&
    (fabs(z-Material.zGGG) <= 0.5*Material.thinGGG)
  ) return IndGGG;
// Get Air
  else 
    return IndAir;
};
