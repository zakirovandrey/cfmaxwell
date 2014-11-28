#ifndef R0_STRUCT_HPP
#define R0_STRUCT_HPP
#include <math.h>

#include "Coff.hpp"

struct AuxField {
  MainType Ei[3];
  MainType Eim[3];
  MainType Ji[3][Md]; MainType Jim[3][Md];
  MainType rot[3];
  int I[3];
};

struct fieldsDDD{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz

  indx ind;
};
struct fieldsSDD{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_1; MainType Di_pml1_1; MainType Ei_pml0_1; MainType Ei_pml1_1;  // Ei_pmla_b = E_{bc}; c=PMLfunc(a,b); E_b=\sum_c{...}
  MainType Di_pml0_2; MainType Di_pml1_2; MainType Ei_pml0_2; MainType Ei_pml1_2; 
  MainType Bi_pml0_1; MainType Bi_pml1_1; MainType Hi_pml0_1; MainType Hi_pml1_1; // так как Hi_pmla_2[pmliz] не нужны
  MainType Bi_pml0_2; MainType Bi_pml1_2; MainType Hi_pml0_2; MainType Hi_pml1_2; 

  indx ind;
  SrcIndx SrcInd;
  short int pmlix;
};
struct fieldsDSD{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_0; MainType Di_pml1_0; MainType Ei_pml0_0; MainType Ei_pml1_0; 
  MainType Di_pml0_2; MainType Di_pml1_2; MainType Ei_pml0_2; MainType Ei_pml1_2; 
  MainType Bi_pml0_0; MainType Bi_pml1_0; MainType Hi_pml0_0; MainType Hi_pml1_0; 
  MainType Bi_pml0_2; MainType Bi_pml1_2; MainType Hi_pml0_2; MainType Hi_pml1_2; 

  indx ind;
  SrcIndx SrcInd;
  short int pmliy;
};
struct fieldsDDS{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_0; MainType Di_pml1_0; MainType Ei_pml0_0; MainType Ei_pml1_0; 
  MainType Di_pml0_1; MainType Di_pml1_1; MainType Ei_pml0_1; MainType Ei_pml1_1; 
  MainType Bi_pml0_0; MainType Bi_pml1_0; MainType Hi_pml0_0; MainType Hi_pml1_0; 
  MainType Bi_pml0_1; MainType Bi_pml1_1; MainType Hi_pml0_1; MainType Hi_pml1_1; 

  indx ind;
  SrcIndx SrcInd;
  short int pmliz;
};
struct fieldsDSS{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_0; MainType Di_pml1_0; MainType Ei_pml0_0; MainType Ei_pml1_0; 
  MainType Di_pml0_1; MainType Di_pml1_1; MainType Ei_pml0_1; MainType Ei_pml1_1; 
  MainType Di_pml0_2; MainType Di_pml1_2; MainType Ei_pml0_2; MainType Ei_pml1_2; 
  MainType Bi_pml0_0; MainType Bi_pml1_0; MainType Hi_pml0_0; MainType Hi_pml1_0; 
  MainType Bi_pml0_1; MainType Bi_pml1_1; MainType Hi_pml0_1; MainType Hi_pml1_1; 
  MainType Bi_pml0_2; MainType Bi_pml1_2; MainType Hi_pml0_2; MainType Hi_pml1_2; 

  indx ind;
  SrcIndx SrcInd;
  short int pmliy, pmliz;
};
struct fieldsSDS{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_0; MainType Di_pml1_0; MainType Ei_pml0_0; MainType Ei_pml1_0;
  MainType Di_pml0_1; MainType Di_pml1_1; MainType Ei_pml0_1; MainType Ei_pml1_1;
  MainType Di_pml0_2; MainType Di_pml1_2; MainType Ei_pml0_2; MainType Ei_pml1_2;
  MainType Bi_pml0_0; MainType Bi_pml1_0; MainType Hi_pml0_0; MainType Hi_pml1_0;
  MainType Bi_pml0_1; MainType Bi_pml1_1; MainType Hi_pml0_1; MainType Hi_pml1_1;
  MainType Bi_pml0_2; MainType Bi_pml1_2; MainType Hi_pml0_2; MainType Hi_pml1_2;

  indx ind;
  SrcIndx SrcInd;
  short int pmlix, pmliz;
};
struct fieldsSSD{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_0; MainType Di_pml1_0; MainType Ei_pml0_0; MainType Ei_pml1_0;
  MainType Di_pml0_1; MainType Di_pml1_1; MainType Ei_pml0_1; MainType Ei_pml1_1;
  MainType Di_pml0_2; MainType Di_pml1_2; MainType Ei_pml0_2; MainType Ei_pml1_2;
  MainType Bi_pml0_0; MainType Bi_pml1_0; MainType Hi_pml0_0; MainType Hi_pml1_0;
  MainType Bi_pml0_1; MainType Bi_pml1_1; MainType Hi_pml0_1; MainType Hi_pml1_1;
  MainType Bi_pml0_2; MainType Bi_pml1_2; MainType Hi_pml0_2; MainType Hi_pml1_2;

  indx ind;
  SrcIndx SrcInd;
  short int pmlix, pmliy;
};
struct fieldsSSS{
  MainType Di[3]; MainType Ei[3]; //Ex,Ey,Ez
  MainType Bi[3]; MainType Hi[3]; //Hx,Hy,Hz
  MainType Di_pml0_0; MainType Di_pml1_0; MainType Ei_pml0_0; MainType Ei_pml1_0;
  MainType Di_pml0_1; MainType Di_pml1_1; MainType Ei_pml0_1; MainType Ei_pml1_1;
  MainType Di_pml0_2; MainType Di_pml1_2; MainType Ei_pml0_2; MainType Ei_pml1_2;
  MainType Bi_pml0_0; MainType Bi_pml1_0; MainType Hi_pml0_0; MainType Hi_pml1_0;
  MainType Bi_pml0_1; MainType Bi_pml1_1; MainType Hi_pml0_1; MainType Hi_pml1_1;
  MainType Bi_pml0_2; MainType Bi_pml1_2; MainType Hi_pml0_2; MainType Hi_pml1_2;

  indx ind;
  SrcIndx SrcInd;
  short int pmlix, pmliy, pmliz;
};
//------------------------------------------------------------------------------------------------------------
struct fieldsXDD: public fieldsDDD{ int iY; int iZ; long int it; bool isSource; };
struct fieldsDXD: public fieldsDDD{ int iX; int iZ; long int it; bool isSource; };
struct fieldsDDX: public fieldsDDD{ int iX; int iY; long int it; bool isSource; };
struct fieldsDXX: public fieldsDDD{ int iX;         long int it; bool isSource; };
struct fieldsXDX: public fieldsDDD{ int iY;         long int it; bool isSource; };
struct fieldsXXD: public fieldsDDD{ int iZ;         long int it; bool isSource; };
struct fieldsXSD: public fieldsDSD{ int iY; int iZ; long int it; bool isSource; };
struct fieldsXDS: public fieldsDDS{ int iY; int iZ; long int it; bool isSource; };
struct fieldsSXD: public fieldsSDD{ int iX; int iZ; long int it; bool isSource; };
struct fieldsDXS: public fieldsDDS{ int iX; int iZ; long int it; bool isSource; };
struct fieldsSDX: public fieldsSDD{ int iX; int iY; long int it; bool isSource; };
struct fieldsDSX: public fieldsDSD{ int iX; int iY; long int it; bool isSource; };
struct fieldsXSS: public fieldsDSS{ int iY; int iZ; long int it; bool isSource; };
struct fieldsSXS: public fieldsSDS{ int iX; int iZ; long int it; bool isSource; };
struct fieldsSSX: public fieldsSSD{ int iX; int iY; long int it; bool isSource; };
struct fieldsSXX: public fieldsSDD{ int iX;         long int it; bool isSource; };
struct fieldsXSX: public fieldsDSD{ int iY;         long int it; bool isSource; };
struct fieldsXXS: public fieldsDDS{ int iZ;         long int it; bool isSource; };
struct fieldsXXX: public fieldsDDD{                 long int it; bool isSource; };

#endif//R0_STRUCT_HPP

