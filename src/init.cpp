#include <stdio.h>
#include <Python.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <xmmintrin.h>
#include "params.hpp"
#include "cubeLR.hpp"
#include "r0strct.hpp"
#include "omp.h"
#include <aivlib/vectorD.hpp>
#include <aivlib/arrayTD.hpp>
#include <aivlib/indexD.hpp>

size_t ProblemSize=0;

#include "init.hpp"
#include "nLAinit.cpp"

FILE* numaLogs[8];
NUMAnodesStruct Bricks;

extern long long INITIALIZED;
extern long long BOUNDARY_CELLS;
extern double STARTTIME;
extern double cctime1;
extern double cctime2;

DropParams DropP;
DropParams* getDropP() { return &DropP; };
NUMAnodesStruct* getBricks() { return &Bricks; };

inline void init_index(indx& Ind, double x, double y, int iz) {
//  for(int iz; iz<2*PMLNzV+NzV; iz++) Ind.I[iz] = IndVac;

//  CoffStruct coffs;
//  ParsStructIso Fast(Rho, cSrcP, cSrcS);
//  coffs.set(Fast);
//  for(int i=0; i<3; i++) {
//    Ind.kS[i][0].sRef(iz) = coffs.KS[i][0]; Ind.kS[i][1].sRef(iz) = coffs.KS[i][1]; Ind.kS[i][2].sRef(iz) = coffs.KS[i][2];
//    Ind.kV[i][0].sRef(iz) = coffs.KV[i][0]; Ind.kV[i][1].sRef(iz) = coffs.KV[i][1]; Ind.kV[i][2].sRef(iz) = coffs.KV[i][2];
//    Ind.kT[i][0].sRef(iz) = coffs.KT[i][0]; Ind.kT[i][1].sRef(iz) = coffs.KT[i][1];
//  }
}
bool isNearFieldReg(indx& Ind, double norm[3]) {
  return false;
/*  // don't forget about norm.abs()=1 !!!
  norm[0]=0.0; norm[1]=0.0; norm[2]=0.0; 
  if(Ind.x>pars.NearXp || Ind.x<pars.NearXm || Ind.y>pars.NearYp || Ind.y<pars.NearYm || Ind.z>pars.NearZp || Ind.z<pars.NearZm) { return false; }
  bool retval=false;
  if(Ind.x==pars.NearXm) { norm[0]=-1.0; retval=true; }
  if(Ind.x==pars.NearXp) { norm[0]= 1.0; retval=true; }
  if(Ind.y==pars.NearYm) { norm[1]=-1.0; retval=true; }
  if(Ind.y==pars.NearYp) { norm[1]= 1.0; retval=true; }
  if(Ind.z==pars.NearZm) { norm[2]=-1.0; retval=true; }
  if(Ind.z==pars.NearZp) { norm[2]= 1.0; retval=true; }
  return retval;*/
}

inline void init_index(indx& Ind, int x, int y, int z, int& M0, int& M1) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  coordinates c(x,y,z);
  vector<coordinates>::iterator it;
  it = find(DropP.Sensors.begin(), DropP.Sensors.end(), c);
  if(it!=DropP.Sensors.end()) Ind.isDrop = distance(DropP.Sensors.begin(), it);
  else Ind.isDrop=-1;
  Ind.isBound=0;
  
  MainType p[3]={x*pars.dx,y*pars.dy,z*pars.dz};
  Ind.I = getMaterial(p);
//  if(DispArr[Ind.I].isSet==1) Ind.isDisp=1;
//  else Ind.isDisp=0;

//  return;
  if(MatArr[Ind.I].isTensor==1) Ind.isBound=1; 
  if (pars.subpixel) {

    float averR[3] = {pars.averRx,pars.averRy,pars.averRz};
    float mCorn[3] = {(x-1)*pars.dx-averR[0],(y-1)*pars.dy-averR[1],(z-1)*pars.dz-averR[2]};
    M0 = getMaterial(mCorn); M1=M0;
    for(int i=1;i<8;i++) {
      float curCorn[3]={mCorn[0]+(2*pars.dx+averR[0])*(i&1), mCorn[1]+(2*pars.dy+averR[1])*((i>>1)&1), mCorn[2]+(2*pars.dz+averR[2])*((i>>2)&1)};
      int M = getMaterial(curCorn);
      if (M!=M0) { Ind.isBound=1; M1=M; if(DispArr[M0].isSet==1 || DispArr[M].isSet==1) break; }
    }
 
    if (Ind.isBound==0) Ind.I = M0;
    else {
//    #pragma omp critical
  	  {
      float cent[3] = {x*pars.dx,y*pars.dy,z*pars.dz};
      Ind.I = MatArr.size();
      MatStruct lastArr;
      CoffStruct cfArr;
      cfArr.set(); lastArr.set(MatArr[M0].eps,MatArr[M0].mu);
      lastArr.count_coffs(cent, averR);
      map<MatStruct, int>::iterator it = MaterialMap.find(lastArr);

      if (it==MaterialMap.end()) {
        MatArr.push_back(lastArr); CoffArr.push_back(cfArr);
        MaterialMap[lastArr] = Ind.I;
      }
      else {
        Ind.I = MaterialMap[lastArr];
      }
	    }
    }
  }
  
//  Ind.x=x; Ind.y=y; Ind.z=z; Ind.time=0;
//  #pragma omp critical 
  {
  if(Ind.isBound) BOUNDARY_CELLS++;
  long long int Xn=(1<<MaxRank)*Nx+(1<<PMLrank)+(1<<PMLrank);
  long long int Yn=(1<<MaxRank)*Ny+(1<<PMLrank)+(1<<PMLrank);
  long long int Zn=(1<<MaxRank)*Nz+(1<<PMLrank)+(1<<PMLrank);
  double time4sp=0.;
  if(INITIALIZED%10000==0 || INITIALIZED+1==Xn*Yn*Zn) {
    time4sp = omp_get_wtime()-STARTTIME;
    double speed = INITIALIZED/time4sp;
    printf("\rcompleted %d\% | index=%d | bd cells=%d | VolMap=%d | %9.8g cells/sec | cc1=%6.5g | cc2=%6.5g | cc3=%6.5g",
              50*(INITIALIZED+1)/(Xn*Yn*Zn),MatArr.size(), BOUNDARY_CELLS, VolumeMap.size(), speed, cctime1, cctime2, cctime3); 
    fflush(stdout);
  }
//  if(time4sp>1.) STARTTIME=omp_get_wtime();
  INITIALIZED++;
  }
}
inline void init_AuxIndex(AuxField* F, int x, int y, int z, int M0, int M1, float f1[3], float f2[3]){
  for(int i: {0,1,2}) {
    MainType crd[3] = {(x+1)*pars.dx-0.5*pars.dx*(i!=0), (y+1)*pars.dy-0.5*pars.dy*(i!=1), (z+1)*pars.dz-0.5*pars.dz*(i!=2)};
    F[0].I[i] = getMaterial(crd);
    F[1].I[i] = (F[0].I[i]==M1)?M0:M1;
  }
  for(int i: {0,1,2}) {
    F[2].I[i] = DispArr.size();
    DispStruct* last = new DispStruct;
    last->set(DispArr[F[0].I[i]], DispArr[F[1].I[i]], f1[i], f2[i]); // set
    int isexists = DispMaterialMap.count(*last);
    if (isexists) {
      F[2].I[i] = DispMaterialMap[*last];
      delete last;
    }
    else {
      DispArr.push_back(*last);
      DispMaterialMap[*last] = F[2].I[i];
    }
  }
}
inline void init_indexSrc(SrcIndx& Ind, int x, int y, int z){
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--; Ny--; Nz--;
  Ind.iX = x; Ind.iY = y;
  Ind.it = 0;
  Ind.isSource = (z==(1<<MaxRank)*Nz);
}
bool isSourceInside(indx& Ind) {
  /*if(Ind.z==pars.Nz*(1<<MaxRank)-2) return true;
  else*/ return false;
}

inline void init(AuxField*& F,int x,int y,int z, int& M0, int& M1, float f1[3], float f2[3]){
/*  F = new AuxField[3]; for(int i:{0,1,2}) for(int ai=0;ai<3;ai++) { F[ai].Ei[i] = 0.; F[ai].Eim[i] = 0.; F[ai].rot[i]=0.; for(int ip=0; ip<Md; ip++) { F[ai].Ji[i][ip] = 0.; F[ai].Jim[i][ip] = 0.; }  }
  init_AuxIndex(F, x,y,z, M0,M1, f1,f2);*/
}

void init(fieldsDDD& F, int x, int y, int z) {
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.;
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.;
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.; }
  int M0,M1; init_index(F.ind, x, y, z, M0,M1); 
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
}

void init(fieldsSDD& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.; 
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_1 = 0.; F.Di_pml1_1 = 0.; F.Ei_pml0_1 = 0.; F.Ei_pml1_1 = 0.;
  F.Di_pml0_2 = 0.; F.Di_pml1_2 = 0.; F.Ei_pml0_2 = 0.; F.Ei_pml1_2 = 0.;
  F.Bi_pml0_1 = 0.; F.Bi_pml1_1 = 0.; F.Hi_pml0_1 = 0.; F.Hi_pml1_1 = 0.;
  F.Bi_pml0_2 = 0.; F.Bi_pml1_2 = 0.; F.Hi_pml0_2 = 0.; F.Hi_pml1_2 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1); 
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (x<0) F.pmlix = x+(1<<PMLrank); else F.pmlix = x+(1<<PMLrank)-(1<<MaxRank)*Nx; 
};
void init(fieldsDSD& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.;
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_0 = 0.; F.Di_pml1_0 = 0.; F.Ei_pml0_0 = 0.; F.Ei_pml1_0 = 0.;
  F.Di_pml0_2 = 0.; F.Di_pml1_2 = 0.; F.Ei_pml0_2 = 0.; F.Ei_pml1_2 = 0.;
  F.Bi_pml0_0 = 0.; F.Bi_pml1_0 = 0.; F.Hi_pml0_0 = 0.; F.Hi_pml1_0 = 0.;
  F.Bi_pml0_2 = 0.; F.Bi_pml1_2 = 0.; F.Hi_pml0_2 = 0.; F.Hi_pml1_2 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1);
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (y<0) F.pmliy = y+(1<<PMLrank); else F.pmliy = y+(1<<PMLrank)-(1<<MaxRank)*Ny; 
};
void init(fieldsDDS& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.;
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_0 = 0.; F.Di_pml1_0 = 0.; F.Ei_pml0_0 = 0.; F.Ei_pml1_0 = 0.;
  F.Di_pml0_1 = 0.; F.Di_pml1_1 = 0.; F.Ei_pml0_1 = 0.; F.Ei_pml1_1 = 0.;
  F.Bi_pml0_0 = 0.; F.Bi_pml1_0 = 0.; F.Hi_pml0_0 = 0.; F.Hi_pml1_0 = 0.;
  F.Bi_pml0_1 = 0.; F.Bi_pml1_1 = 0.; F.Hi_pml0_1 = 0.; F.Hi_pml1_1 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1);
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (z<0) F.pmliz = z+(1<<PMLrank); else F.pmliz = z+(1<<PMLrank)-(1<<MaxRank)*Nz; 
};
void init(fieldsDSS& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.;
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_0 = 0.; F.Di_pml1_0 = 0.; F.Ei_pml0_0 = 0.; F.Ei_pml1_0 = 0.;
  F.Di_pml0_1 = 0.; F.Di_pml1_1 = 0.; F.Ei_pml0_1 = 0.; F.Ei_pml1_1 = 0.;
  F.Di_pml0_2 = 0.; F.Di_pml1_2 = 0.; F.Ei_pml0_2 = 0.; F.Ei_pml1_2 = 0.;
  F.Bi_pml0_0 = 0.; F.Bi_pml1_0 = 0.; F.Hi_pml0_0 = 0.; F.Hi_pml1_0 = 0.;
  F.Bi_pml0_1 = 0.; F.Bi_pml1_1 = 0.; F.Hi_pml0_1 = 0.; F.Hi_pml1_1 = 0.;
  F.Bi_pml0_2 = 0.; F.Bi_pml1_2 = 0.; F.Hi_pml0_2 = 0.; F.Hi_pml1_2 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1);
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (y<0) F.pmliy = y+(1<<PMLrank); else F.pmliy = y+(1<<PMLrank)-(1<<MaxRank)*Ny; 
  if (z<0) F.pmliz = z+(1<<PMLrank); else F.pmliz = z+(1<<PMLrank)-(1<<MaxRank)*Nz; 
};
void init(fieldsSDS& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.;
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_0 = 0.; F.Di_pml1_0 = 0.; F.Ei_pml0_0 = 0.; F.Ei_pml1_0 = 0.;
  F.Di_pml0_1 = 0.; F.Di_pml1_1 = 0.; F.Ei_pml0_1 = 0.; F.Ei_pml1_1 = 0.;
  F.Di_pml0_2 = 0.; F.Di_pml1_2 = 0.; F.Ei_pml0_2 = 0.; F.Ei_pml1_2 = 0.;
  F.Bi_pml0_0 = 0.; F.Bi_pml1_0 = 0.; F.Hi_pml0_0 = 0.; F.Hi_pml1_0 = 0.;
  F.Bi_pml0_1 = 0.; F.Bi_pml1_1 = 0.; F.Hi_pml0_1 = 0.; F.Hi_pml1_1 = 0.;
  F.Bi_pml0_2 = 0.; F.Bi_pml1_2 = 0.; F.Hi_pml0_2 = 0.; F.Hi_pml1_2 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1);
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (x<0) F.pmlix = x+(1<<PMLrank); else F.pmlix = x+(1<<PMLrank)-(1<<MaxRank)*Nx; 
  if (z<0) F.pmliz = z+(1<<PMLrank); else F.pmliz = z+(1<<PMLrank)-(1<<MaxRank)*Nz; 
};
void init(fieldsSSD& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.;
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_0 = 0.; F.Di_pml1_0 = 0.; F.Ei_pml0_0 = 0.; F.Ei_pml1_0 = 0.;
  F.Di_pml0_1 = 0.; F.Di_pml1_1 = 0.; F.Ei_pml0_1 = 0.; F.Ei_pml1_1 = 0.;
  F.Di_pml0_2 = 0.; F.Di_pml1_2 = 0.; F.Ei_pml0_2 = 0.; F.Ei_pml1_2 = 0.;
  F.Bi_pml0_0 = 0.; F.Bi_pml1_0 = 0.; F.Hi_pml0_0 = 0.; F.Hi_pml1_0 = 0.;
  F.Bi_pml0_1 = 0.; F.Bi_pml1_1 = 0.; F.Hi_pml0_1 = 0.; F.Hi_pml1_1 = 0.;
  F.Bi_pml0_2 = 0.; F.Bi_pml1_2 = 0.; F.Hi_pml0_2 = 0.; F.Hi_pml1_2 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1);
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (x<0) F.pmlix = x+(1<<PMLrank); else F.pmlix = x+(1<<PMLrank)-(1<<MaxRank)*Nx; 
  if (y<0) F.pmliy = y+(1<<PMLrank); else F.pmliy = y+(1<<PMLrank)-(1<<MaxRank)*Ny; 
};
void init(fieldsSSS& F, int x, int y, int z) {
  int Nx=Bricks.mapNUMAnodes[0][0].size();
  int Ny=Bricks.mapNUMAnodes[0].size();
  int Nz=Bricks.mapNUMAnodes.size();
  Nx--;Ny--;Nz--;
  for (int i=0;i<3;i++) { F.Di[i] = 0.; F.Ei[i] = 0.; F.Bi[i] = 0.; F.Hi[i] = 0.; }
//                          F.EHi_m[i] = 0.; F.EHi_p[i] = 0.; F.Si[i] = 0.; 
//                          F.avHx[i] = 0.;  F.avHy[i] = 0.;  F.avHz[i] = 0.; 
//                          F.Eim[i] = 0.; for(int ip=0; ip<Md; ip++) { F.Ji[i][ip] = 0.; F.Jim[i][ip] = 0.; } F.rot[i] = 0.;  }
  F.Di_pml0_0 = 0.; F.Di_pml1_0 = 0.; F.Ei_pml0_0 = 0.; F.Ei_pml1_0 = 0.;
  F.Di_pml0_1 = 0.; F.Di_pml1_1 = 0.; F.Ei_pml0_1 = 0.; F.Ei_pml1_1 = 0.;
  F.Di_pml0_2 = 0.; F.Di_pml1_2 = 0.; F.Ei_pml0_2 = 0.; F.Ei_pml1_2 = 0.;
  F.Bi_pml0_0 = 0.; F.Bi_pml1_0 = 0.; F.Hi_pml0_0 = 0.; F.Hi_pml1_0 = 0.;
  F.Bi_pml0_1 = 0.; F.Bi_pml1_1 = 0.; F.Hi_pml0_1 = 0.; F.Hi_pml1_1 = 0.;
  F.Bi_pml0_2 = 0.; F.Bi_pml1_2 = 0.; F.Hi_pml0_2 = 0.; F.Hi_pml1_2 = 0.;
  int M0,M1; init_index(F.ind ,x, y, z, M0,M1);
  int M0_p, M1_p; indx Ind_p; init_index(Ind_p,x+1, y+1, z+1, M0_p,M1_p);
  init_indexSrc(F.SrcInd, x, y, z);
//  if(Ind_p.isBound && Ind_p.isDisp) init(F.auxfld,x,y,z, M0_p,M1_p, MatArr[Ind_p.I].f1, MatArr[Ind_p.I].f2);
  if (x<0) F.pmlix = x+(1<<PMLrank); else F.pmlix = x+(1<<PMLrank)-(1<<MaxRank)*Nx; 
  if (y<0) F.pmliy = y+(1<<PMLrank); else F.pmliy = y+(1<<PMLrank)-(1<<MaxRank)*Ny; 
  if (z<0) F.pmliz = z+(1<<PMLrank); else F.pmliz = z+(1<<PMLrank)-(1<<MaxRank)*Nz; 
};
void init(fieldsDDX& F, int x, int y, int z) { 
  init((fieldsDDD&)F, x, y, z); 
  F.iX = x; F.iY = y; 
  F.it=0; F.isSource=true;
};
void init(fieldsSDX& F, int x, int y, int z) {
  init((fieldsSDD&)F, x, y, z); 
  F.iX = x; F.iY = y; 
  F.it=0; F.isSource=false;
};
void init(fieldsDSX& F, int x, int y, int z) {
  init((fieldsDSD&)F, x, y, z); 
  F.iX = x; F.iY = y; 
  F.it=0; F.isSource=false;
};
void init(fieldsSSX& F, int x, int y, int z) {
  init((fieldsSSD&)F, x, y, z); 
  F.iX = x; F.iY = y; 
  F.it=0; F.isSource=false;
};

//-----------------------------------------------------------------------------------------------------//

template <class T> void drop(T& F, FILE* file, int comp) {
  MainType fld=-1000;
  if(DropP.Fields4Drop[comp]=="Ex") fld=F.Ei[0]; else
  if(DropP.Fields4Drop[comp]=="Ey") fld=F.Ei[1]; else
  if(DropP.Fields4Drop[comp]=="Ez") fld=F.Ei[2]; else
//    if(DropP.Fields4Drop[comp]== "Dx") fld=F.Di[0]; else
//    if(DropP.Fields4Drop[comp]== "Dy") fld=F.Di[1]; else
//    if(DropP.Fields4Drop[comp]== "Dz") fld=F.Di[2]; else
//    if(DropP.Fields4Drop[comp]== "Hx") fld=F.Hi[0]; else
//    if(DropP.Fields4Drop[comp]== "Hy") fld=F.Hi[1]; else
//    if(DropP.Fields4Drop[comp]== "Hz") fld=F.Hi[2]; else
//    if(DropP.Fields4Drop[comp]== "Bx") fld=F.Bi[0]; else
//    if(DropP.Fields4Drop[comp]== "By") fld=F.Bi[1]; else
//    if(DropP.Fields4Drop[comp]== "Bz") fld=F.Bi[2]; else
    if(DropP.Fields4Drop[comp]=="Ef") fld=sqrt(F.Ei[0]*F.Ei[0]+F.Ei[1]*F.Ei[1]+F.Ei[2]*F.Ei[2]); else
//    if(DropP.Fields4Drop[comp]== "Df") fld=sqrt(F.Di[0]*F.Di[0]+F.Di[1]*F.Di[1]+F.Di[2]*F.Di[2]); else
//    if(DropP.Fields4Drop[comp]== "Bf") fld=sqrt(F.Bi[0]*F.Bi[0]+F.Bi[1]*F.Bi[1]+F.Bi[2]*F.Bi[2]); else
//    if(DropP.Fields4Drop[comp]== "Hf") fld=sqrt(F.Hi[0]*F.Hi[0]+F.Hi[1]*F.Hi[1]+F.Hi[2]*F.Hi[2]); else
//    if(DropP.Fields4Drop[comp]== "Sx")  fld=F.Si[0]; else
//    if(DropP.Fields4Drop[comp]== "Sy")  fld=F.Si[1]; else
//    if(DropP.Fields4Drop[comp]== "Sz")  fld=F.Si[2]; else
//    if(DropP.Fields4Drop[comp]== "Mat") fld=F.ind.I; else
    if(DropP.Fields4Drop[comp]== "Mat") fld = F.ind.I; else 
//	{ if (F.ind.isBound) 
//		fld = 1./3.*(1./MatArr[F.ind.I].depsXX + 1./MatArr[F.ind.I].depsYY + 1./MatArr[F.ind.I].depsZZ); 
//	else fld = MatArr[F.ind.I].eps; }
//    if(DropP.Fields4Drop[comp]== "DepsXX") fld = MatArr[F.ind.I].depsXX; else
//    if(DropP.Fields4Drop[comp]== "DepsYY") fld = MatArr[F.ind.I].depsYY; else
//    if(DropP.Fields4Drop[comp]== "DepsZZ") fld = MatArr[F.ind.I].depsZZ; else
//    if(DropP.Fields4Drop[comp]== "DepsXY") fld = 0.5*(MatArr[F.ind.I].depsXY1+MatArr[F.ind.I].depsXY2); else
//    if(DropP.Fields4Drop[comp]== "DepsYX") fld = 0.5*(MatArr[F.ind.I].depsYX1+MatArr[F.ind.I].depsYX2); else
//    if(DropP.Fields4Drop[comp]== "DepsXZ") fld = 0.5*(MatArr[F.ind.I].depsXZ1+MatArr[F.ind.I].depsXZ2); else
//    if(DropP.Fields4Drop[comp]== "DepsZX") fld = 0.5*(MatArr[F.ind.I].depsZX1+MatArr[F.ind.I].depsZX2); else
//    if(DropP.Fields4Drop[comp]== "DepsYZ") fld = 0.5*(MatArr[F.ind.I].depsYZ1+MatArr[F.ind.I].depsYZ2); else
//    if(DropP.Fields4Drop[comp]== "DepsZY") fld = 0.5*(MatArr[F.ind.I].depsZY1+MatArr[F.ind.I].depsZY2); else
  {}
  fwrite(&fld,1*sizeof(MainType),1,file);
}
template<> void drop(fieldsDDX& F, FILE* file, int comp) {};
template<> void drop(fieldsSDX& F, FILE* file, int comp) { };
template<> void drop(fieldsDSX& F, FILE* file, int comp) { };
template<> void drop(fieldsSSX& F, FILE* file, int comp) { };

#include <algorithm>
#include <boost/serialization/map.hpp>

aiv::array<float,3> PPoutXm[500][6]; aiv::array<float,3> dPPoutXm[500][6]; 
aiv::array<float,3> PPoutXp[500][6]; aiv::array<float,3> dPPoutXp[500][6];
aiv::array<float,3> PPoutYm[500][6]; aiv::array<float,3> dPPoutYm[500][6]; 
aiv::array<float,3> PPoutYp[500][6]; aiv::array<float,3> dPPoutYp[500][6];
aiv::array<float,3> PPoutZm[500][6]; aiv::array<float,3> dPPoutZm[500][6]; 
aiv::array<float,3> PPoutZp[500][6]; aiv::array<float,3> dPPoutZp[500][6];
pthread_mutex_t NtF_mutex=PTHREAD_MUTEX_INITIALIZER;

#include "SaveNf.cpp"
int InitializeData(){
  _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON );
  _MM_SET_ROUNDING_MODE (  _MM_ROUND_NEAREST );

  nLAarr=reset_numa(Bricks.mapNUMAnodes, Bricks.numNUMAnodes);

//-----------------------------------Serialization of Volumes----------------------------------------//
  char s[256]; sprintf(s,"dump/VolumeMap-deep%d-res%g-(%gx%gx%g)-DIGEQUIVALENT.dump",pars.deep3d,log(pars.mapResolution),pars.averRx,pars.averRy,pars.averRz);
  ifstream ifs(s);
  if(ifs.good()) {
    printf("Using %s... ",s); fflush(stdout);
    boost::archive::binary_iarchive from(ifs);
    from >> VolumeMap;
    printf("VolumeMap loaded\n");
  }
//  exit(-1);
// ----------------------------Initialising-Data-----------------------------------------------------//
  INITIALIZED=0;   // progress if initializing
  BOUNDARY_CELLS=0;   // the size of boundary area
  STARTTIME=omp_get_wtime();
  cctime1=0;
  cctime2=0;
  cctime3=0;
  int tt=omp_get_wtime();
//  numa_run_on_nodes(&nLAnetArrayNUMA_init, nLAarr, numNUMAnodes);

  nLAnetArrayNUMA<3> ATs[Bricks.numNUMAnodes];

  for(int ithr=0; ithr<Bricks.numNUMAnodes; ithr++) {
    ATs[ithr].prepare(ithr, Bricks.numNUMAnodes);
    ATs[ithr].base_net = nLAarr+ithr;
    ATs[ithr].iNode = ithr;
    nLAnetArrayNUMA_init(ATs+ithr);
//    printf("nLAnetArrayNUMA_init ok\n");
  }
  
  printf("-----------------------------------------------------\n");
  numa_run_on_nodes(&nLAnetArrayNUMA_prepChessFold, nLAarr, Bricks.numNUMAnodes);
  
  ofstream ofs(s);
  boost::archive::binary_oarchive to(ofs);
  to << VolumeMap;
  printf("\nTime of initializing=%9.8gs\n",omp_get_wtime()-tt);
  
  char fname[256];
  vector<string>::iterator sfield_it; vector<coordinates>::iterator sen_it;
  for(int sf =0; sf <DropP.Fields4Sensors.size(); sf++ )
  for(int sen=0; sen<DropP.Sensors.size()       ; sen++) {
    sprintf(fname,"%s/tt%s-x%04dy%04dz%04d.dat",DropP.DropPath.c_str(),DropP.Fields4Sensors[sf].c_str(),DropP.Sensors[sen].x,DropP.Sensors[sen].y,DropP.Sensors[sen].z);
    remove(fname);
  }
//-------------------------------------------------------------------------------------------//
/*  aiv::indx<3> tmpind(0); tmpind[2]=(1<<MaxRank)/DropP.timesep;
  for(int c=0; c<6; c++) {
    tmpind[0]=pars.Ny*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutXm[0][c].init( tmpind ); dPPoutXm[0][c].init( tmpind );
    tmpind[0]=pars.Ny*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutXp[0][c].init( tmpind ); dPPoutXp[0][c].init( tmpind );
    tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutYm[0][c].init( tmpind ); dPPoutYm[0][c].init( tmpind );
    tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutYp[0][c].init( tmpind ); dPPoutYp[0][c].init( tmpind );
    tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Ny*(1<<MaxRank); PPoutZm[0][c].init( tmpind ); dPPoutZm[0][c].init( tmpind );
    tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Ny*(1<<MaxRank); PPoutZp[0][c].init( tmpind ); dPPoutZp[0][c].init( tmpind );
  }*/
//-------------------------------------------------------------------------------------------//
  printf("Size of one LRnLA cell: %d bytes\n",sizeof(fieldsDDD));
  printf("Sizeof PPouts: %d Mbytes\n",(pars.Nx*pars.Ny+pars.Nx*pars.Nz+pars.Ny*pars.Nz)*2*(1<<(3*MaxRank)*2*6*sizeof(float)/1024/1024));
  printf("Size of CoffArr struct: %d Bytes * %d + head %d bytes\n", sizeof(CoffStruct), CoffArr.size(), sizeof(CoffArr));
  printf("Size of MatArr  struct: %d Bytes * %d + head %d bytes\n", sizeof(MatStruct) , MatArr.size() , sizeof(MatArr) );
  printf("Size of DispArr struct: %d Bytes * %d + head %d bytes\n", sizeof(DispStruct), (int)(DispArr.size()), sizeof(DispArr));
  printf("VolumeMap      : (%d Bytes -> %d Bytes )* %d + head %d bytes\n", sizeof(keyVolStr), sizeof(Vols), VolumeMap.size(), sizeof(VolumeMap));
  printf("MaterialMap    : (%d Bytes -> %d Bytes )* %d + head %d bytes\n", sizeof(MatStruct), sizeof(int), MaterialMap.size(), sizeof(MaterialMap));
  printf("DispMaterialMap: (%d Bytes -> %d Bytes )* %d + head %d bytes\n", sizeof(DispStruct), sizeof(int), DispMaterialMap.size(), sizeof(DispMaterialMap));
}
//--------------------------------------------------------------------------------------------------//

int InitializeCoffs(){ 
// CoffArr   = (CoffStruct*)memalign(64, sizeof(CoffStruct)*100);
// MatArr    = ( MatStruct*)memalign(64, sizeof(MatStruct )*100);
  KpmlX = (PMLcoffs*)memalign(64, sizeof(PMLcoffs)*4*PMLN);
  KpmlY = (PMLcoffs*)memalign(64, sizeof(PMLcoffs)*4*PMLN);
  KpmlZ = (PMLcoffs*)memalign(64, sizeof(PMLcoffs)*4*PMLN);
  printf("sizeof CoffStruct %lu\n", sizeof(CoffStruct));
  printf("sizeof MatStruct  %lu\n", sizeof(MatStruct));
  printf("sizeof PMLcoffs   %lu\n", sizeof(PMLcoffs));
  for(int i=0; i<4*PMLN; i++) KpmlX[i].setX(i);
  for(int i=0; i<4*PMLN; i++) KpmlY[i].setY(i);
  for(int i=0; i<4*PMLN; i++) KpmlZ[i].setZ(i);
  //setCoffsIso(IndVac,1.,1.,0.,0.,0.,0.);
  //setCoffsIso(IndMat,3.,1.,0.,0.,0.,0.);
}

//-----------------------------------------------------------------------------------------------------//

set<int> OutsDropped;
set<int> OutsLoaded;

int loadOuts(int it) {
  char fname[256];
    if(OutsDropped.find(it)!=OutsDropped.end()){
      for(int c=0; c<6; c++) {
        sprintf(fname, "%s/PPoutXm%d_%05d.arr", DropP.DropPath.c_str(),c,it);  PPoutXm[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/PPoutXp%d_%05d.arr", DropP.DropPath.c_str(),c,it);  PPoutXp[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/PPoutYm%d_%05d.arr", DropP.DropPath.c_str(),c,it);  PPoutYm[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/PPoutYp%d_%05d.arr", DropP.DropPath.c_str(),c,it);  PPoutYp[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/PPoutZm%d_%05d.arr", DropP.DropPath.c_str(),c,it);  PPoutZm[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/PPoutZp%d_%05d.arr", DropP.DropPath.c_str(),c,it);  PPoutZp[it][c].load( aiv::Ifile(fname).self );
        
        sprintf(fname, "%s/dPPoutXm%d_%05d.arr",DropP.DropPath.c_str(),c,it); dPPoutXm[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/dPPoutXp%d_%05d.arr",DropP.DropPath.c_str(),c,it); dPPoutXp[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/dPPoutYm%d_%05d.arr",DropP.DropPath.c_str(),c,it); dPPoutYm[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/dPPoutYp%d_%05d.arr",DropP.DropPath.c_str(),c,it); dPPoutYp[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/dPPoutZm%d_%05d.arr",DropP.DropPath.c_str(),c,it); dPPoutZm[it][c].load( aiv::Ifile(fname).self );
        sprintf(fname, "%s/dPPoutZp%d_%05d.arr",DropP.DropPath.c_str(),c,it); dPPoutZp[it][c].load( aiv::Ifile(fname).self );
      }
    }
    else {
      aiv::indx<3> tmpind(0); tmpind[2]=(1<<MaxRank)/DropP.timesep;
      for(int c=0; c<6; c++) {
        tmpind[0]=pars.Ny*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutXm[it][c].init( tmpind ); dPPoutXm[it][c].init( tmpind );
        tmpind[0]=pars.Ny*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutXp[it][c].init( tmpind ); dPPoutXp[it][c].init( tmpind );
        tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutYm[it][c].init( tmpind ); dPPoutYm[it][c].init( tmpind );
        tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Nz*(1<<MaxRank); PPoutYp[it][c].init( tmpind ); dPPoutYp[it][c].init( tmpind );
        tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Ny*(1<<MaxRank); PPoutZm[it][c].init( tmpind ); dPPoutZm[it][c].init( tmpind );
        tmpind[0]=pars.Nx*(1<<MaxRank); tmpind[1]=pars.Ny*(1<<MaxRank); PPoutZp[it][c].init( tmpind ); dPPoutZp[it][c].init( tmpind );
      }
    }
    OutsLoaded.insert(it);
}

void dropPP(){
  double tstart = omp_get_wtime();
  char fname[256];
  set<int>::iterator it;
  for(it=OutsLoaded.begin(); it!=OutsLoaded.end(); ++it) {
    aiv::indx<3> tmpind(0); tmpind[0]=0; tmpind[1]=0; tmpind[2]=0;
    for(int i=0; i<6; i++) {
      sprintf(fname, "%s/PPoutXm%d_%05d.arr", DropP.DropPath.c_str(),i,*it);  PPoutXm[*it][i].dump( aiv::Ofile(fname).self );  PPoutXm[*it][i].init(tmpind);
      sprintf(fname, "%s/PPoutXp%d_%05d.arr", DropP.DropPath.c_str(),i,*it);  PPoutXp[*it][i].dump( aiv::Ofile(fname).self );  PPoutXp[*it][i].init(tmpind);
      sprintf(fname, "%s/PPoutYm%d_%05d.arr", DropP.DropPath.c_str(),i,*it);  PPoutYm[*it][i].dump( aiv::Ofile(fname).self );  PPoutYm[*it][i].init(tmpind);
      sprintf(fname, "%s/PPoutYp%d_%05d.arr", DropP.DropPath.c_str(),i,*it);  PPoutYp[*it][i].dump( aiv::Ofile(fname).self );  PPoutYp[*it][i].init(tmpind);
      sprintf(fname, "%s/PPoutZm%d_%05d.arr", DropP.DropPath.c_str(),i,*it);  PPoutZm[*it][i].dump( aiv::Ofile(fname).self );  PPoutZm[*it][i].init(tmpind);
      sprintf(fname, "%s/PPoutZp%d_%05d.arr", DropP.DropPath.c_str(),i,*it);  PPoutZp[*it][i].dump( aiv::Ofile(fname).self );  PPoutZp[*it][i].init(tmpind);
                                                                                                                                                   
      sprintf(fname, "%s/dPPoutXm%d_%05d.arr",DropP.DropPath.c_str(),i,*it); dPPoutXm[*it][i].dump( aiv::Ofile(fname).self ); dPPoutXm[*it][i].init(tmpind);
      sprintf(fname, "%s/dPPoutXp%d_%05d.arr",DropP.DropPath.c_str(),i,*it); dPPoutXp[*it][i].dump( aiv::Ofile(fname).self ); dPPoutXp[*it][i].init(tmpind);
      sprintf(fname, "%s/dPPoutYm%d_%05d.arr",DropP.DropPath.c_str(),i,*it); dPPoutYm[*it][i].dump( aiv::Ofile(fname).self ); dPPoutYm[*it][i].init(tmpind);
      sprintf(fname, "%s/dPPoutYp%d_%05d.arr",DropP.DropPath.c_str(),i,*it); dPPoutYp[*it][i].dump( aiv::Ofile(fname).self ); dPPoutYp[*it][i].init(tmpind);
      sprintf(fname, "%s/dPPoutZm%d_%05d.arr",DropP.DropPath.c_str(),i,*it); dPPoutZm[*it][i].dump( aiv::Ofile(fname).self ); dPPoutZm[*it][i].init(tmpind);
      sprintf(fname, "%s/dPPoutZp%d_%05d.arr",DropP.DropPath.c_str(),i,*it); dPPoutZp[*it][i].dump( aiv::Ofile(fname).self ); dPPoutZp[*it][i].init(tmpind);
    }
    OutsDropped.insert(*it);
    for(int c=0; c<6; c++) {
      PPoutXm[*it][c].init( tmpind ); dPPoutXm[*it][c].init( tmpind );
      PPoutXp[*it][c].init( tmpind ); dPPoutXp[*it][c].init( tmpind );
      PPoutYm[*it][c].init( tmpind ); dPPoutYm[*it][c].init( tmpind );
      PPoutYp[*it][c].init( tmpind ); dPPoutYp[*it][c].init( tmpind );
      PPoutZm[*it][c].init( tmpind ); dPPoutZm[*it][c].init( tmpind );
      PPoutZp[*it][c].init( tmpind ); dPPoutZp[*it][c].init( tmpind );
    }
    OutsLoaded.erase(*it);
  }
  printf("Time for drop PPouts %gs\n",omp_get_wtime()-tstart);
}

