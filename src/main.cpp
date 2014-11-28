#include <stdio.h>
#include <sys/sysinfo.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sched.h>
//#include <numaif.h>
#include <fstream>
#include <iostream>
#include "omp.h"
#include <time.h>
#include <malloc.h>
#include <vector>
#include <deque>
#include <map>
#include <string>
#include <algorithm>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <xmmintrin.h>
#include <aivlib/arrayTD.hpp>
#include <aivlib/indexD.hpp>
#include <aivlib/vectorD.hpp>
#include "params.hpp"

#include "cubeLR.hpp"
#include "Coff.hpp"
#include "r0strct.hpp"
#include "Signal.hpp"

#include "Data_def.hpp"
#include "main.hpp"

//Global params
PMLcoffs* KpmlX=0;
PMLcoffs* KpmlY=0;
PMLcoffs* KpmlZ=0;
vector<CoffStruct> CoffArr;
vector<DispStruct> DispArr;
vector<MatStruct> MatArr;
map<keyVolStr, Vols> VolumeMap;
map<MatStruct, int> MaterialMap;
map<DispStruct, int> DispMaterialMap;
long long int INITIALIZED;   // progress if initializing
long long int BOUNDARY_CELLS;   // the size of boundary area
double STARTTIME;
double cctime1;
double cctime2;
double cctime3;
parsMxw pars;
parsMxw* getPars() { return &pars; }

#include "SaveNf.cpp"
#include "UpdateJ.hpp"
#include "Update.inc.hpp"
#include "CF3Dpml.inc.hpp"

#define DIM_3D
#define LEFT_BC_COMPLEX
#include "nLAnet.hpp"
#include "nLAarr.hpp"
#include "nLAnuma.hpp"

#include "nLArun.cpp"

template<> RegionArray<3> netArray<3, ChessCell<3> >::region=0;
template<> int nLAnetArrayNUMA<3>::thrs4node = 1;
template<> int nLAnetArrayNUMA<3>::MaxStep = 0;

netArray<3, ChessCell<3> >* nLAarr;

using namespace std;
//using namespace aiv;

void setCoffsTensorInv(int ind, double* deps, double* dmu){
  if (ind+1>CoffArr.size()) { CoffArr.resize(ind+1); DispArr.resize(ind+1); MatArr.resize(ind+1); }
  CoffArr[ind].set();
  MatArr[ind].setTensorInv(deps,dmu);
  DispArr[ind].setNondisp(MatArr[ind]);
}

void setCoffsAniso(int ind, MainType* eps, MainType* mu){
  if (ind+1>CoffArr.size()) { CoffArr.resize(ind+1); DispArr.resize(ind+1); MatArr.resize(ind+1); }
  CoffArr[ind].set();
  MatArr[ind].set(eps,mu);
  DispArr[ind].setNondisp(MatArr[ind]);
}
void setCoffsIso(int ind, double eps, double mu){
  if (ind+1>CoffArr.size()) { CoffArr.resize(ind+1); DispArr.resize(ind+1); MatArr.resize(ind+1); }
  CoffArr[ind].set();
  MatArr[ind].set(eps,mu);
  DispArr[ind].setNondisp(MatArr[ind]);
}
void setCoffsDisp(int ind, MainType eps_inf, vector<LorenzDisp> LorenzD, MainType mu, MainType sigma) {
  if (ind+1>CoffArr.size()) { CoffArr.resize(ind+1); DispArr.resize(ind+1); MatArr.resize(ind+1); }
  CoffArr[ind].set();
  DispArr[ind].set(eps_inf, LorenzD, LorenzD.size(), sigma); 
  MatArr[ind].set(eps_inf,mu);
}

void parsMxw::setGrid(double _dx, double _dy, double _dz, double _dt){
    dx=_dx; dy=_dy; dz=_dz; dt=_dt;
    Nx=Bricks.Nx; Ny=Bricks.Ny; Nz=Bricks.Nz;
    Cen[0] = 0.5*(Nx*dx)*(1<<MaxRank);
    Cen[1] = 0.5*(Ny*dy)*(1<<MaxRank);
    Cen[2] = 0.5*(Nz*dz)*(1<<MaxRank);
    NearXm = 5; NearXp = Nx*(1<<MaxRank)-5;
    NearYm = 5; NearYp = Ny*(1<<MaxRank)-5;
    NearZm = 5; NearZp = Nz*(1<<MaxRank)-5;
    sigma_max = 11.*(3./(2*((1<<PMLrank)*dx-0.5*dx)));
    attenuation_factor = 2.; 
    averRx=0.5*dx; averRy=0.5*dy; averRz=0.5*dz; 
    mapResolution=1.e-100;
    deep1d=10; deep3d=10;
    subpixel=false;
    MaxNodes=1;
    Nthreads=8;
    nchip=0;
    BigSteps=Nthreads*1;
    swapDir="";
};

//inline void DDDact ( fieldsDDD* const Fdd, const int ixdd, const int iydd, const int izdd){  DDDactPML(Fdd,ixdd,iydd,izdd);   }
void BBB(void* pFldData, void* pFldSide[6], void* pFldEdge[12], void* pFldCorn[8]){
/*  RRXact(areaCCC, areaPCC, areaCPC, areaPPC, areaCCP, areaPCP, areaCPP, areaPPP);
  LRXact(areaMCC, areaCCC, areaMPC, areaCPC, areaMCP, areaCCP, areaMPP, areaCPP);
  RLXact(areaCMC, areaPMC, areaCCC, areaPCC, areaCMP, areaPMP, areaCCP, areaPCP);
  RRLact(areaCCM, areaPCM, areaCPM, areaPPM, areaCCC, areaPCC, areaCPC, areaPPC);
  LLXact(areaMMC, areaCMC, areaMCC, areaCCC, areaMMP, areaCMP, areaMCP, areaCCP);
  LRLact(areaMCM, areaCCM, areaMPM, areaCPM, areaMCC, areaCCC, areaMPC, areaCPC);
  RLLact(areaCMM, areaPMM, areaCCM, areaPCM, areaCMC, areaPMC, areaCCC, areaPCC);
  LLLact(areaMMM, areaCMM, areaMCM, areaCCM, areaMMC, areaCMC, areaMCC, areaCCC);*/
} 

/*void copyVec(int iz){
#ifdef VECTORIZATION
  for(size_t i=0; i<sizeof(areaCCC_type)/sizeof(fieldsDDD); i++) for(size_t f=0; f<sizeof(fieldsDDD)/sizeof(MainType); f++)
  (((fieldsDDD*)(areaCCC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsDDDvec*)(areaCCCvec(0,0,iz)))+i)[f]))[iz/VecL];
  for(size_t i=0; i<sizeof(areaMCC_type)/sizeof(fieldsSDD); i++) for(size_t f=0; f<sizeof(fieldsSDD)/sizeof(MainType); f++) {
  (((fieldsSDD*)(areaMCC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsSDDvec*)(areaMCCvec(0,0,iz)))+i)[f]))[iz/VecL];
  (((fieldsSDD*)(areaPCC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsSDDvec*)(areaPCCvec(0,0,iz)))+i)[f]))[iz/VecL];
  (((fieldsDSD*)(areaCMC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsDSDvec*)(areaCMCvec(0,0,iz)))+i)[f]))[iz/VecL];
  (((fieldsDSD*)(areaCPC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsDSDvec*)(areaCPCvec(0,0,iz)))+i)[f]))[iz/VecL];                }
  for(size_t i=0; i<sizeof(areaMMC_type)/sizeof(fieldsSSD); i++) for(size_t f=0; f<sizeof(fieldsSSD)/sizeof(MainType); f++) {
  (((fieldsSSD*)(areaMMC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsSSDvec*)(areaMMCvec(0,0,iz)))+i)[f]))[iz/VecL];
  (((fieldsSSD*)(areaPMC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsSSDvec*)(areaPMCvec(0,0,iz)))+i)[f]))[iz/VecL];
  (((fieldsSSD*)(areaMPC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsSSDvec*)(areaMPCvec(0,0,iz)))+i)[f]))[iz/VecL];
  (((fieldsSSD*)(areaPPC(0,0,iz)))+i)[f] = ((MainType*)&((((fieldsSSDvec*)(areaPPCvec(0,0,iz)))+i)[f]))[iz/VecL];                }
#endif
}
void copy(int iz){
#ifdef VECTORIZATION
  for(size_t i=0; i<sizeof(areaCCC_type)/sizeof(fieldsDDD); i++) for(size_t f=0; f<sizeof(fieldsDDD)/sizeof(MainType); f++)
  ((MainType*)&((((fieldsDDDvec*)(areaCCCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsDDD*)(areaCCC(0,0,iz)))+i)[f];
  for(size_t i=0; i<sizeof(areaMCC_type)/sizeof(fieldsSDD); i++) for(size_t f=0; f<sizeof(fieldsSDD)/sizeof(MainType); f++) {
  ((MainType*)&((((fieldsSDDvec*)(areaMCCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsSDD*)(areaMCC(0,0,iz)))+i)[f];
  ((MainType*)&((((fieldsSDDvec*)(areaPCCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsSDD*)(areaPCC(0,0,iz)))+i)[f];
  ((MainType*)&((((fieldsDSDvec*)(areaCMCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsDSD*)(areaCMC(0,0,iz)))+i)[f];
  ((MainType*)&((((fieldsDSDvec*)(areaCPCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsDSD*)(areaCPC(0,0,iz)))+i)[f];                }
  for(size_t i=0; i<sizeof(areaMMC_type)/sizeof(fieldsSSD); i++) for(size_t f=0; f<sizeof(fieldsSSD)/sizeof(MainType); f++) {
  ((MainType*)&((((fieldsSSDvec*)(areaMMCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsSSD*)(areaMMC(0,0,iz)))+i)[f];
  ((MainType*)&((((fieldsSSDvec*)(areaPMCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsSSD*)(areaPMC(0,0,iz)))+i)[f];
  ((MainType*)&((((fieldsSSDvec*)(areaMPCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsSSD*)(areaMPC(0,0,iz)))+i)[f];
  ((MainType*)&((((fieldsSSDvec*)(areaPPCvec(0,0,iz)))+i)[f]))[iz/VecL] = (((fieldsSSD*)(areaPPC(0,0,iz)))+i)[f];                }
#endif
}

struct CFinfo{
  int n; // count
  int lx,ly,lz;  //coordinates of left ground of ConeFold
  int core;
  int nodes[8];
  double counttime;
  int thtime;
  CFinfo(): n(0),thtime(-1) {}
  void setinfo(int in, int layer) {
    n=0;
    int iz =  in/((pars.NbrX+1)*(pars.NbrY+1));
    int iy = (in%((pars.NbrX+1)*(pars.NbrY+1)))/(pars.NbrX+1);
    int ix =  in%( pars.NbrY+1 );
    int cubeNn = ix+iy+iz;
    int br_time_step = layer-cubeNn;
    lx=pars.NbrX-1-ix; ly=pars.NbrY-1-iy; lz=pars.NbrZ-1-iz;
    if(cubeNn%(dim+1)==layer%(dim+1) && cubeNn<=layer && br_time_step<pars.BigSteps) {
     n++;
     int rx=lx+1, ry=ly+1, rz=lz+1;
          if(rx==pars.NbrX && ry==pars.NbrY && rz==pars.NbrZ         ) { thtime=3; } 
     else if(rx==pars.NbrX && ry==pars.NbrY && rz< pars.NbrZ && lz>=0) { thtime=4; }
     else if(lx==-1        && ry==pars.NbrY && rz==pars.NbrZ         ) { thtime=1; }
     else if(lx==-1        && ry==pars.NbrY && rz< pars.NbrZ && lz>=0) { thtime=2; }
     else if(rx==pars.NbrX && ly==-1        && rz==pars.NbrZ         ) { thtime=1; }
     else if(rx==pars.NbrX && ly==-1        && rz< pars.NbrZ && lz>=0) { thtime=2; }
     else if(rx==pars.NbrX && ry==pars.NbrY && lz==-1                ) { thtime=1; }
     else if(lx==-1        && ly==-1        && rz==pars.NbrZ         ) { thtime=1; }
     else if(lx==-1        && ly==-1        && rz< pars.NbrZ && lz>=0) { thtime=4; }
     else if(lx==-1        && ry==pars.NbrY && lz==-1                ) { thtime=1; }
     else if(rx==pars.NbrX && ly==-1        && lz==-1                ) { thtime=1; }
     else if(lx==-1        && ly==-1        && lz==-1                ) { thtime=3; }
    }
  }
  void act(void* pFldData, void* pFldSide[6], void* pFldEdge[12], void* pFldCorn[8]){
    int rx=lx+1, ry=ly+1, rz=lz+1;
    core = sched_getcpu();
         if(rx==pars.NbrX && ry==pars.NbrY && rz==pars.NbrZ         ) { double stime=omp_get_wtime(); defRRRact; counttime=(omp_get_wtime()-stime); } 
    else if(rx==pars.NbrX && ry==pars.NbrY && rz< pars.NbrZ && lz>=0) { double stime=omp_get_wtime(); defRRDact; counttime=(omp_get_wtime()-stime); }
    else if(lx==-1        && ry==pars.NbrY && rz==pars.NbrZ         ) { double stime=omp_get_wtime(); defLRRact; counttime=(omp_get_wtime()-stime); }
    else if(lx==-1        && ry==pars.NbrY && rz< pars.NbrZ && lz>=0) { double stime=omp_get_wtime(); defLRDact; counttime=(omp_get_wtime()-stime); }
    else if(rx==pars.NbrX && ly==-1        && rz==pars.NbrZ         ) { double stime=omp_get_wtime(); defRLRact; counttime=(omp_get_wtime()-stime); }
    else if(rx==pars.NbrX && ly==-1        && rz< pars.NbrZ && lz>=0) { double stime=omp_get_wtime(); defRLDact; counttime=(omp_get_wtime()-stime); }
    else if(rx==pars.NbrX && ry==pars.NbrY && lz==-1                ) { double stime=omp_get_wtime(); defRRLact; counttime=(omp_get_wtime()-stime); }
    else if(lx==-1        && ly==-1        && rz==pars.NbrZ         ) { double stime=omp_get_wtime(); defLLRact; counttime=(omp_get_wtime()-stime); }
    else if(lx==-1        && ly==-1        && rz< pars.NbrZ && lz>=0) { double stime=omp_get_wtime(); defLLDact; counttime=(omp_get_wtime()-stime); }
    else if(lx==-1        && ry==pars.NbrY && lz==-1                ) { double stime=omp_get_wtime(); defLRLact; counttime=(omp_get_wtime()-stime); }
    else if(rx==pars.NbrX && ly==-1        && lz==-1                ) { double stime=omp_get_wtime(); defRLLact; counttime=(omp_get_wtime()-stime); }
    else if(lx==-1        && ly==-1        && lz==-1                ) { double stime=omp_get_wtime(); defLLLact; counttime=(omp_get_wtime()-stime); }
//---------------------
    else if(rx< pars.NbrX && lx>=0 && ry==pars.NbrY          && rz==pars.NbrZ         ) { printf("DRRact\n"); defDRRact_3d; }
    else if(rx==pars.NbrX          && ry< pars.NbrY && ly>=0 && rz==pars.NbrZ         ) { printf("RDRact\n"); defRDRact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ry< pars.NbrY && ly>=0 && rz==pars.NbrZ         ) { printf("DDRact\n"); defDDRact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ry==pars.NbrY          && rz< pars.NbrZ && lz>=0) { printf("DRDact\n"); defDRDact_3d; }
    else if(rx==pars.NbrX          && ry< pars.NbrY && ly>=0 && rz< pars.NbrZ && lz>=0) { printf("RDDact\n"); defRDDact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ry< pars.NbrY && ly>=0 && lz==-1                ) { printf("DDLact\n"); defDDLact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ly==-1                 && rz< pars.NbrZ && lz>=0) { printf("DLDact\n"); defDLDact_3d; }
    else if(lx==-1                 && ry< pars.NbrY && ly>=0 && rz< pars.NbrZ && lz>=0) { printf("LDDact\n"); defLDDact_3d; }
    else if(lx==-1                 && ry< pars.NbrY && ly>=0 && rz==pars.NbrZ         ) { printf("LDRact\n"); defLDRact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ly==-1                 && rz==pars.NbrZ         ) { printf("DLRact\n"); defDLRact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ry==pars.NbrY          && lz==-1                ) { printf("DRLact\n"); defDRLact_3d; }
    else if(rx==pars.NbrX          && ry< pars.NbrY && ly>=0 && lz==-1                ) { printf("RDLact\n"); defRDLact_3d; }
    else if(lx==-1                 && ry< pars.NbrY && ly>=0 && lz==-1                ) { printf("LDLact\n"); defLDLact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ly==-1                 && lz==-1                ) { printf("DLLact\n"); defDLLact_3d; }
    else if(rx< pars.NbrX && lx>=0 && ry< pars.NbrY && ly>=0 && rz< pars.NbrZ && lz>=0) { printf("DDDact\n"); defDDDact_3d; }
  }
  static int thtimesort(const void* a, const void* b){ 
	  return (((CFinfo*)b)->thtime - ((CFinfo*)a)->thtime); 
  }
  void printinfo() {
    if(n>0) printf("|%5.4g(%#2d):%d",counttime,core,thtime); 
    else    printf("|           ");  
  }
  void checkptr(void* p) {
//     int status[1]; status[0]=-1; int ret_code;
//     ret_code=move_pages(0, 1, &p,  NULL, status, 0); 
//     nodes[0] = status[0];
//     if(status[0]==-1) printf("\n");
  }
};*/

/*void BBB(void* pFldData, void* pFldSide[6], void* pFldEdge[12], void* pFldCorn[8]){
  double starttime,endtime,tt,t1=omp_get_wtime(); 
	printf("BBB\n");
  for(int layer=0; layer<pars.NbrX+1+pars.NbrY+1+pars.NbrZ+1+pars.BigSteps+1; layer++) {
    if ((layer+1)%4==0) t1=omp_get_wtime(); 
    starttime = omp_get_wtime();
    int CFnum=(pars.NbrX+1)*(pars.NbrY+1)*(pars.NbrZ+1);
    CFinfo cf[CFnum]; // LL_,RL_,LR_,RR_
    for(int iBr=0; iBr<CFnum; iBr++) cf[iBr].setinfo(iBr, layer);
    qsort(cf,CFnum,sizeof(CFinfo),CFinfo::thtimesort);
    #pragma omp parallel for shared(cf) schedule(dynamic,1)
    for(int iBr=0; iBr<CFnum; iBr++) if(cf[iBr].n>0) cf[iBr].act(pFldData, pFldSide, pFldEdge, pFldCorn);
    endtime = omp_get_wtime();
    tt=1./pars.CPUfreq;
    int totCF=0; for(int iBr=0; iBr<CFnum; iBr++) totCF+=cf[iBr].n;
    if (endtime-starttime>1.e-3) {
      printf("#Layer %i, Time %#06gs, Async ConeFolds=%d,           |Time(core,numa_nodes):theor.ntu|\n",layer,endtime-starttime, totCF);
      printf("___|     L     "); for(int z=1;z<pars.NbrZ;z++) printf("|     D     "); printf("|     R     \n");
      for(int xy=3;xy>=0;xy--) {
        if(xy==0) printf("LL_"); if(xy==1) printf("RL_"); if(xy==2) printf("LR_"); if(xy==3) printf("RR_"); 
        for(int z=0;z<pars.NbrZ+1;z++) for(int i=0;i<CFnum;i++) if(cf[i].lx==xy%2-1 && cf[i].ly==xy/2-1 && cf[i].lz==z-1) cf[i].printinfo(); 
        printf("\n");
      }
    }
    printf("Fulltime=%g\n",omp_get_wtime()-t1); 
    fflush(stdout);
  }
}*/

//void* pFldData;
//void* pFldSide[6]; void* pFldEdge[12]; void* pFldCorn[8];
//void *pbuf;

int Run() {
  printf("\nStarting...\n");
  for (int iT = 0; iT<=70; iT++){
    //for(int i=0;i<12;i++) printf("CoffArr[0][%d]=%g\n",i,((float*)(&CoffArr[0]))[i]);
    //for(int i=0;i<12;i++) printf("CoffArr[1][%d]=%g\n",i,((float*)(&CoffArr[1]))[i]);
    //for(int i=0;i<38;i++) printf("MatArr[0][%d]=%g\n",i,((float*)(&MatArr[0]))[i]);
    //for(int i=0;i<38;i++) printf("MatArr[1][%d]=%g\n",i,((float*)(&MatArr[1]))[i]);
    if(iT>0) { UpdateOneStep(iT);  };
    if((iT%1)==0) { DropAllData(iT); }
  }
  DeleteAllData();
};

int UpdateOneStep(int it=0){
  double starttime, endtime, tt;
  nLAnetArrayNUMA<3>::thrs4node = pars.Nthreads;
  nLAnetArrayNUMA<3>::MaxStep = it;
  starttime = omp_get_wtime();
  numa_run_on_nodes(nLAnetArrayNUMA_runChessFold, nLAarr, Bricks.numNUMAnodes);
/*    for(int i=nLAarr[0].Nnet-1;i>=0;i--) {
     nLAarr[0].netArr[i].check4grow();
     ConeFold<0>(nLAarr[0].netArr+i,0);
     nLAarr[0].netArr[i].finish_grow();
    }*/
  if(it>0) {
    endtime = omp_get_wtime();
    size_t PPsize=0;
    for(int tt=0; tt<300; tt++) for(int cc=0;cc<6;cc++)
    PPsize+= (PPoutXm[tt][cc].N[0]*PPoutXm[tt][cc].N[1]*PPoutXm[tt][cc].N[2]+
              PPoutXp[tt][cc].N[0]*PPoutXp[tt][cc].N[1]*PPoutXp[tt][cc].N[2]+
              PPoutYm[tt][cc].N[0]*PPoutYm[tt][cc].N[1]*PPoutYm[tt][cc].N[2]+
              PPoutYp[tt][cc].N[0]*PPoutYp[tt][cc].N[1]*PPoutYp[tt][cc].N[2]+
              PPoutZm[tt][cc].N[0]*PPoutZm[tt][cc].N[1]*PPoutZm[tt][cc].N[2]+
              PPoutZp[tt][cc].N[0]*PPoutZp[tt][cc].N[1]*PPoutZp[tt][cc].N[2] )*sizeof(float);
    if(it%10==0 && false) dropPP();
    tt=1./pars.CPUfreq;
    int Nx=Bricks.mapNUMAnodes[0][0].size();
    int Ny=Bricks.mapNUMAnodes[0].size();
    int Nz=Bricks.mapNUMAnodes.size();
    Nx--;Ny--;Nz--;
//    int Nbricks=0; for(int ix=0;ix<=Nx;ix++) for(int iy=0;iy<=Ny;iy++) for(int iz=0;iz<=Nz;iz++) if(ix+iy+iz<it) Nbricks++;
    int Nbricks=Nx*Ny*Nz*(it*4>Nx+Ny+Nz);
    printf("#Step %i, Time %gs, Cells/s=%g [ PPsize=%g GB\n",it,endtime-starttime, double(((long long unsigned int)(1))<<(MaxRank*4))*Nbricks/(endtime-starttime), PPsize/1024./1024./1024.);
  }
  fflush(stdout);
//  for(int in=0; in<numNUMAnodes; in++) {
//    nLAarr[in].clear();
//  }
//  delete[] nLAarr;
}
int DropAllData(int it=0) { 
  return 1;
} 
int DeleteAllData(){
  return 1;
}

//#include "pc3d.hpp"
#include "materials.hpp"

int Main() {
  parsMxw* pars = getPars();
  parsSig* parsSignal = getParsSrc();
  parsMaterial* parsMat = getParsMat();
  pars->setGrid(0.02, 0.02, 0.02, 0.01);
  pars->Tcount=100;
  parsSignal->setPars();
  parsMat->setPars();
  pars->CPUfreq=2.3e9;
  pars->Cores=-1;
  
  InitializeCoffs();

//  setCoffsIso(IndVac, 1. ,1.,0.,0.,0.,0.);
//  setCoffsIso(IndMat, 12.,1.,0.,0.,0.,0.);

  InitializeData(); // only after coofs were inited!!!
//  Run(); 
}
int main() { Main(); return 0; }
