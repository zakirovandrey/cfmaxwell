#ifndef COFF_HPP
#define COFF_HPP

#include <vector>
#include <map>
#include "materials.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <omp.h>

using namespace std;

struct CoffStruct {
  MainType kD[3][2]; MainType kB[3][2];// MainTypeV kPD[3][2]; MainTypeV kPB[3][2];
  MainType add[4];   // for sizeof CoffStruct to be good
  public:
   void set() { //Вот это и есть осн-й интерфейс иниц-ции стр-ры пар-ров
    MainType dtdx[]={pars.dt/pars.dx, pars.dt/pars.dy, pars.dt/pars.dz};
    for(int i=0; i<3; i++) {
      int ip=(i<2)?i+1:0, im=(i>0)?i-1:2;
       kB[i][0] = -dtdx[ip];        kB[i][1] =  dtdx[im];
       kD[i][0] =  dtdx[ip];        kD[i][1] = -dtdx[im];
    }
  }
};

struct MatStruct;
struct DispStruct;

extern vector<CoffStruct> CoffArr;
extern vector<DispStruct> DispArr;
extern vector<MatStruct> MatArr;

extern double cctime1;
extern double cctime2;
extern double cctime3;

#define MACRO_str(D) #D
#ifdef MAINTYPEisFLOAT
#define DIGEQUIVALENT int
#define DIGEQUIVALENT_str "int"
#elif MAINTYPEisDOUBLE
#define DIGEQUIVALENT long int
#define DIGEQUIVALENT_str "long int"
#endif

struct keyVolStr{ 
  MainType Plane[3][3]; 
  void recount(MainType pl[3][3]) {for(int i=0;i<3;i++) for(int j=0;j<3;j++) Plane[i][j]=pl[i][j]; }
  keyVolStr(MainType pl[3][3]) {for(int i=0;i<3;i++) for(int j=0;j<3;j++) Plane[i][j]=pl[i][j]; }
  keyVolStr() {}
  bool operator <(const keyVolStr &p) const {
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
      if      ( fabs(Plane[i][j]-p.Plane[i][j]) >= pars.mapResolution && Plane[i][j]<p.Plane[i][j] ) return true;
      else if ( fabs(Plane[i][j]-p.Plane[i][j]) >= pars.mapResolution && Plane[i][j]>p.Plane[i][j] ) return false;
    }
    return false;
  } 
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { 
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) ar & (*((DIGEQUIVALENT*)(&Plane[i][j]))); 
  }
};
struct Vols{ 
  MainType Volume[2]; 
  MainType Normal[3]; 
  Vols(MainType v[2],MainType n[3]) { for(int i=0;i<2;i++) Volume[i]=v[i]; for(int i=0;i<3;i++) Normal[i]=n[i]; }
  Vols() {for(int i=0;i<2;i++) Volume[i]=0.; for(int i=0;i<3;i++) Normal[i]=0.; }
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { 
    for(int i=0;i<2;i++) ar & (*((DIGEQUIVALENT*)(&Volume[i]))); 
    for(int i=0;i<3;i++) ar & (*((DIGEQUIVALENT*)(&Normal[i]))); 
  }
};
extern map<keyVolStr, Vols> VolumeMap;
extern map<MatStruct, int> MaterialMap;
extern map<DispStruct, int> DispMaterialMap;

struct ProjMatrix{
  MainType XX, YY, ZZ; 
  MainType XY1, XY2, XZ1, XZ2; 
  MainType YZ1, YZ2, YX1, YX2; 
  MainType ZX1, ZX2, ZY1, ZY2; 
  bool operator <(const ProjMatrix &p) const {
    MainType members_l[30]={   XX,   YY,   ZZ,   XY1,   XY2,   XZ1,   XZ2,   YZ1,   YZ2,   YX1,   YX2,   ZX1,   ZX2,   ZY1,   ZY2 };
    MainType members_r[30]={ p.XX, p.YY, p.ZZ, p.XY1, p.XY2, p.XZ1, p.XZ2, p.YZ1, p.YZ2, p.YX1, p.YX2, p.ZX1, p.ZX2, p.ZY1, p.ZY2 };
    for(int i=0;i<15;i++) {
      if      ( fabs(members_l[i]-members_r[i]) >= pars.mapResolution && members_l[i]<members_r[i] ) return true;
      else if ( fabs(members_l[i]-members_r[i]) >= pars.mapResolution && members_l[i]>members_r[i] ) return false;
    }
    return false;
  }
};

struct MatStruct {
  MainType eps, mu, deps, dmu;
  MainType deps_n, deps_t, dmu_n, dmu_t;
  MainType depsXX, depsYY, depsZZ; 
  MainType depsXY1, depsXY2, depsXZ1, depsXZ2; 
  MainType depsYZ1, depsYZ2, depsYX1, depsYX2; 
  MainType depsZX1, depsZX2, depsZY1, depsZY2; 
  MainType dmuXX, dmuYY, dmuZZ; 
  MainType dmuXY1, dmuXY2, dmuXZ1, dmuXZ2; 
  MainType dmuYZ1, dmuYZ2, dmuYX1, dmuYX2; 
  MainType dmuZX1, dmuZX2, dmuZY1, dmuZY2;
  ProjMatrix Proj;
  MainType f1[3],f2[3];
  bool isTensor;
  bool operator <(const MatStruct &p) const { 
    MainType members_l[30]={
//      eps, mu, deps, dmu, deps_n, deps_t, dmu_n, dmu_t,
      depsXX, depsYY, depsZZ, depsXY1, depsXY2, depsXZ1, depsXZ2, depsYZ1, depsYZ2, depsYX1, depsYX2, depsZX1, depsZX2, depsZY1, depsZY2,
      dmuXX, dmuYY, dmuZZ, depsXY1, dmuXY2, dmuXZ1, dmuXZ2, dmuYZ1, dmuYZ2, dmuYX1, dmuYX2, dmuZX1, dmuZX2, dmuZY1, dmuZY2 };
    MainType members_r[30]={
//      p.eps, p.mu, p.deps, p.dmu, p.deps_n, p.deps_t, p.dmu_n, p.dmu_t,
      p.depsXX, p.depsYY, p.depsZZ, p.depsXY1, p.depsXY2, p.depsXZ1, p.depsXZ2, p.depsYZ1, p.depsYZ2, p.depsYX1, p.depsYX2, p.depsZX1, p.depsZX2, p.depsZY1, p.depsZY2,
      p.dmuXX, p.dmuYY, p.dmuZZ, p.depsXY1, p.dmuXY2, p.dmuXZ1, p.dmuXZ2, p.dmuYZ1, p.dmuYZ2, p.dmuYX1, p.dmuYX2, p.dmuZX1, p.dmuZX2, p.dmuZY1, p.dmuZY2 };
    for(int i=0;i<30;i++) {
      if      ( fabs(members_l[i]-members_r[i]) >= pars.mapResolution && members_l[i]<members_r[i] ) return true;
      else if ( fabs(members_l[i]-members_r[i]) >= pars.mapResolution && members_l[i]>members_r[i] ) return false;
    }
    if      ( Proj<p.Proj ) return true;
    else if ( p.Proj<Proj ) return false;
    for(int i=0;i<3;i++) {
      if      ( fabs(f1[i]-p.f1[i]) >= pars.mapResolution && f1[i]<p.f1[i] ) return true;
      else if ( fabs(f1[i]-p.f1[i]) >= pars.mapResolution && f1[i]>p.f1[i] ) return false;
      if      ( fabs(f2[i]-p.f2[i]) >= pars.mapResolution && f2[i]<p.f2[i] ) return true;
      else if ( fabs(f2[i]-p.f2[i]) >= pars.mapResolution && f2[i]>p.f2[i] ) return false;
    }
    return false;
  }
  void setTensorInv(double de[9], double dm[9]){
    isTensor=1;
    depsXX            =     de[0];  dmuXX           =     dm[0];
    depsXY1 = depsXY2 = 0.5*de[1];  dmuXY1 = dmuXY2 = 0.5*dm[1];
    depsXZ1 = depsXZ2 = 0.5*de[2];  dmuXZ1 = dmuXZ2 = 0.5*dm[2];
    depsYX1 = depsYX2 = 0.5*de[3];  dmuYX1 = dmuYX2 = 0.5*dm[3];
    depsYY            =     de[4];  dmuYY           =     dm[4];
    depsYZ1 = depsYZ2 = 0.5*de[5];  dmuYZ1 = dmuYZ2 = 0.5*dm[5];
    depsZX1 = depsZX2 = 0.5*de[6];  dmuZX1 = dmuZX2 = 0.5*dm[6];
    depsZY1 = depsZY2 = 0.5*de[7];  dmuZY1 = dmuZY2 = 0.5*dm[7];
    depsZZ            =     de[8];  dmuZZ           =     dm[8];
    eps = 0; mu=0; deps=0; dmu=0;
  }
  void set(float e[3], float m[3]) { 
    eps = e[0]   ; mu = m[0];
    deps = 1./eps; dmu = 1./mu;
  }
  void set(float e, float m) { 
    eps = e   ; mu = m;
    deps = 1./eps; dmu = 1./mu;
    isTensor=0;
  }
  inline void getBetween(float point1[3], float point2[3], float pointC[3], float mesh_c[3], int deep=0) {
     MainType point_real1[3], point_real2[3], point_realC[3];
     for (int i=0;i<3;i++) pointC[i] = 0.5*(point1[i]+point2[i]);
     for(int icr=0;icr<3;icr++) point_real1[icr]=point1[icr]+mesh_c[icr]; 
     for(int icr=0;icr<3;icr++) point_real2[icr]=point2[icr]+mesh_c[icr]; 
     for(int icr=0;icr<3;icr++) point_realC[icr]=pointC[icr]+mesh_c[icr]; 
     if (deep>pars.deep1d) return;
     else deep++;
     if      (getMaterial(point_real1)!=getMaterial(point_realC)) getBetween(point1, pointC, pointC, mesh_c, deep);
     else if (getMaterial(point_real2)!=getMaterial(point_realC)) getBetween(point2, pointC, pointC, mesh_c, deep);
   }
   inline bool getNormal(float Plane[3][3], float N[3]) {
     float v1[3], v2[3]; for(int i=0;i<3;i++) { v1[i] = Plane[1][i]-Plane[0][i]; v2[i] = Plane[2][i]-Plane[1][i]; }
     for (int i=0;i<3;i++) N[i] = v1[(i+1)%3]*v2[(i+2)%3]-v1[(i+2)%3]*v2[(i+1)%3]; 
     if (N[0]==0. && N[1]==0. && N[2]==0.) return false; 
     for (int i=0;i<3;i++) N[i]/= sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
     return true;
   }
   inline void countVolume(float corns[8][3], float mesh_c[3], float V[2], int Mat_base, float averRadius[3], int deep=0) {
     float cords_real[3];
     if (deep>pars.deep3d) return;
     for(int icr=0;icr<3;icr++) cords_real[icr]=corns[0][icr]+mesh_c[icr]; 
     int Mat0=getMaterial(cords_real); int Mat1=Mat0;
     for (int i=1;i<8;i++) { for(int icr=0;icr<3;icr++) cords_real[icr]=corns[i][icr]+mesh_c[icr]; int M=getMaterial(cords_real); if(M!=Mat0) { Mat1=M; break; } }
     float sizeX = 2*averRadius[0]/(1<<deep); 
     float sizeY = 2*averRadius[1]/(1<<deep); 
     float sizeZ = 2*averRadius[2]/(1<<deep); 
     if (Mat1==Mat0) { if (Mat1==Mat_base) V[0]+=sizeX*sizeY*sizeZ; else V[1]+=sizeX*sizeY*sizeZ; }
     else {
       for(int is=0;is<8;is++) {
         float corns_split[8][3];
         for(int i=0;i<8;i++) { 
           corns_split[i][0] = corns[0][0] + sizeX*0.5*(  (i    &1) + ( is    &1) ); 
           corns_split[i][1] = corns[0][1] + sizeY*0.5*( ((i>>1)&1) + ((is>>1)&1) ); 
           corns_split[i][2] = corns[0][2] + sizeZ*0.5*( ((i>>2)&1) + ((is>>2)&1) ); 
         }
         countVolume(corns_split, mesh_c, V, Mat_base, averRadius, deep+1);
       }
     }
   }
   void printLog(float center[3], float corns[8][3], float corns_path[16][3], int Mat0, int Mat1, float Plane[3][3], float Normal[3], int n) {
     printf("averaging in (%g,%g,%g)...\n",center[0]/pars.dx,center[1]/pars.dy,center[2]/pars.dz);
     printf("Corners:\n")     ; for (int i=0;i<8 ;i++) printf("(%g,%g,%g)\n",corns[i][0]/pars.dx,corns[i][1]/pars.dy,corns[i][2]/pars.dz);
     printf("Corners Path:\n"); for (int i=0;i<16;i++) printf("(%g,%g,%g)\n",corns_path[i][0]/pars.dx,corns_path[i][1]/pars.dy,corns_path[i][2]/pars.dz);
     printf("Material0: %d\n",Mat0);
     printf("Material1: %d\n",Mat1);
     printf("Materials in corners:\n"); for (int i=0;i<16;i++) printf("%d\n",getMaterial(corns_path[i]));
     printf("n=%d\n",n);
     printf("Plane[0]: (%g,%g,%g)\n",Plane[0][0]/pars.dx,Plane[0][1]/pars.dy,Plane[0][2]/pars.dz);
     printf("Plane[1]: (%g,%g,%g)\n",Plane[1][0]/pars.dx,Plane[1][1]/pars.dy,Plane[1][2]/pars.dz);
     printf("Plane[2]: (%g,%g,%g)\n",Plane[2][0]/pars.dx,Plane[2][1]/pars.dy,Plane[2][2]/pars.dz);
   }

   void averaging(float center[3], float mesh_c[3], float averaging_radius[3], float Normal[3], float& calcVolume0, float& calcVolume1, int& calcMat) {
     float cords_real[3];
     float ar[3]; for(int i=0;i<3;i++) ar[i]=averaging_radius[i]; 
     // counting the normal
     float corns[8][3];
     for (int i=0;i<8;i++) { corns[i][0]=center[0]-ar[0]+2*ar[0]*(i&1); corns[i][1]=center[1]-ar[1]+2*ar[1]*((i>>1)&1); corns[i][2]=center[2]-ar[2]+2*ar[2]*((i>>2)&1); }
     float corns_path[16][3];  // углы по которым обходится куб в правильном порядке
     const int p[16]={0,1,3,7,5,4,6,2,0,4,1,5,3,2,7,6};
     for (int i=0;i<3;i++) for (int ip=0;ip<16;ip++) corns_path[ip][i]=corns[p[ip]][i];
     float Plane[3][3];
     for (int i=0;i<3;i++) for (int j=0 ;j<3  ;j++ ) Plane[i][j]=0.;  // точки по которым строится плоскость границы материалов
     int n=0;
     MainType Volume[2] = {0.,0.};        // finding parts of Volumes
     double tt3=omp_get_wtime(); 
     keyVolStr key(Plane); int isCounted=VolumeMap.count(key);   // check if already counted
     cctime3+=omp_get_wtime()-tt3;
     double tt1=omp_get_wtime(); 
     for(int icr=0;icr<3;icr++) cords_real[icr]=corns_path[0][icr]+mesh_c[icr]; 
     int Mat0 = getMaterial(cords_real); int Mat1=Mat0;
     for (int i=0;i<15;i++) { 
       for(int icr=0;icr<3;icr++) cords_real[icr]=corns_path[i  ][icr]+mesh_c[icr]; int M0 = getMaterial(cords_real); 
       for(int icr=0;icr<3;icr++) cords_real[icr]=corns_path[i+1][icr]+mesh_c[icr]; int M1 = getMaterial(cords_real);
       if (M0!=M1) { if (Mat1==Mat0) Mat1=M1; getBetween(corns_path[i],corns_path[i+1],Plane[n],mesh_c); n++; }
       if (n==3) {
         key.recount(Plane); 
         tt3=omp_get_wtime(); isCounted = VolumeMap.count(key); cctime3+=omp_get_wtime()-tt3;
         if(isCounted>0) break;
         else { if(getNormal(Plane, Normal)) break; else n=2; }
       }
     }
     if (n<3 && Mat0!=Mat1) {
       printf ("Cannot find normal! (center cell in (%g,%g,%g))\n", center[0]/pars.dx,center[1]/pars.dy,center[2]/pars.dz);
       printLog(center,corns,corns_path,Mat0,Mat1, Plane,Normal,n);
     }
     cctime1+=omp_get_wtime()-tt1;
     if (isCounted>0) {
       tt3=omp_get_wtime(); 
       Vols vol(VolumeMap[key]);
       Volume[0] = vol.Volume[0]; Volume[1] = vol.Volume[1];
       for(int i=0;i<3;i++) Normal[i] = vol.Normal[i];
       cctime3+=omp_get_wtime()-tt3;
     }
     else {
       double tt2=omp_get_wtime();
       if (Mat0!=Mat1) countVolume(corns, mesh_c, Volume, Mat0, ar);
       else { Volume[0] = ar[0]*ar[1]*ar[2]; Volume[1] = 0.; }
       Vols vol(Volume,Normal); tt3=omp_get_wtime(); VolumeMap[key] = vol; cctime3+=omp_get_wtime()-tt3;
       cctime2+=omp_get_wtime()-tt2;
     }
     if (Mat0!=Mat1) {
       // counting aver deps and dmu parallel to the normal (i.e. <eps^-1> )
       deps_n = (MatArr[Mat0].deps*Volume[0]+MatArr[Mat1].deps*Volume[1])/(Volume[0]+Volume[1]);
       dmu_n  = (MatArr[Mat0].dmu *Volume[0]+MatArr[Mat1].dmu *Volume[1])/(Volume[0]+Volume[1]);
       // counting aver deps and dmu ortogonal to the normal (i.e. <eps>^-1 )
       deps_t = (Volume[0]+Volume[1])/(MatArr[Mat0].eps*Volume[0]+MatArr[Mat1].eps*Volume[1]);
       dmu_t  = (Volume[0]+Volume[1])/(MatArr[Mat0].mu *Volume[0]+MatArr[Mat1].mu *Volume[1]);
     }
     else { 
       Volume[0] = ar[0]*ar[1]*ar[2]; Volume[1] = 0.;
       deps_n = MatArr[Mat0].deps; dmu_n  = MatArr[Mat0].dmu; 
       deps_t = deps_n; dmu_t=dmu_n; 
       for(int i=0;i<3;i++) Normal[i] = 0.;
     }
     calcVolume0 = Volume[0]; calcVolume1 = Volume[1]; calcMat = Mat0;
   }
   void count_coffs(float mesh_c[3],float aver_radius[3]) {
     float c[3]; // center_new;
     float Norm[3];    // Нормаль
     int mat[3]; float nv0,nv1; int nv2;
     // rounding to get tensor deps and dmu in true coordinat system
     c[0]=0.            ; c[1]=0.-0.5*pars.dy; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, f1[0], f2[0], mat[0]); Proj.XX  = Norm[0]*Norm[0]; depsXX  = Proj.XX*(deps_n-deps_t) + deps_t;
     c[0]=0.-0.5*pars.dx; c[1]=0.            ; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, f1[1], f2[1], mat[1]); Proj.YY  = Norm[1]*Norm[1]; depsYY  = Proj.YY*(deps_n-deps_t) + deps_t;
     c[0]=0.-0.5*pars.dx; c[1]=0.-0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, f1[2], f2[2], mat[2]); Proj.ZZ  = Norm[2]*Norm[2]; depsZZ  = Proj.ZZ*(deps_n-deps_t) + deps_t;         
     c[0]=0.            ; c[1]=0.-    pars.dy; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.XY1 = Norm[0]*Norm[1]; depsXY1 = 0.5*Proj.XY1*(deps_n-deps_t);
     c[0]=0.            ; c[1]=0.-0.5*pars.dy; c[2]=0.-    pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.XZ1 = Norm[0]*Norm[2]; depsXZ1 = 0.5*Proj.XZ1*(deps_n-deps_t);
     c[0]=0.-0.5*pars.dx; c[1]=0.            ; c[2]=0.-    pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.YZ1 = Norm[1]*Norm[2]; depsYZ1 = 0.5*Proj.YZ1*(deps_n-deps_t);
     c[0]=0.-    pars.dx; c[1]=0.            ; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.YX1 = Norm[1]*Norm[0]; depsYX1 = 0.5*Proj.YX1*(deps_n-deps_t);
     c[0]=0.-    pars.dx; c[1]=0.-0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.ZX1 = Norm[2]*Norm[0]; depsZX1 = 0.5*Proj.ZX1*(deps_n-deps_t);
     c[0]=0.-0.5*pars.dx; c[1]=0.-    pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.ZY1 = Norm[2]*Norm[1]; depsZY1 = 0.5*Proj.ZY1*(deps_n-deps_t);
     c[0]=0.            ; c[1]=0.            ; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.XY2 = Norm[0]*Norm[1]; depsXY2 = 0.5*Proj.XY2*(deps_n-deps_t); Proj.YX2 = Proj.XY2; depsYX2 = depsXY2; 
     c[0]=0.            ; c[1]=0.-0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.XZ2 = Norm[0]*Norm[2]; depsXZ2 = 0.5*Proj.XZ2*(deps_n-deps_t); Proj.ZX2 = Proj.XZ2; depsZX2 = depsXZ2; 
     c[0]=0.-0.5*pars.dx; c[1]=0.            ; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); Proj.YZ2 = Norm[1]*Norm[2]; depsYZ2 = 0.5*Proj.YZ2*(deps_n-deps_t); Proj.ZY2 = Proj.YZ2; depsZY2 = depsYZ2; 

     c[0]=0.+0.5*pars.dx; c[1]=0.            ; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuXX  = Norm[0]*Norm[0]*(dmu_n-dmu_t) + dmu_t;
     c[0]=0.            ; c[1]=0.+0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuYY  = Norm[1]*Norm[1]*(dmu_n-dmu_t) + dmu_t;
     c[0]=0.            ; c[1]=0.            ; c[2]=0.+0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuZZ  = Norm[2]*Norm[2]*(dmu_n-dmu_t) + dmu_t;          
     c[0]=0.+0.5*pars.dx; c[1]=0.-0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuXY1 = 0.5*Norm[0]*Norm[1]*(dmu_n-dmu_t);
     c[0]=0.+0.5*pars.dx; c[1]=0.            ; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuXZ1 = 0.5*Norm[0]*Norm[2]*(dmu_n-dmu_t);
     c[0]=0.            ; c[1]=0.+0.5*pars.dy; c[2]=0.-0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuYZ1 = 0.5*Norm[1]*Norm[2]*(dmu_n-dmu_t);
     c[0]=0.-0.5*pars.dx; c[1]=0.+0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuYX1 = 0.5*Norm[1]*Norm[0]*(dmu_n-dmu_t);
     c[0]=0.-0.5*pars.dx; c[1]=0.            ; c[2]=0.+0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuZX1 = 0.5*Norm[2]*Norm[0]*(dmu_n-dmu_t);
     c[0]=0.            ; c[1]=0.-0.5*pars.dy; c[2]=0.+0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuZY1 = 0.5*Norm[2]*Norm[1]*(dmu_n-dmu_t);
     c[0]=0.+0.5*pars.dx; c[1]=0.+0.5*pars.dy; c[2]=0.            ; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuXY2 = 0.5*Norm[0]*Norm[1]*(dmu_n-dmu_t); dmuYX2 = dmuXY2;
     c[0]=0.+0.5*pars.dx; c[1]=0.            ; c[2]=0.+0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuXZ2 = 0.5*Norm[0]*Norm[2]*(dmu_n-dmu_t); dmuZX2 = dmuXZ2;
     c[0]=0.            ; c[1]=0.+0.5*pars.dy; c[2]=0.+0.5*pars.dz; averaging(c, mesh_c,aver_radius, Norm, nv0,nv1,nv2); dmuYZ2 = 0.5*Norm[1]*Norm[2]*(dmu_n-dmu_t); dmuZY2 = dmuYZ2;

     for(int i:{0,1,2}) {
       MainType crd[3] = {mesh_c[0]-0.5*pars.dx*(i!=0), mesh_c[1]-0.5*pars.dy*(i!=1), mesh_c[2]-0.5*pars.dz*(i!=2)};
       MainType fsum = f1[i]+f2[i]; f1[i]/=fsum; f2[i]/=fsum; 
       if ( mat[i] != getMaterial( crd ) ) swap(f1[i],f2[i]);
//       if ( mat[i] != getMaterial(mesh_c) ) swap(f1[i],f2[i]);
     }
     
     if(mesh_c[0]<1.01*pars.dx) { depsYX1 = 0.; depsZX1 = 0.; }
     if(mesh_c[1]<1.01*pars.dy) { depsXY1 = 0.; depsZY1 = 0.; }
     if(mesh_c[2]<1.01*pars.dz) { depsXZ1 = 0.; depsYZ1 = 0.; }
   }
}; 

const int Md=1;
struct LorenzDisp {
  double a0,a1,b0,b1,b2;
  LorenzDisp() : a0(0.),a1(0.),b0(0.),b1(0.),b2(0.) {};
  void set(double Deps=1., double wp=0., double gp=0., double gp_=0.) {
    a0 = Deps*wp*wp;
    a1 = Deps*gp_;
    b0 = wp*wp;
    b1 = 2*gp;
    b2 = 1;
  }
  void setDrude(double Deps=1., double wp=0., double gp=0., double gp_=0.) {
    a0 = Deps*wp*wp;
    a1 = Deps*gp_;
    b0 = 0.0;
    b1 = 2*gp;
    b2 = 1;
  }
};

struct DispStruct {
#ifndef SWIG
  MainType kEJ[3][Md][2]; 
  MainType kE[3][3];
  MainType kJ[3][Md][5];
#endif
  int isSet;
  double epsinf, sigma;
  vector<LorenzDisp> Ld;

  DispStruct(): isSet(0) {}
  void setNondisp(MatStruct& mat) { set(mat.eps, vector<LorenzDisp>(), 0, 0.); isSet=0; }
  void set(DispStruct& disp1, DispStruct& disp2, float f1, float f2){
    vector<LorenzDisp> Ld_t = disp1.Ld;
    Ld_t.insert( Ld_t.end(), disp2.Ld.begin(), disp2.Ld.end() );
    for(int i=0; i<Ld_t.size(); i++) { float f; if(i<disp1.Ld.size()) f=f1; else f=f2; Ld_t[i].a0*=f; Ld_t[i].a1*=f; }
    set(disp1.epsinf*f1+disp2.epsinf*f2, Ld_t, Ld_t.size(), disp1.sigma*f1+disp2.sigma*f2);
  }
  void set(double einf, vector<LorenzDisp> _Ld, int mD, double sm) {
    isSet = 1;

    epsinf=einf; sigma=sm; 

    Ld = _Ld;

    double alpha[Md];
    double ksi[Md]; 
    double zetp[Md];
    double zetm[Md];
    double zet[Md];
    double zetm_sum=0., zetp_sum=0., zet_sum=0.;
    for(int i=0;i<mD;i++) {
      alpha[i] = (4-2*Ld[i].b0*pars.dt*pars.dt)/(Ld[i].b1*pars.dt+2);
      ksi[i]   = (Ld[i].b1*pars.dt-2)/(Ld[i].b1*pars.dt+2);
      zetp[i]  = (+Ld[i].a0*pars.dt*pars.dt+2*Ld[i].a1*pars.dt)/(Ld[i].b1*pars.dt+2);
      zetm[i]  = (-Ld[i].a0*pars.dt*pars.dt+2*Ld[i].a1*pars.dt)/(Ld[i].b1*pars.dt+2);
      zet[i]   = -(4*Ld[i].a1*pars.dt)/(Ld[i].b1*pars.dt+2);
      
      zetm_sum+=zetm[i]; zetp_sum+=zetp[i]; zet_sum+=zet[i];
    }
    for(int i=mD;i<Md;i++) {
      alpha[i] = 2;
      ksi[i]   = -1;
      zetp[i]  = 0.0;
      zetm[i]  = 0.0;
      zet[i]   = 0;
      zetm_sum+=zetm[i]; zetp_sum+=zetp[i]; zet_sum+=zet[i];
    }
 
    double C1 = ( -(zetm_sum)                          ) / ( 2*epsinf + sigma*pars.dt + (zetp_sum) );
    double C2 = ( 2*epsinf - sigma*pars.dt - (zet_sum) ) / ( 2*epsinf + sigma*pars.dt + (zetp_sum) );
    double C3 = ( 2*pars.dt                            ) / ( 2*epsinf + sigma*pars.dt + (zetp_sum) );
    for(int xc=0; xc<3; xc++) {
      for(int i=0; i<Md; i++) {
        kJ[xc][i][0] = alpha[i]; kJ[xc][i][1] = ksi[i]; kJ[xc][i][2] = zetp[i]/pars.dt; kJ[xc][i][3] = zetm[i]/pars.dt; kJ[xc][i][4] = zet[i]/pars.dt; 
        kEJ[xc][i][0] = C3*0.5*(1+alpha[i]); kEJ[xc][i][1] = C3*0.5*ksi[i];
      }
      kE[xc][0] = C1; kE[xc][1] = C2; kE[xc][2] = C3/pars.dt;
    }
  }
//  void setSilicon() { set(1., 0., 8.93, 1.855, 3.42*2*M_PI, 2.72*2*M_PI, 0.425*2*M_PI, 0.123*2*M_PI, 0.087*2*M_PI, 2.678*2*M_PI); }
  //void setSilicon() { set(1., 0., 8.93, 1.855, 3.42, 2.72, 0.0, 0.0, 0.0, 0.0); }
  bool operator <(const DispStruct &p) const { 
     for (int xc=0; xc<3; xc++) for (int i=0; i<Md; i++) for (int j=0; j<2; j++) {
    if ( fabs(kEJ[xc][i][j]-p.kEJ[xc][i][j]) >= pars.mapResolution && kEJ[xc][i][j]<p.kEJ[xc][i][j] ) return true;
    if ( fabs(kEJ[xc][i][j]-p.kEJ[xc][i][j]) >= pars.mapResolution && kEJ[xc][i][j]>p.kEJ[xc][i][j] ) return false;
  }  for (int xc=0; xc<3; xc++) for (int j=0; j<3; j++) {
    if ( fabs(kE[xc][j]-p.kE[xc][j]) >= pars.mapResolution && kE[xc][j]<p.kE[xc][j] ) return true;
    if ( fabs(kE[xc][j]-p.kE[xc][j]) >= pars.mapResolution && kE[xc][j]>p.kE[xc][j] ) return false;
  }  for (int xc=0; xc<3; xc++) for (int i=0; i<Md; i++) for (int j=0; j<5; j++) {
    if ( fabs(kJ[xc][i][j]-p.kJ[xc][i][j]) >= pars.mapResolution && kJ[xc][i][j]<p.kJ[xc][i][j] ) return true;
    if ( fabs(kJ[xc][i][j]-p.kJ[xc][i][j]) >= pars.mapResolution && kJ[xc][i][j]>p.kJ[xc][i][j] ) return false;
  }
    if ( fabs(epsinf-p.epsinf) >= pars.mapResolution && epsinf<p.epsinf ) return true;
    if ( fabs(epsinf-p.epsinf) >= pars.mapResolution && epsinf>p.epsinf ) return false;
    if ( fabs(sigma-p.sigma) >= pars.mapResolution && sigma<p.sigma ) return true;
    if ( fabs(sigma-p.sigma) >= pars.mapResolution && sigma>p.sigma ) return false;
    return false;
  }
};


struct indx {
  int isBound;// int isDisp; 
  int I;
  int isDrop;
//  int x,y,z;
  int time;
  //int I[2*PMLNzV+NzV];
};
struct SrcIndx {
  int iX, iY;
  long int it;
  bool isSource;
};

struct PMLcoffs {
  MainType k1_E, k2_E, k1_H, k2_H;   // Gammas for primary waves ( V -- gamma, S (S or T) -- gamma* )
  void setX(short int ii) {
    int NN = PMLNx*2-1;
    k2_E = 1.0/(1.0+0.5*pars.dt*gamma_funcX_E(NN-abs(ii-NN),NN));
    k2_H = 1.0/(1.0+0.5*pars.dt*gamma_funcX_H(NN-abs(ii-NN),NN));
    k1_E = 2.0*k2_E-1.0; k1_H = 2.0*k2_H-1.0; 
  }
  void setY(short int ii) {
    int NN = PMLNy*2-1;
    k2_E = 1.0/(1.0+0.5*pars.dt*gamma_funcY_E(NN-abs(ii-NN),NN));
    k2_H = 1.0/(1.0+0.5*pars.dt*gamma_funcY_H(NN-abs(ii-NN),NN));
    k1_E = 2.0*k2_E-1.0; k1_H = 2.0*k2_H-1.0; 
  }
  void setZ(short int ii) {
    int NN = PMLNz*2-1;
    k2_E = 1.0/(1.0+0.5*pars.dt*gamma_funcZ_E(NN-abs(ii-NN),NN));
    k2_H = 1.0/(1.0+0.5*pars.dt*gamma_funcZ_H(NN-abs(ii-NN),NN));
    k1_E = 2.0*k2_E-1.0; k1_H = 2.0*k2_H-1.0; 
  }
};

extern PMLcoffs* KpmlX;
extern PMLcoffs* KpmlY;
extern PMLcoffs* KpmlZ;

#endif//COFF_HPP
