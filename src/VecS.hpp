#ifndef PML_VECTOR_HPP
#define PML_VECTOR_HPP
#include "VecB.hpp"

template <class T> struct vecSrotL;
template <class T> struct vecSrotR;
template <class T, int NvC, int Npml> struct vecS;

#define PML_VECTOR_DEF(type,typeV,nV)\
template <int tNvC, int tNpml> struct vecS<typeV, tNvC, tNpml> {\
  const static int NvV=(tNvC+tNpml*2)/nV;\
  const static int NvCV=tNvC/nV, NvSV=(tNpml*2)/nV;\
  typeV v[NvV];\
  inline typeV& operator [](const int i) { return v[i]; }\
  inline typeV  operator ()(const int i) { return v[i]; }\
  inline type& sRef(const int i) { if(i<2*tNpml) return sSSref(i); else return sCCref(i-2*tNpml); }\
  inline type sVal(const int i) { return sRef(i); }\
  inline type& sCCref(const int i) { return ((type*)&(v[NvSV+(i%(NvCV))]))[i/(NvCV)]; } inline type  sCC(const int i) { return sCCref(i); }\
  inline type& sSSref(const int i) { return ((type*)&(v[i%NvSV]))[i/NvSV]; } inline type  sSS(const int i) { return sSSref(i); }\
  inline typeV L(const int i) { if(i<0) return LS(i); else return LC(i-NvSV); }\
  inline typeV R(const int i) { if(i<NvCV+NvSV) return RS(i); else return RC(i-NvSV); }\
  inline typeV CC(const int i) { return CCref(i); } inline typeV& CCref(const int i) { return v[i+NvSV]; }\
  inline typeV SS(const int i) { return SSref(i); } inline typeV& SSref(const int i) { return v[i]; }\
  inline typeV LS(const int i) { type* pv=(type*)&v[NvSV+i], fN=sCC(tNvC+i); return ROT_R2(fN,pv); }\
  inline typeV RS(const int i) { type* pv=(type*)&v[i-NvSV], f0=sCC(i-NvSV); return ROT_L2(pv,f0); }\
  inline typeV LC(const int i) { type fL=sSS(tNpml+i),* pv=(type*)&v[NvV+i]; return ROT_R(fL,pv); }\
  inline typeV RC(const int i) { type* pv=(type*)&v[NvSV-NvCV+i], fR=sSS(i-NvCV+tNpml); return ROT_L(pv,fR); }\
};

#define PML_ONLY_VECTOR_DEF(type,typeV,nV)\
template <int tNpml> struct vecS<typeV, 0, tNpml> {\
const static int NvV=(tNpml*2)/nV;\
typeV v[NvV];\
inline typeV& operator [](const int i) { return v[i]; }\
inline typeV  operator ()(const int i) { return v[i]; }\
inline type& sRef(const int i) { return sSSref(i); }\
inline type sVal(const int i) { return sRef(i); }\
inline type& sSSref(const int i) { return ((type*)&(v[i%NvV]))[i/NvV]; } inline type  sSS(const int i) { return sSSref(i); }\
inline typeV SS(const int i) { return SSref(i); } inline typeV& SSref(const int i) { return v[i]; }\
};

#define PML_noVECTOR_DEF(type)\
template <int tNvC, int tNpml> struct vecS<type, tNvC, tNpml> {\
  type v[tNvC+2*tNpml];\
  inline type& operator [](const int i) { return v[i]; }\
  inline type  operator ()(const int i) { return v[i]; }\
  inline type& sRef(const int i) { return v[i]; }\
  inline type  sVal(const int i) { return v[i]; }\
  inline type CC(const int i) { return v[2*tNpml+i]; }\
  inline type SS(const int i) { return v[i]; }\
  inline type L(const int i) { if(i<0) return 0.0; else return LC(i-2*tNpml); }\
  inline type R(const int i) { if(i<tNvC+2*tNpml) return RC(i); else return 0.0; }\
  inline type LC(const int i) { return v[i+tNpml]; }\
  inline type RC(const int i) { return v[i+(tNpml-tNvC)]; }\
  inline type LS(const int i) { return 0.0; }\
  inline type RS(const int i) { return 0.0; }\
};

#if defined(BASIC_VECTOR_AVX)//======BASIC_VECTOR_(AVX/SSE/noVEC)========================
#define ROT_L(pv,fR) (doubleV4) { pv[1],pv[2],pv[3],fR }
#define ROT_R(fL,pv) (doubleV4) { fL,pv[0],pv[1],pv[2] }
#define ROT_L2(pv,f0) (doubleV4) { pv[1],f0,pv[3],0.0 }
#define ROT_R2(fN,pv) (doubleV4) { 0.0,pv[0],fN,pv[2] }
PML_VECTOR_DEF(double,doubleV4, 4);
PML_ONLY_VECTOR_DEF(double,doubleV4, 4);
#undef ROT_L
#undef ROT_R
#undef ROT_L2
#undef ROT_R2
#define ROT_L(pv,fR) (floatV8) { pv[1],pv[2],pv[3],pv[4],pv[5],pv[6],pv[7],fR }
#define ROT_R(fL,pv) (floatV8) { fL,pv[0],pv[1],pv[2],pv[3],pv[4],pv[5],pv[6] }
#define ROT_L2(pv,f0) (floatV8) { pv[1],pv[2],pv[3],f0,pv[5],pv[6],pv[7],0.0f }
#define ROT_R2(fN,pv) (floatV8) { 0.0f,pv[0],pv[1],pv[2],fN,pv[4],pv[5],pv[6] }
PML_VECTOR_DEF(float,floatV8, 8);
PML_ONLY_VECTOR_DEF(float,floatV8, 8);
#undef ROT_L
#undef ROT_R
#undef ROT_L2
#undef ROT_R2
#elif defined(BASIC_VECTOR_SSE)//======BASIC_VECTOR_(AVX/SSE/noVEC)========================
#define ROT_L(pv,fR) (doubleV2) { pv[1],fR }
#define ROT_R(fL,pv) (doubleV2) { fL,pv[0] }
#define ROT_L2(pv,f0) (doubleV2) { f0,0.0 }
#define ROT_R2(fN,pv) (doubleV2) { 0.0,fN }
PML_VECTOR_DEF(double,doubleV2, 2);
PML_ONLY_VECTOR_DEF(double,doubleV2, 2);
#undef ROT_L
#undef ROT_R
#undef ROT_L2
#undef ROT_R2
#define ROT_L(pv,fR) (floatV4) { pv[1],pv[2],pv[3],fR }
#define ROT_R(fL,pv) (floatV4) { fL,pv[0],pv[1],pv[2] }
#define ROT_L2(pv,f0) (floatV4) { pv[1],f0,pv[3],0.0f }
#define ROT_R2(fN,pv) (floatV4) { 0.0f,pv[0],fN,pv[2] }
PML_VECTOR_DEF(float,floatV4, 4);
PML_ONLY_VECTOR_DEF(float,floatV4, 4);
#undef ROT_L
#undef ROT_R
#undef ROT_L2
#undef ROT_R2
#elif defined(BASIC_VECTOR_noVEC)//=====BASIC_VECTOR_(AVX/SSE/noVEC)=========================

#endif//====BASIC_VECTOR_(AVX/SSE/noVEC)==========================
PML_noVECTOR_DEF(double);
PML_noVECTOR_DEF(float);
#undef PML_VECTOR_DEF
#undef PML_noVECTOR_DEF
#endif//BASIC_VECTOR_HPP