#ifndef BASIC_VECTOR_HPP
#define BASIC_VECTOR_HPP
#include <math.h>

#define VEC_DEC(type,typeV,nv,nvV) type __attribute__ (( vector_size(nv*sizeof(type)))) __attribute__ ((aligned(nvV)))
template <class T, int Nv> struct vecB;
#define BASIC_VECTOR_DEF(type, typeV, nV)\
template <int tNv> struct vecB<typeV, tNv> {\
  const static int NvV=tNv/nV;\
  typeV v[NvV];\
  inline typeV& operator [](const int i) { return v[i]; }\
  inline typeV  operator ()(const int i) { return v[i]; }\
  inline type& sRef(const int i) { return ((type*)&(v[i%(NvV)]))[i/(NvV)]; }\
  inline type  sVal(const int i) { return sRef(i); }\
};
#define BASIC_noVECTOR_DEF(type)\
template <int Nv> struct vecB<type, Nv> {\
  type v[Nv];\
  inline type& operator [](const int i) { return v[i]; }\
  inline type  operator ()(const int i) { return v[i]; }\
  inline type& sRef(const int i) { return v[i]; }\
  inline type  sVal(const int i) { return v[i]; }\
};
template <class T> struct anyV;

#if defined(__AVX__)//======BASIC_VECTOR_(AVX/SSE/noVEC)========================
#define BASIC_VECTOR_AVX
#warning AVX vectorization
#include <immintrin.h>
#define F2V(v) _mm_set1_pd(v)
typedef VEC_DEC(float,__m256,8,32) floatV8;
typedef VEC_DEC(double,__m256d,4,32) doubleV4;
inline const int getNV(floatV8&) { return 8; }
inline const int getNV(doubleV4&) { return 4; }

template <> struct anyV<doubleV4> { static inline doubleV4 s2v(const double v) { return _mm256_set1_pd(v); } };
template <> struct anyV<floatV8> { static inline floatV8 s2v(const float v) { return _mm256_set1_ps(v); } };

inline doubleV4 s2v(double const* a) { return _mm256_broadcast_sd(a); }
inline floatV8 s2v(float const* a) { return _mm256_broadcast_ss(a); }
inline doubleV4 s2v(double const& a) { return _mm256_set1_pd(a); }
inline floatV8 s2v(float const& a) { return _mm256_set1_ps(a); }
inline double v2s(doubleV4& a, int i) { return ((double*)&a)[i]; }
inline float v2s(floatV8& a, int i) { return ((float*)&a)[i]; }
inline double& v2sRef(doubleV4& a, int i) { return ((double*)&a)[i]; }
inline float& v2sRef(floatV8& a, int i) { return ((float*)&a)[i]; }
const union {int i[8]; __m256 m; } fabs_mask_cheat_float = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
const union {long long int li[4];__m256d m;} fabs_mask_cheat_double = {0x7fffffffffffffff, 0x7fffffffffffffff, 0x7fffffffffffffff, 0x7fffffffffffffff};
namespace floatVconst {
#ifndef SWIG
  const int NV=8;
  const floatV8 zero=s2v(0.0f);
  const floatV8 half=s2v(0.5f);
  const floatV8 one=s2v(1.0f);
  const floatV8 two=s2v(2.0f);
#endif
};
namespace doubleVconst {
  const int NV=4;
  const doubleV4 zero=s2v(0.0);
  const doubleV4 half=s2v(0.5);
  const doubleV4 one=s2v(1.0);
  const doubleV4 two=s2v(2.0);
};
inline floatV8 fabs(floatV8 v) { return _mm256_and_ps(v, fabs_mask_cheat_float.m); }
inline doubleV4 fabs(doubleV4 v) { return _mm256_and_pd(v,fabs_mask_cheat_double.m); }
//{ return (doubleV4){fabs(v2s(v,0)), fabs(v2s(v,1)), fabs(v2s(v,2)), fabs(v2s(v,3))}; }
//floatV8 fabs(floatV8 v) { return (floatV8){fabsf(v2s(v,0)), fabsf(v2s(v,1)), fabsf(v2s(v,2)), fabsf(v2s(v,3)), fabsf(v2s(v,4)), fabsf(v2s(v,5)), fabsf(v2s(v,6)), fabsf(v2s(v,7))}; }
//doubleV4 fabs(doubleV4 v) { return (doubleV4){fabs(v2s(v,0)), fabs(v2s(v,1)), fabs(v2s(v,2)), fabs(v2s(v,3))}; }
inline floatV8 Max(floatV8 v1,floatV8 v2) { return _mm256_max_ps(v1,v2); }
inline floatV8 Min(floatV8 v1,floatV8 v2) { return _mm256_min_ps(v1,v2); }
inline doubleV4 Max(doubleV4 v1,doubleV4 v2) { return _mm256_max_pd(v1,v2); }
inline doubleV4 Min(doubleV4 v1,doubleV4 v2) { return _mm256_min_pd(v1,v2); }
inline floatV8 rcp(floatV8 const a) { return _mm256_rcp_ps(a); }
inline floatV8 rcp_sqrt(floatV8 const a) { return _mm256_rsqrt_ps(a); }
inline floatV8 sqrt(floatV8 const a) { return _mm256_sqrt_ps(a); }
inline doubleV4 sqrt(doubleV4 const a) { return _mm256_sqrt_pd(a); }
inline doubleV4 rcp(doubleV4 const a) { return _mm256_div_pd(doubleVconst::one, a); }
inline doubleV4 rcp_sqrt(doubleV4 const a) { return _mm256_div_pd(doubleVconst::one, _mm256_sqrt_pd(a)); }

typedef doubleV4 doubleV;
typedef floatV8 floatV;
//typedef float floatT;
//typedef doubleV4 floatV;
//typedef double floatT;
//typedef double __attribute__ ((vector_size(4*sizeof(double)))) __attribute__ ((aligned(4*sizeof(type)))) floatV;
//#define Iv2V(N,i) __extension__ (doubleV4) {iFUNC((i)), iFUNC((i+N)), iFUNC((i+2*N)), iFUNC((i+3*N))}
#define Iv2V(N,i) (doubleV4) {iFUNC((i)), iFUNC((i+N)), iFUNC((i+2*N)), iFUNC((i+3*N))}
#define Src2V(N,t,x,i) (doubleV4) {srcFUNC(t,x,(i)), srcFUNC(t,x,(i+N)), srcFUNC(t,x,(i+2*N)), srcFUNC(t,x,(i+3*N))}
//#define Iv2V(N,i) (floatV8) {iFUNC((i)), iFUNC((i+N)), iFUNC((i+2*N)), iFUNC((i+3*N)), iFUNC((i+4*N)), iFUNC((i+5*N)), iFUNC((i+6*N)), iFUNC((i+7*N))}

BASIC_VECTOR_DEF(double,doubleV4, 4);
BASIC_VECTOR_DEF(float,floatV8, 8);
#elif defined(__SSE2__) //======BASIC_VECTOR_(AVX/SSE/noVEC)========================
#define BASIC_VECTOR_SSE
#warning SSE vectorization
#include <xmmintrin.h>
typedef VEC_DEC(float,__m128,4,16) floatV4;
typedef VEC_DEC(double,__m128d,2,16) doubleV2;
inline const int getNV( floatV4&) { return 4; }
inline const int getNV(doubleV2&) { return 2; }

template <> struct anyV<doubleV2> { static inline doubleV2 s2v(const double v) { return _mm_set1_pd(v); } };
template <> struct anyV<floatV4> { static inline floatV4 s2v(const float v) { return _mm_set1_ps(v); } };

inline doubleV2 s2v(double const a) { return _mm_set1_pd(a); }
inline floatV4 s2v(float const a) { return _mm_set1_ps(a); }
inline double v2s(doubleV2& a, int i) { return ((double*)&a)[i]; }
inline float v2s(floatV4& a, int i) { return ((float*)&a)[i]; }
inline double& v2sRef(doubleV2& a, int i) { return ((double*)&a)[i]; }
inline float& v2sRef(floatV4& a, int i) { return ((float*)&a)[i]; }
const union {int i[4];__m128 m;} fabs_mask_cheat_float = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};
const union {long long int li[2];__m128d m;} fabs_mask_cheat_double = {0x7fffffffffffffff, 0x7fffffffffffffff};

namespace floatVconst {
#ifndef SWIG
  const int NV=4;
  const floatV4 zero=s2v(0.0f);
  const floatV4 half=s2v(0.5f);
  const floatV4 one=s2v(1.0f);
  const floatV4 two=s2v(2.0f);
#endif
};
namespace doubleVconst {
  const int NV=2;
  const doubleV2 zero=s2v(0.0);
  const doubleV2 half=s2v(0.5);
  const doubleV2 one=s2v(1.0);
  const doubleV2 two=s2v(2.0);
};

inline floatV4 fabs(const floatV4 v) { return _mm_and_ps(v, fabs_mask_cheat_float.m); }
inline doubleV2 fabs(const doubleV2 v) { return _mm_and_pd(v,fabs_mask_cheat_double.m); }
inline floatV4 Max(floatV4 v1,floatV4 v2) { return _mm_max_ps(v1,v2); }
inline floatV4 Min(floatV4 v1,floatV4 v2) { return _mm_min_ps(v1,v2); }
inline doubleV2 Max(doubleV2 v1,doubleV2 v2) { return _mm_max_pd(v1,v2); }
inline doubleV2 Min(doubleV2 v1,doubleV2 v2) { return _mm_min_pd(v1,v2); }
inline floatV4 rcp(floatV4 const a) { return _mm_rcp_ps(a); }
inline floatV4 rcp_sqrt(floatV4 const a) { return _mm_rsqrt_ps(a); }
inline floatV4 sqrt(floatV4 const a) { return _mm_sqrt_ps(a); }
inline doubleV2 sqrt(doubleV2 const a) { return _mm_sqrt_pd(a); }
inline doubleV2 rcp(doubleV2 const a) { return _mm_div_pd(doubleVconst::one, a); }
inline doubleV2 rcp_sqrt(doubleV2 const a) { return _mm_div_pd(doubleVconst::one, _mm_sqrt_pd(a)); }
typedef doubleV2 doubleV;
typedef floatV4 floatV;
//typedef float floatT;
//typedef doubleV2 floatV;
//typedef double floatT;

#define Iv2V(N,i) (doubleV2) {iFUNC((i)), iFUNC((i+N))}
//#define Iv2V(N,i) (floatV4) {iFUNC((i)), iFUNC((i+N)), iFUNC((i+2*N)), iFUNC((i+3*N))}
BASIC_VECTOR_DEF(double,doubleV2, 2);
BASIC_VECTOR_DEF(float,floatV4, 4);

#else //======BASIC_VECTOR_(AVX/SSE/noVEC)========================
#define BASIC_VECTOR_noVEC
#warning without SSE/AVX
const int NV=1;
typedef double doubleV;
//typedef double floatT;
typedef float floatV;
//typedef float floatT;
#define Iv2V(N,i) iFUNC((i))

template <> struct anyV<double> { static inline double s2v(const double v) { return v; } };
template <> struct anyV<float> { static inline float s2v(const float v) { return v; } };

inline double s2v(double const* a) { return (*a); }
inline float s2v(float const* a) { return (*a); }
inline double s2v(double const& a) { return (a); }
inline float s2v(float const& a) { return (a); }
inline double v2s(double const& a, int i) { return a; }
inline float v2s(float const& a, int i) { return a; }
inline double v2sRef(double const& a, int i) { return a; }
inline float v2sRef(float const& a, int i) { return a; }
inline float rcp(float const a) { return (1.0f/a); }
inline float rcp_sqrt(float const a) { return (1.0f/sqrtf(a)); }
inline float sqrt(float const a) { return sqrtf(a); }
inline double sqrt(double const a) { return sqrt(a); }


namespace floatVconst {
  const float zero=0.0f;
  const float half=0.5f;
  const float one=1.0f;
  const float two=2.0f;
};
namespace doubleVconst {
  const double zero=0.0;
  const double half=0.5;
  const double one =1.0;
  const double two =2.0;
};
#endif//======BASIC_VECTOR_(AVX/SSE/noVEC)========================

BASIC_noVECTOR_DEF(double);
BASIC_noVECTOR_DEF(float);
#undef BASIC_VECTOR_DEF
#undef BASIC_noVECTOR_DEF
#endif//BASIC_VECTOR_HPP
