#ifndef LR_CUBE_HPP
#define LR_CUBE_HPP
#include <cstring>
//=======================
template <int dim, class T, int rank> struct cubeLR {
  cubeLR<dim, T, rank-1> datas[1<<dim];
  T& operator[](size_t i) { return ((T*)datas)[i]; }
  size_t getSize() { return ((size_t)1)<<(rank*dim); }
};
template <int dim, class T> struct cubeLR<dim, T, 1> {
  T datas[1<<dim];
  T& operator[](size_t i) { return datas[i]; }
  size_t getSize() { return 1<<dim; }
};
inline int getLRshift(int x, int y, int rank) {
  int Xpos = x/(1<<(rank-1));
  int Ypos = y/(1<<(rank-1));
  int Xsh  = x%(1<<(rank-1));
  int Ysh  = y%(1<<(rank-1));
  int Bsh = (Xpos+2*Ypos)*(1<<(2*(rank-1)));
  if (rank==1) return Bsh;
  else return Bsh+getLRshift(Xsh,Ysh,rank-1);
}
#endif//LR_CUBE_HPP

#ifndef LR_ZIPPER_HPP
#define LR_ZIPPER_HPP
//=======================
namespace LRindex {
const unsigned int mask[]={0, 0, 0x55555555, 011111111111};
const unsigned long long int maskL[]={0L, 0L, 0x5555555555555555L, 0111111111111111111111L};

template <int dim, int rank> struct knot {
  static const unsigned int rMask=(1<<(dim*rank))-1;
  const unsigned int index;
  inline unsigned int inc(int s) { return rMask&((((index|~(mask[dim]<<s))+(1<<s))&(mask[dim]<<s))|(index&~(mask[dim]<<s))); }
  inline unsigned int dec(int s) { return rMask&((((index|~(mask[dim]<<s))-(1<<s))&(mask[dim]<<s))|(index&~(mask[dim]<<s))); }
  inline unsigned int cyclicRrot() { return ((index<<(dim-1))&mask[dim])|((index>>1)&~mask[dim]); }
  inline unsigned int cyclicLrot() { return ((index>>(dim-1))&mask[dim])|((index<<1)&~mask[dim]); }
  bool isL(int s) { return (rMask&(index & mask[dim]<<s)) == 0; }
  bool isR(int s) { return (rMask&(index|~(mask[dim]<<s))) == rMask; }
};
template <int dim, int rank> struct knotL {
  static const unsigned int rMask=(1L<<(dim*rank))-1;
  const unsigned long long int index;
  inline unsigned long long int inc(int s) { return (((index|~(maskL[dim]<<s))+(1<<s))&(maskL[dim]<<s))|(index&~(maskL[dim]<<s)); }
  inline unsigned long long int dec(int s) { return (((index|~(maskL[dim]<<s))-(1<<s))&(maskL[dim]<<s))|(index&~(maskL[dim]<<s)); }
  inline unsigned long long int cyclicRrot() { return ((index<<(dim-1))&maskL[dim])|((index>>1)&~maskL[dim]); }
  inline unsigned long long int cyclicLrot() { return ((index>>(dim-1))&maskL[dim])|((index<<1)&~maskL[dim]); }
  bool isL(int s) { return (rMask&(index & maskL[dim]<<s)) == 0; }
  bool isR(int s) { return (rMask&(index|~(maskL[dim]<<s))) == rMask; }
};

template <int dim> const unsigned int inc(const unsigned int index, int s)
{ return (((index|~(mask[dim]<<s))+(1<<s))&(mask[dim]<<s))|(index&~(mask[dim]<<s)); }
template <int dim> const unsigned int dec(const unsigned int index, int s)
{ return (((index|~(mask[dim]<<s))-(1<<s))&(mask[dim]<<s))|(index&~(mask[dim]<<s)); }

template <int dim> const unsigned long long int inc(const unsigned long long int index, int s)
{ return (((index|~(maskL[dim]<<s))+(1<<s))&(maskL[dim]<<s))|(index&~(maskL[dim]<<s)); }
template <int dim> const unsigned long long int dec(const unsigned long long int index, int s)
{ return (((index|~(maskL[dim]<<s))-(1<<s))&(maskL[dim]<<s))|(index&~(maskL[dim]<<s)); }

template <int dim> const unsigned int cyclicLrot(const unsigned int index)
{ return ((index>>2)&mask[dim])|((index<<1)&~mask[dim]); }
template <int dim> const unsigned int apply_mask(const unsigned int index, int s)
{ return index&(mask[dim]<<s); }
template <int dim> const unsigned int applyUmask(const unsigned int index, int s)
{ return index|~(mask[dim]<<s); }

//============================
template <int rank, int s> inline unsigned int ind3_2(unsigned int I0) {
  unsigned int I;
  //const unsigned int M0=0xDB6DB6DB, M[]={0x030C30C3,0x030C30C3<<3,0x0F00F00F,0x0F00F00F<<6,0x0F0000FF,0x000FF000,0x0000FFFF,0x0F000000};
  const unsigned int M0=033333333333, M[]={030303030303,0x0F00F00F,0x0F0000FF,0x0000FFFF};
  if(s==0) I=I0>>1 & M0;
  else if(s==1) I=I0>>2 & M0;
  else if(s==2) I=I0&M0;
  if(rank>1) I = I&M[0] | (I&M[0]<<3)>>1;
  if(rank>2) I = I&M[1] | (I&M[1]<<6)>>2;
  if(rank>4) I = I&M[2] | (I&M[2]<<12)>>4;
  if(rank>8) I = I&M[3] | (I&M[3]<<24)>>8;
  if(s==1) I=I<<1 | I0&1;
  return I;
}
template <int rank, int s> inline unsigned int ind3_1(unsigned int I0) {
  unsigned int I;
  //const unsigned int M0=0x49249249, M[]={0x41041041,0x28282828,0x03003003,0x000C00C0,0x0300000F,0x0000F000,0x000000FF,0x03000000};
  const unsigned int M0=011111111111, M[]={010101010101,0x03003003,0x0300000F,0x000000FF};
  if(s>0) I=I0>>s & M0;
  else I=I0&M0;
  if(rank>1) I = I&M[0] | (I&M[0]<<3)>>2;
  if(rank>2) I = I&M[1] | (I&M[1]<<6)>>4;
  if(rank>4) I = I&M[2] | (I&M[2]<<12)>>8;
  if(rank>8) I = I&M[3] | (I&M[3]<<24)>>16;
  return I;
}

template <int dim, int rank> struct zipper {
  unsigned int i[1];
  inline unsigned int& operator [](int s) { return i[0]; }
  inline void unzip(unsigned int I0) { i[0]=I0; }
  inline unsigned int zip() { return i[0]; }
};
struct zipper2D {
  unsigned int i[2];
  inline unsigned int& operator [](int s) { return i[s]; }
  inline void unzip(unsigned int I0) {
    const unsigned int M[]={0x11111111,0x03030303,0x000F000F,0x000000FF};
    unsigned int I1=I0>>1;
    I0 = (I0&M[0]<< 2)>>1 | I0&M[0]; I1 = (I1&M[0]<< 2)>>1 | I1&M[0];
    I0 = (I0&M[1]<< 4)>>2 | I0&M[1]; I1 = (I1&M[1]<< 4)>>2 | I1&M[1];
    I0 = (I0&M[2]<< 8)>>4 | I0&M[2]; I1 = (I1&M[2]<< 8)>>4 | I1&M[2];
    I0 = (I0&M[3]<<16)>>8 | I0&M[3]; I1 = (I1&M[3]<<16)>>8 | I1&M[3];
    i[0]=I0; i[1]=I1;
  }
  inline unsigned int zip() {
    const unsigned int M[]={0x11111111,0x03030303,0x000F000F,0x000000FF};
    unsigned int I0=i[0], I1=i[1];
    I0 = (I0&M[3]<<8)<<8 | I0&M[3]; I1 = (I1&M[3]<<8)<<8 | I1&M[3];
    I0 = (I0&M[2]<<4)<<4 | I0&M[2]; I1 = (I1&M[2]<<4)<<4 | I1&M[2];
    I0 = (I0&M[1]<<2)<<2 | I0&M[1]; I1 = (I1&M[1]<<2)<<2 | I1&M[1];
    I0 = (I0&M[0]<<1)<<1 | I0&M[0]; I1 = (I1&M[0]<<1)<<1 | I1&M[0];
    return I0 | I1<<1;
  }
};
template <int rank> struct zipper<2,rank>: public zipper2D {
  inline void unzip(unsigned int I0) {
    const unsigned int M[]={0x11111111,0x03030303,0x000F000F,0x000000FF};
    unsigned int I1=I0>>1;
    if(rank>1) { I0 = (I0&M[0]<< 2)>>1 | I0&M[0]; I1 = (I1&M[0]<< 2)>>1 | I1&M[0]; }
    if(rank>2) { I0 = (I0&M[1]<< 4)>>2 | I0&M[1]; I1 = (I1&M[1]<< 4)>>2 | I1&M[1]; }
    if(rank>4) { I0 = (I0&M[2]<< 8)>>4 | I0&M[2]; I1 = (I1&M[2]<< 8)>>4 | I1&M[2]; }
    if(rank>8) { I0 = (I0&M[3]<<16)>>8 | I0&M[3]; I1 = (I1&M[3]<<16)>>8 | I1&M[3]; }
    //is = I0 | (unsigned long long int)I1<<32;
    i[0]=I0; i[1]=I1;
  }
  inline unsigned int zip() {
    const unsigned int M[]={0x11111111,0x03030303,0x000F000F,0x000000FF};
    unsigned int I0=i[0], I1=i[1];
    //unsigned long long int IL=is;
    //unsigned int I0=IL&0xFFFFFFFF, I1=IL>>32&0xFFFFFFFF;
    if(rank>8) { I0 = (I0&M[3]<<8)<<8 | I0&M[3]; I1 = (I1&M[3]<<8)<<8 | I1&M[3]; }
    if(rank>4) { I0 = (I0&M[2]<<4)<<4 | I0&M[2]; I1 = (I1&M[2]<<4)<<4 | I1&M[2]; }
    if(rank>2) { I0 = (I0&M[1]<<2)<<2 | I0&M[1]; I1 = (I1&M[1]<<2)<<2 | I1&M[1]; }
    if(rank>1) { I0 = (I0&M[0]<<1)<<1 | I0&M[0]; I1 = (I1&M[0]<<1)<<1 | I1&M[0]; }
    return I0 | I1<<1;
  }
};
template <> struct zipper<2,1>: public zipper2D {
  inline void unzip(unsigned int I0) { i[0]=I0&1; i[1]=I0>>1&1; }
  inline unsigned int zip() { return i[0] | i[1]<<1; }
};
struct zipper3D {
  unsigned short int i[4];
  inline unsigned short int& operator [](int s) { return i[s]; }
  inline void unzip(unsigned int Is) {
    unsigned long long int I0=Is, I=I0 | I0<<32;
    const unsigned long long int M[]={0xAAAAAAAAAAAAAAAA,0xCCCCCCCCCCCCCCCC,0xF0F0F0F0F0F0F0F0,0xFF00FF00FF00FF00};
    I = I&M[0]>>1 | (I&M[0])>>2;
    I = I&M[1]>>2 | (I&M[1])>>4;
    I = I&M[2]>>4 | (I&M[2])>>8;
    I = I&M[3]>>8 | (I&M[3])>>16;
    const unsigned int M0=0x3FF;
    i[0] = I&M0; i[1] = I>>11 & M0; i[2] = I>>22 & M0;
  }
  inline unsigned int zip() {
    unsigned long long int i2=i[2];
    unsigned long long int I=i[0] | i[1]<<11 | i2<<22;
    const unsigned long long int M[]={0xAAAAAAAAAAAAAAAA,0xCCCCCCCCCCCCCCCC,0xF0F0F0F0F0F0F0F0,0xFF00FF00FF00FF00};
    I = I&M[3]>>8 | (I&M[3])<<16;
    I = I&M[2]>>4 | (I&M[2])<< 8;
    I = I&M[1]>>2 | (I&M[1])<< 4;
    I = I&M[0]>>1 | (I&M[0])<< 2;
    const unsigned int M0=0xFFFFFFFF, I1=I&M0 | I>>32&M0;
    return I1;
  }
};
template <int rank> struct zipper<3,rank>: public zipper3D {
  inline void unzip(unsigned int I0) { ((zipper<3,rank+1>*)this)->unzip(I0); }
  inline unsigned int zip() { return ((zipper<3,rank+1>*)this)->zip(); }
};
template <> struct zipper<1,1>: public zipper3D {
  inline void unzip(unsigned int I) { i[0] = I&1; i[1] = I>>1 & 1; i[2] = I>>2 & 1; }
  inline unsigned int zip() { return i[0] | i[1]<<1 | i[2]<<2; }
};
template <> struct zipper<3,3>: public zipper3D {
  inline void unzip(unsigned int I0) {
    unsigned int I=I0 | (I0&0x26)<<8;
    const unsigned int M[]={0xAAAA,0xCCCC};
    I = I&M[0]>>1 | (I&M[0])>>2;
    I = I&M[1]>>2 | (I&M[1])>>4;
    const unsigned int M0=7;
    i[0] = I&M0; i[1] = I>>3 & M0; i[2] = I>>6 & M0;
  }
  inline unsigned int zip() {
    unsigned int I=i[0] | i[1]<<3 | i[2]<<6;
    const unsigned int M[]={0xAAAA,0xCCCC};
    I = I&M[1]>>2 | (I&M[1])<< 4;
    I = I&M[0]>>1 | (I&M[0])<< 2;
    return I&0x1FF | I>>8&0x26;
  }
};
template <> struct zipper<3,5>: public zipper3D {
  inline void unzip(unsigned int I0) {
    unsigned int I=I0 | I0<<16;
    const unsigned int M[]={0xAAAAAAAA,0xCCCCCCCC,0xF0F0F0F0};
    I = I&M[0]>>1 | (I&M[0])>>2;
    I = I&M[1]>>2 | (I&M[1])>>4;
    I = I&M[2]>>4 | (I&M[2])>>8;
    const unsigned int M0=0x1F;
    i[0] = I&M0; i[1] = I>>11 & M0; i[2] = I>>6 & M0;
  }
  inline unsigned int zip() {
    unsigned int I=i[0] | i[1]<<11 | i[2]<<6;
    const unsigned int M[]={0xAAAAAAAA,0xCCCCCCCC,0xF0F0F0F0};
    I = I&M[2]>>4 | (I&M[2])<< 8;
    I = I&M[1]>>2 | (I&M[1])<< 4;
    I = I&M[0]>>1 | (I&M[0])<< 2;
    const unsigned int M0=0xFFFF;
    return I&M0 | I>>16&M0;
  }
};
template <> struct zipper<3,10>: public zipper3D {
  inline void unzip(unsigned int Is) {
    unsigned long long int I0=Is, I=I0 | I0<<32;
    const unsigned long long int M[]={0xAAAAAAAAAAAAAAAA,0xCCCCCCCCCCCCCCCC,0xF0F0F0F0F0F0F0F0,0xFF00FF00FF00FF00};
    I = I&M[0]>>1 | (I&M[0])>>2;
    I = I&M[1]>>2 | (I&M[1])>>4;
    I = I&M[2]>>4 | (I&M[2])>>8;
    I = I&M[3]>>8 | (I&M[3])>>16;
    const unsigned int M0=0x3FF;
    i[0] = I&M0; i[1] = I>>11 & M0; i[2] = I>>22 & M0;
  }
  inline unsigned int zip() {
    unsigned long long int i2=i[2];
    unsigned long long int I=i[0] | i[1]<<11 | i2<<22;
    const unsigned long long int M[]={0xAAAAAAAAAAAAAAAA,0xCCCCCCCCCCCCCCCC,0xF0F0F0F0F0F0F0F0,0xFF00FF00FF00FF00};
    I = I&M[3]>>8 | (I&M[3])<<16;
    I = I&M[2]>>4 | (I&M[2])<< 8;
    I = I&M[1]>>2 | (I&M[1])<< 4;
    I = I&M[0]>>1 | (I&M[0])<< 2;
    const unsigned int M0=0xFFFFFFFF, I1=I&M0 | I>>32&M0;
    return I1;
  }
};
};
/*
namespace LRind2D {
  const unsigned int M=0x55555555;
  unsigned int inc_x(unsigned int I) { return (I|~M)+1&M|I&~M; }
  unsigned int inc_y(unsigned int I) { return (I|M)+2&~M|I&M; }
  unsigned int inc_xy(unsigned int I) { return (I|~M)+1&M|(I|M)+2&~M; }
}
namespace LRind3D {
  const unsigned int M=0x49249249;
  unsigned int inc_x(unsigned int I) { return (I|~M)+1&M|I&~M; }
  unsigned int inc_y(unsigned int I) { return (I|~(M<<1))+2&(M<<1)|I&~(M<<1); }
  unsigned int inc_z(unsigned int I) { return (I|~(M<<2))+4&(M<<2)|I&~(M<<2); }
  unsigned int inc_xyz(unsigned int I) { return (I|~M)+1&M|(I|~(M<<1))+2&(M<<1)|(I|~(M<<2))+4&(M<<2); }
};
*/
#endif//LR_ZIPPER_HPP
