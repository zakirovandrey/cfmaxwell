
inline void err_linkType(char* linkType, int s) { printf("Illegal linkType char (%s)[%d]\n",linkType,s); }

#ifdef DIM_3D
using LRindex::ind3_2;
using LRindex::ind3_1;

#include "nLAdat.inc.hpp"
template <> void ConeFold<nLArank>(nLAnet<3>* net, unsigned int I) {
  using namespace LRindex;
  const int dim=3, iB=1, rank=nLArank; // iB=1 ConeFold "падает направо"; iB=0 "налево"
  const unsigned int M=mask[3]&((1<<dim*rank)-1);
  //brick_header<dim>* d=net->datas[7];
  //printf("%s CF: (%d%d%d) I=%d\n", net->linkType, d->index[0],d->index[1],d->index[2], I);
           int IM0=I&M , IS0=IM0-1, I0=IS0&M , Sh0=IM0-I0, iX =IS0<0?1:0;
  int M1=M<<1, IM1=I&M1, IS1=IM1-2, I1=IS1&M1, Sh1=IM1-I1; iX+=IS1<0?2:0;
  int M2=M<<2, IM2=I&M2, IS2=IM2-4, I2=IS2&M2, Sh2=IM2-I2; iX+=IS2<0?4:0;
#include "nLAnet.inc.hpp"
}
#endif //DIM_3D

#ifdef DIM_2D

template <int rank,int ind> cubeLR<2,dataDD,MaxRank-rank>* get_dd(nLAnet<2>* net, unsigned int I) { return (((multi_brick<2,2,dataDD,MaxRank-rank>*)net->datas[ind])->data)+((1<<2*rank)-1-I); }
template <int rank,int ind> cubeLR<1,dataRD,MaxRank-rank>* get_xd(nLAnet<2>* net, unsigned int I) { return (((multi_brick<2,1,dataRD,MaxRank-rank>*)net->datas[ind])->data)+((1<<rank)-1-ind2_1<rank,1>(I)); }
template <int rank,int ind> cubeLR<1,dataDR,MaxRank-rank>* get_dx(nLAnet<2>* net, unsigned int I) { return (((multi_brick<2,1,dataDR,MaxRank-rank>*)net->datas[ind])->data)+((1<<rank)-1-ind2_1<rank,0>(I)); }
template <int rank,int ind> cubeLR<0,dataRR,MaxRank-rank>* get_xx(nLAnet<2>* net, unsigned int I) { return (((multi_brick<2,0,dataRR,MaxRank-rank>*)net->datas[ind])->data); }

template <> void ConeFold<nLArank>(nLAnet<2>* net, unsigned int I) {
  const int dim=2, iB=1, rank=nLArank, datRank=MaxRank-nLArank; // iB=1 ConeFold "падает направо"; iB=0 "налево"
  const unsigned int M=0x55555555&((1<<(dim*rank))-1);
  //brick_header<dim>* d=net->datas[3];
  //printf("%s CF: (%d%d) I=%d\n", net->linkType, d->index[0],d->index[1], I);
           int IM0=I&M , IS0=IM0-1, I0=IS0&M , Sh0=IM0-I0, iX =IS0<0?1:0;
  int M1=M<<1, IM1=I&M1, IS1=IM1-2, I1=IS1&M1, Sh1=IM1-I1; iX+=IS1<0?2:0;
#include "nLAnet2D.inc.hpp"
}
#endif //DIM_2D

void* nLAnetArrayNUMA_run_step(void* t) { nLAnetArrayNUMA<3>* pt=(nLAnetArrayNUMA<3>*)t; pt->run_step(); return t; }
void* nLAnetArrayNUMA_runChessFold(void* t) { nLAnetArrayNUMA<3>* pt=(nLAnetArrayNUMA<3>*)t; pt->runChessFold(); return t; }
