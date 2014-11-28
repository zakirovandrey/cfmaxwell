#ifndef ACT_TORRE_FOLD_HEADER
#define ACT_TORRE_FOLD_HEADER
#include "ActThr.hpp"
#include "SemArr.hpp"

//ConeFold ранга rankTF разбивается на ConeTorre, которые состоят из ConeFold ранга rankCT
template <int dim, class T0, int rankTF, int rankCT> struct TorreFold {
  const static int NCT=1<<(rankTF-rankCT), NarrI=(NCT<<1)-1;
  SemCube<dim, rankTF-rankCT>& semC;
  cubeLR<dim,T0,rankTF>* data;
  int shifts[NarrI][dim];//массив сдвигов (относительных смещений по осям xs,s=0..dim-1) соседних элементов cubeLR<dim,T,rankCT>, участвующих в выполнении данного алгоритма
  inline int getSh(int iax, int ish) { return shifts[ish][iax]; }
  inline void setShN(int Ish[dim], int ish=NCT-1) { for(int s=0; s<dim; s++) shifts[ish][s] = Ish[s]; }
  inline void setShS(int Ish, int ish) { for(int s=0; s<dim; s++) shifts[ish][s] = Ish<<s; }
  TorreFold(SemCube<dim,rankTF-rankCT>& _semC, cubeLR<dim,T0,rankTF>* _data): semC(_semC), data(_data) {
    int Ish=1;
    for(int sh=2; sh<NarrI; sh<<=1) {
      for(int i=sh/2-1; i<NarrI; i+=sh) setShS(Ish, i);
      Ish=(Ish<<dim)-1;
    }
  }
  inline void initDd(int IshA[dim]) {
    for(int sh=2; sh<NarrI; sh<<=1) {
      for(int s=0; s<dim; s++) IshA[s]=(IshA[s]<<dim)-(1<<s);
    } setShN(IshA);
  }
};
//2**(dim*(rankH-rankJump)) ConeTorre, сложенные ConeFold ранга rankCT высотой 2**rankH
template <int dim, class T0, int rankCT, int rankH> struct ConeTorre: public ActionThreadStruct {
  const static int NCT=1<<rankH;
  TorreFold<dim,T0,rankCT+rankH,rankCT>* baseTF;
  cubeLR<dim,T0,rankCT>* data;
  int curIc[dim+1];
  inline int getIc(int iax) { return curIc[iax]; }
  inline int getSh(int iax, int ish=0) { return baseTF->getSh(iax,curIc[iax]+ish); }
  inline void IcLR2arr(int IcLR) {
    for(int s=0; s<=dim; s++) curIc[s] = 0;
    for(int r=0, ix=IcLR; r<rankH; r++) for(int s=0; s<dim; s++) {
      curIc[s] += (ix&1)<<r; ix >>= 1;
    }
    data = ((cubeLR<dim,T0,rankCT>*)baseTF->data)+IcLR;
  }
  inline bool next_shift() {
    if(++curIc[dim] >= NCT) return false;
    int sh=0;
    for(int s=0; s<dim; s++) { sh += getSh(s); curIc[s]++; }
    data += sh;
    return true;
  }
  void CTact() {
    const int IcLRmax=(1<<(dim*rankH))-1;
    for(int IcLR=IcLRmax-iThr; IcLR>=0; IcLR-=nThr) {
      IcLR2arr(IcLR);
      do {
        //printf("CTact thr[%d] %d: (%d,%d,%d)\n", iThr,IcLR, curIc[0],curIc[1],curIc[2]);
        baseTF->semC.wait(IcLRmax-IcLR);
        DDDact_(data, getSh(0), getSh(1), getSh(2));
        baseTF->semC.post(IcLRmax-IcLR);
      } while(next_shift());
    }
  }
};
#endif //ACT_TORRE_FOLD_HEADER

void* void_CTact_func(void* vact) {
  ActionThreadStruct* act=(ActionThreadStruct*)vact;
  //set_chip_affinity(nchip);
  act->set_core_affinity();
//  act->print_CPU_affinity();
  act->CTact();
  return 0;
}

#include <pthread.h>
#include <malloc.h>

//const int rankTF=4;//, dim=3;

#define DefAct(rankTF, rankCS)\
template <class T0> inline void DDDact(cubeLR<dim,T0,rankTF>* const data0, const int _Ipdd, const int _Idpd, const int _Iddp) {\
  const int rankCT=rankTF-rankCS, Nthr=Cores_perChip;\
  typedef ConeTorre<dim,T0,rankCT,rankCS> typeConeTorre;\
  typedef TorreFold<dim,T0,rankTF,rankCT> typeTorreFold;\
  SemCube<dim,rankCS> sems; sems.init();\
  int Ish[3]; Ish[0] = _Ipdd; Ish[1] = _Idpd; Ish[2] = _Iddp;\
  typeTorreFold TF(sems, data0); TF.initDd(Ish);\
  typeConeTorre CTs[Nthr];\
  for(int ithr=0; ithr<Nthr; ithr++) {\
    CTs[ithr].prepare(ithr, Nthr);\
    CTs[ithr].baseTF = &TF;\
    pthread_create(&CTs[ithr].ThreadID, NULL, &void_CTact_func, CTs+ithr);\
  }\
  void* res=0;\
  for(int ithr=0; ithr<Nthr; ithr++) pthread_join(CTs[ithr].ThreadID,&res);\
  sems.destroy();\
}
//DefAct(1,1)
//DefAct(2,2)
DefAct(3,2)
DefAct(4,2)
DefAct(5,3)
DefAct(6,3)
DefAct(7,4)
DefAct(8,4)
DefAct(9,4)
/*
template <class T0> inline void ActDDD(cubeLR<dim,T0,rankTF>* const data0, const int _Ipdd, const int _Idpd, const int _Iddp) {
  const int rankCS=3, rankCT=rankTF-rankCS, Nthr=Cores_perChip;
  typedef ConeTorre<dim,T0,rankCT,rankCS> typeConeTorre;
  typedef TorreFold<dim,T0,rankTF,rankCT> typeTorreFold;
  SemCube<dim,rankCS> sems; sems.init();
  int Ish[3]; Ish[0] = _Ipdd; Ish[1] = _Idpd; Ish[2] = _Iddp;
  typeTorreFold TF(sems, data0); TF.initDd(Ish);
  typeConeTorre CTs[Nthr];
  for(int ithr=0; ithr<Nthr; ithr++) {
    CTs[ithr].prepare(ithr, Nthr);
    CTs[ithr].baseTF = &TF;
    pthread_create(&CTs[ithr].ThreadID, NULL, &void_CTact_func, CTs+ithr);
  }
  void* res=0;
  for(int ithr=0; ithr<Nthr; ithr++) pthread_join(CTs[ithr].ThreadID,&res);
  sems.destroy();
}*/
