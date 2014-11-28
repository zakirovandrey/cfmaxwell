#ifndef nLA_DATA_HPP
#define nLA_DATA_HPP
#include <fcntl.h>
#include <sys/mman.h>

template <int dim> struct BrickPos;

//--------------------------
template <int dim> struct brick_header {
  int netID, XDind, data_dim, iNUMAnode;
  static int Pid;
  virtual size_t get_size()=0;
  brick_header(int ID, int inode=-1): netID(ID), iNUMAnode(inode), XDind(0) {
    //for(int s=0; s<dim; s++) index[s] = ind[s];
  }
  void operator delete(void* p);
  void setXDind(int xd) {
    XDind = xd; data_dim = dim;
    for(int s=0; s<dim; s++) if(xd&(1<<s)) data_dim--;
  }
  virtual void init(BrickPos<dim>& pos) = 0;
  virtual void drop(BrickPos<dim>& pos) = 0;
};
template<int dim> int brick_header<dim>::Pid=0;
template <int dim, int data_dim, class Tdat, int rank> struct brick: public brick_header<dim> {
  cubeLR<data_dim,Tdat,rank> __attribute__ ((aligned(64))) data;
  virtual size_t get_size() { return sizeof(*this); }
  void* operator new(size_t sz, int inode=-1);
  brick(int ID, int inode=-1): brick_header<dim>(ID, inode) {}
  virtual void init(BrickPos<dim>& pos);
  virtual void drop(BrickPos<dim>& pos);
};

template <int dim, int data_dim, class Tdat, int rank> struct brickPML: public brick_header<dim> {
  cubeLR<data_dim,Tdat,rank> __attribute__ ((aligned(64))) data;
  virtual size_t get_size() { return sizeof(*this); }
  void* operator new(size_t sz, int inode=-1);
  brickPML(int ID, int inode=-1): brick_header<dim>(ID, inode) {}
  virtual void init(BrickPos<dim>& pos);
  virtual void drop(BrickPos<dim>& pos);
};
template <int dim, int data_dim, class Tdat, int rank> struct multi_brick: public brick_header<dim> {
  cubeLR<data_dim,Tdat,rank>* data;
  int shift[data_dim];
  virtual size_t get_size() { return sizeof(*this); }
  //void* operator new(size_t sz, cubeLR<data_dim,Tdat,rank>* data, int inode=-1);
  multi_brick(int ind[dim], int* sh, cubeLR<data_dim,Tdat,rank>* dat, int inode=-1):
    brick_header<dim>(ind, inode), data(dat) { for(int s=0; s<data_dim; s++) shift[s] = sh[s]; }
  virtual void init(BrickPos<dim>& pos);
  virtual void drop(BrickPos<dim>& pos);
  //virtual void drop(int it, char*, char*);
};
#ifdef NUMA
#include <numa.h>
#endif
extern size_t ProblemSize;

template <int dim>
void brick_header<dim>::operator delete(void* p) {
  brick_header<dim>* br=(brick_header<dim>*)p;
  //printf("%p, %d ind=(%d,%d,%d)\n", p, br->iNUMAnode, br->index[0], br->index[1], br->index[2]);
  #ifdef NUMA
  if(br->iNUMAnode>=0) numa_free(p, br->get_size());
  else
  #endif
  free(p);
}

template <int dim, int data_dim, class Tdat, int rank>
void* brickPML<dim, data_dim,Tdat,rank>::operator new(size_t sz, int inode) {
  ProblemSize += sz;
//  printf("new brickPML<dim=%d, data_dim=%d, Tdat, rank=%d> Inode=%d, Size: %ld of %.3gG\n", dim, data_dim, rank, inode, sz, ProblemSize/(1024*1024*1024.));
  #ifdef NUMA
  fflush(stdout);
  if(inode>=0) return numa_alloc_onnode(sz, inode);
  else
  #endif
  #ifndef USE_SWAP
  return valloc(sz);
  #endif
  void* pf=0;
  char fname[256]; sprintf(fname,"%s/SwapPML%d.swp",pars.swapDir.c_str(),brickPML::Pid);
  remove(fname);
  int swp_file = open(fname, O_CREAT|O_RDWR, 0666);
//  posix_fadvise(swp_file, 0, sz, POSIX_FADV_DONOTNEED);
  lseek(swp_file, sz, SEEK_SET);
  write(swp_file, "", 1);
  lseek(swp_file, 0, SEEK_SET);
  if ( (pf = mmap(NULL, sz, PROT_READ|PROT_WRITE, MAP_SHARED, swp_file, 0)) == MAP_FAILED ) perror("Cannot mmap!\n");
  close(swp_file);
  printf("swap in PML %p\n",pf);
  brickPML::Pid++;
  return pf;
}
template <int dim, int data_dim, class Tdat, int rank>
void* brick<dim, data_dim,Tdat,rank>::operator new(size_t sz, int inode) {
  ProblemSize += sz;
//  printf("new brick<dim=%d, data_dim=%d, Tdat, rank=%d> Inode=%d, Size: %ld of %.3gG\n", dim, data_dim, rank, inode, sz, ProblemSize/(1024*1024*1024.));
  #ifdef NUMA
  if(inode>=0) return numa_alloc_onnode(sz, inode);
  else 
  #endif
  #ifndef USE_SWAP
  return valloc(sz);
  #endif
  void* pf=0;
  char fname[256]; sprintf(fname,"%s/Swap%d.swp",pars.swapDir.c_str(),brick::Pid);
  remove(fname);
  int swp_file = open(fname, O_CREAT|O_RDWR, 0666);
//  posix_fadvise(swp_file, 0, sz, POSIX_FADV_DONOTNEED);
  lseek(swp_file, sz, SEEK_SET);
  write(swp_file, "", 1);
  lseek(swp_file, 0, SEEK_SET);
  if ( (pf = mmap(NULL, sz, PROT_READ|PROT_WRITE, MAP_SHARED, swp_file, 0)) == MAP_FAILED ) perror("Cannot mmap!\n");
  close(swp_file);
  printf("swap in %p\n",pf);
  brick::Pid++;
  return pf;
}

//===============================
template <int dim> struct nLAnet {
  //static const unsigned int iB=0, iF=1;//ConeFold "падает налево",  индекс I растёт слева неправо
  static const unsigned int iB=1, iF=0;//ConeFold "падает направо", индекс I растёт справа налево
  nLAnet<dim>* net[dim+1][2];
  nLAnet<dim>* netX[1<<dim];
  brick_header<dim>* datas[1<<dim];
  unsigned int indL, indR;
  int ID, iStep, iDGtier;
  char linkType[dim+1];
  double Trun;
  void clear() {
    iStep = 0; iDGtier = 0; Trun = 0.0;
    for(int s=0; s<=dim; s++) net[s][0] = net[s][1] = 0;
    for(int i=0; i<(1<<dim); i++) netX[i] = 0;
    for(int i=0; i<(1<<dim); i++) datas[i] = 0;
    linkType[0]='\0'; linkType[dim]='\0';
    indL = indR = 0;
  }
  void setLinkL(int s, nLAnet<dim>& n) { net[s][0] = &n; n.net[s][1] = this; }
  void setLinkR(int s, nLAnet<dim>& n) { net[s][1] = &n; n.net[s][0] = this; }
  void copyNBdatasL() {//копирует указатели на данные ячеек, на которые уже ссылаются "левые" от текущего узлы
    for(int s=0; s<dim; s++) if(net[s][0]) { // перебираем все оси s=0..d-1, и для каждого "левого" узла по этой оси, если он существует, выполняем:
      for(int ib=0; ib<(1<<(dim-1)); ib++) { // перебираем 2^{d-1} индексов, по числу ячеек, одинаковых у текущего и левого узлов
        int im=ib, ip=ib<<1; // искомый индекс разбиваем на две части, в бинарной записи: ip[01]im
        for(int i=0; i<=s; i++) ip &= ~(1<<i);
        for(int i=s; i<dim; i++) im &= ~(1<<i);
        //printf("s=%d, ib=%d, im=%d, ip=%d, ind=(%d,%d)\n", s, ib, im, ip, im+ip, im+ip+(1<<s));
        datas[im+ip] = net[s][0]->datas[im+ip+(1<<s)];
      }
    }
  }
  void set_netX(const int iB=1) {
    for(int i=0; i<(1<<dim); i++) {
      nLAnet<dim>* cur=this;
      for(int s=0; cur && s<dim; s++) if(i&(1<<s)) cur = cur->net[s][iB];
      netX[i] = cur;
    }
    //net[dim][iB] = netX[(1<<dim)-1];
    //if(net[dim][iB]) net[dim][iB]->net[dim][1-iB] = this;
  }
  virtual bool check4grow() {
    for(int s=0; s<dim; s++) if(net[s][iB] && iStep >= net[s][iB]->iStep) return false;
    if(net[dim][iF] && iStep > net[dim][iF]->iStep) return false;
    return true;
  }// printf("start grow, %d\n", iDGtier); }
  virtual void finish_grow() { iStep++; }
};

//тут rank растёт от 0 до nLArank, что соответствует уменьшению rank в LR от MaxRank до MaxRank-nLArank
template <int rank> void ConeFold(nLAnet<3>* net, unsigned int I) {
  const int dim=3; const unsigned int M=0x49249249&((1<<(dim*(rank+1)))-1);
  //const unsigned int iB=0, netX=net->indR;//ConeFold "падает налево",  индекс I растёт слева неправо
  const unsigned int iB=1, netX=net->indL;//ConeFold "падает направо", индекс I растёт справа налево
//brick_header<dim>* d=net->datas[7];
//printf("%s CF<%d>(%2d): %d%d%d, netX=%d", net->linkType, rank, I, d->index[0],d->index[1],d->index[2], netX);
  I<<=dim; int IM0=I&M , IS0=IM0-1, I0=IS0&M, iX=IS0<0?1:0;
  int M1=M<<1, IM1=I&M1, IS1=IM1-2, I1=IS1&M1; iX+=IS1<0?2:0;
  int M2=M<<2, IM2=I&M2, IS2=IM2-4, I2=IS2&M2; iX+=IS2<0?4:0;
//printf(" iX=%d, (%d,%d,%d)\n", iX, I0,I1,I2);
  if(iX + netX == 0) {
    ConeFold<rank+1>(net,I  ); ConeFold<rank+1>(net,I0|I1|I2);
    ConeFold<rank+1>(net,I+1); ConeFold<rank+1>(net,IM0|I1|I2);
    ConeFold<rank+1>(net,I+2); ConeFold<rank+1>(net,I0|IM1|I2);
    ConeFold<rank+1>(net,I+3); ConeFold<rank+1>(net,I&~M2|I2);
    ConeFold<rank+1>(net,I+4); ConeFold<rank+1>(net,I0|I1|IM2);
    ConeFold<rank+1>(net,I+5); ConeFold<rank+1>(net,I&~M1|I1);
    ConeFold<rank+1>(net,I+6); ConeFold<rank+1>(net,I&~M |I0);
    ConeFold<rank+1>(net,I+7); ConeFold<rank+1>(net,I);
  } else {
                    ConeFold<rank+1>(net,I  ); if(net->netX[iX  ]) ConeFold<rank+1>(net->netX[iX  ],I0|I1|I2);
    if((netX&1)==0) ConeFold<rank+1>(net,I+1); if(net->netX[iX&6]) ConeFold<rank+1>(net->netX[iX&6],IM0|I1|I2);
    if((netX&2)==0) ConeFold<rank+1>(net,I+2); if(net->netX[iX&5]) ConeFold<rank+1>(net->netX[iX&5],I0|IM1|I2);
    if((netX&4)==0) ConeFold<rank+1>(net,I+4); if(net->netX[iX&3]) ConeFold<rank+1>(net->netX[iX&3],I0|I1|IM2);
    if((netX&3)==0) ConeFold<rank+1>(net,I+3); if(net->netX[iX&4]) ConeFold<rank+1>(net->netX[iX&4],I&~M2|I2);
    if((netX&5)==0) ConeFold<rank+1>(net,I+5); if(net->netX[iX&2]) ConeFold<rank+1>(net->netX[iX&2],I&~M1|I1);
    if((netX&6)==0) ConeFold<rank+1>(net,I+6); if(net->netX[iX&1]) ConeFold<rank+1>(net->netX[iX&1],I&~M |I0);
    if((netX  )==0) ConeFold<rank+1>(net,I+7); ConeFold<rank+1>(net,I);
  }
}

//тут rank растёт от 0 до nLArank, что соответствует уменьшению rank в LR от MaxRank до MaxRank-nLArank
template <int rank> void ConeFold(nLAnet<2>* net, unsigned int I) {
  const int dim=2; const unsigned int M=0x55555555&((1<<(dim*(rank+1)))-1);
  //const unsigned int iB=0, netX=net->indR;//ConeFold "падает налево",  индекс I растёт слева неправо
  const unsigned int iB=1, netX=net->indL;//ConeFold "падает направо", индекс I растёт справа налево
//brick_header<dim>* d=net->datas[(1<<dim)-1];
//printf("%s CF<%d>(%2d): %d%d, netX=%d", net->linkType, rank, I, d->index[0],d->index[1], netX);
  I<<=dim; int IM0=I&M , IS0=IM0-1, I0=IS0&M , iX =IS0<0?1:0;
  int M1=M<<1, IM1=I&M1, IS1=IM1-2, I1=IS1&M1; iX+=IS1<0?2:0;
//printf(" iX=%d, (%d,%d)\n", iX, I0,I1);
  if(iX + netX == 0) {
    ConeFold<rank+1>(net,I  ); ConeFold<rank+1>(net,I0|I1);
    ConeFold<rank+1>(net,I+1); ConeFold<rank+1>(net,IM0|I1);
    ConeFold<rank+1>(net,I+2); ConeFold<rank+1>(net,I0|IM1);
    ConeFold<rank+1>(net,I+3); ConeFold<rank+1>(net,I);
  } else {
                    ConeFold<rank+1>(net,I  ); if(net->netX[iX  ]) ConeFold<rank+1>(net->netX[iX  ],I0|I1);
    if((netX&1)==0) ConeFold<rank+1>(net,I+1); if(net->netX[iX&2]) ConeFold<rank+1>(net->netX[iX&2],IM0|I1);
    if((netX&2)==0) ConeFold<rank+1>(net,I+2); if(net->netX[iX&1]) ConeFold<rank+1>(net->netX[iX&1],I0|IM1);
    if((netX  )==0) ConeFold<rank+1>(net,I+3); ConeFold<rank+1>(net,I);
  }
}

#include "ActThr.hpp"

//============================
template <int dim> struct CTindex {
  nLAnet<dim>* net;
  unsigned int I, it, rank;
  CTindex(nLAnet<dim>* _net, unsigned int _I, int r): net(_net), I(_I), it(0), rank(r) {}
};
template <int dim> inline bool has_cur(CTindex<dim>& ind);
template <int dim> inline bool has_next(CTindex<dim>& ind);
template <> inline bool has_cur<3>(CTindex<3>& ind) {
  if(ind.net->indL==0) return true;
  const unsigned int M=0x49249249;
  if(ind.I&M    && ind.net->indL&1) return false;
  if(ind.I&M<<1 && ind.net->indL&2) return false;
  if(ind.I&M<<2 && ind.net->indL&4) return false;
  return true;
}
template <> inline bool has_next<3>(CTindex<3>& ind) {
  const unsigned int M=0x49249249&((1<<3*ind.rank)-1);
  ind.it++; if(ind.it >= (1<<ind.rank)) return false;
           int IM0=ind.I&M , IS0=IM0-1, I0=IS0&M , iX =IS0<0?1:0;
  int M1=M<<1, IM1=ind.I&M1, IS1=IM1-2, I1=IS1&M1; iX+=IS1<0?2:0;
  int M2=M<<2, IM2=ind.I&M2, IS2=IM2-4, I2=IS2&M2; iX+=IS2<0?4:0;
  ind.I = I0|I1|I2;
  if(iX) {
    ind.net = ind.net->netX[iX];
    if(ind.net == 0) return false;
  }
  return true;
}
template <> inline bool has_cur<2>(CTindex<2>& ind) {
  if(ind.net->indL==0) return true;
  const unsigned int M=0x55555555;
  if(ind.I&M    && ind.net->indL&1) return false;
  if(ind.I&M<<1 && ind.net->indL&2) return false;
  return true;
}
template <> inline bool has_next<2>(CTindex<2>& ind) {
  const unsigned int M=0x55555555&((1<<2*ind.rank)-1);
  ind.it++; if(ind.it >= (1<<ind.rank)) return false;
           int IM0=ind.I&M , IS0=IM0-1, I0=IS0&M , iX =IS0<0?1:0;
  int M1=M<<1, IM1=ind.I&M1, IS1=IM1-2, I1=IS1&M1; iX+=IS1<0?2:0;
  ind.I = I0|I1;
  if(iX) {
    ind.net = ind.net->netX[iX];
    if(ind.net == 0) return false;
  }
  return true;
}

#include <sched.h>

//2**(dim*(rankH-rankJump)) ConeTorre, сложенные ConeFold ранга rankCT высотой 2**rankH
//============================
template <int dim, int rank> struct ConeTorre: public ActionThreadStruct {
  SemCube<dim,rank>* semC;
  nLAnet<dim>* base_net;
  void* CTact() {
    //set_chip_affinity(nchip);
    //set_core_affinity();
    //print_CPU_affinity();
//cpu_set_t cpu_set;
//for(int i=0; i<6; i++) {
//int ic=(i*2+i/3)%6;
//if((iThr%6)/2 == i/2) CPU_SET(ic, &cpu_set);
//else CPU_CLR(ic, &cpu_set);
//}
//sched_setaffinity(0,6,&cpu_set);
    for(unsigned int I=iThr; I<1<<dim*rank; I+=nThr) {
      CTindex<dim> ind(base_net, I, rank);
      do {
        semC->wait(I);
        if(has_cur(ind)) ConeFold<rank>(ind.net, ind.I);
        semC->post(I);
      } while(has_next(ind));
      for(int it=ind.it; it<1<<rank; it++) { semC->wait(I); semC->post(I); }
    } return this;
  }
};

#include <pthread.h>
//============================
template <int dim, int rank> void TorreFold(nLAnet<dim>* net, int nCores) {
  if(!net->check4grow()) { printf("net is not checked 4 grow\n"); return; }
  const int Nthr=nCores;
  void* (*act_func)(void*)=(void* (*)(void*))&ConeTorre<dim,rank>::CTact;
  SemCube<dim,rank> sems; sems.init();
  //ConeTorre<dim,rank> CTs[Nthr];
  ConeTorre<dim,rank>* CTs=new ConeTorre<dim,rank>[Nthr];
  for(int ithr=0; ithr<Nthr; ithr++) {
    CTs[ithr].prepare(ithr, Nthr);
    CTs[ithr].semC = &sems;
    CTs[ithr].base_net = net;
    pthread_create(&CTs[ithr].ThreadID, NULL, act_func, CTs+ithr);
  }
  void* res=0;
  for(int ithr=0; ithr<Nthr; ithr++) pthread_join(CTs[ithr].ThreadID,&res);
  delete [] CTs;
  sems.destroy();
  net->finish_grow();
}

#endif//nLA_DATA_HPP
