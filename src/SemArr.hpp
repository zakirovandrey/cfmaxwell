#ifndef SEM_ARRAY_HEADER
#define SEM_ARRAY_HEADER

#include <semaphore.h>
#include <errno.h>

struct SemAny {
  void error_print(const char* func) {
    printf("%s Error: ", func);
    if(errno==EINVAL) printf("Sem is not a valid semaphore"); 
    else if(errno==EOVERFLOW) printf("value exceeds SEM_VALUE_MAX"); 
    else if(errno==EINTR) printf("Semaphore was restarted due to interrupt by a signal handler");
    else if(errno==ENOSYS) printf("pshared is nonzero, but the system does not support process-shared semaphores"); 
    else if(errno==EAGAIN) printf("The operation could not be performed without blocking (i.e., the semaphore currently has the value zero)");
    else printf("unknown error %d", errno);
    printf("\n");
  }
  int getVal_any(sem_t* sem) { int rs=0; while(sem_getvalue(sem, &rs)==-1) error_print("GetVal"); return rs; }
  void wait_any(sem_t* sem) { while(sem_wait(sem)==-1) error_print("Wait"); }
  bool trywait_any(sem_t* sem) {
    if(sem_trywait(sem)==0) return true;
    if(errno != EAGAIN) error_print("TryWait");
    return false;
  }
  void post_any(sem_t* sem) { if(sem_post(sem) !=0) error_print("Post"); }
  void init_any(sem_t* sem, int val) { if(sem_init(sem, 0, val) !=0) error_print("Init"); }
  void destroy_any(sem_t* sem) { if(sem_destroy(sem) !=0) error_print("Destroy"); }
};
/*
template <int dim> struct LRindex {
};

template <> struct LRindex<3> {
  static unsigned const int mask=011111111111;
  unsigned const int index;
  inline unsigned int apply_mask(int s) { return index&(mask<<s); }
  inline unsigned int applyUmask(int s) { return index|~(mask<<s); }
  inline unsigned int inc(int s) { return (((index|~(mask<<s))+(1<<s))&(mask<<s))|(index&~(mask<<s)); }
  inline unsigned int dec(int s) { return (((index|~(mask<<s))-(1<<s))&(mask<<s))|(index&~(mask<<s)); }
  inline unsigned int cyclicLrot() { return ((index>>2)&mask)|((index<<1)&~mask); }
};
template <> struct LRindex<2> {
  static unsigned const int mask=0x55555555;
  unsigned const int index;
  inline unsigned int apply_mask(int s) { return index&(mask<<s); }
  inline unsigned int applyUmask(int s) { return index|~(mask<<s); }
  inline unsigned int inc(int s) { return (((index|~(mask<<s))+(1<<s))&(mask<<s))|(index&~(mask<<s)); }
  inline unsigned int dec(int s) { return (((index|~(mask<<s))-(1<<s))&(mask<<s))|(index&~(mask<<s)); }
  //inline unsigned int cyclicLrot() { return ((index>>1)&mask)|((index<<1)&~mask); }
};
*/

#include "cubeLR.hpp"

template <int dim, int rank> struct SemCube: public SemAny {
  static const unsigned int Ns=1<<(dim*rank), Nsem=dim*Ns;//, rMask=Ns-1;
  sem_t sems[Nsem];
  //void error_print(const char* func) {
  //  printf("SemCube<dim=%d,rank=%d> ", dim,rank);
  //  this->SemAny::error_print(func);
  //}
  //sem_t* get_sem(int i=0) { if(i>=Nsem) return 0; return sems+i; }
  void init() { for(int i=0; i<Nsem; i++) init_any(sems+i, 0); }
  void destroy() { for(int i=0; i<Nsem; i++) destroy_any(sems+i); }
  void wait(int ix) {
    int ixD=ix*dim; LRindex::knot<dim,rank>& ixLR=(LRindex::knot<dim,rank>&)ix;
    for(int s=0; s<dim; s++) { if(!ixLR.isL(s)) wait_any(sems+(ixD+s)); }
  }
  bool trywait(int ix) {
    int ixD=ix*dim; LRindex::knot<dim,rank>& ixLR=(LRindex::knot<dim,rank>&)ix;
    for(int s=0; s<dim; s++) {
      if(!ixLR.isL(s) && !trywait_any(sems+(ixD+s))) {
        for(int bs=0; bs<s; bs++) { if(!ixLR.isL(bs)) post_any(sems+(ixD+bs)); }
        return false;
      }
    } return true;
  }
  void post(int ix) {
    LRindex::knot<dim,rank>& ixLR=(LRindex::knot<dim,rank>&)ix;
    for(int s=0; s<dim; s++) { if(!ixLR.isR(s)) post_any(sems+(ixLR.inc(s)*dim+s)); }
  }
};

struct SemArray: public SemAny {
  sem_t* sems;
  int Nsem, Nx, Ny, Nz;
  int Nxy, Ndiag;
 public:
  SemArray(): sems(0), Nsem(0) {}
  ~SemArray() { destroy(); }

  //void error_print(const char* func) {
  //  printf("SemArray[%d,%d,%d] ", Nx,Ny,Nz);
  //  this->SemAny::error_print(func);
  //}
  bool init(int nx, int ny=1, int nz=1) {
    if(sems) destroy();
    Nx = nx; Ny = ny; Nz = nz;
    Nxy = Nx*Ny; Ndiag = 1+(Ny>1?Nx:0)+(Nz>1?Nxy:0);
    Nsem = Nxy*nz; sems = new sem_t[Nsem];
    for(int i=0; i<Nsem; i++) init_any(sems+i, 1);
    return true;
  }
  void destroy() {
    if(sems==0) return;
    for(int i=0; i<Nsem; i++) destroy_any(sems+i);
    delete[] sems;
    sems = 0;
  }
  bool wait(int ix, int iy=0, int iz=0) {
    //printf("wait for %d,%d,%d...\n", ix,iy,iz);
    sem_t* sem_i=sems+(iz*Nxy+iy*Nx+ix);
    if(ix>0) wait_any(sem_i);
    if(iy>0) wait_any(sem_i);
    if(iz>0) wait_any(sem_i);
    wait_any(sem_i);
    return true;
  }
  bool trywait(int ix, int iy=0, int iz=0) {
    //printf("wait for %d,%d,%d...\n", ix,iy,iz);
    sem_t* sem_i=sems+(iz*Nxy+iy*Nx+ix);
    if(!trywait_any(sem_i)) return false;
    int sem_val=0, val=0, res=sem_getvalue(sem_i, &sem_val);
    if(ix>0) val++;
    if(iy>0) val++;
    if(iz>0) val++;
    if(sem_val<val) { post_any(sem_i); return false; }
    for(int i=0; i<val; i++) wait_any(sem_i);
    return true;
  }
  bool post(int ix, int iy=0, int iz=0) {
    sem_t* sem_i=sems+(iz*Nxy+iy*Nx+ix);
    if(ix<Nx-1) post_any(sem_i+1);
    if(iy<Ny-1) post_any(sem_i+Nx);
    if(iz<Nz-1) post_any(sem_i+Nxy);
    if(ix>0 && (Ny>1 || iy>0) && (Nz>1 || iz>0)) post_any(sem_i-Ndiag);
    return true;
  }
};
#endif//SEM_ARRAY_HEADER
