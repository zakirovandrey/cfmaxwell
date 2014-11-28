#ifndef ACT_THREADS_HEADER
#define ACT_THREADS_HEADER

#include "SemArr.hpp"

struct __attribute__ ((aligned(64))) ActionThreadStruct {
public:
  pthread_t ThreadID;
  int thread_num;//, acts_finished;
  int iThr, nThr;
  //int type;
  //ActionThreadStructDiagnoz Diagnoz;
  //enum { DD, XD, DX, XX } actTypes;
  void prepare(int _iThr, int _nThr) { iThr = _iThr; nThr = _nThr; }

  void bind2core();
  void bind2chip();
  void bind2node();
  void set_CPU_affinity(int nthr, int Nthr);
  bool check_CPU_affinity(int nthr, int Nthr);
  void print_CPU_affinity();
  void get_CPU_affinity(char* cpus);
};
//------------
namespace HWthreadsStruct {
  const int lim_CPUs_number=128;
  const int max_CPUs_number=48;

  const int Threads_perCore=1;
  const int Cores_perChip=6;
  const int CPUs_perSMPnode=1;
  const int NUMAnodes=8;

//#include "Source/basic.MT.inc"
};

#endif//ACT_THREADS_HEADER
