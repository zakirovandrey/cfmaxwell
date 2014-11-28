#ifndef SPH_MANAGER_HPP
#define SPH_MANAGER_HPP
#include <deque>
#include <aivlib/sphereT.hpp>

struct _P_ {
  MainType _p_[6];
};

const int SphereRank=2;

extern struct SphManager {
  deque< aiv::sphere<_P_> > sph;
  double ctime;
  int N;
  int time0, timeL, timeMin;
  deque< aiv::sphere<_P_> >::iterator it;
  void drop() {
/*    char fname[256];
    string Comps[] = {"Ex","Ey","Ez","Bx","By","Bz"};
    for( deque< aiv::sphere<_P_> >::iterator sn = sph.begin()+timeMin; sn!=sph.end(); ++sn ) {
      aiv::sphere<MainType> sph4save[6]; for(int _c=0; _c<6; _c++) sph4save[_c].init(SphereRank);
      for(int id=0; id<sph4save[0].N; id++) for(int _c=0; _c<6; _c++) { sph4save[_c][id] = (*sn)[id]._p_[_c]; }
      for(int _c=0; _c<6; _c++) {
        sprintf(fname,"%s/dropSphere%s%05d.sph",DropP.DropPath.c_str(),Comps[_c].c_str(),sn-sph.begin());
        aiv::Ofile file(fname);
        sph4save[_c].dump(file);
      }
    }
    timeMin=timeL;
    printf("dropped up to timeL=%d\n",timeL);
*/  }
  void add(int timeN, int i_dir, int c, double val){
    if(timeN<time0) {
//      printf("timeN<time0: %d<%d\n",timeN,time0);
      for(int i=0; i<time0-timeN; i++) {
        aiv::sphere<_P_>* s = new aiv::sphere<_P_>(SphereRank); for(int id=0; id<s->N; id++) for(int _c=0; _c<6; _c++) (*s)[id]._p_[_c] = 0.;
        sph.push_front(*s);
      }
      sph.front()[i_dir]._p_[c]+= val;
      time0=timeN; timeMin=time0;
      printf("new min time0=%d\n",time0);
    }
    else if(timeN>=timeL) {
//      printf("timeN>=timeL: %d>=%d\n",timeN,timeL);
      for(int i=0; i<=timeN-timeL; i++) {
        aiv::sphere<_P_>* s = new aiv::sphere<_P_>(SphereRank); for(int id=0; id<s->N; id++) for(int _c=0; _c<6; _c++) (*s)[id]._p_[_c] = 0.;
        sph.push_back(*s);
      }
      sph.back()[i_dir]._p_[c]+= val;
      timeL=timeN+1;
      printf("new timeL=%d, calc_time=%gs\n",timeL,omp_get_wtime()-ctime);
      ctime=omp_get_wtime();
    }
    else { 
//      printf("timeN: %d<%d<=%d\n",time0,timeN,timeL);
      it = sph.begin()+timeN-time0;
      (*it)[i_dir]._p_[c]+= val;
      timeMin=timeN;
    }
  }
  aiv::vctr<3> get_cell_center(int i) { return sph.front().get_cell_center(i); }
  SphManager(): time0(0), timeL(1), timeMin(0) {
    aiv::sphere<_P_>* s = new aiv::sphere<_P_>(SphereRank); for(int id=0; id<s->N; id++) for(int _c=0; _c<6; _c++) (*s)[id]._p_[_c] = 0.;
    sph.push_front(*s);
    N = s->N;
  }
} SphMan;

#endif
