#ifndef nLA_NET_HPP
#define nLA_NET_HPP

#include "cubeLR.hpp"
#include "nLAnet.hpp"

//==================
template <int dim> struct MultiArray {};

#include <stdarg.h>

#define ARGS2ARR(Tint, i, i0) Tint i[dim]; i[0] = i0; {\
  va_list inds; va_start(inds, i0);\
  for(int s=1; s<dim; s++) i[s] = va_arg(inds, Tint);\
  va_end(inds);}
//----------------------------
template <int dim> struct RegionArray {
  int NcBr[dim], netSh[dim+1];
  int Nnet, Nbricks;
  nLAnet<dim>** netArr;
 public:
  RegionArray(int d): netArr(0) {}
  //RegionArray(int d[dim]) {}
  int inet(int s, int netID) { return (netID/netSh[s])%(NcBr[s]+1); }
  int get_ind(int i0, ...) { ARGS2ARR(int,i,i0); return get_ind(i); }
  int get_ind(int i[dim]) { int ind=0; for(int s=0; s<dim; s++) ind += i[s]*netSh[s]; return ind; }
  void resetRegionSize(int i0, ...) { ARGS2ARR(int,i,i0); resetRegionSize(i); }
  void resetRegionSize(int i[dim]) {
    delete[] netArr;
    Nnet = Nbricks = 1; netSh[dim] = 0;
    for(int s=0; s<dim; s++) { NcBr[s] = i[s]; netSh[s] = Nnet; netSh[dim] += Nnet; Nbricks *= NcBr[s]; Nnet *= NcBr[s]+1; }
    netArr = new nLAnet<dim>*[Nnet];
    for(int netID=0; netID<Nnet; netID++) netArr[netID] = 0;
  }
  bool isL(int s, int netID) { return inet(s,netID)==0; }
  bool isR(int s, int netID) { return inet(s,netID)==NcBr[s]; }
};
#undef ARGS2ARR
template<int dim> RegionArray<dim>* getModelRegion();

//----------------------------
template <int dim, class Tnet> struct netArray {
  static RegionArray<dim> region;
  int Nnet;
 //private:
  Tnet* netArr;
  Tnet* getLnet(int s, int in) { return (Tnet*)region.netArr[netArr[in].ID-region.netSh[s]]; }
  Tnet* getRnet(int s, int in) { return (Tnet*)region.netArr[netArr[in].ID+region.netSh[s]]; }
 public:
  netArray(): netArr(0), Nnet(-1) {}
  //Tnet* get_net(int i) { if(i>=0 && i<Nnet) return &netArr[i]; else return 0; }
  MultiArray<dim>* reset_datas4multiArray();
  int inet(int s, int in) { return (0<=in && in<Nnet)?region.inet(s, netArr[in].ID):-1; }
  void reset_datas4RegionSize();
  void set_net() { netArr=new Tnet[Nnet]; for(int in=0; in<Nnet; in++) netArr[in].clear(); }
  void setID(int in, int id) { netArr[in].ID = id; }
  void net2region() { for(int in=0; in<Nnet; in++) region.netArr[netArr[in].ID] = &netArr[in]; }
  void set_link_net();
  void set_linkXnet(const int iB=1);
  void set_datas();
  void init();
  void run_step();
  void runTorreFold(int nCores);
  void runChessFold(int nCores, int nDrops=1);
  //void print_net() { for(int in=0; in<Nnet; in++) print(net, in); }
  void clear();
};

#ifdef DIM_3D
//================== 3-dim specific
void set_data(nLAnet<3>& net, const int inode=-1);
void print(RegionArray<3>& reg, int i);

#endif //DIM_3D

#ifdef DIM_2D
//===================2-dim specific
template <> struct MultiArray<2> {
  int shift[2];
  cubeLR<2,dataDD,MaxRank>* DD;
  cubeLR<1,dataRD,MaxRank>* RD;
  cubeLR<1,dataDR,MaxRank>* DR;
  cubeLR<0,dataRR,MaxRank>* RR;
};

void set_data(nLAnet<2>& net, const int inode=-1);
void set_multi_data(MultiArray<2>& ma, nLAnet<2>& net, int netID, const int inode=-1);
void print(nLAnet<2>* netArr, int in);
void set_data(MultiArray<2>& ma, RegionArray<2>& net);
#endif //DIM_2D

//----------------------------
/*
template <int dim, class Tnet> struct net_iter4CF {
  netArray<dim, Tnet>* netArr;
  int ind;
 public:
  net_iter4CF(netArray<dim, Tnet>* _netArr): netArr(_netArr), ind(-1) { start(); }
  bool isOK() { return netArr && ind>=0 && ind<netArr->Nnet; }
  Tnet* start() { if(netArr) ind = netArr->Nnet-1; return cur(); }
  Tnet* next() { if(ind>=0) ind--; return cur(); }
  Tnet* cur() { if(isOK()) return netArr->get_net(ind); else return 0; }
};*/

//----------------------------
template <int dim, class Tnet> void netArray<dim, Tnet>::set_link_net() {
  //nLAnet structure, base data new
  for(int in=0; in<Nnet; in++) {
    Tnet& net=netArr[in];
    net.indL = net.indR = 0;
    bool diagL=true, diagR=true;
    for(int s=0; s<dim; s++) {
      if(inet(s,in)==0) {
        diagL = false;
        net.linkType[s]='L';
        net.indL += 1<<s;
      } else {
        net.setLinkL(s, *getLnet(s, in));//Надо переделать
        if(inet(s,in)==region.NcBr[s]) {
          diagR = false;
          net.linkType[s]='R';
          net.indR += 1<<s;
        } else {
          net.setLinkR(s, *getRnet(s, in));//Надо переделать
          net.linkType[s]='C';
        }
      }
    }
    if(diagL) net.setLinkL(dim, *getLnet(dim, in));//Надо переделать
    if(diagR) net.setLinkR(dim, *getRnet(dim, in));//Надо переделать
  }
}
template <int dim, class Tnet> void netArray<dim, Tnet>::set_linkXnet(const int iB) {
  for(int in=0; in<Nnet; in++) netArr[in].set_netX(iB);
}
template <int dim, class Tnet> void netArray<dim, Tnet>::clear() {
  for(int in=0; in<Nnet; in++) {
    Tnet& net=netArr[in];
    for(int i=0; i<(1<<dim); i++) if(net.datas[i] && net.ID == net.datas[i]->netID) {
      delete net.datas[i]; net.datas[i] = 0;
    } net.clear(); }// print(net, in); }
  //printf("clear Ok\n"); for(in=0; in<Nnet; in++) print(net, in);
  delete[] netArr; netArr = 0;
}

//----------------------------
template <int dim, class Tnet> MultiArray<dim>* netArray<dim,Tnet>::reset_datas4multiArray() {
  Nnet = region.Nnet;
  set_net(); net2region();
  set_link_net();
  set_linkXnet();
  ProblemSize=0;
  MultiArray<dim>* ma=new MultiArray<dim>(); set_data(*ma, *this);
  for(int in=0; in<Nnet; in++) {
    printf("reset %s: %d\n", netArr[in].linkType, netArr[in].indR);
    set_multi_data(*ma, netArr[in]);
  }
  return ma;
}
template <int dim, class Tnet> void netArray<dim,Tnet>::reset_datas4RegionSize() {
  Nnet = region.Nnet;
  set_net(); net2region();
  set_link_net();
  set_linkXnet();
  set_datas();
}
template <int dim, class Tnet> void netArray<dim,Tnet>::set_datas() {
  ProblemSize=0;
  for(int in=0; in<Nnet; in++) set_data(netArr[in]);
}
template <class numaArr> int getNode(RegionArray<3>& r, numaArr& map, int ind)
{ return map[r.inet(2,ind)][r.inet(1,ind)][r.inet(0,ind)];  }
inline int imax(int a, int b) { return a>b?a:b; }

//----------------------------
template <int dim> struct BrickPos {
  nLAnet<dim>* net;
  int netIndex[dim], coor[dim], DataIndex, BCind[2];
  long long int BaseIndex[dim], ix[dim];
  BrickPos(RegionArray<dim>& reg, int ID): net(reg.netArr[ID]) { for(int s=0; s<dim; s++) netIndex[s] = reg.inet(s,ID); }
  BrickPos(nLAnet<dim>& rnet): net(&rnet) { for(int s=0; s<dim; s++) netIndex[s] = getModelRegion<dim>()->inet(s,net->ID); }
  void set_rank(int rank) { for(int s=0; s<dim; s++) BaseIndex[s] = netIndex[s]*(1<<rank); }
  void set_rankNsubCube(int rank, int subrank=0) {
    for(int s=0; s<dim; s++) BaseIndex[s] = netIndex[s]*(1<<(rank+subrank)) - (1<<subrank)*(BCind[0]>>s & 1);
  }
  void set_shift(int* sh, int rank=0) { for(int s=0; s<dim; s++) { ix[s] = BaseIndex[s]; if(coor[s]>=0) ix[s] += sh[coor[s]]<<rank; } }
  void set_shift(short int* sh, int rank=0) { for(int s=0; s<dim; s++) { ix[s] = BaseIndex[s]; if(coor[s]>=0) ix[s] += (int(sh[coor[s]]))<<rank; } }
  void set_shift(unsigned int* sh, int rank=0) { for(int s=0; s<dim; s++) { ix[s] = BaseIndex[s]; if(coor[s]>=0) ix[s] += sh[coor[s]]<<rank; } }
  void set_shift(unsigned short int* sh, int rank=0) { for(int s=0; s<dim; s++) { ix[s] = BaseIndex[s]; if(coor[s]>=0) ix[s] += (int(sh[coor[s]]))<<rank; } }
  void setBrick(int dataI) {
    BCind[0] = BCind[1] = 0;
    int ic=0;
    for(int s=0; s<dim; s++) {
      int ib=dataI>>s & 1;
      if(net->net[s][ib] == 0) { BCind[ib] += 1<<s; coor[s] = -1; }
      else { coor[s] = ic; ic++; }
    }
  }
};

//----------------------------
template <int dim, class Tnet> void netArray<dim,Tnet>::init() {
  for(int in=0; in<Nnet; in++) {
    Tnet& net=netArr[in]; BrickPos<dim> pos(region, net.ID);
    for(int i=0; i<(1<<dim); i++) if(net.datas[i] && net.ID == net.datas[i]->netID) {
      //printf(net.linkType);
      pos.setBrick(i);
      net.datas[i]->init(pos);
    }
  }
}

//----------------------------
template <int dim, class Tnet> void netArray<dim,Tnet>::run_step() {
  for(int in=Nnet-1; in>=0; in--) {
    Tnet& net=netArr[in];
    if(net.check4grow()) { ConeFold<0>(&net,0); net.finish_grow(); }
  }
}
template <int dim, class Tnet> void netArray<dim,Tnet>::runTorreFold(int nCores) {
  for(int in=Nnet-1; in>=0; in--) TorreFold<dim, nLArank>(netArr+in,nCores);
}

#endif//nLA_NET_HPP
