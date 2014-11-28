#ifndef nLA_NUMA_HPP
#define nLA_NUMA_HPP

#include "nLAnet.hpp"
#include "ChessFold.hpp"
#include <unistd.h>
//============================
template <int dim> struct nodes_map {
  
};
#define DEBUG_NUMA

#ifdef DEBUG_NUMA
#include <omp.h>
extern FILE* numaLogs[8];
#endif //DEBUG_NUMA

//============================
template <int dim> struct nLAnetArrayNUMA: public ActionThreadStruct {
  int iNode;
  static int thrs4node, MaxStep;
  netArray<dim,ChessCell<dim> >* base_net;
  void* init() {
    #ifdef NUMA
    numa_run_on_node(iNode);
    #endif
    base_net->init();
    return this;
  }
  void* run_step() {
    #ifdef NUMA
    numa_run_on_node(iNode);
    #endif
    base_net->run_step();
    return this;
  }
  void* prepChessFold() {
    #ifdef NUMA
    numa_run_on_node(iNode);
    #endif
#ifdef DEBUG_NUMA
    char fname[50]; sprintf(fname,"%s/res%d.log",DropP.LogPath.c_str(),iNode);
    numaLogs[iNode] = fopen(fname,"w");
    fprintf(numaLogs[iNode], "\n#Log file for %d numa\n",iNode);fflush(numaLogs[iNode]);
    for(int in=0; in<base_net->Nnet; in++) {
      ChessCell<dim>& net=base_net->netArr[in];
      net.prepare();
      fprintf(numaLogs[iNode], "%d[%3d]: %d%d%d(iDGtier=%d, sem_post_num=%d iStep=%d\n",in,net.ID,base_net->inet(0,in),base_net->inet(1,in),base_net->inet(2,in),net.iDGtier, net.sem_post_num,net.iStep);
    }
    fprintf(numaLogs[iNode], ";;;");
    for(int in=0; in<base_net->Nnet; in++) fprintf(numaLogs[iNode], "\t%s", base_net->netArr[in].linkType);
    fprintf(numaLogs[iNode], "\n");fflush(numaLogs[iNode]);
#else //DEBUG_NUMA
    for(int in=0; in<base_net->Nnet; in++) base_net->netArr[in].prepare();
#endif //DEBUG_NUMA
    return this;
  }
  void* runChessFold() {
    #ifdef NUMA
    numa_run_on_node(iNode);
    #endif
    int net4wait=0, wait_sec=0, jobs=0;
#ifdef DEBUG_NUMA
    double starttime, endtime;
    starttime = omp_get_wtime();
    fprintf(numaLogs[iNode], "MaxStep = %d, %d Nnet; ",MaxStep,base_net->Nnet);
    do {
      int net4run=0; net4wait = 0;
      for(int in=0; in<base_net->Nnet; in++) {
        ChessCell<dim>& net=base_net->netArr[in];
        fprintf(numaLogs[iNode], "%d..",in);
        if(net.iStep+net.iDGtier/(dim+1) >= MaxStep) continue;
        fprintf(numaLogs[iNode], "[");
        for(int s=0; s<dim; s++) fprintf(numaLogs[iNode], "%d",base_net->inet(s,in));
        fprintf(numaLogs[iNode], "/");
        if(net.checkNlock()) {
          fprintf(numaLogs[iNode], "%d->",net.iStep); fflush(numaLogs[iNode]);
          double begtime=omp_get_wtime();
          net4run++; TorreFold<dim, nLArank>(&net,thrs4node);
          net.Trun = omp_get_wtime()-begtime;
        } else {
          fprintf(numaLogs[iNode], "{%d}",net.SemValue); fflush(numaLogs[iNode]);
        }
        if(net.iStep+net.iDGtier/(dim+1) < MaxStep) net4wait++;
        fprintf(numaLogs[iNode], "%d#%d]",net.iStep,net4run);
      }
      if(net4run) {
        jobs+=net4run;
        if(net4wait) fprintf(numaLogs[iNode], " %d wait,", net4wait);
        fprintf(numaLogs[iNode], " %d run\n", net4run); fflush(numaLogs[iNode]);
      }
      if(net4wait>0 && net4run==0) { sleep(1); wait_sec++; }
    } while(net4wait);
    endtime = omp_get_wtime();
    fprintf(numaLogs[iNode], "=== %d jobs; %d wait_sec; %g sec\n", jobs, wait_sec, endtime-starttime);fflush(numaLogs[iNode]);
    fprintf(numaLogs[iNode], ";;;");
    for(int in=0; in<base_net->Nnet; in++) fprintf(numaLogs[iNode], "\t%.3g", base_net->netArr[in].Trun);
    fprintf(numaLogs[iNode], "\n");fflush(numaLogs[iNode]);
#else //DEBUG_NUMA
    do {
      int net4run=0; net4wait = 0;
      for(int in=0; in<base_net->Nnet; in++) {
        ChessCell<dim>& net=base_net->netArr[in];
        if(net.iStep+net.iDGtier/(dim+1) >= MaxStep) continue;
        if(net.checkNlock()) { net4run++; TorreFold<dim, nLArank>(&net,thrs4node); }
        if(net.iStep+net.iDGtier/(dim+1) < MaxStep) net4wait++;
      }
      if(net4wait>0 && net4run==0) { sleep(1); wait_sec++; }
    } while(net4wait);
#endif //DEBUG_NUMA
    return this;
  }
};

template <int dim> void numa_run_on_nodes(void* (*act_func)(void*), netArray<dim,ChessCell<dim> >* net, int nNodes) {
  nLAnetArrayNUMA<dim> ATs[nNodes];
  for(int ithr=0; ithr<nNodes; ithr++) {
    ATs[ithr].prepare(ithr, nNodes);
    ATs[ithr].base_net = net+ithr;
    ATs[ithr].iNode = ithr;
    pthread_create(&ATs[ithr].ThreadID, NULL, act_func, ATs+ithr);
  }
  void* res=0;
  for(int ithr=0; ithr<nNodes; ithr++) pthread_join(ATs[ithr].ThreadID,&res);
}
template <int dim> void print_net(int in) {
  nLAnet<dim>& net=*getModelRegion<dim>()->netArr[in];
  printf("I am %d, %p\n",net.ID, &net);
  for(int s=0; s<=dim; s++) {
    printf("%d:", s);
    for(int d=0; d<=1; d++)
    if(net.net[s][d]) printf(" [%d]%p", net.net[s][d]->ID, net.net[s][d]);
    else printf(" [--]%p", net.net[s][d]);
    printf("\n");
  }
}
//============================
template <class numaArr> netArray<3,ChessCell<3> >* reset_numa(numaArr& mapNUMAnodes, int nNodes=1) {
  const int dim=3;
  int Nx=mapNUMAnodes[0][0].size();
  int Ny=mapNUMAnodes[0].size();
  int Nz=mapNUMAnodes.size();
  printf("Nx*Ny*Nz = %dx%dx%d\n", Nx-1,Ny-1,Nz-1);
  RegionArray<3>& reg=netArray<3,ChessCell<3> >::region;
  reg.resetRegionSize(Nx-1,Ny-1,Nz-1);
  int iNodeMax=0;
  for(int i=0; i<reg.Nnet; i++) iNodeMax = imax(getNode(reg,mapNUMAnodes, i), iNodeMax);
  if(iNodeMax+1 != nNodes) printf("!!! Error in reset_numa: Datas in %d Nodes, Max is %d\n", nNodes, iNodeMax);
  netArray<3,ChessCell<dim> >* arr=new netArray<3,ChessCell<dim> >[nNodes];
  int indNode[nNodes]; for(int in=0; in<nNodes; in++) indNode[in] = arr[in].Nnet = 0;
  for(int i=0; i<reg.Nnet; i++) arr[getNode(reg,mapNUMAnodes, i)].Nnet++;
//   arr.Nbricks = arr.region.Nbricks;
  for(int in=0; in<nNodes; in++) arr[in].set_net();
  for(int i=0; i<reg.Nnet; i++) {
    int id=getNode(reg,mapNUMAnodes, i);
    arr[id].setID(indNode[id], i);
    indNode[id]++;
  }
  //print_net<dim>(13);
  for(int in=0; in<nNodes; in++) arr[in].net2region();
  for(int in=0; in<nNodes; in++) arr[in].set_link_net();
  for(int in=0; in<nNodes; in++) arr[in].set_linkXnet();
  ProblemSize=0;
  for(int id=0; id<reg.Nnet; id++) {
    int in=getNode(reg, mapNUMAnodes, id);
    //int inet[]={reg.inet(0,i),reg.inet(1,i),reg.inet(2,i)};
    set_data(*reg.netArr[id], in);
  }
//  for(int i=0; i<reg.Nnet; i++) print(reg, i);
//  for(int in=0; in<nNodes; in++) arr[in].set_datas();
  return arr;
}

#endif//nLA_NUMA_HPP
