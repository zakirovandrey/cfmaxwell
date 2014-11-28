#ifdef DIM_3D

//================== 3-dim specific
void set_data(nLAnet<3>& net, const int inode) {
  const int xd=net.indR;
  const int dim=3, iR=(1<<dim)-1;
  //net.iNUMAnode = inode;
  net.copyNBdatasL();
//  printf("%s: ",net.linkType);
#ifdef LEFT_BC_COMPLEX
  if(net.datas[0]==0) net.datas[0] = new(inode) brickPML<3,0,dataLLL,MaxRank-PMLrank>(net.ID, inode);
  if(net.datas[1]==0) {
    if(xd&1) net.datas[1] = new(inode) brickPML<3,0,dataRLL,MaxRank-PMLrank>(net.ID, inode);
    else     net.datas[1] = new(inode) brickPML<3,1,dataDLL,MaxRank-PMLrank>(net.ID, inode);
  }
  if(net.datas[2]==0) {
    if(xd&2) net.datas[2] = new(inode) brickPML<3,0,dataLRL,MaxRank-PMLrank>(net.ID, inode);
    else     net.datas[2] = new(inode) brickPML<3,1,dataLDL,MaxRank-PMLrank>(net.ID, inode);
  }
  if(net.datas[4]==0) {
    if(xd&4) net.datas[4] = new(inode) brickPML<3,0,dataLLR,MaxRank-PMLrank>(net.ID, inode);
    else     net.datas[4] = new(inode) brickPML<3,1,dataLLD,MaxRank-PMLrank>(net.ID, inode);
  }
  if(net.datas[6]==0) {
    size_t ldSz=0; switch(xd&6) {
      case 0: net.datas[6] = new(inode) brickPML<3,2,dataLDD,MaxRank-PMLrank>(net.ID, inode); break;
      case 2: net.datas[6] = new(inode) brickPML<3,1,dataLRD,MaxRank-PMLrank>(net.ID, inode); break;
      case 4: net.datas[6] = new(inode) brickPML<3,1,dataLDR,MaxRank-PMLrank>(net.ID, inode); break;
      case 6: net.datas[6] = new(inode) brickPML<3,0,dataLRR,MaxRank-PMLrank>(net.ID, inode); break;
    }
  }
  if(net.datas[5]==0) {
    size_t ldSz=0; switch(xd&5) {
      case 0: net.datas[5] = new(inode) brickPML<3,2,dataDLD,MaxRank-PMLrank>(net.ID, inode); break;
      case 1: net.datas[5] = new(inode) brickPML<3,1,dataRLD,MaxRank-PMLrank>(net.ID, inode); break;
      case 4: net.datas[5] = new(inode) brickPML<3,1,dataDLR,MaxRank-PMLrank>(net.ID, inode); break;
      case 5: net.datas[5] = new(inode) brickPML<3,0,dataRLR,MaxRank-PMLrank>(net.ID, inode); break;
    }
  }
  if(net.datas[3]==0) {
    size_t ldSz=0; switch(xd&3) {
      case 0: net.datas[3] = new(inode) brickPML<3,2,dataDDL,MaxRank-PMLrank>(net.ID, inode); break;
      case 1: net.datas[3] = new(inode) brickPML<3,1,dataRDL,MaxRank-PMLrank>(net.ID, inode); break;
      case 2: net.datas[3] = new(inode) brickPML<3,1,dataDRL,MaxRank-PMLrank>(net.ID, inode); break;
      case 3: net.datas[3] = new(inode) brickPML<3,0,dataRRL,MaxRank-PMLrank>(net.ID, inode); break;
    }
  }
  switch(xd) {
    case 0: net.datas[iR] = new(inode) brick<dim,3,dataDDD,MaxRank>(net.ID, inode); break;
    case 1: net.datas[iR] = new(inode) brickPML<dim,2,dataRDD,MaxRank-PMLrank>(net.ID, inode); break;
    case 2: net.datas[iR] = new(inode) brickPML<dim,2,dataDRD,MaxRank-PMLrank>(net.ID, inode); break;
    case 4: net.datas[iR] = new(inode) brickPML<dim,2,dataDDR,MaxRank-PMLrank>(net.ID, inode); break;
    case 6: net.datas[iR] = new(inode) brickPML<dim,1,dataDRR,MaxRank-PMLrank>(net.ID, inode); break;
    case 5: net.datas[iR] = new(inode) brickPML<dim,1,dataRDR,MaxRank-PMLrank>(net.ID, inode); break;
    case 3: net.datas[iR] = new(inode) brickPML<dim,1,dataRRD,MaxRank-PMLrank>(net.ID, inode); break;
    case 7: net.datas[iR] = new(inode) brickPML<dim,0,dataRRR,MaxRank-PMLrank>(net.ID, inode); break;
  }
#else //LEFT_BC_COMPLEX
  switch(xd) {
    case 0: net.datas[iR] = new(inode) brick<dim,3,dataDDD,MaxRank>(net.ID, inode); break;
    case 1: net.datas[iR] = new(inode) brick<dim,2,dataRDD,MaxRank>(net.ID, inode); break;
    case 2: net.datas[iR] = new(inode) brick<dim,2,dataDRD,MaxRank>(net.ID, inode); break;
    case 4: net.datas[iR] = new(inode) brick<dim,2,dataDDR,MaxRank>(net.ID, inode); break;
    case 6: net.datas[iR] = new(inode) brick<dim,1,dataDRR,MaxRank>(net.ID, inode); break;
    case 5: net.datas[iR] = new(inode) brick<dim,1,dataRDR,MaxRank>(net.ID, inode); break;
    case 3: net.datas[iR] = new(inode) brick<dim,1,dataRRD,MaxRank>(net.ID, inode); break;
    case 7: net.datas[iR] = new(inode) brick<dim,0,dataRRR,MaxRank>(net.ID, inode); break;
  }
#endif //LEFT_BC_COMPLEX
  net.datas[iR]->setXDind(xd);
}
void print(RegionArray<3>& reg, int i) {
  const int dim=3;
  printf("%3d[", i); for(int s=0; s<dim; s++) printf("%d", reg.inet(s,i));
  printf("] %s: %p, net: ", reg.netArr[i]->linkType, reg.netArr[i]);
  for(int s=0; s<dim; s++) {
    nLAnet<dim>** LR=reg.netArr[i]->net[s];
    printf("(");
    if(LR[0]==0) printf("-");
    else if(reg.isL(s,i)) printf("LLLL");
    else if(LR[0] != reg.netArr[i-reg.netSh[s]]) printf("EEEE");
    else printf("%d", reg.inet(s,i)-1);
    printf(",");
    if(LR[1]==0) printf("-");
    else if(reg.isR(s,i)) printf("RRRR");
    else if(LR[1] != reg.netArr[i+reg.netSh[s]]) printf("EEEE");
    else printf("%d", reg.inet(s,i)+1);
    printf(")");
  }
  printf("; datas: [");
  for(int id=0; id<(1<<dim); id++) {
    if(id) printf(" ");
    if(reg.netArr[i]->datas[id]==0) printf("---"); else {
      int netID=reg.netArr[i]->datas[id]->netID;
      for(int s=0; s<dim; s++) printf("%d", reg.inet(s,netID));
    }
  }
  printf("]\n");
}
#endif //DIM_3D

#ifdef DIM_2D
//===================2-dim specific
void set_data(nLAnet<2>& net, const int inode=-1) {
  const int dim=2, iR=(1<<dim)-1;
  const int xd=net.indR;
  net.iNUMAnode = inode;
  net.copyNBdatasL();
  #ifdef LEFT_BC_COMPLEX
  if(net.datas[0]==0) net.datas[0] = new(inode) brick<2,0,dataLL,MaxRank>(net.ID, inode);
  if(net.datas[1]==0) {
    if(xd&1) net.datas[1] = new(inode) brick<2,0,dataRL,MaxRank>(net.ID, inode);
    else     net.datas[1] = new(inode) brick<2,1,dataDL,MaxRank>(net.ID, inode);
  }
  if(net.datas[2]==0) {
    if(xd&2) net.datas[2] = new(inode) brick<2,0,dataLR,MaxRank>(net.ID, inode);
    else     net.datas[2] = new(inode) brick<2,1,dataLD,MaxRank>(net.ID, inode);
  }
  #endif //LEFT_BC_COMPLEX
  switch(xd) {
    case 0: net.datas[iR] = new(inode) brick<dim,2,dataDD,MaxRank>(net.ID, inode); break;
    case 1: net.datas[iR] = new(inode) brick<dim,1,dataRD,MaxRank>(net.ID, inode); break;
    case 2: net.datas[iR] = new(inode) brick<dim,1,dataDR,MaxRank>(net.ID, inode); break;
    case 3: net.datas[iR] = new(inode) brick<dim,0,dataRR,MaxRank>(net.ID, inode); break;
  }
}

void set_multi_data(MultiArray<2>& ma, nLAnet<2>& net, const int inode=-1) {
  const int dim=2, iR=(1<<dim)-1;
  const int xd=net.indR;
  net.iNUMAnode = inode;
  net.copyNBdatasL();
  switch(xd) {
    case 0: net.datas[iR] = new multi_brick<dim,2,dataDD,MaxRank>(net.ID, ma.shift, ma.DD, inode); break;
    case 1: net.datas[iR] = new multi_brick<dim,1,dataRD,MaxRank>(net.ID, ma.shift+1, ma.RD, inode); break;
    case 2: net.datas[iR] = new multi_brick<dim,1,dataDR,MaxRank>(net.ID, ma.shift, ma.DR, inode); break;
    case 3: net.datas[iR] = new multi_brick<dim,0,dataRR,MaxRank>(net.ID, ma.shift, ma.RR, inode); break;
  }
}

void print(nLAnet<2>* netArr, int in) {
  const int dim=2;
  printf("%2d %s: %p, net: ", in, netArr[in].linkType, &netArr[in]);
  int sz=int(sizeof(nLAnet<2>));
  if(sz&0xF) sz = (sz+0xF)&(~0xF);
  long long unsigned pi=(long long unsigned)netArr;
  for(int s=0; s<=dim; s++){
    long long unsigned pim=(long long unsigned)netArr[in].net[s][0], pip=(long long unsigned)netArr[in].net[s][1];
    if(pim) printf("(%2d:", (pim-pi)/sz);
    else printf("(--:");
    if(pip) printf("%2d)", (pip-pi)/sz);
    else printf("--)");
    //printf("(%p,%p)", netArr[in].net[s][0], netArr[in].net[s][1]);
  }
  printf("; datas (%d): [", netArr[in].iNUMAnode);
  //for(int s=0; s<(1<<dim); s++) printf(" %p", netArr[in].datas[s]);
  for(int i=0; i<(1<<dim); i++) {
    //if(netArr[in].datas[i]) printf(" %p", ((long long unsigned*)netArr[in].datas[i])[-1]);
    if(i) printf(" ");
    if(netArr[in].datas[i]) {
      for(int s=0; s<dim; s++) printf("%d", netArr[in].datas[i]->index[s]);
    }
    else printf("--");
  }
  printf("]\n");
}

void set_data(MultiArray<2>& ma, RegionArray<2>& net) {
  ma.DD = new cubeLR<2,dataDD,MaxRank>[net.Nbricks];
  ma.DR = new cubeLR<1,dataRD,MaxRank>[net.NcBr[0]];
  ma.RD = new cubeLR<1,dataDR,MaxRank>[net.NcBr[1]];
  ma.RR = new cubeLR<0,dataRR,MaxRank>[1];
}
#endif //DIM_2D

//============================
template <int dim, int data_dim, class Tdat, int rank> void brick<dim, data_dim, Tdat, rank>::init(BrickPos<dim>& pos) {
  Tdat* cell=(Tdat*) &data; pos.set_rank(rank);
  //printf(" Brick L/R: %d/%d, axis: %d %d %d, Base: %ld %ld %ld\n", pos.BCind[0],pos.BCind[1], pos.coor[0], pos.coor[1], pos.coor[2], pos.BaseIndex[0], pos.BaseIndex[1], pos.BaseIndex[2]);
  LRindex::zipper<data_dim, rank> zI;
  for(int i=0; i<data.getSize(); i++) {
    zI.unzip(i); pos.set_shift(zI.i);
    //printf("%d: %ld %ld %ld: ", i, pos.ix[0], pos.ix[1], pos.ix[2]);
    ::init(cell[i], pos.ix[0], pos.ix[1], pos.ix[2]);
  }
}
void my_dump_head_multilrc(FILE* file, int Num){
  int zero=0, twelve = 12, DIM = -3, sizeofT = -1*int(sizeof(MainType));
  fwrite(&zero, sizeof(int),1,file);
  fwrite(&DIM, sizeof(int),1,file);//dim =
  fwrite(&sizeofT, sizeof(int),1,file);
  fwrite(&Num, sizeof(int),1,file);
};
void my_dump_head_lrc(FILE* file, int rank, int x,int y,int z){
  fwrite(&x, sizeof(int),1,file);
  fwrite(&y, sizeof(int),1,file);
  fwrite(&z, sizeof(int),1,file);
  int zero=0, DIM = -3, sizeofT = 1*int(sizeof(MainType));
  fwrite(&zero, sizeof(int),1,file);
  fwrite(&DIM, sizeof(int),1,file);
  fwrite(&sizeofT, sizeof(int),1,file);
  fwrite(&rank, sizeof(int),1,file);
};
pthread_mutex_t drop_mutex=PTHREAD_MUTEX_INITIALIZER;

int max_time=0;
template <int dim, int data_dim, class Tdat, int rank> void brick<dim, data_dim, Tdat, rank>::drop(BrickPos<dim>& pos) { //return;
  pthread_mutex_lock(&drop_mutex);
  Tdat* cell=(Tdat*) &data; pos.set_rank(rank);
  LRindex::zipper<data_dim, rank> zI;
  zI.unzip(0); pos.set_shift(zI.i);
  char fname[256];
  if((pos.net->iStep-1)%pars.timeSep==0) {
  if((pos.net->iStep)>max_time) {
    max_time=pos.net->iStep; for (int comp=0; comp<DropP.Fields4Drop.size(); comp++) { 
    sprintf(fname,"%s/drop%s%05d.mlrc",DropP.DropPath.c_str(),DropP.Fields4Drop[comp].c_str(),pos.net->iStep);
    FILE* file; file = fopen(fname,"w");
//    aiv::dump_head( file, aiv::lrCube<MainType[27],3,0>(), MaxRank, "" );
//    int Narrs = (1<<(MaxRank-PMLrank))*(1<<(MaxRank-PMLrank))*5 + (1<<(MaxRank-PMLrank))*8 + 4 + 1; // 5 sides + 8 edges + 4 corners + 1 center
    int Nx=Bricks.Nx; int Ny=Bricks.Ny; int Nz=Bricks.Nz;
    int Narrs = (1<<(MaxRank-PMLrank))*(1<<(MaxRank-PMLrank))*(2*Nx*Ny+2*Nx*Nz+2*Ny*Nz) +
                (1<<(MaxRank-PMLrank))                       *(4*Nx+4*Ny+4*Nz) + 8 + Nx*Ny*Nz; // 6 sides + 12 edges + 8 corners + 1 center
    my_dump_head_multilrc(file,Narrs); fclose(file);  
    }
  }
  for (int comp=0; comp<DropP.Fields4Drop.size(); comp++) { 
    sprintf(fname,"%s/drop%s%05d.mlrc",DropP.DropPath.c_str(),DropP.Fields4Drop[comp].c_str(),pos.net->iStep);
    FILE* file; file = fopen(fname,"a");
    my_dump_head_lrc(file,MaxRank,pos.ix[0]+(1<<PMLrank),pos.ix[1]+(1<<PMLrank),pos.ix[2]+(1<<PMLrank));
    for(int i=0; i<data.getSize(); i++) {
      ::drop(cell[i],file,comp);
    }
    fclose(file);  
  }
  }
  pthread_mutex_unlock(&drop_mutex);
}

template <int dim, int data_dim, class Tdat, int rank> void brickPML<dim, data_dim, Tdat, rank>::init(BrickPos<dim>& pos) {
//  const int pml_dim=dim-(pos.BCind[1]>>2&1); //На дневной границе PML нет
  const int pml_dim=dim;
  //const int rank=MaxRank-PMLrank;
  Tdat* cell=(Tdat*) &data; pos.set_rankNsubCube(rank,PMLrank);
  //printf(" Brick L/R: %d/%d, axis: %d %d %d, Base: %ld %ld %ld\n", pos.BCind[0],pos.BCind[1], pos.coor[0], pos.coor[1], pos.coor[2], pos.BaseIndex[0], pos.BaseIndex[1], pos.BaseIndex[2]);
  LRindex::zipper<data_dim, rank> zI;
  for(int i=0; i<data.getSize(); i++) {
    zI.unzip(i); pos.set_shift(zI.i, PMLrank);
    if(pml_dim == dim) {
      LRindex::zipper<dim, PMLrank> zJ;
      for(int j=0; j<(1<<(dim*PMLrank)); j++) {
        zJ.unzip(j);
        //printf("%d %d -> (%d,%d,%d): ", i,j,zJ[0],zJ[1],zJ[2]);
        ::init(cell[i][j], pos.ix[0]+zJ[0], pos.ix[1]+zJ[1], pos.ix[2]+zJ[2]);
      }
    } else if(pml_dim == dim-1) {
      LRindex::zipper<dim-1, PMLrank> zJ;
      for(int j=0; j<(1<<(pml_dim*PMLrank)); j++) {
        zJ.unzip(j);
        //printf("%d %d -> (%d,%d): ", i,j,zJ[0],zJ[1]);
        ::init(cell[i][j], pos.ix[0]+zJ[0], pos.ix[1]+zJ[1], pos.ix[2]);
      }
    }
  }
}
template <int dim, int data_dim, class Tdat, int rank> void brickPML<dim, data_dim, Tdat, rank>::drop(BrickPos<dim>& pos) { //return;
  pthread_mutex_lock(&drop_mutex);
//  const int pml_dim=dim-(pos.BCind[1]>>2&1); //На дневной границе PML нет
  const int pml_dim=dim;
  //const int rank=MaxRank-PMLrank;
  Tdat* cell=(Tdat*) &data; pos.set_rankNsubCube(rank,PMLrank);
  //printf(" Brick L/R: %d/%d, axis: %d %d %d, Base: %ld %ld %ld\n", pos.BCind[0],pos.BCind[1], pos.coor[0], pos.coor[1], pos.coor[2], pos.BaseIndex[0], pos.BaseIndex[1], pos.BaseIndex[2]);
  LRindex::zipper<data_dim, rank> zI;
  char fname[256];
  if((pos.net->iStep-1)%pars.timeSep==0) {
  if((pos.net->iStep)>max_time) {
    max_time=pos.net->iStep; for (int comp=0; comp<DropP.Fields4Drop.size(); comp++) { 
    sprintf(fname,"%s/drop%s%05d.mlrc",DropP.DropPath.c_str(),DropP.Fields4Drop[comp].c_str(),pos.net->iStep);
    FILE* file; file = fopen(fname,"w");
//    aiv::dump_head( file, aiv::lrCube<MainType[27],3,0>(), MaxRank, "" );
//    int Narrs = (1<<(MaxRank-PMLrank))*(1<<(MaxRank-PMLrank))*5 + (1<<(MaxRank-PMLrank))*8 + 4 + 1; // 5 sides + 8 edges + 4 corners + 1 center
    int Nx=Bricks.Nx; int Ny=Bricks.Ny; int Nz=Bricks.Nz;
    int Narrs = (1<<(MaxRank-PMLrank))*(1<<(MaxRank-PMLrank))*(2*Nx*Ny+2*Nx*Nz+2*Ny*Nz) +
                (1<<(MaxRank-PMLrank))                       *(4*Nx+4*Ny+4*Nz) + 8 + Nx*Ny*Nz; // 6 sides + 12 edges + 8 corners + 1 center
    my_dump_head_multilrc(file,Narrs); fclose(file);  
    }
  }
  for(int i=0; i<data.getSize(); i++) {
    zI.unzip(i); pos.set_shift(zI.i, PMLrank);
    for (int comp=0; comp<DropP.Fields4Drop.size(); comp++) { 
      sprintf(fname,"%s/drop%s%05d.mlrc",DropP.DropPath.c_str(),DropP.Fields4Drop[comp].c_str(),pos.net->iStep);
      FILE* file; file = fopen(fname,"a");
//    aiv::dump_head( fileSx, aiv::lrCube<MainType[27],3,0>(), MaxRank, "" );
      if(pml_dim == dim) {
        my_dump_head_lrc(file,PMLrank,pos.ix[0]+(1<<PMLrank),pos.ix[1]+(1<<PMLrank),pos.ix[2]+(1<<PMLrank));
        for(int j=0; j<(1<<(dim*PMLrank)); j++) ::drop(cell[i][j],file,comp);
      } else if(pml_dim == dim-1) {
/*      LRindex::zipper<dim-1, PMLrank> zJ;
      for(int j=0; j<(1<<(pml_dim*PMLrank)); j++) {
        zJ.unzip(j);
        //printf("%d %d -> (%d,%d): ", i,j,zJ[0],zJ[1]);
        ::drop(cell[i][j]);
      }*/
      }
      fclose(file);
    }
  }
  }
  pthread_mutex_unlock(&drop_mutex);
}

template<int dim> RegionArray<dim>* getModelRegion() {
  return &netArray<dim,ChessCell<dim> >::region;
}

void* nLAnetArrayNUMA_init(void* t) { nLAnetArrayNUMA<3>* pt=(nLAnetArrayNUMA<3>*)t; pt->init(); return t; }
void* nLAnetArrayNUMA_prepChessFold(void* t) { nLAnetArrayNUMA<3>* pt=(nLAnetArrayNUMA<3>*)t; pt->prepChessFold(); return t; }
