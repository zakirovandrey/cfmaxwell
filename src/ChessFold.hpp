#ifndef CHESS_FOLD_HEADER
#define CHESS_FOLD_HEADER

//#include "nLAnet.hpp"

template <int dim> struct ChessCell: public nLAnet<dim> {
  int color, sem_post_num;//, intD[dim];
  int SemValue;
  sem_t lock, ready;
//  int iStep, iDGtier;
  //int get_level() { int res=iDGtier; for(int s=0; s<dim; s++) res += indD[s]; return res; }
  virtual bool check4grow() { return nLAnet<dim>::check4grow(); }
  virtual void finish_grow() {
    nLAnet<dim>::finish_grow();
    BrickPos<dim> pos(*this);
    for(int i=0; i<(1<<dim); i++) if(this->datas[i] && this->ID == this->datas[i]->netID) {
      //printf(net.linkType);
      pos.setBrick(i);
      this->datas[i]->drop(pos);
    }
    for(int s=0; s<dim; s++) if(this->net[s][this->iF]) ((ChessCell<dim>*)this->net[s][this->iF])->post();
    if(this->net[dim][this->iB]) ((ChessCell<dim>*)this->net[dim][this->iB])->post();
  }
  void prepare() {
    sem_post_num = 0; this->iDGtier = 0;
    if(this->net[dim][this->iF]) sem_post_num++;
    for(int s=0; s<dim; s++) if(this->net[s][this->iB]) sem_post_num++;
    nLAnet<dim>* cur=this;
    while(cur->net[dim][this->iB]) { cur = cur->net[dim][this->iB]; this->iDGtier += dim; };
    for(int s=0; s<dim; s++) while(cur->net[s][this->iB]) { cur = cur->net[s][this->iB]; this->iDGtier++; };
    color = this->iDGtier%(dim+1);
    SemAny sem; sem.init_any(&ready,0); sem.init_any(&lock,0);
    if(this->net[dim][this->iF]) sem.post_any(&ready);
    sem.post_any(&lock);
  }
  bool post() {
    SemAny sem; sem.wait_any(&lock);
    sem.post_any(&ready);
    sem.post_any(&lock);
  }
  bool checkNlock() {
    SemAny sem; sem.wait_any(&lock); bool res=false;
    SemValue = sem.getVal_any(&ready);
    if(SemValue == sem_post_num) {
      res = true;
      for(int i=0; i<sem_post_num; i++) sem.wait_any(&ready);
    } sem.post_any(&lock);
    return res;
  }
};
/*
template <int dim, class Tnet> void netArray<dim,Tnet>::runChessFold(int nCores, int nDrops) {
  for(bool has_acts=true; has_acts; has_acts=false) {
    
    for(int in=0; in<Nnet)
    for(it.start(); it.isOK(); it.next()) {
      TorreFold<dim, nLArank>(it.cur(),nCores);
    }
  }
}*/

/*struct ChessInd {
  short int it, ix, iy, iz;
};*/
/*
typedef unsigned long long int IndType;
template <class IndT> struct IndClump {
  static const int IndNum=128/8-2;
  IndClump<IndT>* next;
  int num, start;
  IndT IndArr[IndNum];
  ChessSquareClump(): next(0), num(0), start(0) {}
  bool addInd(IndT id) {
    if(num >= IndNum) {
      if(next==0) return false;
      return next->addInd(id);
    }
    if(num+start>=IndNum) IndArr[--start] = id;
    else IndArr[num+start] = id;
    num++; return true;
  }
};//128bytes
*/
/*
template <int dim> struct ChessBoard {
  IndClump<IndType>* ClumpArr,* ClumpFree,** ColorsArr;
  int N[dim], Ncolor, Nclump, BottomLevel;
  ChessBoard(int n[dim]);
  ~ChessBoard();
  bool iterateInd(IndType& id) {
    ChessInd<dim>& ci=(ChessInd<dim>&)id;
    for(int s=0; s<dim; s++) {
      if((++ci.indD[s])<N[s]) return true;
      ci.indD[s] = 0;
    } return false;
  }
  void Ind2level(IndT id, int lv) {
    lv %= Ncolor;
    if(ColorsArr[lv]!=0 && ColorsArr[lv]->addInd(id)) return;
    if(!ClumpFree) throw("There are not free clumps more...");
    IndClump<IndType>* clt=ClumpFree; ClumpFree = clt->next;
    clt->next = ColorsArr[lv]; ColorsArr[lv] = clt;
    ColorsArr[lv]->addInd(id);
  }
  IndClump<IndType>* free(IndClump<IndType>* cl) {
    IndClump<IndType>* clt=cl->next;
    cl->num = cl->start = 0;
    cl->next = ClumpFree; ClumpFree = cl;
    return clt;
  }
};

template <int dim> ChessBoard::ChessBoard(int n[dim]): Nclump(0), Ncolor(1) {
  ClumpArr = ClumpFree = ColorsArr = 0; BottomLevel = 0;
  for(int s=0; s<dim; s++) { N[s] = n[s]+1; Ncolor += n[s]; }
  ColorsArr = new IndClump<IndType>*[Ncolor];
  int* Squares4Col=new int[Ncolor];
  for(int c=0; c<Ncolor; c++) { ColorsArr[c] = 0; Squares4Col[c] = 0; }
  union {
    IndType m;
    ChessInd<dim> o;
  } Ind;
  do { Squares4Col[Ind.o.get_level()]++; } while(iterateInd(Ind.m));
  for(int c=0; c<Ncolor; c++) Nclump += 1 + Squares4Col[c]/ChessInd<dim>::IndNum;
  delete[] Squares4Col;
  ClumpArr = new IndClump<IndType>[Nclump *= 2];
  for(int i=0; i<Nclump; i++) free(ClumpArr+i);
  Ind.o.clear();
  try { do { Ind2level(Ind.m, Ind.o.get_level()); } while(iterateInd(Ind.m));
  } catch(char* mess) { printf("Error: %s\n", mess); }
}
template <int dim> ChessBoard::~ChessBoard() {
  delete[] ClumpArr;
  delete[] ColorsArr;
}
*/
//#include "nLAnet.hpp"
//template <int dim> struct ChessFold {
//  static const int MaxKnots=1024;
//  ChessInd<ind> knots[MaxKnots];
//};

#endif //CHESS_FOLD_HEADER
