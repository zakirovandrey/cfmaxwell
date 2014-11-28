#ifndef INIT_HPP
#define INIT_HPP

#include <string>
#include <vector>

//#include <aivlib/mystream.hpp>

#include "Data_def.hpp"
#define DIM_3D
#define LEFT_BC_COMPLEX
#include "nLAnet.hpp"
#include "nLAarr.hpp"

using namespace std;

extern struct NUMAnodesStruct {
  int numNUMAnodes;
  int Nx, Ny, Nz;
  vector< vector< vector< int > > > mapNUMAnodes;
  NUMAnodesStruct(): numNUMAnodes(1) {
    set(1,1,1);
  }
  void set(int nx, int ny, int nz) {
    Nx = nx; Ny = ny; Nz = nz;
    pars.Nx=nx; pars.Ny=ny; pars.Nz=nz;
    mapNUMAnodes.resize(nz+1);
    for(int iz=0; iz<nz+1; iz++) {
      mapNUMAnodes[iz].resize(ny+1);
      for(int iy=0; iy<ny+1; iy++) {
        mapNUMAnodes[iz][iy].resize(nx+1);
        for(int ix=0; ix<nx+1; ix++) mapNUMAnodes[iz][iy][ix]=0;
      }
    }
  }
} Bricks;

struct coordinates{
  int x,y,z;
  coordinates(int _x, int _y, int _z): x(_x), y(_y), z(_z) {}
  coordinates() {}
  const bool operator==(const coordinates& r) {return x==r.x && y==r.y && z==r.z; }
};
extern struct DropParams {
  string LogPath;
  vector<string> Fields4Drop;
  vector<string> Fields4Sensors;
  vector<coordinates> Sensors;
  string DropPath;
  int timesep;
  DropParams(): LogPath(""), DropPath(""), timesep(1) {}
  void setSensor(int sx, int sy, int sz) {
    int s = Sensors.size();
    Sensors.resize(s+1);
    Sensors.back().x=sx; Sensors.back().y=sy; Sensors.back().z=sz;
  }
} DropP;

DropParams* getDropP();
NUMAnodesStruct* getBricks();

#include "nLAnuma.hpp"

extern netArray<3,ChessCell<3> >* nLAarr;

void init(fieldsDDD& F, int x, int y, int z);
void init(fieldsSDD& F, int x, int y, int z);
void init(fieldsDSD& F, int x, int y, int z);
void init(fieldsDDS& F, int x, int y, int z);
void init(fieldsDSS& F, int x, int y, int z);
void init(fieldsSDS& F, int x, int y, int z);
void init(fieldsSSD& F, int x, int y, int z);
void init(fieldsSSS& F, int x, int y, int z);
void init(fieldsDDX& F, int x, int y, int z);
void init(fieldsSDX& F, int x, int y, int z);
void init(fieldsDSX& F, int x, int y, int z);
void init(fieldsSSX& F, int x, int y, int z);

bool isSourceInside(indx& Ind);

int loadOuts(int);
template<class T> void drop(T& F, FILE*, int);
void dropPP();

#endif //INIT_HPP
