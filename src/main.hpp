#include <stdio.h>
#include "omp.h"
#include <malloc.h>
#include "params.hpp"

#include "cubeLR.hpp"

#include "Coff.hpp"
#include "r0strct.hpp"

#include "init.hpp"

extern parsMxw pars;
extern PMLcoffs* KpmlX;
extern PMLcoffs* KpmlY;
extern PMLcoffs* KpmlZ;
extern vector<CoffStruct> CoffArr;
extern vector<DispStruct> DispArr;
extern vector<MatStruct> MatArr;
extern map<keyVolStr, Vols> VolumeMap;
extern map<MatStruct, int> MaterialMap;
extern map<DispStruct, int> DispMaterialMap;
extern long long int INITIALIZED;   // progress if initializing
extern long long int BOUNDARY_CELLS;   // the size of boundary area

//voidusetCoffs(int ind, MainType* eps, MainType* mu, MainType* Deps, MainType* Dmu, MainType* gammaE, MainType* gammaH); 
void setCoffsIso(int ind, double eps=1, double mu=1);
void setCoffsAniso(int ind, float* eps, float* mu);
void setCoffsTensorInv(int ind, double* eps, double* mu);
const vector<LorenzDisp> vL0;
void setCoffsDisp(int ind, MainType eps_inf=1, vector<LorenzDisp> LorenzD=vL0, MainType mu=1, MainType sigma=0);

void BBB(void* pFldData, void* pFldSide[6], void* pFldEdge[12], void* pFldCorn[8]);

int InitializeData();
int InitializeCoffs();
int Run();
int UpdateOneStep(int);
int DropAllData(int);
int DeleteAllData();
int Main();
