#include <set>

inline bool isNearFieldReg(indx& Ind, int normalDir[3]){
return false;
/*  normalDir[0]=0; normalDir[1]=0; normalDir[2]=0; bool retval=false; 
  if(Ind.x>pars.NearXp || Ind.x<pars.NearXm || Ind.y>pars.NearYp || Ind.y<pars.NearYm || Ind.z>pars.NearZp || Ind.z<pars.NearZm) { return false; }
  if(Ind.x==pars.NearXm) { normalDir[0]=-1; return true; }
  if(Ind.y==pars.NearYm) { normalDir[0]=-2; return true; }
  if(Ind.z==pars.NearZm) { normalDir[0]=-3; return true; }
  if(Ind.x==pars.NearXp) { normalDir[0]= 1; retval=true; }
  if(Ind.y==pars.NearYp) { normalDir[1]= 2; retval=true; }
  if(Ind.z==pars.NearZp) { normalDir[2]= 3; retval=true; }
  return retval;*/
}
extern aiv::array<float,3> PPoutXm[][6]; extern aiv::array<float,3> dPPoutXm[][6]; 
extern aiv::array<float,3> PPoutXp[][6]; extern aiv::array<float,3> dPPoutXp[][6];
extern aiv::array<float,3> PPoutYm[][6]; extern aiv::array<float,3> dPPoutYm[][6]; 
extern aiv::array<float,3> PPoutYp[][6]; extern aiv::array<float,3> dPPoutYp[][6];
extern aiv::array<float,3> PPoutZm[][6]; extern aiv::array<float,3> dPPoutZm[][6]; 
extern aiv::array<float,3> PPoutZp[][6]; extern aiv::array<float,3> dPPoutZp[][6];

extern set<int> OutsDropped;
extern set<int> OutsLoaded;

extern pthread_mutex_t NtF_mutex;

inline void SaveNearField(fieldsDDD& F0, fieldsDDD& Fx, fieldsDDD& Fy, fieldsDDD& Fxy, fieldsDDD& Fz, fieldsDDD& Fxz, fieldsDDD& Fyz, fieldsDDD& Fxyz, indx& Ind){ return;
/*  double ttstart=omp_get_wtime();
  // Ind is taken from F0[ixddd+iyddd+izddd];
  double norm[3]; int normalDir[3];
  if( !isNearFieldReg(Ind, normalDir) ) { return; }

//  if(Ind.x==63 && Ind.y==31) printf("Ok, x,y,z=%d,%d,%d time=%d, norm=%g,%g,%g\n",Ind.x,Ind.y,Ind.z,Ind.time, norm[0],norm[1],norm[2]);
  pthread_mutex_lock(&NtF_mutex);
  for(int nnn=0;nnn<3;nnn++){
    if(0==normalDir[nnn]) continue;
    norm[0]=0.; norm[1]=0.; norm[2]=0.;
    norm[abs(normalDir[nnn])-1]=normalDir[nnn]/abs(normalDir[nnn]);
    MainType P[6];
    MainType dP[6];
    fieldsDDD Fm,Fp; 
    double ddx=1.0/pars.dx;
    double ddy=1.0/pars.dy;
    double ddz=1.0/pars.dz;
    if(1==abs(normalDir[nnn])) { 
       P[0] =  F0.Ei[0]                  ;  P[1] = 0.25*   ( F0.Di[1]+F0.Ei[1]+Fx.Di[1]+Fx.Ei[1]);  P[2] = 0.25*   ( F0.Di[2]+F0.Ei[2]+Fx.Di[2]+Fx.Ei[2]);
      dP[0] = (Fx.Ei[0]-F0.Di[0])*0.5*ddx; dP[1] = 0.5*ddx*(-F0.Di[1]-F0.Ei[1]+Fx.Di[1]+Fx.Ei[1]); dP[2] = 0.5*ddx*(-F0.Di[2]-F0.Ei[2]+Fx.Di[2]+Fx.Ei[2]);

       P[3] = 0.25*   ( F0.avHx[1]+Fy.avHx[1]+Fx.avHx[1]+Fxy.avHx[1]); 
      dP[3] = 0.5*ddx*(-F0.avHx[1]-Fy.avHx[1]+Fx.avHx[1]+Fxy.avHx[1]); 
       P[4] = 0.5*(Fx.Hi[1]+Fxz.Hi[1]); dP[4] = 0.5*ddx*(-F0.avHy[2]-Fz.avHy[2]+Fx.avHy[2]+Fxz.avHy[2]); 
       P[5] = 0.5*(Fx.Hi[2]+Fxy.Hi[2]); dP[5] = 0.5*ddx*(-F0.avHz[1]-Fy.avHz[1]+Fx.avHz[1]+Fxy.avHz[1]);
    }
    if(2==abs(normalDir[nnn])) { 
       P[1] =  F0.Ei[1]                  ;  P[2] = 0.25*   ( F0.Di[2]+F0.Ei[2]+Fy.Di[2]+Fy.Ei[2]);  P[0] = 0.25*   ( F0.Di[0]+F0.Ei[0]+Fy.Di[0]+Fy.Ei[0]);
      dP[1] = (Fy.Ei[1]-F0.Di[1])*0.5*ddy; dP[2] = 0.5*ddy*(-F0.Di[2]-F0.Ei[2]+Fy.Di[2]+Fy.Ei[2]); dP[0] = 0.5*ddy*(-F0.Di[0]-F0.Ei[0]+Fy.Di[0]+Fy.Ei[0]);

       P[4] = 0.25*   ( F0.avHy[2]+Fz.avHy[2]+Fy.avHy[2]+Fyz.avHy[2]); 
      dP[4] = 0.5*ddy*(-F0.avHy[2]-Fz.avHy[2]+Fy.avHy[2]+Fyz.avHy[2]); 
       P[5] = 0.5*(Fy.Hi[2]+Fxy.Hi[2]); dP[5] = 0.5*ddy*(-F0.avHz[0]-Fx.avHz[0]+Fy.avHz[0]+Fxy.avHz[0]); 
       P[3] = 0.5*(Fy.Hi[0]+Fyz.Hi[0]); dP[3] = 0.5*ddy*(-F0.avHx[2]-Fz.avHx[2]+Fy.avHx[2]+Fyz.avHx[2]);
    }
    if(3==abs(normalDir[nnn])) { 
       P[2] =  F0.Ei[2]                  ;  P[0] = 0.25*   ( F0.Di[0]+F0.Ei[0]+Fz.Di[0]+Fz.Ei[0]);  P[1] = 0.25*   ( F0.Di[1]+F0.Ei[1]+Fz.Di[1]+Fz.Ei[1]);
      dP[2] = (Fz.Ei[2]-F0.Di[2])*0.5*ddz; dP[0] = 0.5*ddz*(-F0.Di[0]-F0.Ei[0]+Fz.Di[0]+Fz.Ei[0]); dP[1] = 0.5*ddz*(-F0.Di[1]-F0.Ei[1]+Fz.Di[1]+Fz.Ei[1]);

       P[5] = 0.25*   ( F0.avHz[0]+Fx.avHz[0]+Fz.avHz[0]+Fxz.avHz[0]); 
      dP[5] = 0.5*ddz*(-F0.avHz[0]-Fx.avHz[0]+Fz.avHz[0]+Fxz.avHz[0]); 
       P[3] = 0.5*(Fz.Hi[0]+Fyz.Hi[0]); dP[3] = 0.5*ddz*(-F0.avHx[1]-Fy.avHx[1]+Fz.avHx[1]+Fyz.avHx[1]); 
       P[4] = 0.5*(Fz.Hi[1]+Fxz.Hi[1]); dP[4] = 0.5*ddz*(-F0.avHy[0]-Fx.avHy[0]+Fz.avHy[0]+Fxz.avHy[0]);
    }
    if(normalDir[nnn]<0) for(int ip=0;ip<6;ip++) dP[ip]*=-1;

    if(OutsLoaded.find(Ind.time/(1<<MaxRank))==OutsLoaded.end()) loadOuts(Ind.time/(1<<MaxRank));
    aiv::indx<3> tmpind(0); tmpind[0]=Ind.y; tmpind[1]=Ind.z; tmpind[2]=Ind.time%(1<<MaxRank);
    for(int c=0; c<6; c++) {
      if(-1==normalDir[nnn]) { tmpind[0]=Ind.y; tmpind[1]=Ind.z; PPoutXm[Ind.time/(1<<MaxRank)][c][tmpind] = P[c]; dPPoutXm[Ind.time/(1<<MaxRank)][c][tmpind] = dP[c]; }
      if( 1==normalDir[nnn]) { tmpind[0]=Ind.y; tmpind[1]=Ind.z; PPoutXp[Ind.time/(1<<MaxRank)][c][tmpind] = P[c]; dPPoutXp[Ind.time/(1<<MaxRank)][c][tmpind] = dP[c]; }
      if(-2==normalDir[nnn]) { tmpind[0]=Ind.x; tmpind[1]=Ind.z; PPoutYm[Ind.time/(1<<MaxRank)][c][tmpind] = P[c]; dPPoutYm[Ind.time/(1<<MaxRank)][c][tmpind] = dP[c]; }
      if( 2==normalDir[nnn]) { tmpind[0]=Ind.x; tmpind[1]=Ind.z; PPoutYp[Ind.time/(1<<MaxRank)][c][tmpind] = P[c]; dPPoutYp[Ind.time/(1<<MaxRank)][c][tmpind] = dP[c]; }
      if(-3==normalDir[nnn]) { tmpind[0]=Ind.x; tmpind[1]=Ind.z; PPoutZm[Ind.time/(1<<MaxRank)][c][tmpind] = P[c]; dPPoutZm[Ind.time/(1<<MaxRank)][c][tmpind] = dP[c]; }
      if( 3==normalDir[nnn]) { tmpind[0]=Ind.x; tmpind[1]=Ind.y; PPoutZp[Ind.time/(1<<MaxRank)][c][tmpind] = P[c]; dPPoutZp[Ind.time/(1<<MaxRank)][c][tmpind] = dP[c]; }
    }
//  if(norm[2]>0. && Ind.x==63 && Ind.y==31) printf("P0=%g, x,y,z=%d,%d,%d time=%d\n",P[0],Ind.x,Ind.y,Ind.z,Ind.time);

  }
  pthread_mutex_unlock(&NtF_mutex);
//  printf("time4count one point %gs\n",omp_get_wtime()-ttstart);
*/
}
