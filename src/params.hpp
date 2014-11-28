#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <math.h>
#include <string>

using namespace std;
//#define USE_SWAP

//--------typedefs------------------------------//
#define MAINTYPEisFLOAT
typedef float MainType;
//----------------------------------------------//

//-----------------Grid----------------------//
//==============MaxRank-nLArank > PMLrank===============//
const int MaxRank = 5;
const int nLArank = 1;
const int PMLrank = 3;
const int PMLN=1<<PMLrank;
const int PMLNx=PMLN;
const int PMLNy=PMLN;
const int PMLNz=PMLN;
const int dim=3;
//-------------------------------------------//

struct parsMxw{
//-----------------Numeric Params------------------------//
  double dx,dy,dz;
  double dt;
  int Nx, Ny, Nz;
  int Tcount,timeSep;
  void setGrid(double dx=0.02, double dy=0.02, double dz=0.02, double dt=0.01);
  parsMxw() { timeSep=1; }
  //---------------Near_to_far_field_tarnsformation params-----------//
  double Cen[3];
  int NearXm,NearXp, NearYm,NearYp, NearZm,NearZp;

  //---------------Counting Area params----------------------//
  int MaxNodes; 
  int Nthreads;
  int BigSteps;

//------------PML params-----------------------------------------------------------//
//      dv/dt + gamma   * v = k_1*du/dx        ////  u is S or T 
//      du/dt + gamma^* * u = k_2*dv/dx        ////               gamma задаются не как в статье!!!
//      (gamma^*)/gamma = 1
//float g_max = 6.* (3.*cSrcP/(2.*PMLN*(accuracyZ-1)*dz));  // при r.*(...) затухание в e^r раз.
  double sigma_max;
  double attenuation_factor;

//--------------Coeffisients params------------------------------------------------//
  double mapResolution;
  int deep1d;
  int deep3d;
  MainType averRx, averRy, averRz;
  bool subpixel;

//-----------------System params----------------------------------------------------//
  double CPUfreq;
  int Cores;
  int nchip;

  string swapDir;

};
extern parsMxw pars;
parsMxw* getPars();

inline MainType gamma_func(int i, int Nn) { 
    MainType x_max=pow(pars.sigma_max,(1./pars.attenuation_factor));
//    int _i=(i<1)?0:(i-1);
    MainType x=x_max-i*(x_max/Nn); // 0<=i<=PMLN 
//    if ( i/(2*(accuracyZ-1))==PMLN-1 || i/(2*(accuracyZ-1))==PMLN ) return 0.;    // left boundary main-PML in first cell is gamma=0 
//    if ( i/(2*(accuracyZ-1))==PMLN ) return 0.;    // left boundary main-PML in first cell is gamma=0 
    return pow(x,(pars.attenuation_factor));
}   
inline MainType gamma_funcX(int i, int Nn) { return gamma_func(i,Nn); }
inline MainType gamma_funcY(int i, int Nn) { return gamma_func(i,Nn); }
inline MainType gamma_funcZ(int i, int Nn) { return gamma_func(i,Nn); }
//const int k_gammas_p = 1./pow(Rho*cSrcP,2);   //  (gamma*)/gamma
//const int k_gammas_s = 1./pow(Rho*cSrcS,2);   //
const int k_gammas = 1.;   //  (gamma*)/gamma
inline MainType gamma_funcX_E(int i, int Nn) { return gamma_funcX(i, Nn);            }
inline MainType gamma_funcX_H(int i, int Nn) { return k_gammas*gamma_funcX_E(i, Nn); }
inline MainType gamma_funcY_E(int i, int Nn) { return gamma_funcY(i, Nn);            }
inline MainType gamma_funcY_H(int i, int Nn) { return k_gammas*gamma_funcY_E(i, Nn); }
inline MainType gamma_funcZ_E(int i, int Nn) { return gamma_funcZ(i, Nn);            }
inline MainType gamma_funcZ_H(int i, int Nn) { return k_gammas*gamma_funcZ_E(i, Nn); }
//---------------------------------------------------------------------------------//

//----------------------------reflect class for counting PML reflection------------------------------------------------//
/*
struct reflect {
  public:
    double deltax,deltat,Courant,k1,k2,speed;           // k1 и k2 --- коэффициенты в уравнениях, speed --- скорость возмущений, speed=1/sqrt(k1*k2)
                                                // в максвелле k1=eps/c, k2=mu/c
                                                // в упругости k1=rho, k2=1/(rho*V*V)
    double omega,k,N,Rth,R,phaseI,phaseR;
    complex<double> M;
    reflect(double dx_, double dt_, double k1_, double k2_): deltax(dx_), deltat(dt_), k1(k1_), k2(k2_) {
        speed = 1./sqrt(k1*k2);
        Courant=dx_/(speed*dt_);
        if(Courant<1.) { perror("Incorrect Courant number > 1.0"); exit(-1); }
    }
    void count_R();
    void count_Rth();
    double getgamma(int);    // gamma function
//  private:
    double A(int);
    double B(int);
    void count_dispersion();
    void count_phase();
    void count_Gform();
};
//sqrt(k1/k2)=Vp
#ifndef SWIG
const complex<double> j(0.,1.);
#endif
inline void reflect::count_Rth() {
    double factor = 0.; for(int i=0;i<=2*N;i++) factor+= 2.*sqrt(k2/k1)*getgamma(i)*(deltax*0.5);
    Rth = exp(-factor);
}
inline double reflect::getgamma(int i) { if(i>=2*N) return 0.; return gamma_func(i,2*N); }
//inline double reflect::getgamma(int i) { if(i>=2*N) return 0.; return gamma_func(i/2,3*PMLN); }
inline double reflect::A(int i){ return (2.-getgamma(i)*deltat)/(2.+getgamma(i)*deltat); }
inline double reflect::B(int i){ return 2./(2.+getgamma(i)*deltat);          } 
inline void reflect::count_dispersion() { k=2./deltax*asin(Courant*sin(omega*deltat*0.5)); }
inline void reflect::count_phase() {}
inline void reflect::count_Gform() {
    complex<double> Wold(1.);                  //   W_0 (on the bound)
    double k12 = (0%2==0)?(k1):(k2); complex<double> g(deltax*k12/deltat*(complex<double>((1-A(0))*cos(omega*deltat*0.5), (1+A(0))*sin(omega*deltat*0.5)))/(-B(0)));
    complex<double> W(complex<double>(-0.5)*g*Wold);
    complex<double> Wm(W);
    for(int i=1; i<2*N+1; i++){
        Wm=W;
        k12 = (i%2==0)?(k1):(k2);
        complex<double> g(deltax*k12/deltat*(complex<double>((1-A(i))*cos(omega*deltat*0.5), (1+A(i))*sin(omega*deltat*0.5)))/(-B(i)));
        W=Wold-g*W;     // W_(i+1) = W_(i-1) - g_i*W_i
        Wold=Wm;
    }
    M = W/Wm;
}
inline void reflect::count_R() {
    count_dispersion();
    count_Gform();
    R = abs((sqrt(k1/k2)*complex<double>(cos(k*deltax*0.5),sin(k*deltax*0.5))-M)/(sqrt(k1/k2)*complex<double>(cos(k*deltax*0.5),-sin(k*deltax*0.5))+M));
}
*/
#endif//PARAMS_HPP


