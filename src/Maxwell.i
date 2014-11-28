%module Maxwell

%include "carrays.i"
%include "cpointer.i"

%{
#include "cubeLR.hpp"
#include "params.hpp"
#include "r0strct.hpp"
#include "Coff.hpp"
#include "main.hpp"
#include "init.hpp"
#include "Data_def.hpp"
#include "materials.hpp"
#include "Signal.hpp"
%}

%feature("kwargs");
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%apply const std::string& {std::string* swapDir};
%apply const std::string& {std::string* DropPath};
%apply const std::string& {std::string* LogPath};
%array_class(double, Array)
namespace std {
  %template(vectorOfCoffStruct) vector<CoffStruct>;
  %template(vectorOfDispStruct) vector<DispStruct>;
  %template(vectorOfMatStruct) vector<MatStruct>;
  %template(mapOfMatstruct) map<MatStruct, int>;
  %template(mapOfVolumes) map<keyVolStr, Vols>;
  %template(Components) vector<string>;
  %template(Dpoints) vector<coordinates>;
  %template(Line) vector<int>;
  %template(Square) vector< vector <int> >;
  %template(mapArray) vector< vector< vector <int> > >;
}
struct LorenzDisp {
  double a0,a1,b0,b1,b2;
  LorenzDisp();
  void set(double Deps, double wp, double gp, double gp_=0.); 
  void setDrude(double Deps, double wp, double gp, double gp_=0.);
};
namespace std {
  %template(DispParams) vector<LorenzDisp>;
}
%include "params.hpp"
%include "init.hpp"
%include "cubeLR.hpp"
//%include "pc3d.hpp"
//%include "r0strct.hpp"
%include "Coff.hpp"
%include "main.hpp"
%include "Data_def.hpp"
%include "materials.hpp"
%include "Signal.hpp"
