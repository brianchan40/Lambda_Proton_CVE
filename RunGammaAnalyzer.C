#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
// #include </star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/structs.h>
using namespace std;

//_________________
void RunGammaAnalyzer(const int cen = 1, const int opt_weight = 1, const int sys_err_opt = 0, const Char_t *inFile = "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/file.list", const TString JobIDName = "1234", const TString lambdatype = "lam") {
  // Next line is not needed if you are not running in a standalone mode
//   gROOT->ProcessLine("#define _VANILLA_ROOT_");
  gSystem->Load("./StRoot/particle_cxx.so");
  gSystem->Load("./StRoot/particle_all_cxx.so");
  gSystem->Load("./StRoot/particle_dau_cxx.so");
  gSystem->Load("./StRoot/event_cxx.so");
  gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/namespaces/gv_gamma_cpp.so");
//   gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/libStPicoDst.so");
  TString str;
  str = ".x Gamma_112_module.C+(";
  str += cen;
  str += ",";
  str += opt_weight;
  str += ",";
  str += sys_err_opt;
  str += ", \"";
  str += inFile;
  str += "\", \"";
  str += JobIDName;
  str += "\")";
  // std::cout<<str.Data()<<std::endl;
  gROOT->ProcessLine( str.Data() );
  // Next line should be commented if you run in a batch mode
  gROOT->ProcessLine(".!rm -f Gamma_112_module_C* ");
}
