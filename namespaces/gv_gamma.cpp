/// C++ headers
#include <iostream>

/// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
//#include "TVector3.hh"
#include <assert.h>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include "TVector3.h"

/// PicoDst headers
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoDstReader.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoDst.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoEvent.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoTrack.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTofHit.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTowHit.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoHelix.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoPhysicalHelix.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTofPidTraits.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/StRefMultCorr.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/CentralityMaker.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/particle.h"
// #include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/event.h"

/// Analysis headers
//#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StV0Type.h"
//#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StDcaService.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/structs.h"
#include "./gv_gamma.h"

namespace gv_gamma{
	std::vector<int> iTrack;

    TVector3 *trk_mom = new TVector3();
    TVector3 *trk_mom2 = new TVector3();
    TVector3 *trk_mom_temp = new TVector3();
	
}
