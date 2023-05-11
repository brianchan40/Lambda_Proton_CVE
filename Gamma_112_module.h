using namespace std;

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TF1.h"
#include "./namespaces/gv_gamma.h"

#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/particle.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/particle_all.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/particle_dau.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/event.h"

const float PI = TMath::Pi();
const float MM = 2 / PI;
const int opt_useBBC = 1;
const int Phibin = 80;
const int order = 10;
const int nHar = 2; // 2nd order EP
const int run_sta = 71000;
const int run_end = 168000;
const float pt_trig_up = 2;
const float pt_trig_lo = .2;
const float pt_asso_up = 2.0;  // 1;
const float pt_asso_lo = 0.15; // 0.15;
const float EtaCut = 1.0;
const float DcaCut = 2;
const float Eta_EP_Cut = 0.1;  // particles in the event plane
const Float_t Vz_cut = 40;     // 40 for 39 Gev, 75 for 11 GeV and 70 for 7GeV
const float Eta_ESE_Cut = 0.5; // particles for event shape engineering in using TPC EP
const int opt_useEPD = 1;      // 0 is TPC EP, 1 is EPD, 11 is 1st-order EPD, 2 is BBC, 3 is ZDC

const double Lambda_Mass = 1.115684;

// new efficiency
float efficiency[9][22], efficiency_anti[9][22], efficiency_err[9][22], efficiency_anti_err[9][22];

// Q2 correction parameter
//  const float v2_averaged[9] = {0.0673874,0.0849566,0.0785422,0.0815793,0.0779556,0.0702452,0.056361,0.0383315,0.0286762};//This is mean TPC averaged
//  const float v2_averaged_EPD[9] = {0.0673874,0.0849566,0.0785422,0.0815793,0.0779556,0.0702452,0.056361,0.0383315,0.0286762};//This is mean EPD averaged
//  const float v2_parent_averaged[9]={0.077609,0.0640669,0.0651383,0.0661118,0.0633571,0.0553656,0.0415035,0.0285145,0.0168886};//This is parent v2
const float v2_parent_averaged_rot[9] = {0.0835153, 0.0663472, 0.0648103, 0.0678791, 0.0642486, 0.0573776, 0.0439953, 0.0307474, 0.0189867}; // This is parent v2
const float v2_parent_averaged_TPC[9] = {0.0642794, 0.0515398, 0.0444409, 0.0419725, 0.03933, 0.0334889, 0.0246591, 0.0163502, 0.0105014};
const float v2_parent_averaged_EPD[9] = {0.110697, 0.0674147, 0.0408674, 0.0408868, 0.0388992, 0.0326713, 0.0245368, 0.0145434, 0.0117897};
const float v2_parent_averaged_EPD1[9] = {-0.112763, 0.0335038, 0.0330448, 0.0378166, 0.0385301, 0.0316791, 0.0239076, 0.0188769, -0.00641116};

const float v2_averaged_TPC[9] = {0.0442338, 0.0482253, 0.0532989, 0.0567928, 0.0564132, 0.0510684, 0.040288, 0.0287219, 0.019306};
const float v2_averaged_EPD[9] = {0.0326293, 0.0431983, 0.0522669, 0.0559102, 0.0551005, 0.0501335, 0.03935, 0.0275008, 0.0189943};
const float v2_averaged_EPD1[9] = {0.0150635, 0.0334534, 0.0431322, 0.0495804, 0.0509755, 0.0468642, 0.036826, 0.025391, 0.013376};

// const float v2_averaged_pairpion_TPC[9] = {0.0519615, 0.0526463, 0.0541021, 0.054624, 0.0514534, 0.0431231, 1, 1, 1};
// const float v2_averaged_pairpion_EPD[9] = {0.0519615, 0.0526463, 0.0541021, 0.054624, 0.0514534, 0.0431231, 1, 1, 1};
// const float v2_averaged_pairpion_EPD1[9] = {0.0519615, 0.0526463, 0.0541021, 0.054624, 0.0514534, 0.0431231, 1, 1, 1};

// const float v2_averaged_pairpion_TPC[9] = {0.520766, 0.526748, 0.541054, 0.546308, 0.514655, 0.431335, 0.273008, 0.254984, 0.236959};
const float v2_averaged_pairpion_TPC[9] = {0.056187, 0.060232, 0.061147, 0.061668, 0.059219, 0.052820, 0.0420832, 0.0282789, 0.0220476};
const float v2_averaged_pairpion_EPD[9] = {0.056187, 0.060232, 0.061147, 0.061668, 0.059219, 0.052820, 0.0420832, 0.0282789, 0.0220476};
const float v2_averaged_pairpion_EPD1[9] = {0.056187, 0.060232, 0.061147, 0.061668, 0.059219, 0.052820, 0.0420832, 0.0282789, 0.0220476};

// defining histograms
TH2D *Ref_TOF = new TH2D("Ref_TOF", "Ref_TOF", 500, 0.5, 500.5, 5000, 0.5, 5000.5);
TProfile *Ref_Day3 = new TProfile("Ref_Day3", "RefMult vs Run", run_end - run_sta, run_sta, run_end, 0, 999);
TProfile *TOF_Day3 = new TProfile("TOF_Day3", "TOFMult vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *NPT_Day3 = new TProfile("NPT_Day3", "NPTracks vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *NPA_Day3 = new TProfile("NPA_Day3", "NPAsso vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *TPC_Day3_cos2 = new TProfile("TPC_Day3_cos2", "cos(2*psi) vs Run", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *TPC_Day3_sin2 = new TProfile("TPC_Day3_sin2", "sin(2*psi) vs Run", run_end - run_sta, run_sta, run_end, -1, 1);

TH1D *hTally = new TH1D("hTally", "hTally", 10, 0.5, 10.5);
TH1D *hTall = new TH1D("hTall ", "hTall ", 11, 0.5, 11.5);
TH1D *hZDCcoin = new TH1D("hZDCcoin", "hZDCcoin", 1000, 0, 100000);
TH1D *hTrigger = new TH1D("hTrigger", "hTrigger", 200, 0.5, 200.5);
TH1D *hCentrality = new TH1D("hCentrality", "hCentrality", 11, -1.5, 9.5);
TH2D *hVertexXY = new TH2D("hVertexXY", "hVertexXY", 300, -3, 3, 300, -3, 3);
TH1D *hVertexZ = new TH1D("hVertexZ", "hVertexZ", 100, -100, 100);
TH1D *hVzDiff = new TH1D("hVzDiff", "hVzDiff", 120, -30, 30);
TH2D *hMult_Vz = new TH2D("hMult_Vz", "hMult_Vz", 500, 0.5, 500.5, 5000, 0.5, 5000.5);
TH2D *hMult_Vz_new = new TH2D("hMult_Vz_new", "hMult_Vz_new", 1000, -0.5, 999.5, 100, -100, 100);

TH1D *num_lam_tree = new TH1D("num_lam_tree", "Number of Lambdas in an Event from Tree", 20, 0, 20);
TH1D *num_antilam_tree = new TH1D("num_antilam_tree", "Number of AntiLambdas in an Event from Tree", 20, 0, 20);
TH1D *num_lam_rot_tree = new TH1D("num_lam_rot_tree", "Number of Lambdas (BKG) in an Event from Tree", 20, 0, 20);
TH1D *num_antilam_rot_tree = new TH1D("num_antilam_rot_tree", "Number of AntiLambdas (BKG) in an Event from Tree", 20, 0, 20);
TH1D *num_proton_tree = new TH1D("num_proton_tree", "Number of Protons in an Event from Tree", 1000, 0, 1000);
TH1D *num_antiproton_tree = new TH1D("num_antiproton_tree", "Number of AntiProtons in an Event from Tree", 1000, 0, 1000);
TH1D *num_lam_final = new TH1D("num_lam_final", "Number of Lambdas in an Event Actually Used", 20, 0, 20);
TH1D *num_proton_final = new TH1D("num_proton_final", "Number of Protons+AntiProtons in an Event Actually Used", 1000, 0, 1000);
TH1D *num_gamma_final = new TH1D("num_gamma_final", "Number of Gamma Entries in an Event", 1000, 0, 1000);

TH1D *proton_overlap_ratio = new TH1D("proton_overlap_ratio", "Number of Protons rejected from EP Calculation Due to Daughter Criteria", 100, 0, 1);

TH1D *num_proton_used = new TH1D("num_proton_used", "Number of Protons in an Event Actually Used Multiplied by Times Used", 1000, 0, 1000);
TH1D *num_antiproton_used = new TH1D("num_antiproton_used", "Number of AntiProtons in an Event Actually Used Multiplied by Times Used", 1000, 0, 1000);

TH1D *Hist_check_trksplitting_mag_ss = new TH1D("Hist_check_trksplitting_mag_ss", "Hist_check_trksplitting_mag_ss", 2000, 0, 10);
TH1D *Hist_check_trksplitting_mag_os = new TH1D("Hist_check_trksplitting_mag_os", "Hist_check_trksplitting_mag_os", 2000, 0, 10);
TH1D *Hist_check_trksplitting_mag_peak_ss = new TH1D("Hist_check_trksplitting_mag_peak_ss", "Hist_check_trksplitting_mag_peak_ss", 2000, 0, 10);
TH1D *Hist_check_trksplitting_mag_peak_os = new TH1D("Hist_check_trksplitting_mag_peak_os", "Hist_check_trksplitting_mag_peak_os", 2000, 0, 10);

TH1D *Hist_check_trksplitting_phi_ss = new TH1D("Hist_check_trksplitting_phi_ss", "Hist_check_trksplitting_phi_ss", 800, -4, 4);
TH1D *Hist_check_trksplitting_phi_os = new TH1D("Hist_check_trksplitting_phi_os", "Hist_check_trksplitting_phi_os", 800, -4, 4);

TH1D *Hist_check_trksplitting_phi_all_ss = new TH1D("Hist_check_trksplitting_phi_all_ss", "Hist_check_trksplitting_phi_all_ss", 800, -4, 4);
TH1D *Hist_check_trksplitting_phi_all_os = new TH1D("Hist_check_trksplitting_phi_all_os", "Hist_check_trksplitting_phi_all_os", 800, -4, 4);

TH1D *Hist_check_trksplitting_phi_p = new TH1D("Hist_check_trksplitting_phi_p", "Hist_check_trksplitting_phi_p", 800, -0.4, 0.4);
TH1D *Hist_check_trksplitting_phi_ap = new TH1D("Hist_check_trksplitting_phi_ap", "Hist_check_trksplitting_phi_ap", 800, -4, 4);

TH2D *Hist_check_trksplitting_pmom_ss = new TH2D("Hist_check_trksplitting_pmom_ss", "Hist_check_trksplitting_pmom_ss", 200, -0.1, 0.1, 200, -0.1, 0.1);
TH2D *Hist_check_trksplitting_pmom_os = new TH2D("Hist_check_trksplitting_pmom_os", "Hist_check_trksplitting_pmom_os", 200, -0.1, 0.1, 200, -0.1, 0.1);

TH1D *Hist_check_trksplitting_ptrkid_ss = new TH1D("Hist_check_trksplitting_ptrkid_ss", "Hist_check_trksplitting_ptrkid_ss", 400, -20, 20);
TH1D *Hist_check_trksplitting_ptrkid_os = new TH1D("Hist_check_trksplitting_ptrkid_os", "Hist_check_trksplitting_ptrkid_os", 400, -20, 20);

TH1D *Hist_check_trksplitting_ptrkid2_ss = new TH1D("Hist_check_trksplitting_ptrkid2_ss", "Hist_check_trksplitting_ptrkid2_ss", 400, -20, 20);
TH1D *Hist_check_trksplitting_ptrkid2_os = new TH1D("Hist_check_trksplitting_ptrkid2_os", "Hist_check_trksplitting_ptrkid2_os", 400, -20, 20);

TH1D *nHits_ratio_under_peak = new TH1D("nHits_ratio_under_peak", "nHits_ratio_under_peak", 300, 0, 1.5);

TH1D *Hist_phi_below1_ss = new TH1D("Hist_phi_below1_ss", "Hist_phi_below1_ss", 800, -4, 4);
TH1D *Hist_phi_after1_ss = new TH1D("Hist_phi_after1_ss", "Hist_phi_after1_ss", 800, -4, 4);
TH1D *Hist_phi_both1_ss = new TH1D("Hist_phi_both1_ss", "Hist_phi_both1_ss", 800, -4, 4);
TH1D *Hist_phi_below1_os = new TH1D("Hist_phi_below1_os", "Hist_phi_below1_os", 800, -4, 4);
TH1D *Hist_phi_after1_os = new TH1D("Hist_phi_after1_os", "Hist_phi_after1_os", 800, -4, 4);
TH1D *Hist_phi_both1_os = new TH1D("Hist_phi_both1_os", "Hist_phi_both1_os", 800, -4, 4);

TH1D *Hist_positive = new TH1D("Hist_positive", "Hist_positive", 500, -0.5, 499.5);
TH1D *Hist_negative = new TH1D("Hist_negative", "Hist_negative", 500, -0.5, 499.5);
TH1D *Hist_Ch = new TH1D("Hist_Ch", "Hist_Ch", 1000, -0.5, 999.5);
TH1D *Hist_netCh = new TH1D("Hist_netCh", "Hist_netCh", 999, -499.5, 499.5);
TH1D *Hist_netChAsym = new TH1D("Hist_netChAsym", "Hist_netChAsym", 600, -1.5 + 0.0025, 1.5 + 0.0025);
TH2D *Hist_netChAsym_bin = new TH2D("Hist_netChAsym_bin", "Hist_netChAsym_bin", 5, 0.5, 5.5, 600, -1.5 + 0.0025, 1.5 + 0.0025);
TProfile *p_netChAsym_RefMult = new TProfile("p_netChAsym_RefMult", "p_netChAsym_RefMult", 300, -1.5 + 0.0025, 1.5 + 0.0025, 0., 900);
TProfile *p_netChAsym_cos = new TProfile("p_netChAsym_cos", "p_netChAsym_cos", 300, -1.5 + 0.0025, 1.5 + 0.0025, -1, 1);
TProfile *Hist_cos_Ach = new TProfile("Hist_cos_Ach", "Hist_cos_Ach", 5, 0.5, 5.5, -1, 1);
TProfile *p_v2_Ach = new TProfile("p_v2_Ach", "p_v2_Ach", 300, -1.5 + 0.0025, 1.5 + 0.0025, -100, 100);
TH2D *Hist_pt_pos_Ach = new TH2D("Hist_pt_pos_Ach", "Hist_pt_pos_Ach", 5, 0.5, 5.5, 300, 0, 15);
TH2D *Hist_pt_neg_Ach = new TH2D("Hist_pt_neg_Ach", "Hist_pt_neg_Ach", 5, 0.5, 5.5, 300, 0, 15);
TH1D *Hist_pt_Ach[5];
TProfile2D *p_v2_pt_pos_Ach = new TProfile2D("p_v2_pt_pos_Ach", "p_v2_pt_pos_Ach", 5, 0.5, 5.5, 300, 0, 15, -100, 100);
TProfile2D *p_v2_pt_neg_Ach = new TProfile2D("p_v2_pt_neg_Ach", "p_v2_pt_neg_Ach", 5, 0.5, 5.5, 300, 0, 15, -100, 100);
TH1D *Hist_inv_Mass = new TH1D("Hist_inv_Mass", "Hist_inv_Mass", 600, -0.5, 5.5);

TProfile2D *pTPCmeanPhi_FF_1 = new TProfile2D("TPCmeanPhi_FF_1", "TPCmeanPhi_FF_1",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_1_p = new TProfile2D("TPCmeanPhi_FF_1_p", "TPCmeanPhi_FF_1_p",
                                                8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_1_rot = new TProfile2D("TPCmeanPhi_FF_1_rot", "TPCmeanPhi_FF_1_rot",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_1 = new TProfile2D("TPCmeanPhi_RF_1", "TPCmeanPhi_RF_1",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_1_p = new TProfile2D("TPCmeanPhi_RF_1_p", "TPCmeanPhi_RF_1_p",
                                                8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_1_rot = new TProfile2D("TPCmeanPhi_RF_1_rot", "TPCmeanPhi_RF_1_rot",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_FF_1 = new TProfile2D("TPCmeanPhiAsso_FF_1", "TPCmeanPhiAsso_FF_1",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_RF_1 = new TProfile2D("TPCmeanPhiAsso_RF_1", "TPCmeanPhiAsso_RF_1",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_2 = new TProfile2D("TPCmeanPhi_FF_2", "TPCmeanPhi_FF_2",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_2_p = new TProfile2D("TPCmeanPhi_FF_2_p", "TPCmeanPhi_FF_2_p",
                                                8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_2_rot = new TProfile2D("TPCmeanPhi_FF_2_rot", "TPCmeanPhi_FF_2_rot",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_2 = new TProfile2D("TPCmeanPhi_RF_2", "TPCmeanPhi_RF_2",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_2_p = new TProfile2D("TPCmeanPhi_RF_2_p", "TPCmeanPhi_RF_2_p",
                                                8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_2_rot = new TProfile2D("TPCmeanPhi_RF_2_rot", "TPCmeanPhi_RF_2_rot",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_FF_2 = new TProfile2D("TPCmeanPhiAsso_FF_2", "TPCmeanPhiAsso_FF_2",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_RF_2 = new TProfile2D("TPCmeanPhiAsso_RF_2", "TPCmeanPhiAsso_RF_2",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_3 = new TProfile2D("TPCmeanPhi_FF_3", "TPCmeanPhi_FF_3",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_3_p = new TProfile2D("TPCmeanPhi_FF_3_p", "TPCmeanPhi_FF_3_p",
                                                8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_3_rot = new TProfile2D("TPCmeanPhi_FF_3_rot", "TPCmeanPhi_FF_3_rot",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_3 = new TProfile2D("TPCmeanPhi_RF_3", "TPCmeanPhi_RF_3",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_3_p = new TProfile2D("TPCmeanPhi_RF_3_p", "TPCmeanPhi_RF_3_p",
                                                8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_3_rot = new TProfile2D("TPCmeanPhi_RF_3_rot", "TPCmeanPhi_RF_3_rot",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_FF_3 = new TProfile2D("TPCmeanPhiAsso_FF_3", "TPCmeanPhiAsso_FF_3",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_RF_3 = new TProfile2D("TPCmeanPhiAsso_RF_3", "TPCmeanPhiAsso_RF_3",
                                                  8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TH2D *Hist_TPC_EP_east = new TH2D("Hist_TPC_EP_east", "Hist_TPC_EP_east", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_west = new TH2D("Hist_TPC_EP_west", "Hist_TPC_EP_west", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_full = new TH2D("Hist_TPC_EP_full", "Hist_TPC_EP_full", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for = new TH2D("Hist_TPC_EP_for", "Hist_TPC_EP_for", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac = new TH2D("Hist_TPC_EP_bac", "Hist_TPC_EP_bac", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_pos = new TH2D("Hist_TPC_EP_for_pos", "Hist_TPC_EP_for_pos", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_neg = new TH2D("Hist_TPC_EP_for_neg", "Hist_TPC_EP_for_neg", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_pos = new TH2D("Hist_TPC_EP_bac_pos", "Hist_TPC_EP_bac_pos", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_neg = new TH2D("Hist_TPC_EP_bac_neg", "Hist_TPC_EP_bac_neg", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_east_flat = new TH2D("Hist_TPC_EP_east_flat", "Hist_TPC_EP_east_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_west_flat = new TH2D("Hist_TPC_EP_west_flat", "Hist_TPC_EP_west_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_full_flat = new TH2D("Hist_TPC_EP_full_flat", "Hist_TPC_EP_full_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_flat = new TH2D("Hist_TPC_EP_for_flat", "Hist_TPC_EP_for_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_flat = new TH2D("Hist_TPC_EP_bac_flat", "Hist_TPC_EP_bac_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_pos_flat = new TH2D("Hist_TPC_EP_for_pos_flat", "Hist_TPC_EP_for_pos_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_neg_flat = new TH2D("Hist_TPC_EP_for_neg_flat", "Hist_TPC_EP_for_neg_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_pos_flat = new TH2D("Hist_TPC_EP_bac_pos_flat", "Hist_TPC_EP_bac_pos_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_neg_flat = new TH2D("Hist_TPC_EP_bac_neg_flat", "Hist_TPC_EP_bac_neg_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH1D *Hist_TPC_EP_full_m1 = new TH1D("Hist_TPC_EP_full_m1", "Hist_TPC_EP_full_m1", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m1_rot = new TH1D("Hist_TPC_EP_full_m1_rot", "Hist_TPC_EP_full_m1_rot", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2 = new TH1D("Hist_TPC_EP_full_m2", "Hist_TPC_EP_full_m2", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2_rot = new TH1D("Hist_TPC_EP_full_m2_rot", "Hist_TPC_EP_full_m2_rot", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m1_flat = new TH1D("Hist_TPC_EP_full_m1_flat", "Hist_TPC_EP_full_m1_flat", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m1_flat_rot = new TH1D("Hist_TPC_EP_full_m1_flat_rot", "Hist_TPC_EP_full_m1_flat_rot", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2_flat = new TH1D("Hist_TPC_EP_full_m2_flat", "Hist_TPC_EP_full_m2_flat", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2_flat_rot = new TH1D("Hist_TPC_EP_full_m2_flat_rot", "Hist_TPC_EP_full_m2_flat_rot", 36, 0, PI);
TProfile2D *pTPC_EP_east = new TProfile2D("pTPC_EP_east", "pTPC_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_west = new TProfile2D("pTPC_EP_west", "pTPC_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_full = new TProfile2D("pTPC_EP_full", "pTPC_EP_full", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_for = new TProfile2D("pTPC_EP_for", "pTPC_EP_for", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_bac = new TProfile2D("pTPC_EP_bac", "pTPC_EP_bac", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_for_pos = new TProfile2D("pTPC_EP_for_pos", "pTPC_EP_for_pos", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_for_neg = new TProfile2D("pTPC_EP_for_neg", "pTPC_EP_for_neg", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_bac_pos = new TProfile2D("pTPC_EP_bac_pos", "pTPC_EP_bac_pos", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_bac_neg = new TProfile2D("pTPC_EP_bac_neg", "pTPC_EP_bac_neg", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile *Hist_cos = new TProfile("Hist_cos", "Hist_cos", 6, 0.5, 6.5, -1, 1, "");
TProfile *Hist_cos_BBC = new TProfile("Hist_cos_BBC", "Hist_cos_BBC", 4, 0.5, 4.5, -1, 1, "");
// TProfile *Hist_cos_EPD = new TProfile("Hist_cos_EPD", "Hist_cos_EPD", 4, 0.5, 4.5, -1, 1, "");
TProfile *Hist_cos_ZDC = new TProfile("Hist_cos_ZDC", "Hist_cos_ZDC", 4, 0.5, 4.5, -1, 1, "");

TH1D *Hist_DCA = new TH1D("Hist_DCA", "Hist_DCA", 100, 0, 10);
TH2D *hEtaPtDist = new TH2D("EtaPtDist", "EtaPtDist", 40, -2, 2, 300, 0, 15);
TH2D *hEtaPtDist_anti = new TH2D("EtaPtDist_anti", "EtaPtDist_anti", 40, -2, 2, 300, 0, 15);
TH2D *hEtaPt_rot_Dist = new TH2D("EtaPt_rot_Dist", "EtaPt_rot_Dist", 40, -2, 2, 300, 0, 15);
TH2D *hEtaPhiDist = new TH2D("hEtaPhiDist", "hEtaPhiDist", 40, -2, 2, Phibin, -PI, PI);
TH2D *hPhiPtDist = new TH2D("PhiPtDist", "PhiPtDist", Phibin, -PI, PI, 300, 0, 15);
TH1D *Hist_Pt = new TH1D("Hist_Pt", "Hist_Pt", 300, 0, 15);
TH1D *Hist_Pt_rot = new TH1D("Hist_Pt_rot", "Hist_Pt_rot", 300, 0, 15);
TH1D *Hist_Pt2 = new TH1D("Hist_Pt2", "Hist_Pt2", 300, 0, 15);
TH1D *Hist_Pt2_rot = new TH1D("Hist_Pt2_rot", "Hist_Pt2_rot", 300, 0, 15);
TH1D *Hist_Pt2_anti = new TH1D("Hist_Pt2_anti", "Hist_Pt2_anti", 300, 0, 15);
TH1D *Hist_Pt2_part = new TH1D("Hist_Pt2_part", "Hist_Pt2_part", 300, 0, 15);
TH1D *Hist_Pt_TOF = new TH1D("Hist_Pt_TOF", "Hist_Pt_TOF", 300, 0, 15);
TH1D *rc;
// TH2D *wt = new TH2D("Order1etaWeight", "Order1etaWeight", 100, 1.5, 6.5, 9, 0, 9);
// TH2D *wt2 = new TH2D("Order2etaWeight", "Order2etaWeight", 100, 1.5, 6.5, 9, 0, 9);
TH2D *Hist_Phi = new TH2D("Hist_Phi", "Hist_Phi", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_1 = new TH2D("Hist_Phi_FF_1", "Hist_Phi_FF_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_1_p = new TH2D("Hist_Phi_FF_1_p", "Hist_Phi_FF_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_1_rot = new TH2D("Hist_Phi_FF_1_rot", "Hist_Phi_FF_1_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_1 = new TH2D("Hist_Phi_RF_1", "Hist_Phi_RF_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_1_p = new TH2D("Hist_Phi_RF_1_p", "Hist_Phi_RF_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_1_rot = new TH2D("Hist_Phi_RF_1_rot", "Hist_Phi_RF_1_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_1 = new TH2D("Hist_Phi_FF_new_1", "Hist_Phi_FF_new_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_1 = new TH2D("Hist_Phi_RF_new_1", "Hist_Phi_RF_new_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_1_p = new TH2D("Hist_Phi_FF_new_1_p", "Hist_Phi_FF_new_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_1_p = new TH2D("Hist_Phi_RF_new_1_p", "Hist_Phi_RF_new_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_1_rot = new TH2D("Hist_Phi_FF_new_1_rot", "Hist_Phi_FF_new_1_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_1_rot = new TH2D("Hist_Phi_RF_new_1_rot", "Hist_Phi_RF_new_1_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_2 = new TH2D("Hist_Phi_FF_2", "Hist_Phi_FF_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_2_p = new TH2D("Hist_Phi_FF_2_p", "Hist_Phi_FF_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_2_rot = new TH2D("Hist_Phi_FF_2_rot", "Hist_Phi_FF_2_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_2 = new TH2D("Hist_Phi_RF_2", "Hist_Phi_RF_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_2_p = new TH2D("Hist_Phi_RF_2_p", "Hist_Phi_RF_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_2_rot = new TH2D("Hist_Phi_RF_2_rot", "Hist_Phi_RF_2_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_2 = new TH2D("Hist_Phi_FF_new_2", "Hist_Phi_FF_new_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_2 = new TH2D("Hist_Phi_RF_new_2", "Hist_Phi_RF_new_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_2_p = new TH2D("Hist_Phi_FF_new_2_p", "Hist_Phi_FF_new_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_2_p = new TH2D("Hist_Phi_RF_new_2_p", "Hist_Phi_RF_new_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_2_rot = new TH2D("Hist_Phi_FF_new_2_rot", "Hist_Phi_FF_new_2_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_2_rot = new TH2D("Hist_Phi_RF_new_2_rot", "Hist_Phi_RF_new_2_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_3 = new TH2D("Hist_Phi_FF_3", "Hist_Phi_FF_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_3_p = new TH2D("Hist_Phi_FF_3_p", "Hist_Phi_FF_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_3_rot = new TH2D("Hist_Phi_FF_3_rot", "Hist_Phi_FF_3_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_3 = new TH2D("Hist_Phi_RF_3", "Hist_Phi_RF_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_3_p = new TH2D("Hist_Phi_RF_3_p", "Hist_Phi_RF_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_3_rot = new TH2D("Hist_Phi_RF_3_rot", "Hist_Phi_RF_3_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_3 = new TH2D("Hist_Phi_FF_new_3", "Hist_Phi_FF_new_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_3 = new TH2D("Hist_Phi_RF_new_3", "Hist_Phi_RF_new_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_3_p = new TH2D("Hist_Phi_FF_new_3_p", "Hist_Phi_FF_new_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_3_p = new TH2D("Hist_Phi_RF_new_3_p", "Hist_Phi_RF_new_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_3_rot = new TH2D("Hist_Phi_FF_new_3_rot", "Hist_Phi_FF_new_3_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_3_rot = new TH2D("Hist_Phi_RF_new_3_rot", "Hist_Phi_RF_new_3_rot", Phibin, -PI, PI, 4, 0.5, 4.5);
TH1D *hDpt = new TH1D("hDpt", "hDpt", 200, 0, 2);
TH1D *hDpt_rot = new TH1D("hDpt_rot", "hDpt_rot", 200, 0, 2);

TProfile *pParity_int_obs1 = new TProfile("Parity_int_obs1", "Parity_int_obs1", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs1_rot = new TProfile("Parity_int_obs1_rot", "Parity_int_obs1_rot", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3 = new TProfile("Parity_int_obs3", "Parity_int_obs3", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_rot = new TProfile("Parity_int_obs3_rot", "Parity_int_obs3_rot", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs1_anti = new TProfile("Parity_int_obs1_anti", "Parity_int_obs1_anti", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs1_anti_rot = new TProfile("Parity_int_obs1_anti_rot", "Parity_int_obs1_anti_rot", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_anti = new TProfile("Parity_int_obs3_anti", "Parity_int_obs3_anti", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_anti_rot = new TProfile("Parity_int_obs3_anti_rot", "Parity_int_obs3_anti_rot", 24, 0.5, 24.5, -100, 100, "");
// TProfile *pParity_int_obs1_splitpt[3][15];
// TProfile *pParity_int_obs1_splitpt_rot[3][15];
// TProfile *pParity_int_obs1_splitpt_anti[3][15];
// TProfile *pParity_int_obs1_splitpt_anti_rot[3][15];
TProfile *pParity_int_obs3_splitpt[6][15][8];
// TProfile *pParity_int_obs3_splitpt_rot[6][15];
// TProfile *pParity_int_obs3_splitpt_anti[6][15];
TProfile *pParity_int_obs3_splitpt_anti_rot[6][15];
TProfile *pParity_int_obs1_QQcut = new TProfile("Parity_int_obs1_QQcut", "Parity_int_obs1_QQcut", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs1_rot_QQcut = new TProfile("Parity_int_obs1_QQcut_rot", "Parity_int_obs1_rot_QQcut", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_QQcut = new TProfile("Parity_int_obs3_QQcut", "Parity_int_obs3_QQcut", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_rot_QQcut = new TProfile("Parity_int_obs3_QQcut_rot", "Parity_int_obs3_rot_QQcut", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs1_QQcut_anti = new TProfile("Parity_int_obs1_QQcut_anti", "Parity_int_obs1_QQcut_anti", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs1_rot_QQcut_anti = new TProfile("Parity_int_obs1_QQcut_anti_rot", "Parity_int_obs1_rot_QQcut_anti", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_QQcut_anti = new TProfile("Parity_int_obs3_QQcut_anti", "Parity_int_obs3_QQcut_anti", 24, 0.5, 24.5, -100, 100, "");
TProfile *pParity_int_obs3_rot_QQcut_anti = new TProfile("Parity_int_obs3_QQcut_anti_rot", "Parity_int_obs3_rot_QQcut_anti", 24, 0.5, 24.5, -100, 100, "");
// TProfile *pParity_int_obs1_splitpt_QQcut[3][15];
// TProfile *pParity_int_obs1_splitpt_rot_QQcut[3][15];
// TProfile *pParity_int_obs1_splitpt_QQcut_anti[3][15];
// TProfile *pParity_int_obs1_splitpt_rot_QQcut_anti[3][15];
// TProfile *pParity_int_obs3_splitpt_QQcut[6][15];
// TProfile *pParity_int_obs3_splitpt_rot_QQcut[6][15];
// TProfile *pParity_int_obs3_splitpt_QQcut_anti[6][15];
TProfile *pParity_int_obs3_splitpt_rot_QQcut_anti[6][15];
TProfile *pParity_int_pt_ss_obs3 = new TProfile("pParity_int_pt_ss_obs3", "Parity_int_pt_ss_obs3", 12, 0, 3, -100, 100, "");
TProfile *pParity_int_pt_os_obs3 = new TProfile("pParity_int_pt_os_obs3", "Parity_int_pt_os_obs3", 12, 0, 3, -100, 100, "");
TProfile *pParity_int_pt_ssb_obs3 = new TProfile("pParity_int_pt_ssb_obs3", "Parity_int_pt_ssb_obs3", 12, 0, 3, -100, 100, "");
TProfile *pParity_int_pt_osb_obs3 = new TProfile("pParity_int_pt_osb_obs3", "Parity_int_pt_osb_obs3", 12, 0, 3, -100, 100, "");
TProfile *pParity_int_ss_obs1 = new TProfile("Parity_int_ss_obs1", "Parity_int_ss_obs1", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs1_rot = new TProfile("Parity_int_ss_obs1_rot", "Parity_int_ss_obs1_rot", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs3 = new TProfile("Parity_int_ss_obs3", "Parity_int_ss_obs3", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_rot = new TProfile("Parity_int_ss_obs3_rot", "Parity_int_ss_obs3_rot", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs1_anti = new TProfile("Parity_int_ss_obs1_anti", "Parity_int_ss_obs1_anti", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs1_anti_rot = new TProfile("Parity_int_ss_obs1_anti_rot", "Parity_int_ss_obs1_anti_rot", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_anti = new TProfile("Parity_int_ss_obs3_anti", "Parity_int_ss_obs3_anti", 12, 0.5, 12.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_anti_rot = new TProfile("Parity_int_ss_obs3_anti_rot", "Parity_int_ss_obs3_anti_rot", 12, 0.5, 12.5, -100, 100, "");
// TProfile *pParity_int_ss_obs1_splitpt[15];
// TProfile *pParity_int_ss_obs1_splitpt_rot[15];
// TProfile *pParity_int_ss_obs1_splitpt_anti[15];
// TProfile *pParity_int_ss_obs1_splitpt_anti_rot[15];
// TProfile *pParity_int_ss_obs3_splitpt[15];
// TProfile *pParity_int_ss_obs3_splitpt_rot[15];
// TProfile *pParity_int_ss_obs3_splitpt_anti[15];
// TProfile *pParity_int_ss_obs3_splitpt_anti_rot[15];
TProfile *pParity_int_ss_obs1_QQcut = new TProfile("Parity_int_ss_obs1_QQcut", "Parity_int_ss_obs1_QQcut", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs1_rot_QQcut = new TProfile("Parity_int_ss_obs1_QQcut_rot", "Parity_int_ss_obs1_rot_QQcut", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_QQcut = new TProfile("Parity_int_ss_obs3_QQcut", "Parity_int_ss_obs3_QQcut", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_rot_QQcut = new TProfile("Parity_int_ss_obs3_QQcut_rot", "Parity_int_ss_obs3_rot_QQcut", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs1_QQcut_anti = new TProfile("Parity_int_ss_obs1_QQcut_anti", "Parity_int_ss_obs1_QQcut_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs1_rot_QQcut_anti = new TProfile("Parity_int_ss_obs1_QQcut_anti_rot", "Parity_int_ss_obs1_rot_QQcut_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_QQcut_anti = new TProfile("Parity_int_ss_obs3_QQcut_anti", "Parity_int_ss_obs3_QQcut_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs3_rot_QQcut_anti = new TProfile("Parity_int_ss_obs3_QQcut_anti_rot", "Parity_int_ss_obs3_rot_QQcut_anti", 4, 0.5, 4.5, -100, 100, "");
// TProfile *pParity_int_ss_obs1_splitpt_QQcut[15];
// TProfile *pParity_int_ss_obs1_splitpt_rot_QQcut[15];
// TProfile *pParity_int_ss_obs1_splitpt_QQcut_anti[15];
// TProfile *pParity_int_ss_obs1_splitpt_rot_QQcut_anti[15];
// TProfile *pParity_int_ss_obs3_splitpt_QQcut[15];
// TProfile *pParity_int_ss_obs3_splitpt_rot_QQcut[15];
// TProfile *pParity_int_ss_obs3_splitpt_QQcut_anti[15];
// TProfile *pParity_int_ss_obs3_splitpt_rot_QQcut_anti[15];
TProfile *pParity_int_ss_same_run = new TProfile("Parity_int_ss_same_run", "Parity_int_ss_same_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile *pParity_int_ss_oppo_run = new TProfile("Parity_int_ss_oppo_run", "Parity_int_ss_oppo_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile2D *pParity_eta_ss_obs1 = new TProfile2D("Parity_eta_ss_obs1", "Parity_eta_ss_obs1", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pParity_eta_ss_obs3 = new TProfile2D("Parity_eta_ss_obs3", "Parity_eta_ss_obs3", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pParity_Deta_ss_obs1 = new TProfile2D("Parity_Deta_ss_obs1", "Parity_Deta_ss_obs1", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pParity_Deta_ss_obs3 = new TProfile2D("Parity_Deta_ss_obs3", "Parity_Deta_ss_obs3", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pParity_pt_ss_obs1 = new TProfile2D("Parity_pt_ss_obs1", "Parity_pt_ss_obs1", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pParity_pt_ss_obs3 = new TProfile2D("Parity_pt_ss_obs3", "Parity_pt_ss_obs3", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pParity_Dpt_ss_obs1 = new TProfile2D("Parity_Dpt_ss_obs1", "Parity_Dpt_ss_obs1", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile2D *pParity_Dpt_ss_obs3 = new TProfile2D("Parity_Dpt_ss_obs3", "Parity_Dpt_ss_obs3", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile *pParity_noHBT_ss_obs1 = new TProfile("Parity_noHBT_ss_obs1", "Parity_noHBT_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_noHBT_ss_obs3 = new TProfile("Parity_noHBT_ss_obs3", "Parity_noHBT_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile2D *pParity_Deta_highDpt_ss_obs1 = new TProfile2D("pParity_Deta_highDpt_ss_obs1", "pParity_Deta_highDpt_ss_obs1", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");
TProfile2D *pParity_Deta_highDpt_ss_obs3 = new TProfile2D("pParity_Deta_highDpt_ss_obs3", "pParity_Deta_highDpt_ss_obs3", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");

TProfile *pDelta_int_ss_obs1 = new TProfile("Delta_int_ss_obs1", "Delta_int_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs1_rot = new TProfile("Delta_int_ss_obs1_rot", "Delta_int_ss_obs1_rot", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3 = new TProfile("Delta_int_ss_obs3", "Delta_int_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_rot = new TProfile("Delta_int_ss_obs3_rot", "Delta_int_ss_obs3_rot", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs1_anti = new TProfile("Delta_int_ss_obs1_anti", "Delta_int_ss_obs1_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs1_anti_rot = new TProfile("Delta_int_ss_obs1_anti_rot", "Delta_int_ss_obs1_anti_rot", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_anti = new TProfile("Delta_int_ss_obs3_anti", "Delta_int_ss_obs3_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_anti_rot = new TProfile("Delta_int_ss_obs3_anti_rot", "Delta_int_ss_obs3_anti_rot", 4, 0.5, 4.5, -100, 100, "");
// TProfile *pDelta_int_ss_obs1_splitpt[15];
// TProfile *pDelta_int_ss_obs1_splitpt_rot[15];
// TProfile *pDelta_int_ss_obs1_splitpt_anti[15];
// TProfile *pDelta_int_ss_obs1_splitpt_anti_rot[15];
// TProfile *pDelta_int_ss_obs3_splitpt[15];
// TProfile *pDelta_int_ss_obs3_splitpt_rot[15];
// TProfile *pDelta_int_ss_obs3_splitpt_anti[15];
// TProfile *pDelta_int_ss_obs3_splitpt_anti_rot[15];
TProfile *pDelta_int_ss_obs1_QQcut = new TProfile("Delta_int_ss_obs1_QQcut", "Delta_int_ss_obs1_QQcut", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs1_rot_QQcut = new TProfile("Delta_int_ss_obs1_QQcut_rot", "Delta_int_ss_obs1_QQcut_rot", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_QQcut = new TProfile("Delta_int_ss_obs3_QQcut", "Delta_int_ss_obs3_QQcut", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_rot_QQcut = new TProfile("Delta_int_ss_obs3_QQcut_rot", "Delta_int_ss_obs3_QQcut_rot", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs1_QQcut_anti = new TProfile("Delta_int_ss_obs1_QQcut_anti", "Delta_int_ss_obs1_QQcut_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs1_rot_QQcut_anti = new TProfile("Delta_int_ss_obs1_QQcut_anti_rot", "Delta_int_ss_obs1_QQcut_anti_rot", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_QQcut_anti = new TProfile("Delta_int_ss_obs3_QQcut_anti", "Delta_int_ss_obs3_QQcut_anti", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3_rot_QQcut_anti = new TProfile("Delta_int_ss_obs3_QQcut_anti_rot", "Delta_int_ss_obs3_QQcut_anti_rot", 4, 0.5, 4.5, -100, 100, "");
// TProfile *pDelta_int_ss_obs1_splitpt_QQcut[15];
// TProfile *pDelta_int_ss_obs1_splitpt_rot_QQcut[15];
// TProfile *pDelta_int_ss_obs1_splitpt_QQcut_anti[15];
// TProfile *pDelta_int_ss_obs1_splitpt_rot_QQcut_anti[15];
// TProfile *pDelta_int_ss_obs3_splitpt_QQcut[15];
// TProfile *pDelta_int_ss_obs3_splitpt_rot_QQcut[15];
// TProfile *pDelta_int_ss_obs3_splitpt_QQcut_anti[15];
// TProfile *pDelta_int_ss_obs3_splitpt_rot_QQcut_anti[15];
TProfile *pDelta_int_ss_same_run = new TProfile("Delta_int_ss_same_run", "Delta_int_ss_same_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile *pDelta_int_ss_oppo_run = new TProfile("Delta_int_ss_oppo_run", "Delta_int_ss_oppo_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile2D *pDelta_eta_ss_obs1 = new TProfile2D("Delta_eta_ss_obs1", "Delta_eta_ss_obs1", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pDelta_eta_ss_obs3 = new TProfile2D("Delta_eta_ss_obs3", "Delta_eta_ss_obs3", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pDelta_Deta_ss_obs1 = new TProfile2D("Delta_Deta_ss_obs1", "Delta_Deta_ss_obs1", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pDelta_Deta_ss_obs3 = new TProfile2D("Delta_Deta_ss_obs3", "Delta_Deta_ss_obs3", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pDelta_pt_ss_obs1 = new TProfile2D("Delta_pt_ss_obs1", "Delta_pt_ss_obs1", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pDelta_pt_ss_obs3 = new TProfile2D("Delta_pt_ss_obs3", "Delta_pt_ss_obs3", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pDelta_Dpt_ss_obs1 = new TProfile2D("Delta_Dpt_ss_obs1", "Delta_Dpt_ss_obs1", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile2D *pDelta_Dpt_ss_obs3 = new TProfile2D("Delta_Dpt_ss_obs3", "Delta_Dpt_ss_obs3", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile *pDelta_noHBT_ss_obs1 = new TProfile("Delta_noHBT_ss_obs1", "Delta_noHBT_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_noHBT_ss_obs3 = new TProfile("Delta_noHBT_ss_obs3", "Delta_noHBT_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile2D *pDelta_Deta_highDpt_ss_obs1 = new TProfile2D("pDelta_Deta_highDpt_ss_obs1", "pDelta_Deta_highDpt_ss_obs1", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");
TProfile2D *pDelta_Deta_highDpt_ss_obs3 = new TProfile2D("pDelta_Deta_highDpt_ss_obs3", "pDelta_Deta_highDpt_ss_obs3", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");

TProfile *Hist_v2_pt_obs1 = new TProfile("Hist_v2_pt_obs1", "Hist_v2_pt_obs1", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs1_rot = new TProfile("Hist_v2_pt_obs1_rot", "Hist_v2_pt_obs1_rot", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs1_p = new TProfile("Hist_v2_pt_obs1_p", "Hist_v2_pt_obs1_p", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs1_ap = new TProfile("Hist_v2_pt_obs1_ap", "Hist_v2_pt_obs1_ap", 300, 0, 15, -100, 100, "");
// TProfile *Hist_v2_pt_obs1_alltrks = new TProfile("Hist_v2_pt_obs1_alltrks", "Hist_v2_pt_obs1_alltrks", 300, 0, 15, -100, 100, "");
// TProfile *Hist_v2_pt_obs1_alltrks_TPC = new TProfile("Hist_v2_pt_obs1_alltrks_TPC", "Hist_v2_pt_obs1_alltrks_TPC", 300, 0, 15, -100, 100, "");
// TProfile *Hist_v2_pt_obs1_alltrks_1stOrder = new TProfile("Hist_v2_pt_obs1_alltrks_1stOrder", "Hist_v2_pt_obs1_alltrks_1stOrder", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs1_alltrks[3];
TProfile *Hist_v2_Ach = new TProfile("Hist_v2_Ach", "V2 against Charge Assymetry", 300, -1.5 + 0.0025, 1.5 + 0.0025, -100, 100);
TProfile *Hist_v2_Ach_p = new TProfile("Hist_v2_Ach_p", "V2 against Charge Assymetry protons", 300, -1.5 + 0.0025, 1.5 + 0.0025, -100, 100);
TProfile *Hist_v2_pt_obs2 = new TProfile("Hist_v2_pt_obs2", "Hist_v2_pt_obs2", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs2_rot = new TProfile("Hist_v2_pt_obs2_rot", "Hist_v2_pt_obs2_rot", 300, 0, 15, -100, 100, "");
TProfile *Hist_112ss_v2_parent_obs2;
TProfile *Hist_112os_v2_parent_obs2;
TProfile *Hist_132ss_v2_parent_obs2;
TProfile *Hist_132os_v2_parent_obs2;
TProfile *Hist_112ss_v2_parent_obs2_low;
TProfile *Hist_112os_v2_parent_obs2_low;
TProfile *Hist_132ss_v2_parent_obs2_low;
TProfile *Hist_132os_v2_parent_obs2_low;
TProfile *Hist_112ss_v2_parent_obs2_high;
TProfile *Hist_112os_v2_parent_obs2_high;
TProfile *Hist_132ss_v2_parent_obs2_high;
TProfile *Hist_132os_v2_parent_obs2_high;

TProfile *Hist_112ss_pt_parent_obs2 = new TProfile("Hist_112ss_pt_parent_obs2", "Hist_112ss_pt_parent_obs2", 100, 0, 5, -100, 100, "");
TProfile *Hist_112os_pt_parent_obs2 = new TProfile("Hist_112os_pt_parent_obs2", "Hist_112os_pt_parent_obs2", 100, 0, 5, -100, 100, "");
TProfile *Hist_132ss_pt_parent_obs2 = new TProfile("Hist_132ss_pt_parent_obs2", "Hist_132ss_pt_parent_obs2", 100, 0, 5, -100, 100, "");
TProfile *Hist_132os_pt_parent_obs2 = new TProfile("Hist_132os_pt_parent_obs2", "Hist_132os_pt_parent_obs2", 100, 0, 5, -100, 100, "");
TProfile *Hist_112ss_m_parent_obs2 = new TProfile("Hist_112ss_m_parent_obs2", "Hist_112ss_m_parent_obs2", 300, -0.5, 5.5, -100, 100, "");
TProfile *Hist_112os_m_parent_obs2 = new TProfile("Hist_112os_m_parent_obs2", "Hist_112os_m_parent_obs2", 300, -0.5, 5.5, -100, 100, "");
TProfile *Hist_v2_parent_lam_pT = new TProfile("Hist_v2_parent_lam_pT", "Hist_v2_parent_lam_pT", 100, 0, 5, -100, 100);
TProfile *Hist_v2_parent_parent_pT = new TProfile("Hist_v2_parent_parent_pT", "Hist_v2_parent_parent_pT", 100, 0, 5, -100, 100);
TProfile *Hist_v2_parent_parent_pT_ss = new TProfile("Hist_v2_parent_parent_pT_ss", "Hist_v2_parent_parent_pT_ss", 100, 0, 5, -100, 100);
TProfile *Hist_v2_parent_parent_pT_os = new TProfile("Hist_v2_parent_parent_pT_os", "Hist_v2_parent_parent_pT_os", 100, 0, 5, -100, 100);
TProfile *Hist_v2_lam_parent_pT = new TProfile("Hist_v2_lam_parent_pT", "Hist_v2_lam_parent_pT", 100, 0, 5, -100, 100);
TProfile *Hist_v2_p_parent_pT = new TProfile("Hist_v2_p_parent_pT", "Hist_v2_p_parent_pT", 100, 0, 5, -100, 100);
TH1D *Hist_parent_phi_low_pT = new TH1D("Hist_parent_phi_low_pT", "Hist_parent_phi_low_pT", 800, -4, 4);
TH1D *Hist_parent_phi_high_pT = new TH1D("Hist_parent_phi_high_pT", "Hist_parent_phi_high_pT", 800, -4, 4);
TH1D *Hist_delta_phi_low_parent_pT = new TH1D("Hist_delta_phi_low_parent_pT", "Hist_delta_phi_low_parent_pT", 800, -4, 4);
TH1D *Hist_delta_phi_high_parent_pT = new TH1D("Hist_delta_phi_high_parent_pT", "Hist_delta_phi_high_parent_pT", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_1 = new TH1D("Hist_delta_phi_parent_pT_1", "Hist_delta_phi_parent_pT_1", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_2 = new TH1D("Hist_delta_phi_parent_pT_2", "Hist_delta_phi_parent_pT_2", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_3 = new TH1D("Hist_delta_phi_parent_pT_3", "Hist_delta_phi_parent_pT_3", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_4 = new TH1D("Hist_delta_phi_parent_pT_4", "Hist_delta_phi_parent_pT_4", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_5 = new TH1D("Hist_delta_phi_parent_pT_5", "Hist_delta_phi_parent_pT_5", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_6 = new TH1D("Hist_delta_phi_parent_pT_6", "Hist_delta_phi_parent_pT_6", 800, -4, 4);
TH1D *Hist_delta_phi_parent_pT_7 = new TH1D("Hist_delta_phi_parent_pT_7", "Hist_delta_phi_parent_pT_7", 800, -4, 4);
TProfile *Hist_gamma112ss_QQ = new TProfile("Hist_gamma112ss_QQ", "Hist_gamma112ss_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112os_QQ = new TProfile("Hist_gamma112os_QQ", "Hist_gamma112os_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112ss_lpt = new TProfile("Hist_gamma112ss_lpt", "Hist_gamma112ss_lpt", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112os_lpt = new TProfile("Hist_gamma112os_lpt", "Hist_gamma112os_lpt", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112ss_ppt = new TProfile("Hist_gamma112ss_ppt", "Hist_gamma112ss_ppt", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112os_ppt = new TProfile("Hist_gamma112os_ppt", "Hist_gamma112os_ppt", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112ss_k = new TProfile("Hist_gamma112ss_k", "Hist_gamma112ss_k", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112os_k = new TProfile("Hist_gamma112os_k", "Hist_gamma112os_k", 100, 0, 5, -100, 100);
TProfile *Hist_gamma132ss_QQ = new TProfile("Hist_gamma132ss_QQ", "Hist_gamma132ss_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma132os_QQ = new TProfile("Hist_gamma132os_QQ", "Hist_gamma132os_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112ss_pi_QQ = new TProfile("Hist_gamma112ss_pi_QQ", "Hist_gamma112ss_pi_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112os_pi_QQ = new TProfile("Hist_gamma112os_pi_QQ", "Hist_gamma112os_pi_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112ss_pibar_QQ = new TProfile("Hist_gamma112ss_pibar_QQ", "Hist_gamma112ss_pibar_QQ", 100, 0, 5, -100, 100);
TProfile *Hist_gamma112os_pibar_QQ = new TProfile("Hist_gamma112os_pibar_QQ", "Hist_gamma112os_pibar_QQ", 100, 0, 5, -100, 100);
TH1D *Hist_QQ_dist_coupling = new TH1D("Hist_QQ_dist_coupling", "Hist_QQ_dist_coupling", 100, 0, 5);
TH1D *Hist_QQ_dist = new TH1D("Hist_QQ_dist", "Hist_QQ_dist", 100, 0, 5);
TH1D *Hist_lam_v2_low_parent_v2 = new TH1D("Hist_lam_v2_low_parent_v2", "Hist_lam_v2_low_parent_v2", 400, -100, 100);
TH1D *Hist_p_v2_low_parent_v2 = new TH1D("Hist_p_v2_low_parent_v2", "Hist_p_v2_low_parent_v2", 400, -100, 100);
TH1D *Hist_lam_pT_low_parent_v2 = new TH1D("Hist_lam_pT_low_parent_v2", "Hist_lam_pT_low_parent_v2", 100, 0, 5);
TH1D *Hist_p_pT_low_parent_v2 = new TH1D("Hist_p_pT_low_parent_v2", "Hist_p_pT_low_parent_v2", 100, 0, 5);
TProfile *Hist_v2_pt_obs2_caysm[5], *Hist_v2_pt_obs2_caysm_os[5];
TProfile *Hist_v2_pt_obs2_p = new TProfile("Hist_v2_pt_obs2_p", "Hist_v2_pt_obs2_p", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs2_ap = new TProfile("Hist_v2_pt_obs2_ap", "Hist_v2_pt_obs2_ap", 300, 0, 15, -100, 100, "");
// TProfile *Hist_v2_pt_obs2_alltrks = new TProfile("Hist_v2_pt_obs2_alltrks", "Hist_v2_pt_obs2_alltrks", 300, 0, 15, -100, 100, "");
// TProfile *Hist_v2_pt_obs2_alltrks_TPC = new TProfile("Hist_v2_pt_obs2_alltrks_TPC", "Hist_v2_pt_obs2_alltrks_TPC", 300, 0, 15, -100, 100, "");
// TProfile *Hist_v2_pt_obs2_alltrks_1stOrder = new TProfile("Hist_v2_pt_obs2_alltrks_1stOrder", "Hist_v2_pt_obs2_alltrks_1stOrder", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs2_alltrks[3];
TProfile *Hist_v2_pt_obs2_p_caysm[5];
TProfile *Hist_112_pt_ss = new TProfile("Hist_112_pt_ss", "Hist_112_pt_ss", 40, 0, 2, -100, 100, "");
TProfile *Hist_112_pt_os = new TProfile("Hist_112_pt_os", "Hist_112_pt_os", 40, 0, 2, -100, 100, "");
TProfile *Hist_132_pt_ss = new TProfile("Hist_132_pt_ss", "Hist_132_pt_ss", 40, 0, 2, -100, 100, "");
TProfile *Hist_132_pt_os = new TProfile("Hist_132_pt_os", "Hist_132_pt_os", 40, 0, 2, -100, 100, "");
TProfile *Hist_112_eta_ss = new TProfile("Hist_112_eta_ss", "Hist_112_eta_ss", 40, 0, 2, -100, 100, "");
TProfile *Hist_112_eta_os = new TProfile("Hist_112_eta_os", "Hist_112_eta_os", 40, 0, 2, -100, 100, "");
TProfile *Hist_132_eta_ss = new TProfile("Hist_132_eta_ss", "Hist_132_eta_ss", 40, 0, 2, -100, 100, "");
TProfile *Hist_132_eta_os = new TProfile("Hist_132_eta_os", "Hist_132_eta_os", 40, 0, 2, -100, 100, "");
TProfile *Hist_112_cc_pt_ss = new TProfile("Hist_112_cc_pt_ss", "Hist_112_cc_pt_ss", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_cc_pt_os = new TProfile("Hist_112_cc_pt_os", "Hist_112_cc_pt_os", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_ss_pt_ss = new TProfile("Hist_112_ss_pt_ss", "Hist_112_ss_pt_ss", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_ss_pt_os = new TProfile("Hist_112_ss_pt_os", "Hist_112_ss_pt_os", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_cc_eta_ss = new TProfile("Hist_112_cc_eta_ss", "Hist_112_cc_eta_ss", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_cc_eta_os = new TProfile("Hist_112_cc_eta_os", "Hist_112_cc_eta_os", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_ss_eta_ss = new TProfile("Hist_112_ss_eta_ss", "Hist_112_ss_eta_ss", 40, 0, 2, -3, 3, "");
TProfile *Hist_112_ss_eta_os = new TProfile("Hist_112_ss_eta_os", "Hist_112_ss_eta_os", 40, 0, 2, -3, 3, "");
TProfile *Hist_v2_eta_obs1 = new TProfile("Hist_v2_eta_obs1", "Hist_v2_eta_obs1", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs2 = new TProfile("Hist_v2_eta_obs2", "Hist_v2_eta_obs2", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs1_rot = new TProfile("Hist_v2_eta_obs1_rot", "Hist_v2_eta_obs1_rot", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs2_rot = new TProfile("Hist_v2_eta_obs2_rot", "Hist_v2_eta_obs2_rot", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs3 = new TProfile("Hist_v2_eta_obs3", "Hist_v2_eta_obs3", 40, -1, 1, -100, 100, "");

TProfile *pTemp_v2 = new TProfile("pTemp_v2", "pTemp_v2", 13, 0.5, 13.5, -100, 100, "");
TProfile *pTemp_v2_chargedhadrons = new TProfile("pTemp_v2_chargedhadrons", "pTemp_v2_chargedhadrons", 18, 0.5, 18.5, -100, 100, "");
TProfile *pTemp_v2_rot = new TProfile("pTemp_v2_rot", "pTemp_v2_rot", 13, 0.5, 13.5, -100, 100, "");
TProfile *pTemp_parity_e = new TProfile("pTemp_parity_e", "pTemp_parity_e", 36, 0.5, 36.5, -100, 100, "");
TProfile *pTemp_parity_w = new TProfile("pTemp_parity_w", "pTemp_parity_w", 36, 0.5, 36.5, -100, 100, "");
TProfile *pTemp_delta = new TProfile("pTemp_delta", "pTemp_delta", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_e_rot = new TProfile("pTemp_parity_e_rot", "pTemp_parity_e_rot", 36, 0.5, 36.5, -100, 100, "");
TProfile *pTemp_parity_w_rot = new TProfile("pTemp_parity_w_rot", "pTemp_parity_w_rot", 36, 0.5, 36.5, -100, 100, "");
TProfile *pTemp_delta_rot = new TProfile("pTemp_delta_rot", "pTemp_delta_rot", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_e_QQcut = new TProfile("pTemp_parity_e_QQcut", "pTemp_parity_e_QQcut", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w_QQcut = new TProfile("pTemp_parity_w_QQcut", "pTemp_parity_w_QQcut", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta_QQcut = new TProfile("pTemp_delta_QQcut", "pTemp_delta_QQcut", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_e_rot_QQcut = new TProfile("pTemp_parity_e_rot_QQcut", "pTemp_parity_e_rot_QQcut", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w_rot_QQcut = new TProfile("pTemp_parity_w_rot_QQcut", "pTemp_parity_w_rot_QQcut", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta_rot_QQcut = new TProfile("pTemp_delta_rot_QQcut", "pTemp_delta_rot_QQcut", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_v2_noHBT = new TProfile("pTemp_v2_noHBT", "pTemp_v2_noHBT", 6, 0.5, 6.5, -100, 100, "");
TProfile *pTemp_parity_e_noHBT = new TProfile("pTemp_parity_e_noHBT", "pTemp_parity_e_noHBT", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w_noHBT = new TProfile("pTemp_parity_w_noHBT", "pTemp_parity_w_noHBT", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta_noHBT = new TProfile("pTemp_delta_noHBT", "pTemp_delta_noHBT", 8, 0.5, 8.5, -100, 100, "");

TH1D *Hist_Q2_EPD = new TH1D("Hist_Q2_EPD", "Hist_Q2_EPD", 2000, 0, 50);
TH1D *Hist_Q2_EPD1 = new TH1D("Hist_Q2_EPD1", "Hist_Q2_EPD1", 2000, 0, 50);
TH1D *Hist_Q2_TPC = new TH1D("Hist_Q2_TPC", "Hist_Q2_TPC", 2000, 0, 50);
TH1D *Hist_Q2_EPD_pion = new TH1D("Hist_Q2_EPD_pion", "Hist_Q2_EPD_pion", 2000, 0, 50);
TH1D *Hist_Q2_EPD1_pion = new TH1D("Hist_Q2_EPD1_pion", "Hist_Q2_EPD1_pion", 2000, 0, 50);
TH1D *Hist_Q2_TPC_pion = new TH1D("Hist_Q2_TPC_pion", "Hist_Q2_TPC_pion", 2000, 0, 50);
TH2D *Hist_RefMult_Q2 = new TH2D("Hist_RefMult_Q2", "Hist_RefMult_Q2", 500, 0, 50, 1000, 0, 1000);
TH2D *Hist_RefMult_Q2_EPD = new TH2D("Hist_RefMult_Q2_EPD", "Hist_RefMult_Q2_EPD", 500, 0, 50, 1000, 0, 1000);
TProfile *p_RefMult_Q2 = new TProfile("p_RefMult_Q2", "p_RefMult_Q2", 250, 0, 25, 0, 1000, "");
TProfile *p_RefMult_Q2_EPD = new TProfile("p_RefMult_Q2_EPD", "p_RefMult_Q2_EPD", 250, 0, 25, 0, 1000, "");
TProfile *p_cos_Q2 = new TProfile("p_cos_Q2", "p_cos_Q2", 250, 0, 25, -1, 1, "");
TProfile *p_cos_Q2_EPD = new TProfile("p_cos_Q2_EPD", "p_cos_Q2_EPD", 250, 0, 25, -1, 1, "");
TProfile *p_cos_Q2_EPD1 = new TProfile("p_cos_Q2_EPD1", "p_cos_Q2_EPD1", 250, 0, 25, -1, 1, "");
TProfile *p_cos_Q2_TPC_pion = new TProfile("p_cos_Q2_pion", "p_cos_Q2_pion", 250, 0, 25, -1, 1, "");
TProfile *p_cos_Q2_EPD_pion = new TProfile("p_cos_Q2_EPD_pion", "p_cos_Q2_EPD_pion", 250, 0, 25, -1, 1, "");
TProfile *p_cos_Q2_EPD1_pion = new TProfile("p_cos_Q2_EPD1_pion", "p_cos_Q2_EPD1_pion", 250, 0, 25, -1, 1, "");
// TProfile *p_v2e_Q2_obs1 = new TProfile("p_v2e_Q2_obs1", "p_v2e_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2_Q2_obs1[6], *p_v2_Q2_obs2[6], *p_v2_Q2_obs1_rot[6], *p_v2_Q2_obs2_rot[6];

// TProfile *p_v2e_Q2_obs2 = new TProfile("p_v2e_Q2_obs2", "p_v2e_Q2_obs2", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_obs1 = new TProfile("p_v2w_Q2_obs1", "p_v2w_Q2_obs1", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_obs2 = new TProfile("p_v2w_Q2_obs2", "p_v2w_Q2_obs2", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_P_obs1 = new TProfile("p_v2e_Q2_P_obs1", "p_v2e_Q2_P_obs1", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_P_obs2 = new TProfile("p_v2e_Q2_P_obs2", "p_v2e_Q2_P_obs2", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_P_obs1 = new TProfile("p_v2w_Q2_P_obs1", "p_v2w_Q2_P_obs1", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_P_obs2 = new TProfile("p_v2w_Q2_P_obs2", "p_v2w_Q2_P_obs2", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_N_obs1 = new TProfile("p_v2e_Q2_N_obs1", "p_v2e_Q2_N_obs1", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_N_obs2 = new TProfile("p_v2e_Q2_N_obs2", "p_v2e_Q2_N_obs2", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_N_obs1 = new TProfile("p_v2w_Q2_N_obs1", "p_v2w_Q2_N_obs1", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_N_obs2 = new TProfile("p_v2w_Q2_N_obs2", "p_v2w_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_obs1_rot = new TProfile("p_v2e_Q2_obs1_rot", "p_v2e_Q2_obs1_rot", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_obs2_rot = new TProfile("p_v2e_Q2_obs2_rot", "p_v2e_Q2_obs2_rot", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_obs1_rot = new TProfile("p_v2w_Q2_obs1_rot", "p_v2w_Q2_obs1_rot", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_obs2_rot = new TProfile("p_v2w_Q2_obs2_rot", "p_v2w_Q2_obs2_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_P_obs1_rot = new TProfile("p_v2e_Q2_P_obs1_rot", "p_v2e_Q2_P_obs1_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_P_obs2_rot = new TProfile("p_v2e_Q2_P_obs2_rot", "p_v2e_Q2_P_obs2_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_P_obs1_rot = new TProfile("p_v2w_Q2_P_obs1_rot", "p_v2w_Q2_P_obs1_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_P_obs2_rot = new TProfile("p_v2w_Q2_P_obs2_rot", "p_v2w_Q2_P_obs2_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_N_obs1_rot = new TProfile("p_v2e_Q2_N_obs1_rot", "p_v2e_Q2_N_obs1_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2e_Q2_N_obs2_rot = new TProfile("p_v2e_Q2_N_obs2_rot", "p_v2e_Q2_N_obs2_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_N_obs1_rot = new TProfile("p_v2w_Q2_N_obs1_rot", "p_v2w_Q2_N_obs1_rot", 250, 0, 25, -100, 100, "");
// TProfile *p_v2w_Q2_N_obs2_rot = new TProfile("p_v2w_Q2_N_obs2_rot", "p_v2w_Q2_N_obs2_rot", 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_Q2_obs1 = new TProfile2D("Parity_e_Q2_obs1", "Parity_e_Q2_obs1", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_Q2_obs2 = new TProfile2D("Parity_e_Q2_obs2", "Parity_e_Q2_obs2", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_Q2_obs1 = new TProfile2D("Parity_w_Q2_obs1", "Parity_w_Q2_obs1", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_Q2_obs2 = new TProfile2D("Parity_w_Q2_obs2", "Parity_w_Q2_obs2", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_Q2_obs1 = new TProfile2D("pDelta_Q2_obs1", "pDelta_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_Q2_obs2 = new TProfile2D("pDelta_Q2_obs2", "pDelta_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_Q2_obs1_rot = new TProfile2D("Parity_e_Q2_obs1_rot", "Parity_e_Q2_obs1_rot", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_Q2_obs2_rot = new TProfile2D("Parity_e_Q2_obs2_rot", "Parity_e_Q2_obs2_rot", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_Q2_obs1_rot = new TProfile2D("Parity_w_Q2_obs1_rot", "Parity_w_Q2_obs1_rot", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_Q2_obs2_rot = new TProfile2D("Parity_w_Q2_obs2_rot", "Parity_w_Q2_obs2_rot", 8, 0.5, 8.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_Q2_obs1_rot = new TProfile2D("pDelta_Q2_obs1_rot", "pDelta_Q2_obs1_rot", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_Q2_obs2_rot = new TProfile2D("pDelta_Q2_obs2_rot", "pDelta_Q2_obs2_rot", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_obs1 = new TProfile("p_v2e_noHBT_Q2_obs1", "p_v2e_noHBT_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_obs2 = new TProfile("p_v2e_noHBT_Q2_obs2", "p_v2e_noHBT_Q2_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_obs1 = new TProfile("p_v2w_noHBT_Q2_obs1", "p_v2w_noHBT_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_obs2 = new TProfile("p_v2w_noHBT_Q2_obs2", "p_v2w_noHBT_Q2_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_P_obs1 = new TProfile("p_v2e_noHBT_Q2_P_obs1", "p_v2e_noHBT_Q2_P_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_P_obs2 = new TProfile("p_v2e_noHBT_Q2_P_obs2", "p_v2e_noHBT_Q2_P_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_P_obs1 = new TProfile("p_v2w_noHBT_Q2_P_obs1", "p_v2w_noHBT_Q2_P_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_P_obs2 = new TProfile("p_v2w_noHBT_Q2_P_obs2", "p_v2w_noHBT_Q2_P_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_N_obs1 = new TProfile("p_v2e_noHBT_Q2_N_obs1", "p_v2e_noHBT_Q2_N_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_N_obs2 = new TProfile("p_v2e_noHBT_Q2_N_obs2", "p_v2e_noHBT_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_N_obs1 = new TProfile("p_v2w_noHBT_Q2_N_obs1", "p_v2w_noHBT_Q2_N_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_N_obs2 = new TProfile("p_v2w_noHBT_Q2_N_obs2", "p_v2w_noHBT_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_noHBT_Q2_obs1 = new TProfile2D("Parity_e_noHBT_Q2_obs1", "Parity_e_noHBT_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_noHBT_Q2_obs2 = new TProfile2D("Parity_e_noHBT_Q2_obs2", "Parity_e_noHBT_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_noHBT_Q2_obs1 = new TProfile2D("Parity_w_noHBT_Q2_obs1", "Parity_w_noHBT_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_noHBT_Q2_obs2 = new TProfile2D("Parity_w_noHBT_Q2_obs2", "Parity_w_noHBT_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_noHBT_Q2_obs1 = new TProfile2D("pDelta_noHBT_Q2_obs1", "pDelta_noHBT_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_noHBT_Q2_obs2 = new TProfile2D("pDelta_noHBT_Q2_obs2", "pDelta_noHBT_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");

TH1D *V0Mass_gamma112 = new TH1D("V0Mass_gamma112", "V0 Inv. Mass in Gamma112", 200, 1.115684 - 0.07, 1.115684 + 0.07);
TH1D *V0Mass_gamma112_rot = new TH1D("V0Mass_gamma112_rot", "V0 Inv. Mass in Gamma112 (BKG)", 200, 1.115684 - 0.07, 1.115684 + 0.07);

TH1D *Hist_charge1_dist = new TH1D("Hist_charge1_dist", "Hist_charge1_dist", 10, -5, 5);
TH1D *Hist_charge2_dist = new TH1D("Hist_charge2_dist", "Hist_charge2_dist", 10, -5, 5);
TH1D *Hist_charge1_ss_dist = new TH1D("Hist_charge1_ss_dist", "Hist_charge1_ss_dist", 10, -5, 5);
TH1D *Hist_charge2_ss_dist = new TH1D("Hist_charge2_ss_dist", "Hist_charge2_ss_dist", 10, -5, 5);
TH1D *Hist_charge1_os_dist = new TH1D("Hist_charge1_os_dist", "Hist_charge1_os_dist", 10, -5, 5);
TH1D *Hist_charge2_os_dist = new TH1D("Hist_charge2_os_dist", "Hist_charge2_os_dist", 10, -5, 5);

TH1I *num_gamma_ss_pre = new TH1I("num_gamma_ss_pre", "num_gamma_ss_pre", 200, 0, 200);
TH1I *num_gamma_os_pre = new TH1I("num_gamma_os_pre", "num_gamma_os_pre", 200, 0, 200);
TH1I *num_gamma_ss_final = new TH1I("num_gamma_ss_final", "num_gamma_ss_final", 200, 0, 200);
TH1I *num_gamma_os_final = new TH1I("num_gamma_os_final", "num_gamma_os_final", 200, 0, 200);

TH1D *nHitsFitQA = new TH1D("nHitsFitQA", "nHitsFitQA", 100, 0, 100);
TH1D *nSigma_dauQA = new TH1D("nSigma_dauQA", "nSigma_dauQA", 100, 0, 10);
TH1D *dca_protonQA = new TH1D("dca_protonQA", "dca_protonQA", 100, 0, 10);
TH1D *dca_pionQA = new TH1D("dca_pionQA", "dca_pionQA", 100, 0, 10);
TH1D *dca_LambdaQA = new TH1D("dca_LambdaQA", "dca_LambdaQA", 100, 0, 10);
TH1D *nSigma_prim_protonQA = new TH1D("nSigma_prim_protonQA", "nSigma_prim_protonQA", 100, 0, 10);

TProfile *pTemp_v2_parent = new TProfile("pTemp_v2_parent", "pTemp_v2_parent", 7, 0.5, 7.5, -100, 100, "");
TProfile *pTemp_v2_parent_rot = new TProfile("pTemp_v2_parent_rot", "pTemp_v2_parent_rot", 7, 0.5, 7.5, -100, 100, "");
TH1D *Hist_Q2_parent = new TH1D("Hist_Q2_parent", "Q2 parent", 10000, 0, 50);
TH1D *Hist_Q2_parent_remove1 = new TH1D("Hist_Q2_parent_remove1", "Q2 parent (remove 1)", 10000, 0, 50);
TH1D *Hist_Q2_parent_pair_ss = new TH1D("Hist_Q2_parent_pair_ss", "Hist_Q2_parent_pair_ss", 2000, 0, 50);
TH1D *Hist_Q2_parent_pair_os = new TH1D("Hist_Q2_parent_pair_os", "Hist_Q2_parent_pair_os", 2000, 0, 50);

TProfile *Hist_Q2_parent_vs_Qcount_parent = new TProfile("Hist_Q2_parent_vs_Qcount_parent", "Hist_Q2_parent_vs_Qcount_parent", 400, 0, 10, 0, 25);

TH1D *Hist_Q2_parent_tp = new TH1D("Hist_Q2_parent_tp", "Q2 parent_tp", 500, 0, 50);
TH1D *Hist_Q2_parent_p = new TH1D("Hist_Q2_parent_p", "Q2 parent_p", 500, 0, 50);
TH1D *Hist_Q2_parent_diff = new TH1D("Hist_Q2_parent_diff", "Q2 parent_diff", 500, 0, 50);
TH1D *Hist_Q2_parent_ap = new TH1D("Hist_Q2_parent_ap", "Q2 parent_ap", 500, 0, 50);
TH1D *Hist_Q2_parent_lam = new TH1D("Hist_Q2_parent_lam", "Q2 parent_lam", 500, 0, 50);
TH1I *Hist_ngamma_wrange = new TH1I("Hist_ngamma_wrange", "Hist_ngamma_wrange", 300, 0, 300);
TH1I *Hist_nlam_wrange = new TH1I("Hist_nlam_wrange", "Hist_nlam_wrange", 300, 0, 300);
TH1I *Hist_np_wrange = new TH1I("Hist_np_wrange", "Hist_np_wrange", 300, 0, 300);
TH1I *Hist_nap_wrange = new TH1I("Hist_nap_wrange", "Hist_nap_wrange", 300, 0, 300);
TH1I *Hist_ntp_wrange = new TH1I("Hist_ntp_wrange", "Hist_ntp_wrange", 300, 0, 300);
TH2D *Hist_check_trksplitting_flowvector_ss = new TH2D("Hist_check_trksplitting_flowvector_ss", "Hist_check_trksplitting_flowvector_ss", 100, -1, 1, 100, -1, 1);
TH2D *Hist_check_trksplitting_flowvector_os = new TH2D("Hist_check_trksplitting_flowvector_os", "Hist_check_trksplitting_flowvector_os", 100, -1, 1, 100, -1, 1);

TH1D *Hist_Q2_parent_rot = new TH1D("Hist_Q2_parent_rot", "Q2 parent_rot", 10000, 0, 50);
TH1D *Hist_Q2_parent_QQcut = new TH1D("Hist_Q2_parent_QQcut", "Q2 parent_QQcut", 10000, 0, 50);
TH1D *Hist_Q2_parent_QQcut_rot = new TH1D("Hist_Q2_parent_QQcut_rot", "Q2 parent_QQcut_rot", 10000, 0, 50);
TH2D *Hist_RefMult_Q2_parent = new TH2D("Hist_RefMult_Q2_parent", "Hist_RefMult_Q2_parent", 1000, 0, 100, 1000, 0, 1000);
TH2D *Hist_RefMult_Q2_parent_rot = new TH2D("Hist_RefMult_Q2_parent_rot", "Hist_RefMult_Q2_parent_rot", 1000, 0, 100, 1000, 0, 1000);
TProfile *p_cos_Q2_parent = new TProfile("p_cos_Q2_parent", "p_cos_Q2_parent", 500, 0, 50, -1, 1, "");
TProfile *p_cos_Q2_parent_rot = new TProfile("p_cos_Q2_parent_rot", "p_cos_Q2_parent_rot", 500, 0, 50, -1, 1, "");

TProfile *p_Parity_v2e_parent_obs1 = new TProfile("p_Parity_v2e_parent_obs1", "Parity3e, v2e parent", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs2 = new TProfile("p_Parity_v2e_parent_obs2", "Parity3e, v2e parent, eff", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs3 = new TProfile("p_Parity_v2e_parent_obs3", "Parity4e, v2e parent", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs4 = new TProfile("p_Parity_v2e_parent_obs4", "Parity4e, v2e parent, eff", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs1 = new TProfile("p_Parity_v2w_parent_obs1", "Parity3w, v2e parent", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs2 = new TProfile("p_Parity_v2w_parent_obs2", "Parity3w, v2e parent, eff", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs3 = new TProfile("p_Parity_v2w_parent_obs3", "Parity4w, v2e parent", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs4 = new TProfile("p_Parity_v2w_parent_obs4", "Parity4w, v2e parent, eff", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs1_rot = new TProfile("p_Parity_v2e_parent_obs1_rot", "Parity3e, v2e parent (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs2_rot = new TProfile("p_Parity_v2e_parent_obs2_rot", "Parity3e, v2e parent, eff (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs3_rot = new TProfile("p_Parity_v2e_parent_obs3_rot", "Parity4e, v2e parent (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_obs4_rot = new TProfile("p_Parity_v2e_parent_obs4_rot", "Parity4e, v2e parent, eff (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs1_rot = new TProfile("p_Parity_v2w_parent_obs1_rot", "Parity3w, v2e parent (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs2_rot = new TProfile("p_Parity_v2w_parent_obs2_rot", "Parity3w, v2e parent, eff (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs3_rot = new TProfile("p_Parity_v2w_parent_obs3_rot", "Parity4w, v2e parent (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_obs4_rot = new TProfile("p_Parity_v2w_parent_obs4_rot", "Parity4w, v2e parent, eff (bkg)", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs1 = new TProfile("p_Parity_v2e_parent_QQcut_obs1", "Parity3e, v2e parent_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs2 = new TProfile("p_Parity_v2e_parent_QQcut_obs2", "Parity3e, v2e parent, eff_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs3 = new TProfile("p_Parity_v2e_parent_QQcut_obs3", "Parity4e, v2e parent_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs4 = new TProfile("p_Parity_v2e_parent_QQcut_obs4", "Parity4e, v2e parent, eff_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs1 = new TProfile("p_Parity_v2w_parent_QQcut_obs1", "Parity3w, v2e parent_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs2 = new TProfile("p_Parity_v2w_parent_QQcut_obs2", "Parity3w, v2e parent, eff_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs3 = new TProfile("p_Parity_v2w_parent_QQcut_obs3", "Parity4w, v2e parent_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs4 = new TProfile("p_Parity_v2w_parent_QQcut_obs4", "Parity4w, v2e parent, eff_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs1_rot = new TProfile("p_Parity_v2e_parent_QQcut_obs1_rot", "Parity3e, v2e parent (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs2_rot = new TProfile("p_Parity_v2e_parent_QQcut_obs2_rot", "Parity3e, v2e parent, eff (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs3_rot = new TProfile("p_Parity_v2e_parent_QQcut_obs3_rot", "Parity4e, v2e parent (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_QQcut_obs4_rot = new TProfile("p_Parity_v2e_parent_QQcut_obs4_rot", "Parity4e, v2e parent, eff (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs1_rot = new TProfile("p_Parity_v2w_parent_QQcut_obs1_rot", "Parity3w, v2e parent (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs2_rot = new TProfile("p_Parity_v2w_parent_QQcut_obs2_rot", "Parity3w, v2e parent, eff (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs3_rot = new TProfile("p_Parity_v2w_parent_QQcut_obs3_rot", "Parity4w, v2e parent (bkg)_QQcut", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_QQcut_obs4_rot = new TProfile("p_Parity_v2w_parent_QQcut_obs4_rot", "Parity4w, v2e parent, eff (bkg)_QQcut", 250, -25, 25, -100, 100, "");

TProfile *p_Parity_v2e_parent_noHBT_obs1 = new TProfile("p_Parity_v2e_parent_noHBT_obs1", "p_Parity_v2e_parent_noHBT_obs1", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_noHBT_obs2 = new TProfile("p_Parity_v2e_parent_noHBT_obs2", "p_Parity_v2e_parent_noHBT_obs2", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_noHBT_obs3 = new TProfile("p_Parity_v2e_parent_noHBT_obs3", "p_Parity_v2e_parent_noHBT_obs3", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2e_parent_noHBT_obs4 = new TProfile("p_Parity_v2e_parent_noHBT_obs4", "p_Parity_v2e_parent_noHBT_obs4", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_noHBT_obs1 = new TProfile("p_Parity_v2w_parent_noHBT_obs1", "p_Parity_v2w_parent_noHBT_obs1", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_noHBT_obs2 = new TProfile("p_Parity_v2w_parent_noHBT_obs2", "p_Parity_v2w_parent_noHBT_obs2", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_noHBT_obs3 = new TProfile("p_Parity_v2w_parent_noHBT_obs3", "p_Parity_v2w_parent_noHBT_obs3", 250, -25, 25, -100, 100, "");
TProfile *p_Parity_v2w_parent_noHBT_obs4 = new TProfile("p_Parity_v2w_parent_noHBT_obs4", "p_Parity_v2w_parent_noHBT_obs4", 250, -25, 25, -100, 100, "");
TProfile *ppv2_parent_CosFW = new TProfile("ppv2_parent_CosFW", "pair by pair ppv2_parent_CosFW", 200, -100, 100, -1, 1, "");
TProfile *pv2parent_CosFW = new TProfile("pv2parent_CosFW", "event by event ppv2_parent_CosFW", 250, -25, 25, -1, 1, "");
TProfile *pv2parent_CosFW_rot = new TProfile("pv2parent_CosFW_rot", "event by event ppv2_parent_CosFW (bkg)", 250, -25, 25, -1, 1, "");
TProfile *p_v2parente_Q2parent_obs1 = new TProfile("p_v2parente_Q2parent_obs1", "p_v2parente_Q2parent_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_obs2 = new TProfile("p_v2parente_Q2parent_obs2", "p_v2parente_Q2parent_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_obs1 = new TProfile("p_v2parentw_Q2parent_obs1", "p_v2parentw_Q2parent_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_obs2 = new TProfile("p_v2parentw_Q2parent_obs2", "p_v2parentw_Q2parent_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2parent_obs1 = new TProfile("p_v2e_Q2parent_obs1", "p_v2e_Q2parent_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2parent_obs2 = new TProfile("p_v2e_Q2parent_obs2", "p_v2e_Q2parent_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2parent_obs1 = new TProfile("p_v2w_Q2parent_obs1", "p_v2w_Q2parent_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2parent_obs2 = new TProfile("p_v2w_Q2parent_obs2", "p_v2w_Q2parent_obs2", 500, 0, 50, -100, 100, "");

TProfile *p_v2e_Q2_obs1 = new TProfile("p_v2e_Q2_obs1", "p_v2e_Q2_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2_obs2 = new TProfile("p_v2e_Q2_obs2", "p_v2e_Q2_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2_obs1 = new TProfile("p_v2w_Q2_obs1", "p_v2w_Q2_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2_obs2 = new TProfile("p_v2w_Q2_obs2", "p_v2w_Q2_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2_EPD_obs1 = new TProfile("p_v2e_Q2_EPD_obs1", "p_v2e_Q2_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2_EPD_obs2 = new TProfile("p_v2e_Q2_EPD_obs2", "p_v2e_Q2_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2_EPD_obs1 = new TProfile("p_v2w_Q2_EPD_obs1", "p_v2w_Q2_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2_EPD_obs2 = new TProfile("p_v2w_Q2_EPD_obs2", "p_v2w_Q2_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2_EPD1_obs1 = new TProfile("p_v2e_Q2_EPD1_obs1", "p_v2e_Q2_EPD1_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2_EPD1_obs2 = new TProfile("p_v2e_Q2_EPD1_obs2", "p_v2e_Q2_EPD1_obs2", 500, 0, 50, -100, 100, "");

TProfile *p_pionv2e_Q2_obs1 = new TProfile("p_pionv2e_Q2_obs1", "p_pionv2e_Q2_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2e_Q2_obs2 = new TProfile("p_pionv2e_Q2_obs2", "p_pionv2e_Q2_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2w_Q2_obs1 = new TProfile("p_pionv2w_Q2_obs1", "p_pionv2w_Q2_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2w_Q2_obs2 = new TProfile("p_pionv2w_Q2_obs2", "p_pionv2w_Q2_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2e_Q2_EPD_obs1 = new TProfile("p_pionv2e_Q2_EPD_obs1", "p_pionv2e_Q2_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2e_Q2_EPD_obs2 = new TProfile("p_pionv2e_Q2_EPD_obs2", "p_pionv2e_Q2_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2w_Q2_EPD_obs1 = new TProfile("p_pionv2w_Q2_EPD_obs1", "p_pionv2w_Q2_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2w_Q2_EPD_obs2 = new TProfile("p_pionv2w_Q2_EPD_obs2", "p_pionv2w_Q2_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2e_Q2_EPD1_obs1 = new TProfile("p_pionv2e_Q2_EPD1_obs1", "p_pionv2e_Q2_EPD1_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_pionv2e_Q2_EPD1_obs2 = new TProfile("p_pionv2e_Q2_EPD1_obs2", "p_pionv2e_Q2_EPD1_obs2", 500, 0, 50, -100, 100, "");

TProfile *p_v2parente_Q2parent_EPD_obs1 = new TProfile("p_v2parente_Q2parent_EPD_obs1", "p_v2parente_Q2parent_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_EPD_obs2 = new TProfile("p_v2parente_Q2parent_EPD_obs2", "p_v2parente_Q2parent_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD_obs1 = new TProfile("p_v2parentw_Q2parent_EPD_obs1", "p_v2parentw_Q2parent_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD_obs2 = new TProfile("p_v2parentw_Q2parent_EPD_obs2", "p_v2parentw_Q2parent_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_EPD1_obs1 = new TProfile("p_v2parente_Q2parent_EPD1_obs1", "p_v2parente_Q2parent_EPD1_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_EPD1_obs2 = new TProfile("p_v2parente_Q2parent_EPD1_obs2", "p_v2parente_Q2parent_EPD1_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD1_obs1 = new TProfile("p_v2parentw_Q2parent_EPD1_obs1", "p_v2parentw_Q2parent_EPD1_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD1_obs2 = new TProfile("p_v2parentw_Q2parent_EPD1_obs2", "p_v2parentw_Q2parent_EPD1_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2parent_EPD_obs1 = new TProfile("p_v2e_Q2parent_EPD_obs1", "p_v2e_Q2parent_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2parent_EPD_obs2 = new TProfile("p_v2e_Q2parent_EPD_obs2", "p_v2e_Q2parent_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2parent_EPD_obs1 = new TProfile("p_v2w_Q2parent_EPD_obs1", "p_v2w_Q2parent_EPD_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2parent_EPD_obs2 = new TProfile("p_v2w_Q2parent_EPD_obs2", "p_v2w_Q2parent_EPD_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2parent_EPD1_obs1 = new TProfile("p_v2e_Q2parent_EPD1_obs1", "p_v2e_Q2parent_EPD1_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2e_Q2parent_EPD1_obs2 = new TProfile("p_v2e_Q2parent_EPD1_obs2", "p_v2e_Q2parent_EPD1_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2parent_EPD1_obs1 = new TProfile("p_v2w_Q2parent_EPD1_obs1", "p_v2w_Q2parent_EPD1_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2w_Q2parent_EPD1_obs2 = new TProfile("p_v2w_Q2parent_EPD1_obs2", "p_v2w_Q2parent_EPD1_obs2", 500, 0, 50, -100, 100, "");

TProfile *p_v2parente_Q2parent_EPD_obs1_rot = new TProfile("p_v2parente_Q2parent_EPD_obs1_rot", "p_v2parente_Q2parent_EPD_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_EPD_obs2_rot = new TProfile("p_v2parente_Q2parent_EPD_obs2_rot", "p_v2parente_Q2parent_EPD_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD_obs1_rot = new TProfile("p_v2parentw_Q2parent_EPD_obs1_rot", "p_v2parentw_Q2parent_EPD_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD_obs2_rot = new TProfile("p_v2parentw_Q2parent_EPD_obs2_rot", "p_v2parentw_Q2parent_EPD_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_EPD1_obs1_rot = new TProfile("p_v2parente_Q2parent_EPD1_obs1_rot", "p_v2parente_Q2parent_EPD1_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_EPD1_obs2_rot = new TProfile("p_v2parente_Q2parent_EPD1_obs2_rot", "p_v2parente_Q2parent_EPD1_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD1_obs1_rot = new TProfile("p_v2parentw_Q2parent_EPD1_obs1_rot", "p_v2parentw_Q2parent_EPD1_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_EPD1_obs2_rot = new TProfile("p_v2parentw_Q2parent_EPD1_obs2_rot", "p_v2parentw_Q2parent_EPD1_obs2_rot", 500, 0, 50, -100, 100, "");

TProfile2D *pParity_e_Q2parent_obs1 = new TProfile2D("pParity_e_Q2parent_obs1", "pParity_e_Q2parent_obs1", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_obs2 = new TProfile2D("pParity_e_Q2parent_obs2", "pParity_e_Q2parent_obs2", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_obs1 = new TProfile2D("pParity_w_Q2parent_obs1", "pParity_w_Q2parent_obs1", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_obs2 = new TProfile2D("pParity_w_Q2parent_obs2", "pParity_w_Q2parent_obs2", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_obs1_rot = new TProfile("p_v2parente_Q2parent_obs1_rot", "p_v2parente_Q2parent_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_obs2_rot = new TProfile("p_v2parente_Q2parent_obs2_rot", "p_v2parente_Q2parent_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_obs1_rot = new TProfile("p_v2parentw_Q2parent_obs1_rot", "p_v2parentw_Q2parent_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_obs2_rot = new TProfile("p_v2parentw_Q2parent_obs2_rot", "p_v2parentw_Q2parent_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_obs1_rot = new TProfile2D("pParity_e_Q2parent_obs1_rot", "pParity_e_Q2parent_obs1_rot", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_obs2_rot = new TProfile2D("pParity_e_Q2parent_obs2_rot", "pParity_e_Q2parent_obs2_rot", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_obs1_rot = new TProfile2D("pParity_w_Q2parent_obs1_rot", "pParity_w_Q2parent_obs1_rot", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_obs2_rot = new TProfile2D("pParity_w_Q2parent_obs2_rot", "pParity_w_Q2parent_obs2_rot", 68, 0.5, 68.5, 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_QQcut_obs1 = new TProfile("p_v2parente_Q2parent_QQcut_obs1", "p_v2parente_Q2parent_QQcut_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_QQcut_obs2 = new TProfile("p_v2parente_Q2parent_QQcut_obs2", "p_v2parente_Q2parent_QQcut_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_QQcut_obs1 = new TProfile("p_v2parentw_Q2parent_QQcut_obs1", "p_v2parentw_Q2parent_QQcut_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_QQcut_obs2 = new TProfile("p_v2parentw_Q2parent_QQcut_obs2", "p_v2parentw_Q2parent_QQcut_obs2", 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_QQcut_obs1 = new TProfile2D("pParity_e_Q2parent_QQcut_obs1", "pParity_e_Q2parent_QQcut_obs1", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_QQcut_obs2 = new TProfile2D("pParity_e_Q2parent_QQcut_obs2", "pParity_e_Q2parent_QQcut_obs2", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_QQcut_obs1 = new TProfile2D("pParity_w_Q2parent_QQcut_obs1", "pParity_w_Q2parent_QQcut_obs1", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_QQcut_obs2 = new TProfile2D("pParity_w_Q2parent_QQcut_obs2", "pParity_w_Q2parent_QQcut_obs2", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_QQcut_obs1_rot = new TProfile("p_v2parente_Q2parent_QQcut_obs1_rot", "p_v2parente_Q2parent_QQcut_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_QQcut_obs2_rot = new TProfile("p_v2parente_Q2parent_QQcut_obs2_rot", "p_v2parente_Q2parent_QQcut_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_QQcut_obs1_rot = new TProfile("p_v2parentw_Q2parent_QQcut_obs1_rot", "p_v2parentw_Q2parent_QQcut_obs1_rot", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_QQcut_obs2_rot = new TProfile("p_v2parentw_Q2parent_QQcut_obs2_rot", "p_v2parentw_Q2parent_QQcut_obs2_rot", 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_QQcut_obs1_rot = new TProfile2D("pParity_e_Q2parent_QQcut_obs1_rot", "pParity_e_Q2parent_QQcut_obs1_rot", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_QQcut_obs2_rot = new TProfile2D("pParity_e_Q2parent_QQcut_obs2_rot", "pParity_e_Q2parent_QQcut_obs2_rot", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_QQcut_obs1_rot = new TProfile2D("pParity_w_Q2parent_QQcut_obs1_rot", "pParity_w_Q2parent_QQcut_obs1_rot", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_QQcut_obs2_rot = new TProfile2D("pParity_w_Q2parent_QQcut_obs2_rot", "pParity_w_Q2parent_QQcut_obs2_rot", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_noHBT_obs1 = new TProfile("p_v2parente_Q2parent_noHBT_obs1", "p_v2parente_Q2parent_noHBT_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parente_Q2parent_noHBT_obs2 = new TProfile("p_v2parente_Q2parent_noHBT_obs2", "p_v2parente_Q2parent_noHBT_obs2", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_noHBT_obs1 = new TProfile("p_v2parentw_Q2parent_noHBT_obs1", "p_v2parentw_Q2parent_noHBT_obs1", 500, 0, 50, -100, 100, "");
TProfile *p_v2parentw_Q2parent_noHBT_obs2 = new TProfile("p_v2parentw_Q2parent_noHBT_obs2", "p_v2parentw_Q2parent_noHBT_obs2", 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_noHBT_obs1 = new TProfile2D("pParity_e_Q2parent_noHBT_obs1", "pParity_e_Q2parent_noHBT_obs1", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_e_Q2parent_noHBT_obs2 = new TProfile2D("pParity_e_Q2parent_noHBT_obs2", "pParity_e_Q2parent_noHBT_obs2", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_noHBT_obs1 = new TProfile2D("pParity_w_Q2parent_noHBT_obs1", "pParity_w_Q2parent_noHBT_obs1", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile2D *pParity_w_Q2parent_noHBT_obs2 = new TProfile2D("pParity_w_Q2parent_noHBT_obs2", "pParity_w_Q2parent_noHBT_obs2", 4, 0.5, 4.5, 500, 0, 50, -100, 100, "");
TProfile *Hist_v2_v2parent_1 = new TProfile("Hist_v2_v2parent_1", "v2 parent, v2 correlation inclusive", 200, -100, 100, -100, 100);
TProfile *Hist_v2_v2parent_2 = new TProfile("Hist_v2_v2parent_2", "v2 parent, v2 correlation ESE cut", 200, -100, 100, -100, 100);
TProfile *Hist_v2_v2parent_3 = new TProfile("Hist_v2_v2parent_3", "v2 parent, v2 correlation ESE event averaged", 200, -100, 100, -100, 100);
TProfile *Hist_v2_v2parent_3_rot = new TProfile("Hist_v2_v2parent_3_rot", "v2 parent, v2 correlation ESE event averaged (bkg)", 200, -100, 100, -100, 100);

TH1D *Hist_parentv2_event = new TH1D("Hist_parentv2_event", "Hist_parentv2_event", 200, -100, 100);
TH1D *Hist_parentv2_event_rot = new TH1D("Hist_parentv2_event_rot", "Hist_parentv2_event_rot", 200, -100, 100);

TProfile *Hist_v2parent_eta_obs5 = new TProfile("Hist_v2parent_eta_obs5", "v2 parent TPC", 200, -10, 10, -100, 100, "");
// TProfile *Hist_v2parent_pt_obs5 = new TProfile("Hist_v2parent_pt_obs5", "v2 parent TPC", 50, 0, 5, -100, 100, "");
TProfile *Hist_v2parent_pt_obs5[3];
TProfile *Hist_v2parent_eta_obs5_rot = new TProfile("Hist_v2parent_eta_obs5_rot", "v2 parent TPC (bkg)", 200, -10, 10, -100, 100, "");
// TProfile *Hist_v2parent_pt_obs5_rot = new TProfile("Hist_v2parent_pt_obs5_rot", "v2 parent TPC (bkg)", 50, 0, 5, -100, 100, "");
TProfile *Hist_v2parent_pt_obs5_rot[3];
TProfile *Hist_v2pion_pt_obs5[3];
TProfile *v2_2_pion = new TProfile("v2_2_pion", "", 1, 0, 1);

TProfile *EPDe_Day3_cos1 = new TProfile("EPDe_Day3_cos1", "EPDe_Day3_cos1", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDe_Day3_sin1 = new TProfile("EPDe_Day3_sin1", "EPDe_Day3_sin1", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDw_Day3_cos1 = new TProfile("EPDw_Day3_cos1", "EPDw_Day3_cos1", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDw_Day3_sin1 = new TProfile("EPDw_Day3_sin1", "EPDw_Day3_sin1", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDe_Day3_cos2 = new TProfile("EPDe_Day3_cos2", "EPDe_Day3_cos2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDe_Day3_sin2 = new TProfile("EPDe_Day3_sin2", "EPDe_Day3_sin2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDw_Day3_cos2 = new TProfile("EPDw_Day3_cos2", "EPDw_Day3_cos2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDw_Day3_sin2 = new TProfile("EPDw_Day3_sin2", "EPDw_Day3_sin2", run_end - run_sta, run_sta, run_end, -1, 1);

TProfile *Hist_v2_pt_EPD1_obs = new TProfile("Hist_v2_pt_EPD1_obs", "Hist_v2_pt_EPD1_obs", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_EPD_obs = new TProfile("Hist_v2_pt_EPD_obs", "Hist_v2_pt_EPD_obs", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_eta_EPD1_obs = new TProfile("Hist_v2_eta_EPD1_obs", "Hist_v2_eta_EPD1_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPD_obs = new TProfile("Hist_v2_eta_EPD_obs", "Hist_v2_eta_EPD_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPDe_obs = new TProfile("Hist_v2_eta_EPDe_obs", "Hist_v2_eta_EPDe_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPDw_obs = new TProfile("Hist_v2_eta_EPDw_obs", "Hist_v2_eta_EPDw_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_pt_EPD1_obs_rot = new TProfile("Hist_v2_pt_EPD1_obs_rot", "Hist_v2_pt_EPD1_obs_rot", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_EPD_obs_rot = new TProfile("Hist_v2_pt_EPD_obs_rot", "Hist_v2_pt_EPD_obs_rot", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_eta_EPD1_obs_rot = new TProfile("Hist_v2_eta_EPD1_obs_rot", "Hist_v2_eta_EPD1_obs_rot", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPD_obs_rot = new TProfile("Hist_v2_eta_EPD_obs_rot", "Hist_v2_eta_EPD_obs_rot", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPDe_obs_rot = new TProfile("Hist_v2_eta_EPDe_obs_rot", "Hist_v2_eta_EPDe_obs_rot", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPDw_obs_rot = new TProfile("Hist_v2_eta_EPDw_obs_rot", "Hist_v2_eta_EPDw_obs_rot", 40, -1, 1, -100, 100, "");

TProfile *Hist_cos_EPD = new TProfile("Hist_cos_EPD", "Hist_cos_EPD", 7, 0.5, 7.5, -1, 1, "");

TProfile2D *Read_BBC_EP_east = 0, *Read_BBC_EP_west = 0, *Read_EPD_EP_east = 0, *Read_EPD_EP_west = 0, *Read_EPD_EP1_east = 0, *Read_EPD_EP1_west = 0;

TH2D *Hist_EPD_EP1_east = new TH2D("Hist_EPD_EP1_east", "Hist_EPD_EP1_east", 72, 0, 2 * PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP1_west = new TH2D("Hist_EPD_EP1_west", "Hist_EPD_EP1_west", 72, 0, 2 * PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_east = new TH2D("Hist_EPD_EP_east", "Hist_EPD_EP_east", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_west = new TH2D("Hist_EPD_EP_west", "Hist_EPD_EP_west", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP1_east_flat = new TH2D("Hist_EPD_EP1_east_flat", "Hist_EPD_EP1_east_flat", 72, -PI, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP1_west_flat = new TH2D("Hist_EPD_EP1_west_flat", "Hist_EPD_EP1_west_flat", 72, -PI, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_east_flat = new TH2D("Hist_EPD_EP_east_flat", "Hist_EPD_EP_east_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_west_flat = new TH2D("Hist_EPD_EP_west_flat", "Hist_EPD_EP_west_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TProfile2D *pEPD_EP1_east = new TProfile2D("pEPD_EP1_east", "pEPD_EP1_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pEPD_EP1_west = new TProfile2D("pEPD_EP1_west", "pEPD_EP1_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pEPD_EP_east = new TProfile2D("pEPD_EP_east", "pEPD_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pEPD_EP_west = new TProfile2D("pEPD_EP_west", "pEPD_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");

TH2D *Hist_EPD_vs_TPC = new TH2D("Hist_EPD_vs_TPC", "Hist_EPD_vs_TPC", 30, 0, PI, 30, 0, PI);

TProfile *EPD1_cor_Day3 = new TProfile("EPD1_cor_Day3", "EPD1_cor_Day3", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPD_cor_Day3 = new TProfile("EPD_cor_Day3", "EPD_cor_Day3", run_end - run_sta, run_sta, run_end, -1, 1);

TProfile *v2_single_lamda = new TProfile("v2_single_lamda", "", 3, 0, 3);
TProfile *v2_single_proton = new TProfile("v2_single_proton", "", 3, 0, 3);
TProfile *v2_single_lamp = new TProfile("v2_single_lamp", "", 3, 0, 3);
TProfile *v2_single_lamda_err = new TProfile("v2_single_lamda_err", "", 3, 0, 3);
TProfile *v2_single_proton_err = new TProfile("v2_single_proton_err", "", 3, 0, 3);
TProfile *v2_single_lamp_err = new TProfile("v2_single_lamp_err", "", 3, 0, 3);

TH1F *Hist_V0_Mass, *Hist_V0_Mass_anti;
// TH1F *V0Mass_Q2[9][20][3], *V0Mass_Q2_anti[9][20][3];
TH2F *V0Mass_Q2[3], *V0Mass_Q2_anti[3], *V0Mass_Q2_pion[3], *V0Mass_Q2_pion_anti[3];

// defining variables
float pVx, pVy, pVz, VPDvz, BBCco, ZDCcoin, net_Nch_Asym, mQx, mQy, mQx1, mQy1, mQx2, mQy2;                                                                                                     // run, event info
int Run, Day, Day2, Day3, Trigger, RefMult, TOFMult, Centrality, NPTracks, n_particle1, n_particle2, n_particle1_anti, n_particle2_anti, n_particle1_rot, n_particle1_anti_rot, Fcount, Scount; //
int Charge, Charge2, ChargeAsso;
int net_charge_asym_bin = -1;
float ndEdx, nSigma_p, nSigma_pi, DCAGlobal, Eta, Theta, Phi, Pt, eff, TOFflag, rapidity;                 // track info
float ndEdx2, nSigma_p2, nSigma_pi2, DCAGlobal2, Eta2, Theta2, Phi2, Pt2, Pt3, eff2, TOFflag2, rapidity2; // 2nd track info
float Eta_rot, Pt_rot, Charge_rot, Phi_rot, Theta_rot, eff_rot;
float DCAGlobalAsso, EtaAsso, PhiAsso, PtAsso, TOFflagAsso;
TVector3 pV;
int Fcount_parent;

float BBC_EP_east = 0, BBC_EP_west = 0, EPD_EP_east = 0, EPD_EP_west = 0, ZDC_EP_east = 0, ZDC_EP_west = 0;
float BBC_EP_east_new = 0, BBC_EP_west_new = 0, EPD_EP_east_new = 0, EPD_EP_west_new = 0, ZDC_EP_east_new = 0, ZDC_EP_west_new = 0;
float EPD_EP1_east = 0, EPD_EP1_west = 0, EPD_EP1_east_new = 0, EPD_EP1_west_new = 0;
float Psi_BBC_E[2 * order] = {0}, Psi_BBC_W[2 * order] = {0}, Psi_EPD_E[2 * order] = {0}, Psi_EPD_W[2 * order] = {0}, Psi1_EPD_E[2 * order] = {0}, Psi1_EPD_W[2 * order] = {0};

// struct StRefMultCorr refmultCorrUtil  = StRefMultCorr("refmult") ;
float Eweight = 1;
int Weight_Read = 0, EPD_Weight_Read = 0;
TH1D *TOF_eff = 0;
TProfile2D *TPCmean_FF_1 = 0, *TPCmean_RF_1 = 0, *TPCmean_FF_1_p = 0, *TPCmean_RF_1_p = 0, *TPCmean_FF_1_rot = 0, *TPCmean_RF_1_rot = 0, *TPCmeanAsso_FF_1 = 0, *TPCmeanAsso_RF_1 = 0;
TProfile2D *TPCmean_FF_2 = 0, *TPCmean_RF_2 = 0, *TPCmean_FF_2_p = 0, *TPCmean_RF_2_p = 0, *TPCmean_FF_2_rot = 0, *TPCmean_RF_2_rot = 0, *TPCmeanAsso_FF_2 = 0, *TPCmeanAsso_RF_2 = 0;
TProfile2D *TPCmean_FF_3 = 0, *TPCmean_RF_3 = 0, *TPCmean_FF_3_p = 0, *TPCmean_RF_3_p = 0, *TPCmean_FF_3_rot = 0, *TPCmean_RF_3_rot = 0, *TPCmeanAsso_FF_3 = 0, *TPCmeanAsso_RF_3 = 0;
float PhiMean_sin[order] = {0}, PhiMean_cos[order] = {0}, PhiMean_sin_p[order] = {0}, PhiMean_cos_p[order] = {0}, PhiMean_sin_rot[order] = {0}, PhiMean_cos_rot[order] = {0}, PhiMeanAsso_sin[order] = {0}, PhiMeanAsso_cos[order] = {0};
vector<float> PhiAsso_new, Phi_new, Phi_new_rot, Phi_pnew;
TProfile2D *Read_TPC_EP_full = 0, *Read_TPC_EP_east = 0, *Read_TPC_EP_west = 0, *Read_TPC_EP_for = 0, *Read_TPC_EP_bac = 0, *Read_TPC_EP_for_pos = 0, *Read_TPC_EP_for_neg = 0, *Read_TPC_EP_bac_pos = 0, *Read_TPC_EP_bac_neg = 0;
float PsiMean_F[2 * order] = {0}, PsiMean_E[2 * order] = {0}, PsiMean_W[2 * order] = {0}, PsiMean_f[2 * order] = {0}, PsiMean_b[2 * order] = {0}, PsiMean_fp[2 * order] = {0}, PsiMean_fn[2 * order] = {0}, PsiMean_bp[2 * order] = {0}, PsiMean_bn[2 * order] = {0};
float MeanNetChargeAsym, RMSNetChargeAsym, StdDevNetChargeAsym;

float TPC_EP_full = 0, TPC_EP_east = 0, TPC_EP_west = 0, TPC_EP_for = 0, TPC_EP_bac = 0, TPC_EP_for_pos = 0, TPC_EP_for_neg = 0, TPC_EP_bac_pos = 0, TPC_EP_bac_neg = 0;
float TPC_EP_full_new = 0, TPC_EP_east_new = 0, TPC_EP_west_new = 0, TPC_EP_for_new = 0, TPC_EP_bac_new = 0, TPC_EP_for_pos_new = 0, TPC_EP_for_neg_new = 0, TPC_EP_bac_pos_new = 0, TPC_EP_bac_neg_new = 0;
// float Q2 = 0, Q2_proper = 0;
float Q2_proper = 0, Q2_proper_EPD = 0, Q2_proper_EPD1 = 0;
float Q2_pion_TPC = 0, Q2_pion_EPD = 0, Q2_pion_EPD1 = 0;
float v2 = 0, v2e = 0, v2w = 0, v2epos = 0, v2eneg = 0, v2wpos = 0, v2wneg = 0, v2_sub = 0, v2_sub_ss = 0, v2_sub_os = 0, v2p = 0, v2pe = 0, v2pw = 0;
float v2_rot = 0, v2e_rot = 0, v2w_rot = 0, v2epos_rot = 0, v2eneg_rot = 0, v2wpos_rot = 0, v2wneg_rot = 0, v2_sub_rot = 0, v2_sub_ss_rot = 0, v2_sub_os_rot = 0, v2p_rot = 0, v2pe_rot = 0, v2pw_rot = 0;
// float correlator0_rot = 0, correlator3_rot = 0, correlator4_rot = 0, correlator4e_rot = 0, correlator4w_rot = 0, correlator0_alt_rot = 0;
int Ntof = 0, Npos = 0, Nneg = 0, n_ss_pairs = 0, n_os_pairs = 0, n_ss_pairs_pre = 0, n_os_pairs_pre = 0;

float mQQx_parent = 0., mQQy_parent = 0., Q2_parent = 0.;
float mQQx_parent_rot = 0., mQQy_parent_rot = 0., Q2_parent_rot = 0.;
float mQQx_parent_QQcut = 0., mQQy_parent_QQcut = 0., Q2_parent_QQcut = 0.;
float mQQx_parent_QQcut_rot = 0., mQQy_parent_QQcut_rot = 0., Q2_parent_QQcut_rot = 0.;
float temp_eff = -1, temp_proton_eff = -1;
float temp_eff_rot = -1, temp_proton_eff_rot = -1;

// float correlator0 = 0, correlator0e = 0, correlator0w = 0, correlator3 = 0, correlator4 = 0, correlator4e = 0, correlator4w = 0, correlator0_alt = 0, correlator0e_alt = 0, correlator0w_alt = 0;
// float correlator4e_EPD = 0, correlator4w_EPD = 0, correlator4_EPD = 0, correlator0e_EPD = 0, correlator0w_EPD = 0, correlator0_EPD = 0, correlator0e_alt_EPD = 0, correlator0w_alt_EPD = 0, correlator0_alt_EPD = 0;
// float correlator4_EPD1 = 0, correlator0_EPD1 = 0, correlator0_alt_EPD1 = 0;
float correlator3 = 0, correlator3_rot = 0;
float correlator0[3] = {0.}, correlator0e[3] = {0.}, correlator0w[3] = {0.}, correlator4[3] = {0.}, correlator4e[3] = {0.}, correlator4w[3] = {0.}, correlator0_alt[3] = {0.}, correlator0e_alt[3] = {0.}, correlator0w_alt[3] = {0.};
float correlator0_rot[3] = {0.}, correlator0e_rot[3] = {0.}, correlator0w_rot[3] = {0.}, correlator4_rot[3] = {0.}, correlator4e_rot[3] = {0.}, correlator4w_rot[3] = {0.}, correlator0_alt_rot[3] = {0.}, correlator0e_alt_rot[3] = {0.}, correlator0w_alt_rot[3] = {0.};

float v2_EPD1 = 0, v2_EPDe = 0, v2_EPDw = 0, v2_EPD1_rot = 0, v2_EPDe_rot = 0, v2_EPDw_rot = 0;

int n_gamma = 0;

double QQ = 0, QQ_rot = 0;

int particle_option1, particle_option2;

bool debug_1, debug_2, debug_3, debug_4;

bool isrotated;

int num_tracks_index;

ofstream debugfile;

TF1 *fit_eff_pT_final[9], *fit_tof_eff_pT_final[9];

TString fname_new;

// StEpdEpFinder *mEpFinder;
// TClonesArray *mEpdHits;

// Structure 2 - particles
std::vector<particle> *particle1 = new vector<particle>;
std::vector<particle_dau> *particle1_dau = new vector<particle_dau>;
std::vector<particle> *particle2 = new vector<particle>;
// std::vector<particle> *particle1_anti = new vector<particle>;
// std::vector<particle_dau> *particle1_anti_dau = new vector<particle_dau>;
// std::vector<particle> *particle2_anti = new vector<particle>;
std::vector<particle> *particle1_rot = new vector<particle>;
std::vector<particle_dau> *particle1_dau_rot = new vector<particle_dau>;
// std::vector<particle> *particle1_anti_rot = new vector<particle>;
// std::vector<particle_dau> *particle1_anti_dau_rot = new vector<particle_dau>;

std::vector<particle> *p_proton_buffer_tmp = new vector<particle>;
std::deque<particle> p_proton_buffer[28][100];
std::vector<TVector3> *p_pion_list = new vector<TVector3>;

std::vector<int> mQ1_list, mQ2_list;

event *details;
std::vector<particle> *p_lambda = new vector<particle>;
std::vector<particle_all> *p_all = new vector<particle_all>;
std::vector<particle_dau> *p_dau = new vector<particle_dau>;
std::vector<particle> *p_proton = new vector<particle>;
std::vector<particle> *p_pion = new vector<particle>;

TChain *corr_tree;

// Event variable to get in tree
int dt_cent, dt_num_trk, dt_nLambda, dt_nLambdaRot, dt_Run, dt_TOFMult, dt_RefMult, dt_n_Proton, dt_n_Pion, dt_EventID;
float dt_VPDvz, dt_PVtxz, dt_PVtxx, dt_PVtxy, dt_Eweight, dt_EPD_EP1_east, dt_EPD_EP1_west, dt_EPD_EP_east, dt_EPD_EP_west, dt_Magn;
int n_lam, n_dau;

// Lambda Particle to get in tree
std::vector<float> *lambda_px = new vector<float>;
std::vector<float> *lambda_py = new vector<float>;
std::vector<float> *lambda_pz = new vector<float>;
std::vector<float> *lambda_Charge = new vector<float>;
std::vector<float> *lambda_dcaglobal = new vector<float>;
std::vector<float> *lambda_nsigma = new vector<float>;
std::vector<float> *lambda_mass = new vector<float>;
std::vector<int> *lambda_trk_id = new vector<int>;
std::vector<float> *lambda_hits_ratio = new vector<float>;
std::vector<float> *lambda_nhitsfit = new vector<float>;
std::vector<float> *lambda_nhitsmax = new vector<float>;

// // Proton Particle to get in tree
// std::vector<float> *proton_px = new vector<float>;
// std::vector<float> *proton_py = new vector<float>;
// std::vector<float> *proton_pz = new vector<float>;
// std::vector<float> *proton_Charge = new vector<float>;
// std::vector<float> *proton_dcaglobal = new vector<float>;
// std::vector<float> *proton_nsigma = new vector<float>;
// std::vector<float> *proton_mass = new vector<float>;
// std::vector<int> *proton_trk_id = new vector<int>;
// std::vector<float> *proton_hits_ratio = new vector<float>;
// std::vector<float> *proton_nhitsfit = new vector<float>;
// std::vector<float> *proton_nhitsmax = new vector<float>;

// // Pion Particle to get in tree
// std::vector<float> *pion_px = new vector<float>;
// std::vector<float> *pion_py = new vector<float>;
// std::vector<float> *pion_pz = new vector<float>;
// std::vector<float> *pion_Charge = new vector<float>;
// std::vector<float> *pion_dcaglobal = new vector<float>;
// std::vector<float> *pion_nsigma = new vector<float>;
// std::vector<float> *pion_mass = new vector<float>;
// std::vector<int> *pion_trk_id = new vector<int>;
// std::vector<float> *pion_hits_ratio = new vector<float>;
// std::vector<float> *pion_nhitsfit = new vector<float>;
// std::vector<float> *pion_nhitsmax = new vector<float>;

// All Particles to get in tree
std::vector<float> *all_px = new vector<float>;
std::vector<float> *all_py = new vector<float>;
std::vector<float> *all_pz = new vector<float>;
std::vector<float> *all_Charge = new vector<float>;
std::vector<float> *all_dcaglobal = new vector<float>;
std::vector<float> *all_nSigmaProton = new vector<float>;
std::vector<float> *all_nSigmaPion = new vector<float>;
std::vector<bool> *all_isTofTrack = new vector<bool>;
std::vector<int> *all_trk_id = new vector<int>;
std::vector<bool> *all_is_pion = new vector<bool>;
std::vector<bool> *all_is_proton = new vector<bool>;
std::vector<bool> *all_is_all = new vector<bool>;
std::vector<float> *all_nhitsmax = new vector<float>;
std::vector<float> *all_nhitsfit = new vector<float>;

// Daughter Particles to get in tree
std::vector<float> *dau_px = new vector<float>;
std::vector<float> *dau_py = new vector<float>;
std::vector<float> *dau_pz = new vector<float>;
std::vector<float> *dau_dcaglobal = new vector<float>;
std::vector<float> *dau_nSigma = new vector<float>;
std::vector<int> *dau_trk_id = new vector<int>;
std::vector<float> *dau_nHitsFit = new vector<float>;
std::vector<float> *dau_nHitsMax = new vector<float>;

// Structure 3 - Event details
//  event *details = new event();

// sub-functions
////////////////////////// Setting Initial Conditions for Analysis //////////////////////////
void init_Gamma112(int cen, TString JobIDName);
void read_eff(const TString lambdatype);
void read_proton_eff();
void set_debug();
int ReadWeight(char *InFileName);    // input from weight file
int ReadEPDWeight(char *InFileName); // input from EPD weight file
void create_hists(int cen);
void getting_info_from_Trees();
////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////// Event-by-Event Analysis //////////////////////////
void construct_needed_variables();
int Make(int cen, int opt_weight, int sys_err_opt, const TString JobIDName /*, TClonesArray *mEpdHits*/);
void initialize_event();
void determine_particlelists(const int sys_err_opt);
bool IsGoodLambda(float p, float rap);
bool IsGoodPrimary(float r, float d);
bool IsGoodEvent(int c);                                        //select good events and fill event level histograms
bool CountCharge();                                 //count good tracks, charge asymmetry
void MakeTPC_EP();                        //reconstruct TPC EPs
bool IsGoodAsso(float p, float e, float d);                     //cuts on EP particles
void FillPhiAsso();                                             //particles for EP, shift parameters to make phi distribution flat
void ShiftPhiAsso(int tr);                                      //flatten the phi distribution
void ShiftPsi();                                                //flatten EPs
void FillPhiPOI();
void ShiftPhiPOI(int tr);                                       //flatten the phi distribution
void FillPhiPOI_rot();                                              //particles of interest, shift parameters to make phi distribution flat
void ShiftPhiPOI_rot(int tr);                                       //flatten the phi distribution
void FillCMW();                                                 //v2 for CMW
void Fillv2();

void WriteHistogram(int c, int o, TString JobIDName); // into result ROOT file
void WriteWeight(TString OutFileName);                // into weight file
void FillEP_resolution();                                       //Fill profiles for EP resolution
void FillPhiPOI_p();                                              //particles of interest, shift parameters to make phi distribution flat
void ShiftPhiPOI_p(int tr);                                       //flatten the phi distribution
void FillGamma(int ord);                                        //correlations
void FillGamma_rot(int ord);                                        //correlations


void WrapUpESE();
void Fillv2_rot();
void WrapUpESE_rot();

void finish_Gamma112(int cen, int opt_weight, TString JobIDName);

// bool IsGoodBBC(StPicoEvent *e);                                 //cuts on BBC ADCs
// bool IsGoodZDC(StPicoEvent *e);                                 //cuts on ZDC ADCs
// bool IsGoodPOI(float p, float e, float d);                      //cuts on particles of interest
// bool IsGoodTrack(StPicoTrack *p);                               //cuts on tracks
//  bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt);         //cuts on pions
//  bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt);         //cuts on kaons
//  bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt);       //cuts on protons
// void MakeBBC_EP(StPicoEvent *ev);                               //reconstruct BBC EPs
// void MakeZDC_EP(StPicoEvent *ev);                               //reconstruct ZDC EPs
// void CountMSC(int *iTr, float Phi_ne);                          //MSC correlator
// void WrapUpMSC();                                               //MSC
// float GetPhiInBBC(int e_w, int bbcN);                           //input est_wst(0,1),BBC PMT#
// float GetXYInZDC(int e_w, int v_h, int zdcN, int opt_raw = 0);  //input ver_hor(0,1),ZDCSMD slat#
