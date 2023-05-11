using namespace std;

/// C++ headers
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
#include <string>
#include <cstring>
#include <TLorentzVector.h>
#include <deque>
#include "TF1.h"
#include <fstream>
#include <iostream>

/// Analysis headers
#include "./namespaces/gv_gamma.h"
#include "./Gamma_112_module.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/particle.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/particle_all.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/particle_dau.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/test_KFParticle_gamma/StRoot/event.h"

////////////////////////////////////////////////////////////////// Setting Initial Conditions for Analysis //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void init_Gamma112(int cen, TString JobIDName)
{
    read_eff("lam");
    read_proton_eff();

    cout << "RUNNING CENT " << cen << endl;

    set_debug(); // setting debug status: print statements or not; (IMPORTANT) setting particle types

    delete gRandom;
    gRandom = new TRandom3(0);

    /////////////////////////////// Setting Files, Reading Inputs and Creating Outputs ///////////////////////////////

    char fname_old[200];
    sprintf(fname_old, "cen%d.weight_112_module.root", cen);
    Weight_Read = ReadWeight(fname_old);

    TString Name = "sched";
    Name.Append(JobIDName);
    Name.Append(TString::Format("cen%d.weight_112_module_new.root", cen));
    std::cout << "Name = " << Name << endl;
    fname_new = Name;
    std::cout << "fname_new = " << fname_new << endl;

    char EPD_fname_old[200], EPD_fname_new[200];
    sprintf(EPD_fname_new, "cen%d.weight_1%d%d_nop_EPD_new.root", cen, nHar - 1, nHar);
    sprintf(EPD_fname_old, "cen%d.weight_1%d%d_nop_EPD.root", cen, nHar - 1, nHar);
    EPD_Weight_Read = ReadEPDWeight(EPD_fname_old);

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // // Prepare EPD
    // // mEpdHits = new TClonesArray("StPicoEpdHit");
    // mEpFinder = new StEpdEpFinder(9, EPD_fname_new, EPD_fname_old);
    // mEpFinder->SetnMipThreshold(0.3); // recommended by EPD group
    // mEpFinder->SetMaxTileWeight(2.0); // recommended by EPD group
    // mEpFinder->SetEpdHitFormat(2);    // 2=pico
    // mEpFinder->SetEtaWeights(1, wt);  // eta weight for 1st-order EP
    // mEpFinder->SetEtaWeights(2, wt2); // eta weight for 2nd-order EP, select different eta range

    create_hists(cen);

    //    debugfile.open(TString::Format("debug_file_%s.txt", JobIDName.Data()));
}

void read_eff(const TString lambdatype)
{
    std::fstream eff_lam_file("./efficiency_final_lam.txt", std::ios_base::in);
    std::fstream eff_antilam_file("./efficiency_final_antilam.txt", std::ios_base::in);

    float a = 0;
    int count = 0;

    while (eff_lam_file >> a)
    {
        if (count % 2 == 0)
            efficiency[count / 44][(count - (count / 44) * 44) / 2] = a;
        else if (count % 2 == 1)
            efficiency_err[count / 44][(count - (count / 44) * 44) / 2] = a;
        count++;
    }

    a = 0;
    count = 0;

    while (eff_antilam_file >> a)
    {
        if (count % 2 == 0)
            efficiency_anti[count / 44][(count - (count / 44) * 44) / 2] = a;
        else if (count % 2 == 1)
            efficiency_anti_err[count / 44][(count - (count / 44) * 44) / 2] = a;
        count++;
    }
}

void read_proton_eff()
{
    std::fstream eff_proton_file("./proton_efficiency_coefficients.txt", std::ios_base::in);

    float a = 0;
    int count = 0;
    float par_sig[5] = {0.};

    while (eff_proton_file >> a)
    {
        par_sig[count % 5] = a;

        if (count % 5 == 4)
        {
            if (debug_1)
                cout << "count / 5 = " << count / 5 << endl;
            fit_eff_pT_final[count / 5] = new TF1("fit_eff_pT_final", "([0]+[3]*x+[4]*x*x)*exp(-pow([1]/x,[2]))", 0, 4);
            fit_eff_pT_final[count / 5]->SetParameters(par_sig[0], par_sig[1], par_sig[2], par_sig[3], par_sig[4]);
        }

        count++;
    }

    std::fstream tofeff_proton_file("./proton_tofefficiency_coefficients.txt", std::ios_base::in);

    a = 0;
    count = 0;

    float par_sig_tof[4] = {0.};

    while (tofeff_proton_file >> a)
    {
        par_sig_tof[count % 4] = a;

        if (count % 4 == 3)
        {
            if (debug_1)
                cout << "count / 4 = " << count / 4 << endl;
            fit_tof_eff_pT_final[count / 4] = new TF1("fit_tof_eff_pT_final", "[0] + [1]*sqrt(x) + [2]*x + [3]*pow(x,2)", 0.2, 2.5);
            fit_tof_eff_pT_final[count / 4]->SetParameters(par_sig_tof[0], par_sig_tof[1], par_sig_tof[2], par_sig_tof[3]);

            if (debug_1)
                cout << "parameters are: " << par_sig_tof[0] << ", " << par_sig_tof[1] << ", " << par_sig_tof[2] << ", " << par_sig_tof[3] << endl;
        }

        count++;
    }
}

void set_debug()
{
    debug_1 = false;
    debug_2 = false;
    debug_3 = false;
    debug_4 = false;

    particle_option1 = 0; // 0 - Lambda, 1 - proton
    particle_option2 = 0; // 0 - proton, 1 - pion
}

int ReadWeight(char *InFileName)
{
    TFile *fWgt = new TFile(InFileName, "READ");
    if (!fWgt->IsOpen())
        return 0;
    if (fWgt->IsOpen())
    {
        TOF_eff = (TH1D *)fWgt->Get("rc");
        if (TOF_eff && TOF_eff->GetEntries())
        {
            float cont = TOF_eff->GetBinContent(20);
            TOF_eff->Scale(1.25 / cont);
        }
        TPCmean_FF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_1");
        TPCmean_RF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_1");
        TPCmean_FF_1_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_1_p");
        TPCmean_RF_1_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_1_p");
        TPCmean_FF_1_rot = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_1_rot");
        TPCmean_RF_1_rot = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_1_rot");
        TPCmeanAsso_FF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_FF_1");
        TPCmeanAsso_RF_1 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_RF_1");
        TPCmean_FF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_2");
        TPCmean_RF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_2");
        TPCmean_FF_2_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_2_p");
        TPCmean_RF_2_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_2_p");
        TPCmean_FF_2_rot = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_2_rot");
        TPCmean_RF_2_rot = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_2_rot");
        TPCmeanAsso_FF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_FF_2");
        TPCmeanAsso_RF_2 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_RF_2");
        TPCmean_FF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_3");
        TPCmean_RF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_3");
        TPCmean_FF_3_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_3_p");
        TPCmean_RF_3_p = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_3_p");
        TPCmean_FF_3_rot = (TProfile2D *)fWgt->Get("TPCmeanPhi_FF_3_rot");
        TPCmean_RF_3_rot = (TProfile2D *)fWgt->Get("TPCmeanPhi_RF_3_rot");
        TPCmeanAsso_FF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_FF_3");
        TPCmeanAsso_RF_3 = (TProfile2D *)fWgt->Get("TPCmeanPhiAsso_RF_3");
        Read_TPC_EP_full = (TProfile2D *)fWgt->Get("pTPC_EP_full");
        Read_TPC_EP_east = (TProfile2D *)fWgt->Get("pTPC_EP_east");
        Read_TPC_EP_west = (TProfile2D *)fWgt->Get("pTPC_EP_west");
        Read_TPC_EP_for = (TProfile2D *)fWgt->Get("pTPC_EP_for");
        Read_TPC_EP_bac = (TProfile2D *)fWgt->Get("pTPC_EP_bac");
        Read_TPC_EP_for_pos = (TProfile2D *)fWgt->Get("pTPC_EP_for_pos");
        Read_TPC_EP_for_neg = (TProfile2D *)fWgt->Get("pTPC_EP_for_neg");
        Read_TPC_EP_bac_pos = (TProfile2D *)fWgt->Get("pTPC_EP_bac_pos");
        Read_TPC_EP_bac_neg = (TProfile2D *)fWgt->Get("pTPC_EP_bac_neg");
        cout << "Loaded: TPC/BBC/EPD/ZDC EP corrections" << endl;
        TH1D *Read_netChAsym = (TH1D *)fWgt->Get("Hist_netChAsym");
        MeanNetChargeAsym = Read_netChAsym->GetMean();
        RMSNetChargeAsym = Read_netChAsym->GetRMS();
        StdDevNetChargeAsym = Read_netChAsym->GetStdDev();
        cout << "Loaded: charge asymmetry " << endl;
        cout << "MeanNetChargeAsym = " << MeanNetChargeAsym << endl;
        cout << "StdDevNetChargeAsym = " << StdDevNetChargeAsym << endl;
    }
    return 1;
}

int ReadEPDWeight(char *InFileName)
{
    TFile *fWgt_EPD = new TFile(InFileName, "READ");
    if (!fWgt_EPD->IsOpen())
        return 0;
    if (fWgt_EPD->IsOpen())
    {
        Read_EPD_EP1_east = (TProfile2D *)fWgt_EPD->Get("pEPD_EP1_east");
        Read_EPD_EP1_west = (TProfile2D *)fWgt_EPD->Get("pEPD_EP1_west");
        Read_EPD_EP_east = (TProfile2D *)fWgt_EPD->Get("pEPD_EP_east");
        Read_EPD_EP_west = (TProfile2D *)fWgt_EPD->Get("pEPD_EP_west");

        // float lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
        // float cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
        // float par1[9] = {474.651, 474.651, 474.651, 474.651, 474.651, 3.27243e+02, 1.72351, 1.72351, 1.72351};
        // float par2[9] = {3.55515, 3.55515, 3.55515, 3.55515, 3.55515, 3.56938, -7.69075, -7.69075, -7.69075};
        // float par3[9] = {1.80162, 1.80162, 1.80162, 1.80162, 1.80162, 1.67113, 4.72770, 4.72770, 4.72770};

        // for (int iy = 1; iy <= 9; iy++)
        // {
        //     for (int ix = 1; ix < 101; ix++)
        //     {
        //         double eta = wt->GetXaxis()->GetBinCenter(ix);
        //         wt->SetBinContent(ix, iy, (fabs(eta) > 3.8) ? lin[iy - 1] * eta + cub[iy - 1] * pow(eta, 3) : 0);
        //         wt2->SetBinContent(ix, iy, (fabs(eta) < 3.4) ? sqrt(1 - 1 / par1[iy - 1] / par1[iy - 1] / cosh(eta) / cosh(eta)) / (1 + exp((abs(eta) - par2[iy - 1]) / par3[iy - 1])) : 0);
        //         //                              wt2->SetBinContent(ix,iy, (fabs(eta)<3.4)? 1:0);
        //     }
        // }
    }
}

void create_hists(int cen)
{
    char fname[200];
    char title[200];

    std::vector<TString> name_options;
    name_options.push_back("Parity_int_obs3");
    name_options.push_back("Parity_int_ss_obs3");
    name_options.push_back("Delta_int_ss_obs3");
    name_options.push_back("Parity_int_obs1");
    name_options.push_back("Parity_int_ss_obs1");
    name_options.push_back("Delta_int_ss_obs1");
    std::vector<TString> name_options2;
    name_options2.push_back("_");
    name_options2.push_back("_rot_");
    name_options2.push_back("_QQcut_");
    name_options2.push_back("_QQcut_rot_");
    name_options2.push_back("_anti_");
    name_options2.push_back("_anti_rot_");
    name_options2.push_back("_QQcut_anti_");
    name_options2.push_back("_QQcut_anti_rot_");
    std::vector<int> bin_options;
    bin_options.push_back(24);
    bin_options.push_back(12);
    bin_options.push_back(4);

    for (int k = 0; k < 5; k++)
    {
        sprintf(fname, "Hist_v2_pt_obs2_caysm_%d", k);
        sprintf(title, "Hist_v2_pt_obs2_caysm %d", k);
        Hist_v2_pt_obs2_caysm[k] = new TProfile(fname, title, 300, 0, 15, -100, 100, "");
        sprintf(fname, "Hist_v2_pt_obs2_caysm_os_%d", k);
        sprintf(title, "Hist_v2_pt_obs2_caysm_os %d", k);
        Hist_v2_pt_obs2_caysm_os[k] = new TProfile(fname, title, 300, 0, 15, -100, 100, "");
        sprintf(fname, "Hist_v2_pt_obs2_p_caysm_%d", k);
        sprintf(title, "Hist_v2_pt_obs2_p_caysm %d", k);
        Hist_v2_pt_obs2_p_caysm[k] = new TProfile(fname, title, 300, 0, 15, -100, 100, "");
        sprintf(fname, "Hist_pt_Ach_%d", k);
        sprintf(title, "Hist_pt_Ach %d", k);
        Hist_pt_Ach[k] = new TH1D(fname, title, 100, 0, 2.5);
    }

    for (int j = 0; j < 6; j++)
    {
        for (int f = 0; f < 8; f++)
        {
            for (int sp_pt = 0; sp_pt < 15; sp_pt++)
            {
                sprintf(fname, "%s_splitpt%s%d", name_options.at(j).Data(), name_options2.at(f).Data(), sp_pt);
                sprintf(title, "%s_splitpt%s %d", name_options.at(j).Data(), name_options2.at(f).Data(), sp_pt);
                pParity_int_obs3_splitpt[j][sp_pt][f] = new TProfile(fname, title, bin_options.at(j % 3), 0.5, bin_options.at(j % 3) + 0.5, -100, 100, "");
            }
        }
    }

    std::vector<TString> name_options3;
    name_options3.push_back("TPC");
    name_options3.push_back("EPD");
    name_options3.push_back("EPD1");
    for (int jk = 0; jk < 3; jk++)
    {
        Hist_v2parent_pt_obs5[jk] = new TProfile(TString::Format("Hist_v2parent_pt_%s_obs5", name_options3.at(jk).Data()), TString::Format("v2 parent %s", name_options3.at(jk).Data()), 300, 0, 15, -100, 100, "");
        Hist_v2parent_pt_obs5_rot[jk] = new TProfile(TString::Format("Hist_v2parent_pt_%s_obs5_rot", name_options3.at(jk).Data()), TString::Format("v2 parent %s _rot", name_options3.at(jk).Data()), 300, 0, 15, -100, 100, "");
        Hist_v2_pt_obs1_alltrks[jk] = new TProfile(TString::Format("Hist_v2_pt_obs1_%s_alltrks", name_options3.at(jk).Data()), TString::Format("Hist_v2_pt_obs1_%s_alltrks", name_options3.at(jk).Data()), 300, 0, 15, -100, 100, "");
        Hist_v2_pt_obs2_alltrks[jk] = new TProfile(TString::Format("Hist_v2_pt_obs2_%s_alltrks", name_options3.at(jk).Data()), TString::Format("Hist_v2_pt_obs2_%s_alltrks", name_options3.at(jk).Data()), 300, 0, 15, -100, 100, "");

        Hist_v2pion_pt_obs5[jk] = new TProfile(TString::Format("Hist_v2pion_pt_%s_obs5", name_options3.at(jk).Data()), TString::Format("v2 pair_pion %s", name_options3.at(jk).Data()), 300, 0, 15, -100, 100, "");
    }

    std::vector<TString> name_options4;
    name_options4.push_back("v2e_Q2_all");
    name_options4.push_back("v2w_Q2_all");
    name_options4.push_back("v2e_Q2_P");
    name_options4.push_back("v2w_Q2_P");
    name_options4.push_back("v2e_Q2_N");
    name_options4.push_back("v2w_Q2_N");

    for (int b = 0; b < 6; b++)
    {
        p_v2_Q2_obs1[b] = new TProfile(TString::Format("p_%s_obs1", name_options4.at(b).Data()), TString::Format("p_%s_obs1", name_options4.at(b).Data()), 250, 0, 25, -100, 100, "");
        p_v2_Q2_obs2[b] = new TProfile(TString::Format("p_%s_obs2", name_options4.at(b).Data()), TString::Format("p_%s_obs2", name_options4.at(b).Data()), 250, 0, 25, -100, 100, "");
        p_v2_Q2_obs1_rot[b] = new TProfile(TString::Format("p_%s_obs1_rot", name_options4.at(b).Data()), TString::Format("p_%s_obs1_rot", name_options4.at(b).Data()), 250, 0, 25, -100, 100, "");
        p_v2_Q2_obs2_rot[b] = new TProfile(TString::Format("p_%s_obs2_rot", name_options4.at(b).Data()), TString::Format("p_%s_obs2_rot", name_options4.at(b).Data()), 250, 0, 25, -100, 100, "");
    }

    if (cen <= 1)
    {
        Hist_112ss_v2_parent_obs2 = new TProfile("Hist_112ss_v2_parent_obs2", "Hist_112ss_v2_parent_obs2", 50, -100, 100, -100, 100, "");
        Hist_112os_v2_parent_obs2 = new TProfile("Hist_112os_v2_parent_obs2", "Hist_112os_v2_parent_obs2", 50, -100, 100, -100, 100, "");
        Hist_132ss_v2_parent_obs2 = new TProfile("Hist_132ss_v2_parent_obs2", "Hist_132ss_v2_parent_obs2", 50, -100, 100, -100, 100, "");
        Hist_132os_v2_parent_obs2 = new TProfile("Hist_132os_v2_parent_obs2", "Hist_132os_v2_parent_obs2", 50, -100, 100, -100, 100, "");
        Hist_112ss_v2_parent_obs2_low = new TProfile("Hist_112ss_v2_parent_obs2_low", "Hist_112ss_v2_parent_obs2_low", 50, -100, 100, -100, 100, "");
        Hist_112os_v2_parent_obs2_low = new TProfile("Hist_112os_v2_parent_obs2_low", "Hist_112os_v2_parent_obs2_low", 50, -100, 100, -100, 100, "");
        Hist_132ss_v2_parent_obs2_low = new TProfile("Hist_132ss_v2_parent_obs2_low", "Hist_132ss_v2_parent_obs2_low", 50, -100, 100, -100, 100, "");
        Hist_132os_v2_parent_obs2_low = new TProfile("Hist_132os_v2_parent_obs2_low", "Hist_132os_v2_parent_obs2_low", 50, -100, 100, -100, 100, "");
        Hist_112ss_v2_parent_obs2_high = new TProfile("Hist_112ss_v2_parent_obs2_high", "Hist_112ss_v2_parent_obs2_high", 50, -100, 100, -100, 100, "");
        Hist_112os_v2_parent_obs2_high = new TProfile("Hist_112os_v2_parent_obs2_high", "Hist_112os_v2_parent_obs2_high", 50, -100, 100, -100, 100, "");
        Hist_132ss_v2_parent_obs2_high = new TProfile("Hist_132ss_v2_parent_obs2_high", "Hist_132ss_v2_parent_obs2_high", 50, -100, 100, -100, 100, "");
        Hist_132os_v2_parent_obs2_high = new TProfile("Hist_132os_v2_parent_obs2_high", "Hist_132os_v2_parent_obs2_high", 50, -100, 100, -100, 100, "");
    }
    else
    {
        Hist_112ss_v2_parent_obs2 = new TProfile("Hist_112ss_v2_parent_obs2", "Hist_112ss_v2_parent_obs2", 100, -100, 100, -100, 100, "");
        Hist_112os_v2_parent_obs2 = new TProfile("Hist_112os_v2_parent_obs2", "Hist_112os_v2_parent_obs2", 100, -100, 100, -100, 100, "");
        Hist_132ss_v2_parent_obs2 = new TProfile("Hist_132ss_v2_parent_obs2", "Hist_132ss_v2_parent_obs2", 100, -100, 100, -100, 100, "");
        Hist_132os_v2_parent_obs2 = new TProfile("Hist_132os_v2_parent_obs2", "Hist_132os_v2_parent_obs2", 100, -100, 100, -100, 100, "");
        Hist_112ss_v2_parent_obs2_low = new TProfile("Hist_112ss_v2_parent_obs2_low", "Hist_112ss_v2_parent_obs2_low", 100, -100, 100, -100, 100, "");
        Hist_112os_v2_parent_obs2_low = new TProfile("Hist_112os_v2_parent_obs2_low", "Hist_112os_v2_parent_obs2_low", 100, -100, 100, -100, 100, "");
        Hist_132ss_v2_parent_obs2_low = new TProfile("Hist_132ss_v2_parent_obs2_low", "Hist_132ss_v2_parent_obs2_low", 100, -100, 100, -100, 100, "");
        Hist_132os_v2_parent_obs2_low = new TProfile("Hist_132os_v2_parent_obs2_low", "Hist_132os_v2_parent_obs2_low", 100, -100, 100, -100, 100, "");
        Hist_112ss_v2_parent_obs2_high = new TProfile("Hist_112ss_v2_parent_obs2_high", "Hist_112ss_v2_parent_obs2_high", 100, -100, 100, -100, 100, "");
        Hist_112os_v2_parent_obs2_high = new TProfile("Hist_112os_v2_parent_obs2_high", "Hist_112os_v2_parent_obs2_high", 100, -100, 100, -100, 100, "");
        Hist_132ss_v2_parent_obs2_high = new TProfile("Hist_132ss_v2_parent_obs2_high", "Hist_132ss_v2_parent_obs2_high", 100, -100, 100, -100, 100, "");
        Hist_132os_v2_parent_obs2_high = new TProfile("Hist_132os_v2_parent_obs2_high", "Hist_132os_v2_parent_obs2_high", 100, -100, 100, -100, 100, "");
    }

    Hist_V0_Mass = new TH1F("Hist_V0_Mass", "V0 Inv. Mass", 200, 1.115684 - 0.07, 1.115684 + 0.07);
    Hist_V0_Mass_anti = new TH1F("Hist_V0_Mass_anti", "Anti V0 Inv. Mass", 200, 1.115684 - 0.07, 1.115684 + 0.07);
    // for (int i = 0; i < 9; i++)
    // {
    //     for (int j = 0; j < 20; j++)
    //     {
    //         for (int k = 0; k < 3; k++)
    //         {
    //             V0Mass_Q2[i][j][k] = new TH1F(TString::Format("V0Mass_%d_%d_%d", i, j, k), TString::Format("V0 Inv. Mass for Centrality-%d Pt Bin - %d %s", i, j, name_options3.at(k).Data()), 200, 1.115684 - 0.07, 1.115684 + 0.07);
    //             V0Mass_Q2_anti[i][j][k] = new TH1F(TString::Format("V0Mass_%d_%d_%d_anti", i, j, k), TString::Format("V0 Inv. Mass (AntiLam) for Centrality-%d Pt Bin - %d %s", i, j, name_options3.at(k).Data()), 200, 1.115684 - 0.07, 1.115684 + 0.07);
    //         }
    //     }
    // }
    for (int k = 0; k < 3; k++)
    {
        V0Mass_Q2[k] = new TH2F(TString::Format("V0Mass_%d", k), TString::Format("V0 Inv. Mass for %s", name_options3.at(k).Data()), 500, 0, 50, 200, 1.115684 - 0.07, 1.115684 + 0.07);
        V0Mass_Q2_anti[k] = new TH2F(TString::Format("V0Mass_%d_anti", k), TString::Format("V0 Inv. Mass (AntiLam) for %s", name_options3.at(k).Data()), 500, 0, 50, 200, 1.115684 - 0.07, 1.115684 + 0.07);
        V0Mass_Q2_pion[k] = new TH2F(TString::Format("V0Mass_%d_pion", k), TString::Format("V0 Inv. Mass for %s (Pion)", name_options3.at(k).Data()), 500, 0, 50, 200, 1.115684 - 0.07, 1.115684 + 0.07);
        V0Mass_Q2_pion_anti[k] = new TH2F(TString::Format("V0Mass_%d_pion_anti", k), TString::Format("V0 Inv. Mass (AntiLam) for %s (Pion)", name_options3.at(k).Data()), 500, 0, 50, 200, 1.115684 - 0.07, 1.115684 + 0.07);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////// Event-by-Event Analysis //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Gamma_112_module(int cen = 1, int opt_weight = 1, int sys_err_opt = 0, const Char_t *inFile = "./st_physics_12160020_raw_1010001.picoDst.root", const TString JobIDName = "1234")
{
    init_Gamma112(cen, JobIDName);

    corr_tree = new TChain("corr_tree");
    std::ifstream fin(inFile);
    string line1;
    char line[200];
    while (getline(fin, line1)) // loop wiill run till end of file
    {
        strcpy(line, line1.c_str());
        corr_tree->AddFile(line);
    }
    fin.close();

    getting_info_from_Trees();

    int n = corr_tree->GetEntries();
    cout << "Number of Entries = " << n << endl;

    for (int i = 0; i < n; ++i)
    {
        if (i % 1000 == 0)
            cerr << i << " / " << n << endl;
        corr_tree->GetEntry(i);

        construct_needed_variables();

        Make(cen, opt_weight, sys_err_opt, JobIDName);
    }

    finish_Gamma112(cen, opt_weight, JobIDName);
}

void getting_info_from_Trees(){
    // ch->SetBranchAddress("details", &details);
    // ch->SetBranchAddress("p_lambda", &p_lambda);
    // ch->SetBranchAddress("p_proton", &p_proton);
    // ch->SetBranchAddress("p_pion", &p_pion);
    // ch->SetBranchAddress("p_all", &p_all);
    // ch->SetBranchAddress("p_dau", &p_dau);

    // Event Variables
    // ignored: dt_nLambda, dt_nLambdaRot
    corr_tree->SetBranchAddress("dt_cent", &dt_cent);
    corr_tree->SetBranchAddress("dt_num_trk", &dt_num_trk);
    corr_tree->SetBranchAddress("n_lam", &n_lam);
    corr_tree->SetBranchAddress("n_dau", &n_dau);
    corr_tree->SetBranchAddress("dt_Run", &dt_Run);
    corr_tree->SetBranchAddress("dt_TOFMult", &dt_TOFMult);
    corr_tree->SetBranchAddress("dt_RefMult", &dt_RefMult);
    corr_tree->SetBranchAddress("dt_n_Proton", &dt_n_Proton);
    corr_tree->SetBranchAddress("dt_n_Pion", &dt_n_Pion);
    corr_tree->SetBranchAddress("dt_EventID", &dt_EventID);

    corr_tree->SetBranchAddress("dt_VPDvz", &dt_VPDvz);
    corr_tree->SetBranchAddress("dt_PVtxz", &dt_PVtxz);
    corr_tree->SetBranchAddress("dt_PVtxx", &dt_PVtxx);
    corr_tree->SetBranchAddress("dt_PVtxy", &dt_PVtxy);
    corr_tree->SetBranchAddress("dt_Eweight", &dt_Eweight);
    corr_tree->SetBranchAddress("dt_EPD_EP1_east", &dt_EPD_EP1_east);
    corr_tree->SetBranchAddress("dt_EPD_EP1_west", &dt_EPD_EP1_west);
    corr_tree->SetBranchAddress("dt_EPD_EP_east", &dt_EPD_EP_east);
    corr_tree->SetBranchAddress("dt_EPD_EP_west", &dt_EPD_EP_west);
    corr_tree->SetBranchAddress("dt_Magn", &dt_Magn);

    // Lambda Particle
    corr_tree->SetBranchAddress("p_lambda_px", &lambda_px);
    corr_tree->SetBranchAddress("p_lambda_py", &lambda_py);
    corr_tree->SetBranchAddress("p_lambda_pz", &lambda_pz);
    corr_tree->SetBranchAddress("p_lambda_Charge", &lambda_Charge);
    corr_tree->SetBranchAddress("p_lambda_dcaglobal", &lambda_dcaglobal);
    corr_tree->SetBranchAddress("p_lambda_nsigma", &lambda_nsigma);
    corr_tree->SetBranchAddress("p_lambda_mass", &lambda_mass);
    corr_tree->SetBranchAddress("p_lambda_trk_id", &lambda_trk_id);
    corr_tree->SetBranchAddress("p_lambda_hits_ratio", &lambda_hits_ratio);
    corr_tree->SetBranchAddress("p_lambda_nhitsfit", &lambda_nhitsfit);
    corr_tree->SetBranchAddress("p_lambda_nhitsmax", &lambda_nhitsmax);

    // // Proton Particle
    // corr_tree->SetBranchAddress("p_proton_px", &proton_px);
    // corr_tree->SetBranchAddress("p_proton_py", &proton_py);
    // corr_tree->SetBranchAddress("p_proton_pz", &proton_pz);
    // corr_tree->SetBranchAddress("p_proton_Charge", &proton_Charge);
    // corr_tree->SetBranchAddress("p_proton_dcaglobal", &proton_dcaglobal);
    // corr_tree->SetBranchAddress("p_proton_nsigma", &proton_nsigma);
    // corr_tree->SetBranchAddress("p_proton_mass", &proton_mass);
    // corr_tree->SetBranchAddress("p_proton_trk_id", &proton_trk_id);
    // corr_tree->SetBranchAddress("p_proton_hits_ratio", &proton_hits_ratio);
    // corr_tree->SetBranchAddress("p_proton_nhitsfit", &proton_nhitsfit);
    // corr_tree->SetBranchAddress("p_proton_nhitsmax", &proton_nhitsmax);

    // // Pion Particle
    // corr_tree->SetBranchAddress("p_pion_px", &pion_px);
    // corr_tree->SetBranchAddress("p_pion_py", &pion_py);
    // corr_tree->SetBranchAddress("p_pion_pz", &pion_pz);
    // corr_tree->SetBranchAddress("p_pion_Charge", &pion_Charge);
    // corr_tree->SetBranchAddress("p_pion_dcaglobal", &pion_dcaglobal);
    // corr_tree->SetBranchAddress("p_pion_nsigma", &pion_nsigma);
    // corr_tree->SetBranchAddress("p_pion_mass", &pion_mass);
    // corr_tree->SetBranchAddress("p_pion_trk_id", &pion_trk_id);
    // corr_tree->SetBranchAddress("p_pion_hits_ratio", &pion_hits_ratio);
    // corr_tree->SetBranchAddress("p_pion_nhitsfit", &pion_nhitsfit);
    // corr_tree->SetBranchAddress("p_pion_nhitsmax", &pion_nhitsmax);

    // All Particles
    corr_tree->SetBranchAddress("p_all_px", &all_px);
    corr_tree->SetBranchAddress("p_all_py", &all_py);
    corr_tree->SetBranchAddress("p_all_pz", &all_pz);
    corr_tree->SetBranchAddress("p_all_Charge", &all_Charge);
    corr_tree->SetBranchAddress("p_all_dcaglobal", &all_dcaglobal);
    corr_tree->SetBranchAddress("p_all_nSigmaProton", &all_nSigmaProton);
    corr_tree->SetBranchAddress("p_all_nSigmaPion", &all_nSigmaPion);
    corr_tree->SetBranchAddress("p_all_isTofTrack", &all_isTofTrack);
    corr_tree->SetBranchAddress("p_all_trk_id", &all_trk_id);
    corr_tree->SetBranchAddress("p_all_is_pion", &all_is_pion);
    corr_tree->SetBranchAddress("p_all_is_proton", &all_is_proton);
    corr_tree->SetBranchAddress("p_all_is_all", &all_is_all);
    corr_tree->SetBranchAddress("p_all_nhitsmax", &all_nhitsmax);
    corr_tree->SetBranchAddress("p_all_nhitsfit", &all_nhitsfit);

    // Daughter Particles
    corr_tree->SetBranchAddress("p_dau_px", &dau_px);
    corr_tree->SetBranchAddress("p_dau_py", &dau_py);
    corr_tree->SetBranchAddress("p_dau_pz", &dau_pz);
    corr_tree->SetBranchAddress("p_dau_dcaglobal", &dau_dcaglobal);
    corr_tree->SetBranchAddress("p_dau_nSigma", &dau_nSigma);
    corr_tree->SetBranchAddress("p_dau_trk_id", &dau_trk_id);
    corr_tree->SetBranchAddress("p_dau_nHitsFit", &dau_nHitsFit);
    corr_tree->SetBranchAddress("p_dau_nHitsMax", &dau_nHitsMax);
}

void construct_needed_variables(){
    details = new event(dt_cent, dt_num_trk, n_dau, 0, dt_Run, dt_TOFMult, dt_RefMult, dt_n_Proton, dt_n_Pion, dt_VPDvz, dt_PVtxz, dt_PVtxx, dt_PVtxy, dt_Eweight, dt_Magn, dt_EventID, dt_EPD_EP1_east, dt_EPD_EP1_west, dt_EPD_EP_east, dt_EPD_EP_west);
}

int Make(int cen, int opt_weight, int sys_err_opt, const TString JobIDName)
{
    // if (debug_1)
        // std::cout << "NEW EVENT!!" << endl;
    ///////////////////////// Transfer Important Details About Event /////////////////////////
    Centrality = details->cent;
    if (Centrality != cen) return 0;

    Run = details->Run;
    pVz = details->PVtxz;
    pVx = details->PVtxx;
    pVy = details->PVtxy;
    TVector3 pV(pVx, pVy, pVz);
    VPDvz = details->VPDvz;
    RefMult = details->RefMult;
    TOFMult = details->TOFMult;
    NPTracks = details->num_trk;
    Day = (int)((Run % 1000000) / 1000);
    Day2 = (int)((Run % 1000000) / 10);
    Day3 = (int)((Run % 1000000) / 1);
    Eweight = details->Eweight;
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    if (debug_1)
        std::cout << "pVz = " << pVz << endl;
    if ((sys_err_opt == 1) && (pVz >= 0))
        return 1; // Systematic Error Check 1

    initialize_event();                   // Setting Quantities Used to 0 for the Event
    determine_particlelists(sys_err_opt); // Determining what the particles of interest are
    if (debug_2)
        std::cout << "particle1->size() = " << particle1->size() << endl;
    if (debug_2)
        std::cout << "particle2->size() = " << particle2->size() << endl;
    if (debug_2)
        std::cout << "particle1_rot->size() = " << particle1_rot->size() << endl;
    n_particle1 = particle1->size();
    n_particle2 = particle2->size();
    n_particle1_rot = particle1_rot->size();

    if (debug_1)
        std::cout << "Centrality = " << Centrality << endl;

    if (!IsGoodEvent(cen))
        return 1;
    if (debug_1)
        std::cout << "debug_1 1" << endl;

    if (!CountCharge())
        return 1;
    if (debug_1)
        std::cout << "debug_1 2" << endl;

    // shuffle tracks for random EPs
    gv_gamma::iTrack.clear();
    Scount = Fcount / 3;
    for (int q = 0; q < Fcount; q++)
        gv_gamma::iTrack.push_back(q);
    random_shuffle(gv_gamma::iTrack.begin(), gv_gamma::iTrack.end());
    if (debug_1)
        std::cout << "debug_1 3" << endl;

    // TPC EP reconstruction
    MakeTPC_EP();
    if (debug_1)
        std::cout << "debug_1 3a" << endl;

    // // EPD EP
    // StEpdEpInfo result = mEpFinder->Results(mEpdHits, pV, Centrality);
    EPD_EP1_east = details->EPD_EP1_east;
    EPD_EP1_west = details->EPD_EP1_west;
    EPD_EP_east = details->EPD_EP_east;
    EPD_EP_west = details->EPD_EP_west;
    // if (debug_1)
    //     std::cout << "EPD_EP1_east = " << EPD_EP1_east << endl;
    // if (debug_1)
    //     std::cout << "EPD_EP1_west = " << EPD_EP1_west << endl;
    // if (debug_1)
    //     std::cout << "EPD_EP_east = " << EPD_EP_east << endl;
    // if (debug_1)
    //     std::cout << "EPD_EP_west = " << EPD_EP_west << endl;
    // if ((opt_useEPD == 1 || opt_useEPD == 11) && (EPD_EP_east == EPD_EP_west || (EPD_EP1_east > 0.0264 && EPD_EP1_east < 0.0265) || (EPD_EP1_west > 0.0264 && EPD_EP1_west < 0.0265)))
    //     return 1;
    EPDe_Day3_cos1->Fill(Day3, cos(EPD_EP1_east));
    EPDe_Day3_sin1->Fill(Day3, sin(EPD_EP1_east));
    EPDw_Day3_cos1->Fill(Day3, cos(EPD_EP1_west));
    EPDw_Day3_sin1->Fill(Day3, sin(EPD_EP1_west));
    EPDe_Day3_cos2->Fill(Day3, cos(nHar * EPD_EP_east));
    EPDe_Day3_sin2->Fill(Day3, sin(nHar * EPD_EP_east));
    EPDw_Day3_cos2->Fill(Day3, cos(nHar * EPD_EP_west));
    EPDw_Day3_sin2->Fill(Day3, sin(nHar * EPD_EP_west));
    if (debug_1)
        std::cout << "debug_1 4" << endl;

    // flatten EPs
    ShiftPsi();
    if (debug_1)
        std::cout << "debug_1 7" << endl;
    FillEP_resolution();

    // Store the flattened phi for POI (Lambdas)
    if (debug_1 && (n_particle1 != 0))
        std::cout << "n_particle1 1 = " << n_particle1 << endl;
    int number_oflambdas = 0, number_ofantilambdas = 0, number_oflambdas_rot = 0, number_ofantilambdas_rot = 0, number_ofprotons = 0, number_ofantiprotons = 0;

    Phi_new.resize(n_particle1);
    for (int trki = 0; trki < n_particle1; trki++)
    {
        if (debug_1)
            std::cout << "trki nor = " << trki << endl;
        gv_gamma::trk_mom->SetXYZ(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz);
        Eta = gv_gamma::trk_mom->Eta();
        Pt = gv_gamma::trk_mom->Pt();
        Charge = particle1->at(trki).Charge; // 1 for lambda, -1 for antilambdas
        Phi = gv_gamma::trk_mom->Phi();

        if (Charge == 1)
            number_oflambdas++;
        else if (Charge == -1)
            number_ofantilambdas++;
        if (debug_1)
            std::cout << "debug_1 5" << endl;

        FillPhiPOI(); // Charge is needed here
        if (debug_1)
            std::cout << "debug_1 6" << endl;
        ShiftPhiPOI(trki);
        if (debug_1)
            std::cout << "debug_1 7" << endl;
    }
    num_lam_tree->Fill(number_oflambdas);
    num_antilam_tree->Fill(number_ofantilambdas);

    // Store the flattened phi for POI (BKG Lambdas)
    if (debug_1 && (n_particle1_rot != 0))
        std::cout << "n_particle1_rot 1 = " << n_particle1_rot << endl;

    Phi_new_rot.resize(n_particle1_rot);
    for (int trki = 0; trki < n_particle1_rot; trki++)
    {
        if (debug_1)
            std::cout << "trki rot = " << trki << endl;
        gv_gamma::trk_mom->SetXYZ(particle1_rot->at(trki).px, particle1_rot->at(trki).py, particle1_rot->at(trki).pz);
        Eta_rot = gv_gamma::trk_mom->Eta();
        Pt_rot = gv_gamma::trk_mom->Pt();
        Charge_rot = particle1_rot->at(trki).Charge; // 1 for lambda, -1 for antilambdas
        Phi_rot = gv_gamma::trk_mom->Phi();

        if (Charge_rot == 1)
            number_oflambdas_rot++;
        else if (Charge_rot == -1)
            number_ofantilambdas_rot++;
        if (debug_1)
            std::cout << "debug_1 5 rot" << endl;

        FillPhiPOI_rot(); // Charge is needed here
        if (debug_1)
            std::cout << "debug_1 6 rot" << endl;
        ShiftPhiPOI_rot(trki);
        if (debug_1)
            std::cout << "debug_1 7 rot" << endl;
    }
    num_lam_rot_tree->Fill(number_oflambdas_rot);
    num_antilam_rot_tree->Fill(number_ofantilambdas_rot);

    // Store the flattened phi for POI (Protons)
    if (debug_1 && (n_particle2 != 0))
        std::cout << "n_particle2 1 = " << n_particle2 << endl;

    Phi_pnew.resize(n_particle2);
    for (int trki = 0; trki < n_particle2; trki++)
    {
        if (debug_1)
            std::cout << "trki proton = " << trki << endl;
        gv_gamma::trk_mom->SetXYZ(particle2->at(trki).px, particle2->at(trki).py, particle2->at(trki).pz);
        Eta2 = gv_gamma::trk_mom->Eta();
        Pt2 = gv_gamma::trk_mom->Pt();
        Charge2 = particle2->at(trki).Charge;
        DCAGlobal2 = particle2->at(trki).dcaglobal;
        Phi2 = gv_gamma::trk_mom->Phi();

        if (Charge2 == 1)
            number_ofprotons++;
        else if (Charge2 == -1)
            number_ofantiprotons++;
        if (debug_1)
            std::cout << "debug_1 5 proton" << endl;

        FillPhiPOI_p(); // Charge is needed here
        if (debug_1)
            std::cout << "debug_1 6 proton" << endl;
        ShiftPhiPOI_p(trki);
        if (debug_1)
            std::cout << "debug_1 7 proton" << endl;
    }
    num_proton_tree->Fill(number_ofprotons);
    num_antiproton_tree->Fill(number_ofantiprotons);

    if (debug_1)
        std::cout << "Done with weights" << endl;

    //////////Real analysis begins here//////////////////////////////////////////////////////////////////////////////////////
    int n_lam_used = 0, n_proton_used = 0, n_p_used = 0, n_ap_used = 0, n_tp_used = 0, Qcount_parent = 0, Qcount_parent_QQcut = 0;
    n_gamma = 0, mQQx_parent = 0., mQQy_parent = 0., Q2_parent = 0., mQQx_parent_QQcut = 0., mQQy_parent_QQcut = 0., Q2_parent_QQcut = 0.;
    double p_px[2] = {0.}, p_py[2] = {0.}, p_pz[2] = {0.}, mQQx_p[2] = {0.}, mQQy_p[2] = {0.}, p_phi[2] = {0.}, p_trkid[2] = {0.}, p_nHitsRatio[2] = {0.}, p_nhitsfit[2] = {0.}, p_mom[2] = {0.}, parent_phi[2] = {0.};
    double lam_px = 0., lam_py = 0., lam_pz = 0., lam_phi = 0., lam_id = 0., lam_pt = 0., lam_phinew = 0.;
    double lam_dau_px[2] = {0.}, lam_dau_py[2] = {0.}, lam_dau_pz[2] = {0.}, lam_dau_phi[2] = {0.}, lam_dau_id[2] = {0.};
    int p_charge[2] = {0};

    // if(opt_useEPD == 1 || opt_useEPD == 2) Eweight *= 2;

    // loop for the real analysis
    if (debug_1)
        std::cout << "Before Real Analysis Loop" << endl;
    for (int trki = 0; trki < n_particle1; trki++)
    {
        if (debug_1)
            std::cout << "Real Analysis Loop?" << endl;

        gv_gamma::trk_mom->SetXYZ(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz);
        Eta = gv_gamma::trk_mom->Eta();
        Pt = gv_gamma::trk_mom->Pt();
        Charge = particle1->at(trki).Charge; // 1 for lambda, -1 for antilambdas
        Phi = gv_gamma::trk_mom->Phi();
        Theta = 2. * atan(exp(-Eta));
        eff = 1; // Efficiency information later
        float eff_tof = 1;
        eff *= eff_tof;

        ////// for checks on strange features in Q2 distribution //////
        lam_px = particle1->at(trki).px;
        lam_py = particle1->at(trki).py;
        lam_pz = particle1->at(trki).pz;
        lam_phi = Phi;
        lam_id = particle1->at(trki).trk_id;
        lam_pt = Pt;
        lam_phinew = Phi_new[trki];
        lam_dau_px[0] = particle1_dau->at(2 * trki).px;
        lam_dau_py[0] = particle1_dau->at(2 * trki).py;
        lam_dau_pz[0] = particle1_dau->at(2 * trki).pz;
        TVector3 tempdau(particle1_dau->at(2 * trki).px, particle1_dau->at(2 * trki).py, particle1_dau->at(2 * trki).pz);
        lam_dau_phi[0] = tempdau.Phi();
        lam_dau_id[0] = particle1_dau->at(2 * trki).trk_id;
        lam_dau_px[1] = particle1_dau->at(2 * trki + 1).px;
        lam_dau_py[1] = particle1_dau->at(2 * trki + 1).py;
        lam_dau_pz[1] = particle1_dau->at(2 * trki + 1).pz;
        TVector3 tempdau2(particle1_dau->at(2 * trki + 1).px, particle1_dau->at(2 * trki + 1).py, particle1_dau->at(2 * trki + 1).pz);
        lam_dau_phi[1] = tempdau2.Phi();
        lam_dau_id[1] = particle1_dau->at(2 * trki + 1).trk_id;
        //////////////////////////////////////////////////////////////////

        ////// for checks systematic error cuts //////
        nHitsFitQA->Fill(particle1_dau->at(2 * trki).nHitsFit);
        nHitsFitQA->Fill(particle1_dau->at(2 * trki + 1).nHitsFit);
        nSigma_dauQA->Fill(particle1_dau->at(2 * trki).nSigma);
        nSigma_dauQA->Fill(particle1_dau->at(2 * trki + 1).nSigma);
        dca_protonQA->Fill(particle1_dau->at(2 * trki).dcaglobal);
        dca_pionQA->Fill(particle1_dau->at(2 * trki + 1).dcaglobal);
        dca_LambdaQA->Fill(particle1->at(trki).dcaglobal);
        //////////////////////////////////////////////////////////////////

        ////// obtaining efficiency information for Lambdas //////
        if ((particle_option1 == 0) && (Pt < 2.2))
        {
            if (debug_1)
                std::cout << "Reading eff " << endl;
            if (Charge == 1)
            {
                temp_eff = efficiency[Centrality][(int)floor((Pt) / 0.1)];

                if (debug_1)
                    std::cout << "efficiency[" << Centrality << "][" << (int)floor((Pt) / 0.1) << "] = " << efficiency[Centrality][(int)floor((Pt) / 0.1)] << endl;
                if (debug_1)
                    std::cout << "temp eff = " << temp_eff << endl;
            }
            else if (Charge == -1)
            {
                temp_eff = efficiency_anti[Centrality][(int)floor((Pt) / 0.1)];

                if (debug_1)
                    std::cout << "efficiency_anti[" << Centrality << "][" << (int)floor((Pt) / 0.1) << "] = " << efficiency_anti[Centrality][(int)floor((Pt) / 0.1)] << endl;
                if (debug_1)
                    std::cout << "temp eff = " << temp_eff << endl;
            }

            if (debug_1)
                std::cout << "Charge = " << Charge << endl;
        }

        if (debug_1)
            std::cout << "temp eff = " << temp_eff << endl;

        if (temp_eff < 0)
            continue;
        //////////////////////////////////////////////////////////////////

        int lam_index_tmp = 0;
        if (Pt < 2.0)
            lam_index_tmp = (int)floor((Pt - 0.5) / 0.1);
        else if (Pt == 2.0)
            lam_index_tmp = 14;

        hEtaPhiDist->Fill(Eta, Phi, Eweight);
        hPhiPtDist->Fill(Phi, Pt, Eweight);

        double mass_temp = 0;
        if ((particle_option1 == 0))
            mass_temp = 1.115684;
        else if (particle_option1 == 1)
            mass_temp = 0.93827;
        TLorentzVector trk_lorentz(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz, sqrt(pow(gv_gamma::trk_mom->Mag(), 2) + pow(mass_temp, 2)));
        TLorentzVector lA;
        lA.SetPxPyPzE(particle1->at(trki).px, particle1->at(trki).py, particle1->at(trki).pz, sqrt(mass_temp * mass_temp + Pt * Pt * cosh(Eta) * cosh(Eta)));

        V0Mass_gamma112->Fill(particle1->at(trki).mass);
        if (Charge == 1)
            hEtaPtDist->Fill(Eta, Pt, Eweight);
        else if (Charge == -1)
            hEtaPtDist_anti->Fill(Eta, Pt, Eweight);
        Hist_Pt->Fill(Pt, Eweight);

        float mQx_i = mQx, mQy_i = mQy;
        TVector2 mQ_i(mQx_i, mQy_i);
        float psi_F = mQ_i.Phi() / nHar;
        Hist_TPC_EP_full_m1->Fill(psi_F);
        float psi_F_new = psi_F;
        for (int jj = 0; jj < order; jj++)
            psi_F_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F) / nHar / (jj + 1);
        Hist_TPC_EP_full_m1_flat->Fill(psi_F_new);

        ////// obtaining v2 for Lambdas //////
        v2_sub = (Eta > 0) ? cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_new) * 100 : cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_new) * 100;
        FillCMW();
        v2 = cos(nHar * Phi_new[trki] - nHar * psi_F_new) * 100;
        v2e = cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_new) * 100;
        v2w = cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_new) * 100;
        v2epos = cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_pos_new) * 100;
        v2eneg = cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_neg_new) * 100;
        v2wpos = cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_pos_new) * 100;
        v2wneg = cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_neg_new) * 100;
        v2_EPD1 = cos(2. * Phi_new[trki] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
        v2_EPDe = cos(nHar * Phi_new[trki] - nHar * EPD_EP_east_new) * 100;
        v2_EPDw = cos(nHar * Phi_new[trki] - nHar * EPD_EP_west_new) * 100;

        Fillv2();

        float v2_final = 0;
        if (Eta > 0)
        {
            v2_final = v2w;

            if (Charge > 0)
                v2_sub_ss = v2wpos;
            else if (Charge < 0)
                v2_sub_os = v2wneg;
        }
        else if (Eta < 0)
        {
            v2_final = v2e;

            if (Charge > 0)
                v2_sub_ss = v2epos;
            else if (Charge < 0)
                v2_sub_os = v2eneg;
        }

        v2_single_lamda->Fill(0.5, v2, Eweight / temp_eff);
        v2_single_lamda_err->Fill(0.5, v2, Eweight);
        v2_single_lamda->Fill(1.5, (v2_EPDe + v2_EPDw) / 2.0, Eweight / temp_eff);
        v2_single_lamda_err->Fill(1.5, (v2_EPDe + v2_EPDw) / 2.0, Eweight);
        v2_single_lamda->Fill(2.5, v2_EPD1, Eweight / temp_eff);
        v2_single_lamda_err->Fill(2.5, v2_EPD1, Eweight);
        //////////////////////////////////////////////////////////////////

        ////// obtaining charge assymetry v2 plots for Lambdas //////
        int Npos_tmp = Npos, Nneg_tmp = Nneg;
        float net_Nch_Asym_tmp = 0;
        int net_charge_asym_bin_tmp = -1;

        if (particle_option1 == 0)
        {
            for (int trk = 0; trk < NPTracks; trk++)
            {
                if ((particle1_dau->at(2 * trki).trk_id == all_trk_id->at(trk)) || (particle1_dau->at(2 * trki + 1).trk_id == all_trk_id->at(trk)))
                {
                    TVector3 trk_mom_tmp(all_px->at(trk), all_py->at(trk), all_pz->at(trk));
                    float EtaAsso_tmp = trk_mom_tmp.Eta();
                    float PtAsso_tmp = trk_mom_tmp.Pt();
                    float ChargeAsso_tmp = all_Charge->at(trk);
                    float DCAGlobalAsso_tmp = all_dcaglobal->at(trk);
                    float nSigma_p_tmp = all_nSigmaProton->at(trk);
                    if (fabs(EtaAsso_tmp) < 1 && PtAsso_tmp > 0.15 && DCAGlobalAsso_tmp < 1 && !(fabs(nSigma_p_tmp) < 3 && PtAsso_tmp < 0.4))
                    {
                        if (IsGoodAsso(PtAsso_tmp, EtaAsso_tmp, DCAGlobalAsso_tmp))
                        {
                            if (ChargeAsso_tmp > 0)
                                Npos_tmp = Npos_tmp - 1;
                            else if (ChargeAsso_tmp < 0)
                                Nneg_tmp = Nneg_tmp - 1;
                        }
                    }
                }
            }

            if ((Npos_tmp + Nneg_tmp) > 0)
                net_Nch_Asym_tmp = (Npos_tmp - Nneg_tmp) / float(Npos_tmp + Nneg_tmp);
            if (net_Nch_Asym_tmp > -99)
            {
                if (net_Nch_Asym_tmp < (MeanNetChargeAsym - 0.8 * StdDevNetChargeAsym))
                    net_charge_asym_bin_tmp = 0;
                else if (net_Nch_Asym_tmp < (MeanNetChargeAsym - (0.3 * StdDevNetChargeAsym)))
                    net_charge_asym_bin_tmp = 1;
                else if (net_Nch_Asym_tmp < (MeanNetChargeAsym + (0.2 * StdDevNetChargeAsym)))
                    net_charge_asym_bin_tmp = 2;
                else if (net_Nch_Asym_tmp < (MeanNetChargeAsym + 0.8 * StdDevNetChargeAsym))
                    net_charge_asym_bin_tmp = 3;
                else
                    net_charge_asym_bin_tmp = 4;
            }

            Hist_v2_pt_obs1->Fill(Pt, v2_final);
            if (net_Nch_Asym > -99)
                Hist_v2_Ach->Fill(net_Nch_Asym, v2_final);
            Hist_v2_pt_obs2->Fill(Pt, v2_final, Eweight);
        }
        else if (particle_option1 == 1)
        {
            if (Charge == 1)
            {
                Hist_v2_pt_obs1->Fill(Pt, v2_final);
                if (net_Nch_Asym > -99)
                    Hist_v2_Ach->Fill(net_Nch_Asym, v2_final);
                Hist_v2_pt_obs2->Fill(Pt, v2_final, Eweight);
            }
        }
        if (net_charge_asym_bin_tmp >= 0)
        {
            Hist_v2_pt_obs2_caysm[net_charge_asym_bin_tmp]->Fill(Pt, v2_sub_ss, Eweight);
            Hist_v2_pt_obs2_caysm_os[net_charge_asym_bin_tmp]->Fill(Pt, v2_sub_os, Eweight);
        }
        //////////////////////////////////////////////////////////////////

        Hist_v2_eta_obs1->Fill(Eta, v2_final, 1.);
        Hist_v2_eta_obs2->Fill(Eta, v2_final, Eweight);
        Hist_v2_pt_EPD1_obs->Fill(Pt, v2_EPD1, Eweight);
        Hist_v2_pt_EPD_obs->Fill(Pt, 0.5 * (v2_EPDe + v2_EPDw), Eweight);
        Hist_v2_eta_EPD1_obs->Fill(Eta, v2_EPD1, Eweight / eff);
        Hist_v2_eta_EPD_obs->Fill(Eta, 0.5 * (v2_EPDe + v2_EPDw), Eweight / eff);
        Hist_v2_eta_EPDe_obs->Fill(Eta, v2_EPDe, Eweight / eff);
        Hist_v2_eta_EPDw_obs->Fill(Eta, v2_EPDw, Eweight / eff);

        n_lam_used++;

        if (opt_weight == 1)
            continue;
        for (int trkj = 0; trkj < n_particle2; trkj++)
        {
            if (debug_2)
                std::cout << "Proton Loop" << endl;

            Charge2 = particle2->at(trkj).Charge; // 1 for protons, -1 for antiprotons
            gv_gamma::trk_mom2->SetXYZ(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz);
            Eta2 = gv_gamma::trk_mom2->Eta();
            Pt2 = gv_gamma::trk_mom2->Pt();
            DCAGlobal2 = particle2->at(trkj).dcaglobal;
            Phi2 = gv_gamma::trk_mom2->Phi();
            eff2 = 1; // efficiency information acquired later
            float eff_tof2 = 1;
            eff2 *= eff_tof2;
            if (debug_2)
                std::cout << "Starting Track 2" << endl;


            if (trki == 0){
                
                if (debug_2)
                    std::cout << "Starting Proton Track Only" << endl;
                
                float temp_temp_proton_eff = 1;
                if (particle_option2 == 0)
                {
                    temp_temp_proton_eff = fit_eff_pT_final[Centrality]->Eval(Pt2);
                    temp_temp_proton_eff *= fit_tof_eff_pT_final[Centrality]->Eval(Pt2);
                }

                float mQx_j = mQx, mQy_j = mQy;

                TVector2 mQ_j(mQx_j, mQy_j);
                float psi_F2 = mQ_j.Phi() / nHar;
                Hist_TPC_EP_full_m2->Fill(psi_F2);
                float psi_F2_new = psi_F2;
                for (int jj = 0; jj < order; jj++)
                    psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
                if (psi_F2_new < 0)
                    psi_F2_new += PI;
                if (psi_F2_new > PI)
                    psi_F2_new -= PI;
                Hist_TPC_EP_full_m2_flat->Fill(psi_F2_new);

                n_proton_used++;

                v2p = cos(nHar * Phi_pnew[trkj] - nHar * psi_F2_new) * 100;
                // Hist_v2_pt_obs1_p->Fill(Pt2, v2p);
                // if(net_Nch_Asym > -99) Hist_v2_Ach_p->Fill(net_Nch_Asym, v2p);
                // Hist_v2_pt_obs2_p->Fill(Pt2, v2p, Eweight);
                v2pe = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_for_new) * 100;
                v2pw = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_bac_new) * 100;
                float v2p_EPDe = cos(nHar * Phi_pnew[trkj] - nHar * EPD_EP_east_new) * 100;
                float v2p_EPDw = cos(nHar * Phi_pnew[trkj] - nHar * EPD_EP_west_new) * 100;
                float v2p_EPD1 = cos(nHar * Phi_pnew[trkj] - EPD_EP1_west_new - EPD_EP1_east_new) * 100;

                float v2_final = 0;

                if (Eta2 > 0)
                    v2_final = v2pw;
                else if (Eta2 < 0)
                    v2_final = v2pe;

                float v2_final_p[3] = {0.};
                v2_final_p[0] = v2p;
                v2_final_p[1] = (v2p_EPDe + v2p_EPDw) / 0.5;
                v2_final_p[2] = v2p_EPD1;

                for (int v2_p_i = 0; v2_p_i < 3; v2_p_i++)
                {
                    v2_single_proton->Fill((v2_p_i + 0.5), v2_final_p[v2_p_i], Eweight / temp_temp_proton_eff);
                    v2_single_proton_err->Fill((v2_p_i + 0.5), v2_final_p[v2_p_i], Eweight);
                }

                if (((Charge2 == 1) && (particle_option1 == 0)))
                {
                    Hist_v2_pt_obs1_p->Fill(Pt2, v2_final);
                    if (net_Nch_Asym > -99)
                        Hist_v2_Ach_p->Fill(net_Nch_Asym, v2_final);
                    Hist_v2_pt_obs2_p->Fill(Pt2, v2_final, Eweight);

                    Hist_v2_pt_obs2_p_caysm[net_charge_asym_bin]->Fill(Pt2, v2_final, Eweight);
                }
                else if (Charge2 == -1)
                {
                    Hist_v2_pt_obs1_ap->Fill(Pt2, v2_final);
                    Hist_v2_pt_obs2_ap->Fill(Pt2, v2_final, Eweight);
                }

                if (debug_2)
                    std::cout << "Ending Track 2 Only" << endl;
            }

            if ((particle_option1 == 0) && ((particle1_dau->at(2 * trki).trk_id == particle2->at(trkj).trk_id) || (particle1_dau->at(2 * trki + 1).trk_id == particle2->at(trkj).trk_id)))
                continue;

            nSigma_prim_protonQA->Fill(particle2->at(trkj).nsigma);

            if ((particle_option1 == 1) && (particle_option2 == 0) && (particle1->at(trki).trk_id == particle2->at(trkj).trk_id))
                continue;

            TLorentzVector trk_lorentz2(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz, sqrt(pow(gv_gamma::trk_mom2->Mag(), 2) + pow(particle2->at(trkj).mass, 2)));
            rapidity2 = trk_lorentz2.Rapidity();
            TLorentzVector lB;
            lB.SetPxPyPzE(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz, sqrt(particle2->at(trkj).mass * particle2->at(trkj).mass + Pt2 * Pt2 * cosh(Eta2) * cosh(Eta2)));
            if (debug_2)
                std::cout << "rapidity2 = " << rapidity2 << endl;

            ///////////// Lambda Proton Invariant Mass Calculations!! /////////////
            TLorentzVector trk_lorentz3;
            trk_lorentz3 = lA + lB;

            Hist_inv_Mass->Fill(trk_lorentz3.M(), Eweight / efficiency[cen][lam_index_tmp]);
            //////////////////////////////////////////////////////////////////////////////
            if (debug_2)
                std::cout << "Pt2 = " << Pt2 << endl;
            hDpt->Fill(fabs(Pt - Pt2), Eweight);

            float mQx_j = mQx_i, mQy_j = mQy_i;
            TVector2 mQ_j(mQx_j, mQy_j);
            float psi_F2 = mQ_j.Phi() / nHar;
            Hist_TPC_EP_full_m2->Fill(psi_F2);
            float psi_F2_new = psi_F2;
            for (int jj = 0; jj < order; jj++)
                psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
            if (psi_F2_new < 0)
                psi_F2_new += PI;
            if (psi_F2_new > PI)
                psi_F2_new -= PI;
            Hist_TPC_EP_full_m2_flat->Fill(psi_F2_new);

            if (debug_2)
                std::cout << "Before Efficiency " << endl;

            if (particle_option2 == 0)
            {
                temp_proton_eff = fit_eff_pT_final[Centrality]->Eval(Pt2);
                temp_proton_eff *= fit_tof_eff_pT_final[Centrality]->Eval(Pt2);
            }
            else
            {
                temp_proton_eff = 1;
            }
            if (debug_2)
                std::cout << "After Efficiency " << endl;

            // correlations
            correlator3 = cos(Phi_new[trki] - Phi_pnew[trkj]) * 100;

            // TPC
            //  correlator4e = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * TPC_EP_for_new) * 100;
            //  correlator4w = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * TPC_EP_bac_new) * 100;
            correlator4e[0] = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * TPC_EP_east_new) * 100;
            correlator4w[0] = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * TPC_EP_west_new) * 100;
            correlator4[0] = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * psi_F2_new) * 100;
            correlator0[0] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + 2 * psi_F2_new) * 100;
            correlator0e[0] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + 2 * TPC_EP_east_new) * 100;
            correlator0w[0] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + 2 * TPC_EP_west_new) * 100;
            correlator0_alt[0] = cos(Phi_pnew[trkj] - 3 * Phi_new[trki] + 2 * psi_F2_new) * 100;
            correlator0e_alt[0] = cos(Phi_pnew[trkj] - 3 * Phi_new[trki] + 2 * TPC_EP_east_new) * 100;
            correlator0w_alt[0] = cos(Phi_pnew[trkj] - 3 * Phi_new[trki] + 2 * TPC_EP_west_new) * 100;

            // EPD
            correlator4e[1] = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * EPD_EP_east_new) * 100;
            correlator4w[1] = cos(Phi_new[trki] + Phi_pnew[trkj] - 2 * EPD_EP_west_new) * 100;
            correlator4[1] = 0.5 * (correlator4e[1] + correlator4w[1]);
            correlator0e[1] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + 2 * EPD_EP_east_new) * 100;
            correlator0w[1] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + 2 * EPD_EP_west_new) * 100;
            correlator0[1] = 0.5 * (correlator0e[1] + correlator0w[1]);
            correlator0e_alt[1] = cos(Phi_pnew[trkj] - 3 * Phi_new[trki] + 2 * EPD_EP_east_new) * 100;
            correlator0w_alt[1] = cos(Phi_pnew[trkj] - 3 * Phi_new[trki] + 2 * EPD_EP_west_new) * 100;
            correlator0_alt[1] = 0.5 * (correlator0e_alt[1] + correlator0w_alt[1]);

            // EPD 1st Order
            correlator4[2] = cos(Phi_new[trki] + Phi_pnew[trkj] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
            correlator4e[2] = cos(Phi_new[trki] + Phi_pnew[trkj] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
            correlator4w[2] = cos(Phi_new[trki] + Phi_pnew[trkj] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
            correlator0[2] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0e[2] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0w[2] = cos(Phi_new[trki] - 3 * Phi_pnew[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0_alt[2] = cos(Phi_pnew[trki] - 3 * Phi_new[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0e_alt[2] = cos(Phi_pnew[trki] - 3 * Phi_new[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0w_alt[2] = cos(Phi_pnew[trki] - 3 * Phi_new[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;

            float coscos = cos(Phi_new[trki] - psi_F2_new) * cos(Phi_pnew[trkj] - psi_F2_new);
            float sinsin = sin(Phi_new[trki] - psi_F2_new) * sin(Phi_pnew[trkj] - psi_F2_new);

            if (n_gamma <= 1)
            {
                p_px[n_gamma] = particle2->at(trkj).px;
                p_py[n_gamma] = particle2->at(trkj).py;
                p_pz[n_gamma] = particle2->at(trkj).pz;
                p_charge[n_gamma] = particle2->at(trkj).Charge;
                p_phi[n_gamma] = Phi_pnew[trkj];
                p_trkid[n_gamma] = particle2->at(trkj).trk_id;
                p_mom[n_gamma] = Pt2;
                p_nHitsRatio[n_gamma] = particle2->at(trkj).hits_ratio;
                p_nhitsfit[n_gamma] = particle2->at(trkj).nhitsfit;
            }

            if (debug_2)
                cout << "Gamma112 = (correlator4) " << correlator4 << endl;
            if (debug_2)
                cout << "Gamma112 = (correlator4e) " << correlator4e << endl;
            if (debug_2)
                cout << "Gamma112 = (correlator4w) " << correlator4w << endl;

            ///////////// Lambda-Proton Parent v2 calculations for ESE /////////////
            double v2e_parent = cos(nHar * trk_lorentz3.Phi() - nHar * TPC_EP_for_new) * 100;
            double v2w_parent = cos(nHar * trk_lorentz3.Phi() - nHar * TPC_EP_bac_new) * 100;
            double v2e_parent_EPD = cos(nHar * trk_lorentz3.Phi() - nHar * EPD_EP_east_new) * 100;
            double v2w_parent_EPD = cos(nHar * trk_lorentz3.Phi() - nHar * EPD_EP_west_new) * 100;
            double v2_parent_EPD1 = cos(nHar * trk_lorentz3.Phi() - EPD_EP1_west_new - EPD_EP1_east_new) * 100;

            double v2_final_parent[3] = {0.};
            double phi_final_parent = 0;

            if (trk_lorentz3.Eta() > 0)
            {
                v2_final_parent[0] = v2w_parent;
                phi_final_parent = trk_lorentz3.Phi() - TPC_EP_bac_new;
            }
            else
            {
                v2_final_parent[0] = v2e_parent;
                phi_final_parent = trk_lorentz3.Phi() - TPC_EP_for_new;
            }
            v2_final_parent[1] = 0.5 * (v2e_parent_EPD + v2w_parent_EPD);
            v2_final_parent[2] = v2_parent_EPD1;

            if (phi_final_parent < -PI)
                phi_final_parent = phi_final_parent + 2 * PI;

            Hist_v2_parent_lam_pT->Fill(Pt, v2_final_parent[0], Eweight);
            Hist_v2_parent_parent_pT->Fill(trk_lorentz3.Perp(), v2_final_parent[0], Eweight);
            Hist_v2_lam_parent_pT->Fill(trk_lorentz3.Perp(), v2_final, Eweight);

            for (int v2_i = 0; v2_i < 3; v2_i++)
            {
                v2_single_lamp->Fill((v2_i + 0.5), v2_final_parent[v2_i], Eweight / (temp_eff * temp_proton_eff));
                v2_single_lamp_err->Fill((v2_i + 0.5), v2_final_parent[v2_i], Eweight);
            }
            ///////////////////////////////////////////////////////////////////////////////////////////

            TVector3 kstar((trk_lorentz2.Px() - trk_lorentz.Px()) / 2.0, (trk_lorentz2.Py() - trk_lorentz.Py()) / 2.0, (trk_lorentz2.Pz() - trk_lorentz.Pz()) / 2.0);
            QQ = pow(trk_lorentz2.E() - trk_lorentz.E(), 2) - pow(trk_lorentz2.Px() - trk_lorentz.Px(), 2) - pow(trk_lorentz2.Py() - trk_lorentz.Py(), 2) - pow(trk_lorentz2.Pz() - trk_lorentz.Pz(), 2);
            if (QQ <= 0)
                QQ = sqrt(-QQ);
            else
                QQ = -99;

            Hist_charge1_dist->Fill(Charge);
            Hist_charge2_dist->Fill(Charge2);

            ////////////////////////// Parent Q2 Analysis //////////////////////////
            mQQx_parent += cos(trk_lorentz3.Phi() * nHar);
            mQQy_parent += sin(trk_lorentz3.Phi() * nHar);
            Qcount_parent++;

            if (n_gamma <= 1)
            {
                mQQx_p[n_gamma] = cos(trk_lorentz3.Phi() * nHar);
                mQQy_p[n_gamma] = sin(trk_lorentz3.Phi() * nHar);
                parent_phi[n_gamma] = trk_lorentz3.Phi();
            }

            if (QQ >= 0.8)
            {
                mQQx_parent_QQcut += cos(trk_lorentz3.Phi() * nHar);
                mQQy_parent_QQcut += sin(trk_lorentz3.Phi() * nHar);
                Qcount_parent_QQcut++;
            }

            pTemp_v2_parent->Fill(1, v2e_parent, 1. / temp_eff / temp_proton_eff);
            pTemp_v2_parent->Fill(2, v2w_parent, 1. / temp_eff / temp_proton_eff);

            if (fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
            {
                pTemp_v2_parent->Fill(3, v2e_parent, 1. / temp_eff / temp_proton_eff);
                pTemp_v2_parent->Fill(4, v2w_parent, 1. / temp_eff / temp_proton_eff);
            }

            pTemp_v2_parent->Fill(5, v2e_parent_EPD, 1. / temp_eff / temp_proton_eff);
            pTemp_v2_parent->Fill(6, v2w_parent_EPD, 1. / temp_eff / temp_proton_eff);
            pTemp_v2_parent->Fill(7, v2_parent_EPD1, 1. / temp_eff / temp_proton_eff);

            Hist_v2parent_eta_obs5->Fill(trk_lorentz3.Eta(), v2_final_parent[0], Eweight / temp_eff / temp_proton_eff);
            for (int kl = 0; kl < 3; kl++)
                Hist_v2parent_pt_obs5[kl]->Fill(trk_lorentz3.Pt(), v2_final_parent[kl], Eweight / temp_eff / temp_proton_eff);

            ////////////////////////// End of Parent Q2 Analysis //////////////////////////

            n_gamma++;

            ////////////////////////// Checks to Understand Gamma132 > Gamma112 //////////////////////////
            if (Charge * Charge2 > 0)
            {
                Hist_112_pt_ss->Fill(fabs(Pt - Pt2), correlator4[0], Eweight);
                Hist_132_pt_ss->Fill(fabs(Pt - Pt2), correlator0[0], Eweight);
                Hist_112_cc_pt_ss->Fill(fabs(Pt - Pt2), coscos, Eweight);
                Hist_112_ss_pt_ss->Fill(fabs(Pt - Pt2), sinsin, Eweight);
                Hist_112_eta_ss->Fill(fabs(Eta - Eta2), correlator4[0], Eweight);
                Hist_132_eta_ss->Fill(fabs(Eta - Eta2), correlator0[0], Eweight);
                Hist_112_cc_eta_ss->Fill(fabs(Eta - Eta2), coscos, Eweight);
                Hist_112_ss_eta_ss->Fill(fabs(Eta - Eta2), sinsin, Eweight);

                Hist_112ss_v2_parent_obs2->Fill(v2_final_parent[0], correlator4[0], Eweight);
                Hist_112ss_pt_parent_obs2->Fill(trk_lorentz3.Perp(), correlator4[0], Eweight);
                Hist_132ss_v2_parent_obs2->Fill(v2_final_parent[0], (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                Hist_132ss_pt_parent_obs2->Fill(trk_lorentz3.Perp(), (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                Hist_gamma112ss_QQ->Fill(QQ, correlator4[0], Eweight);
                Hist_gamma112ss_k->Fill(kstar.Mag(), correlator4[0], Eweight);
                Hist_gamma112ss_lpt->Fill(Pt, correlator4[0], Eweight);
                Hist_gamma112ss_ppt->Fill(Pt2, correlator4[0], Eweight);
                Hist_gamma132ss_QQ->Fill(QQ, (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                Hist_112ss_m_parent_obs2->Fill(trk_lorentz3.M(), correlator4[0], Eweight);
                Hist_v2_parent_parent_pT_ss->Fill(trk_lorentz3.Perp(), v2_final_parent[0], Eweight);

                Hist_charge1_ss_dist->Fill(Charge);
                Hist_charge2_ss_dist->Fill(Charge2);

                if (Pt < 1.0)
                {
                    Hist_112ss_v2_parent_obs2_low->Fill(v2_final_parent[0], correlator4[0], Eweight);
                    Hist_132ss_v2_parent_obs2_low->Fill(v2_final_parent[0], (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                }
                else
                {
                    Hist_112ss_v2_parent_obs2_high->Fill(v2_final_parent[0], correlator4[0], Eweight);
                    Hist_132ss_v2_parent_obs2_high->Fill(v2_final_parent[0], (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                }

                Hist_QQ_dist_coupling->Fill(QQ, Eweight / efficiency[cen][lam_index_tmp]);

                n_ss_pairs_pre++;
            }
            else if (Charge * Charge2 < 0)
            {
                Hist_112_pt_os->Fill(fabs(Pt - Pt2), correlator4[0], Eweight);
                Hist_132_pt_os->Fill(fabs(Pt - Pt2), correlator0[0], Eweight);
                Hist_112_cc_pt_os->Fill(fabs(Pt - Pt2), coscos, Eweight);
                Hist_112_ss_pt_os->Fill(fabs(Pt - Pt2), sinsin, Eweight);
                Hist_112_eta_os->Fill(fabs(Eta - Eta2), correlator4[0], Eweight);
                Hist_132_eta_os->Fill(fabs(Eta - Eta2), correlator0[0], Eweight);
                Hist_112_cc_eta_os->Fill(fabs(Eta - Eta2), coscos, Eweight);
                Hist_112_ss_eta_os->Fill(fabs(Eta - Eta2), sinsin, Eweight);

                Hist_112os_v2_parent_obs2->Fill(v2_final_parent[0], correlator4[0], Eweight);
                Hist_112os_pt_parent_obs2->Fill(trk_lorentz3.Perp(), correlator4[0], Eweight);
                Hist_132os_v2_parent_obs2->Fill(v2_final_parent[0], (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                Hist_132os_pt_parent_obs2->Fill(trk_lorentz3.Perp(), (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                Hist_gamma112os_QQ->Fill(QQ, correlator4[0], Eweight);
                Hist_gamma112os_k->Fill(kstar.Mag(), correlator4[0], Eweight);
                Hist_gamma112os_lpt->Fill(Pt, correlator4[0], Eweight);
                Hist_gamma112os_ppt->Fill(Pt2, correlator4[0], Eweight);
                Hist_gamma132os_QQ->Fill(QQ, (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                Hist_112os_m_parent_obs2->Fill(trk_lorentz3.M(), correlator4[0], Eweight);
                Hist_v2_parent_parent_pT_os->Fill(trk_lorentz3.Perp(), v2_final_parent[0], Eweight);

                Hist_charge1_os_dist->Fill(Charge);
                Hist_charge2_os_dist->Fill(Charge2);

                if (trk_lorentz3.Perp() < 1.0)
                {
                    Hist_112os_v2_parent_obs2_low->Fill(v2_final_parent[0], correlator4[0], Eweight);
                    Hist_132os_v2_parent_obs2_low->Fill(v2_final_parent[0], (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                }
                else
                {
                    Hist_112os_v2_parent_obs2_high->Fill(v2_final_parent[0], correlator4[0], Eweight);
                    Hist_132os_v2_parent_obs2_high->Fill(v2_final_parent[0], (correlator0[0] + correlator0_alt[0]) / 2.0, Eweight);
                }

                n_os_pairs_pre++;
            }

            if (particle_option2 == 1)
            {
                if ((Charge2 == 1) && (Charge == 1))
                    Hist_gamma112ss_pi_QQ->Fill(QQ, correlator4[0], Eweight);
                else if ((Charge2 == 1) && (Charge == -1))
                    Hist_gamma112os_pi_QQ->Fill(QQ, correlator4[0], Eweight);
                else if ((Charge2 == -1) && (Charge == 1))
                    Hist_gamma112os_pibar_QQ->Fill(QQ, correlator4[0], Eweight);
                else if ((Charge2 == -1) && (Charge == -1))
                    Hist_gamma112ss_pibar_QQ->Fill(QQ, correlator4[0], Eweight);
            }

            Hist_QQ_dist->Fill(QQ, Eweight);

            double delta_phi = Phi_new[trki] - Phi_pnew[trkj];
            if (delta_phi < -PI)
                delta_phi = delta_phi + 2 * PI;
            if (delta_phi > PI)
                delta_phi = delta_phi - 2 * PI;

            double v2pe_tmp = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_for_new) * 100;
            double v2pw_tmp = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_bac_new) * 100;

            float v2_final2 = 0;

            if (Eta2 > 0)
                v2_final2 = v2pw_tmp;
            else if (Eta2 < 0)
                v2_final2 = v2pe_tmp;

            Hist_v2_p_parent_pT->Fill(trk_lorentz3.Perp(), v2_final2, Eweight);

            if (trk_lorentz3.Perp() <= 1.0)
            {
                Hist_parent_phi_low_pT->Fill(phi_final_parent);
                Hist_lam_v2_low_parent_v2->Fill(v2_final, Eweight);
                Hist_lam_pT_low_parent_v2->Fill(Pt, Eweight);

                Hist_delta_phi_low_parent_pT->Fill(delta_phi);

                Hist_p_v2_low_parent_v2->Fill(v2_final2, Eweight);
                Hist_p_pT_low_parent_v2->Fill(Pt2, Eweight);
            }
            else
            {
                Hist_parent_phi_high_pT->Fill(phi_final_parent);
                Hist_delta_phi_high_parent_pT->Fill(delta_phi);
            }

            if (trk_lorentz3.Perp() <= 0.3)
                Hist_delta_phi_parent_pT_1->Fill(delta_phi);
            else if (trk_lorentz3.Perp() <= 0.6)
                Hist_delta_phi_parent_pT_2->Fill(delta_phi);
            else if (trk_lorentz3.Perp() <= 0.8)
                Hist_delta_phi_parent_pT_3->Fill(delta_phi);
            else if (trk_lorentz3.Perp() <= 1.0)
                Hist_delta_phi_parent_pT_4->Fill(delta_phi);
            else if (trk_lorentz3.Perp() <= 1.2)
                Hist_delta_phi_parent_pT_5->Fill(delta_phi);
            else if (trk_lorentz3.Perp() <= 1.5)
                Hist_delta_phi_parent_pT_6->Fill(delta_phi);
            else
                Hist_delta_phi_parent_pT_7->Fill(delta_phi);

            if (trki == 0)
            {
                if (Charge2 == 1)
                    n_p_used++;
                else if (Charge2 == -1)
                    n_ap_used++;
                n_tp_used++;

                Hist_Pt2->Fill(Pt2, Eweight);
                if (Charge2 == -1)
                    Hist_Pt2_anti->Fill(Pt2, Eweight);
                else if (Charge2 == 1)
                    Hist_Pt2_part->Fill(Pt2, Eweight);
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            if (Charge > 0 && Charge2 > 0)
                FillGamma(1);
            if (Charge < 0 && Charge2 < 0)
                FillGamma(2);
            if (Charge * Charge2 > 0)
            {
                FillGamma(3);
                n_ss_pairs++;
            }
            if (Charge * Charge2 < 0)
            {
                FillGamma(4);
                n_os_pairs++;
            }

            trk_lorentz2.Clear();

            if (debug_2)
                cout << "Here?" << endl;

            if (fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
            {
                Hist_v2_eta_obs3->Fill(Eta, v2, Eweight / eff);
                if (fabs(Eta) < 0.5 && fabs(Eta2) < 0.5)
                {
                    pTemp_v2_noHBT->Fill(1, v2e, 1. / eff);
                    pTemp_v2_noHBT->Fill(2, v2w, 1. / eff);
                    if (Charge > 0)
                    {
                        pTemp_v2_noHBT->Fill(3, v2e, 1. / eff);
                        pTemp_v2_noHBT->Fill(4, v2w, 1. / eff);
                    }
                    if (Charge < 0)
                    {
                        pTemp_v2_noHBT->Fill(5, v2e, 1. / eff);
                        pTemp_v2_noHBT->Fill(6, v2w, 1. / eff);
                    }
                }
            }
            if (debug_2)
                std::cout << "Ending Track 2" << endl;
        } // 2nd track
        trk_lorentz.Clear();
    } // 1st Track

    // cout << "after first analysis loop" << endl;

    Q2_parent = double(mQQx_parent * mQQx_parent + mQQy_parent * mQQy_parent) / 10000.0 / (double(Qcount_parent) / 10000.0 + pow((double(Qcount_parent) / 100.0 * v2_parent_averaged_TPC[Centrality]), 2));
    Q2_parent_QQcut = double(mQQx_parent_QQcut * mQQx_parent_QQcut + mQQy_parent_QQcut * mQQy_parent_QQcut) / 100.0 / (double(Qcount_parent_QQcut) / 100.0 + pow((double(Qcount_parent_QQcut) / 10.0 * v2_parent_averaged_TPC[Centrality]), 2));
    Hist_Q2_parent->Fill(Q2_parent);

    ////////////////////////// QA to understand strange structures in parent Q2 distribution //////////////////////////
    if (n_gamma > 1)
    {
        Hist_Q2_parent_remove1->Fill(Q2_parent);
        Hist_Q2_parent_QQcut->Fill(Q2_parent_QQcut);
        Hist_Q2_parent_vs_Qcount_parent->Fill(Q2_parent, Qcount_parent);

        if ((Q2_parent <= 2) && (Q2_parent >= 1.8))
        {
            debugfile << "mQQx_parent = " << mQQx_parent << "\n";
            debugfile << "mQQy_parent = " << mQQy_parent << "\n";
            debugfile << "double(mQQx_parent * mQQx_parent + mQQy_parent * mQQy_parent) = " << double(mQQx_parent * mQQx_parent + mQQy_parent * mQQy_parent) << "\n";
            debugfile << "Qcount_parent = " << Qcount_parent << "\n";
            debugfile << "v2_parent_averaged_TPC[Centrality] = " << v2_parent_averaged_TPC[Centrality] << "\n";
            debugfile << "(double(Qcount_parent)/100.0 + pow((double(Qcount_parent)/10.0 * v2_parent_averaged_TPC[Centrality]), 2)) = " << (double(Qcount_parent) / 100.0 + pow((double(Qcount_parent) / 10.0 * v2_parent_averaged_TPC[Centrality]), 2)) << "\n";
            debugfile << "Q2_parent = " << Q2_parent << "\n";
            debugfile << "n_lam_used = " << n_lam_used << "\n";
            debugfile << "\n";
        }
    }

    if (n_gamma == 2)
    {
        if (n_tp_used == 2)
        {
            // debugfile << "-----------------------------------------------------------------" << "\n";
            // debugfile << "Lambda Information:" << "\n";
            // debugfile << "momentum vector = [" << lam_px << ", " << lam_py << ", " << lam_pz << "]" << "\n";
            // debugfile << "phi = " << lam_phi << "\n";
            // debugfile << "Daughter 1 Information:" << "\n";
            // debugfile << "momentum vector = [" << lam_dau_px[0] << ", " << lam_dau_py[0] << ", " << lam_dau_pz[0] << "]" << "\n";
            // debugfile << "phi = " << lam_dau_phi[0] << "\n";
            // debugfile << "trk id = " << lam_dau_id[0] << "\n";
            // debugfile << "Daughter 2 Information:" << "\n";
            // debugfile << "momentum vector = [" << lam_dau_px[1] << ", " << lam_dau_py[1] << ", " << lam_dau_pz[1] << "]" << "\n";
            // debugfile << "phi = " << lam_dau_phi[1] << "\n";
            // debugfile << "trk id = " << lam_dau_id[1] << "\n";
            // debugfile << "Proton 1 Information:" << "\n";
            // debugfile << "momentum vector = [" << p_px[0] << ", " << p_py[0] << ", " << p_pz[0] << "]" << "\n";
            // debugfile << "phi = " << p_phi[0] << "\n";
            // debugfile << "trk id = " << p_trkid[0] << "\n";
            // debugfile << "Proton 2 Information:" << "\n";
            // debugfile << "momentum vector = [" << p_px[1] << ", " << p_py[1] << ", " << p_pz[1] << "]" << "\n";
            // debugfile << "phi = " << p_phi[1] << "\n";
            // debugfile << "trk id = " << p_trkid[1] << "\n";
            // debugfile << "\n";
            // debugfile << "phi_pair[0] = " << parent_phi[0] << "\n";
            // debugfile << "phi_pair[1] = " << parent_phi[1] << "\n";
            // debugfile << "cos[0] = " << mQQx_p[0] << "\n";
            // debugfile << "cos[1] = " << mQQx_p[1] << "\n";
            // debugfile << "sin[0] = " << mQQy_p[0] << "\n";
            // debugfile << "sin[1] = " << mQQy_p[1] << "\n";
            // debugfile << "cos^2 = " << mQQx_parent * mQQx_parent << "\n";
            // debugfile << "sin^2 = " << mQQy_parent * mQQy_parent << "\n";
            // debugfile << "v2_pair^2 = " << v2_parent_averaged[Centrality] * v2_parent_averaged[Centrality] << "\n";
            // debugfile << "-----------------------------------------------------------------" << "\n";
            // debugfile << "\n";
            

            // if((Q2_parent <= 1.985) && (Q2_parent >= 1.98))
            // {
            // debugfile << "new event" << "\n";
            // cout << lam_px << " " << lam_py << " " << lam_pz << " " << lam_pt << " " << lam_phinew << " " << p_px[0] << " " << p_py[0] << " " << p_pz[0] << " " << p_px[1] << " " << p_py[1] << " " << p_pz[1] << " " << p_phi[0] << " " << p_phi[1] << " " << p_mom[0] << " " << p_mom[1] << " " << Q2_parent << " " << Qcount_parent << " " << p_charge[0] << " " << p_charge[1] << " " << p_nhitsfit[0] << " " << p_nhitsfit[1] << " " << p_nHitsRatio[0] << " " << p_nHitsRatio[1] << " " << pVz << endl;
            // }

            if (p_charge[0] == p_charge[1])
            {
                Hist_check_trksplitting_mag_ss->Fill(sqrt(pow((p_px[1] - p_px[0]), 2) + pow((p_py[1] - p_py[0]), 2) + pow((p_pz[1] - p_pz[0]), 2)));
                Hist_check_trksplitting_flowvector_ss->Fill((mQQx_p[1] - mQQx_p[0]), (mQQy_p[1] - mQQy_p[0]));

                double tmp_delta_phi = p_phi[1] - p_phi[0];
                if (tmp_delta_phi < -PI)
                    tmp_delta_phi = tmp_delta_phi + 2 * PI;
                if (tmp_delta_phi > PI)
                    tmp_delta_phi = tmp_delta_phi - 2 * PI;

                if ((Q2_parent <= 2) && (Q2_parent >= 1.9))
                {
                    Hist_check_trksplitting_phi_ss->Fill(tmp_delta_phi);

                    if ((p_mom[0] < 1) && (p_mom[1] < 1))
                        Hist_phi_below1_ss->Fill(tmp_delta_phi);
                    else if ((p_mom[0] >= 1) && (p_mom[1] >= 1))
                        Hist_phi_after1_ss->Fill(tmp_delta_phi);
                    else
                        Hist_phi_both1_ss->Fill(tmp_delta_phi);

                    Hist_check_trksplitting_mag_peak_ss->Fill((pow((p_px[1] - p_px[0]), 2) + pow((p_py[1] - p_py[0]), 2) + pow((p_pz[1] - p_pz[0]), 2)));
                }

                Hist_check_trksplitting_phi_all_ss->Fill(tmp_delta_phi);
                if (n_p_used == 2)
                    Hist_check_trksplitting_phi_p->Fill(tmp_delta_phi);
                else if (n_ap_used == 2)
                    Hist_check_trksplitting_phi_ap->Fill(tmp_delta_phi);

                if ((fabs(mQQx_p[1] - mQQx_p[0]) <= 0.02) && (fabs(mQQy_p[1] - mQQy_p[0]) <= 0.02))
                    Hist_check_trksplitting_pmom_ss->Fill((p_px[1] - p_px[0]), (p_py[1] - p_py[0]));
                if (fabs(tmp_delta_phi) <= 0.01)
                    Hist_check_trksplitting_ptrkid_ss->Fill(p_trkid[1] - p_trkid[0]);
                if ((pow((p_px[1] - p_px[0]), 2) + pow((p_py[1] - p_py[0]), 2) + pow((p_pz[1] - p_pz[0]), 2)) <= 0.02)
                    Hist_check_trksplitting_ptrkid2_ss->Fill(p_trkid[1] - p_trkid[0]);

                Hist_Q2_parent_pair_ss->Fill(Q2_parent);
            }
            else
            {
                Hist_check_trksplitting_mag_os->Fill((pow((p_px[1] - p_px[0]), 2) + pow((p_py[1] - p_py[0]), 2) + pow((p_pz[1] - p_pz[0]), 2)));
                Hist_check_trksplitting_flowvector_os->Fill((mQQx_p[1] - mQQx_p[0]), (mQQy_p[1] - mQQy_p[0]));

                double tmp_delta_phi = p_phi[1] - p_phi[0];
                if (tmp_delta_phi < -PI)
                    tmp_delta_phi = tmp_delta_phi + 2 * PI;
                if (tmp_delta_phi > PI)
                    tmp_delta_phi = tmp_delta_phi - 2 * PI;
                if ((Q2_parent <= 2) && (Q2_parent >= 1.9))
                {
                    Hist_check_trksplitting_phi_os->Fill(tmp_delta_phi);

                    if ((p_mom[0] < 1) && (p_mom[1] < 1))
                        Hist_phi_below1_os->Fill(tmp_delta_phi);
                    else if ((p_mom[0] >= 1) && (p_mom[1] >= 1))
                        Hist_phi_after1_os->Fill(tmp_delta_phi);
                    else
                        Hist_phi_both1_os->Fill(tmp_delta_phi);

                    Hist_check_trksplitting_mag_peak_os->Fill((pow((p_px[1] - p_px[0]), 2) + pow((p_py[1] - p_py[0]), 2) + pow((p_pz[1] - p_pz[0]), 2)));
                }

                Hist_check_trksplitting_phi_all_os->Fill(tmp_delta_phi);

                if ((fabs(mQQx_p[1] - mQQx_p[0]) <= 0.02) && (fabs(mQQy_p[1] - mQQy_p[0]) <= 0.02))
                    Hist_check_trksplitting_pmom_os->Fill((p_px[1] - p_px[0]), (p_py[1] - p_py[0]));
                if (tmp_delta_phi == 0)
                    Hist_check_trksplitting_ptrkid_os->Fill(p_trkid[1] - p_trkid[0]);
                if ((pow((p_px[1] - p_px[0]), 2) + pow((p_py[1] - p_py[0]), 2) + pow((p_pz[1] - p_pz[0]), 2)) == 0)
                    Hist_check_trksplitting_ptrkid2_os->Fill(p_trkid[1] - p_trkid[0]);

                Hist_Q2_parent_pair_os->Fill(Q2_parent);
            }

            Hist_Q2_parent_tp->Fill(Q2_parent);

            if (n_p_used == 2)
                Hist_Q2_parent_p->Fill(Q2_parent);
            if (n_ap_used == 2)
                Hist_Q2_parent_ap->Fill(Q2_parent);
            if ((n_p_used != 2) && (n_ap_used != 2))
                Hist_Q2_parent_diff->Fill(Q2_parent);
            if (n_lam_used == 2)
                Hist_Q2_parent_lam->Fill(Q2_parent);
            if ((Q2_parent <= 2) && (Q2_parent >= 1.9))
            {
                nHits_ratio_under_peak->Fill(p_nHitsRatio[0]);
                nHits_ratio_under_peak->Fill(p_nHitsRatio[1]);
            }
        }
    }

    if ((Q2_parent <= 2) && (Q2_parent >= 1.9))
    {
        Hist_ngamma_wrange->Fill(n_gamma);
        Hist_nlam_wrange->Fill(n_lam_used);
        Hist_np_wrange->Fill(n_p_used);
        Hist_nap_wrange->Fill(n_ap_used);
        Hist_ntp_wrange->Fill(n_tp_used);
    }

    Hist_RefMult_Q2_parent->Fill(Q2_parent, Fcount_parent);

    ////////////////////////// QA to understand strange structures in parent Q2 distribution //////////////////////////

    if (opt_weight != 1)
    {
        /*for (int trk = 0; trk < NPTracks; trk++)
        {
            gv_gamma::trk_mom_temp->SetXYZ(all_px->at(trk), all_py->at(trk), all_pz->at(trk));
            float EtaAsso_tmpv2 = gv_gamma::trk_mom_temp->Eta();
            float PtAsso_tmpv2 = gv_gamma::trk_mom_temp->Pt();
            float DCAGlobalAsso_tmpv2 = all_dcaglobal->at(trk);

            if (!IsGoodAsso(PtAsso_tmpv2, EtaAsso_tmpv2, DCAGlobalAsso_tmpv2))
                continue;

            float mQx_j = mQx, mQy_j = mQy;

            TVector2 mQ_j(mQx_j, mQy_j);
            float psi_F2 = mQ_j.Phi() / nHar;
            float psi_F2_new = psi_F2;
            for (int jj = 0; jj < order; jj++)
                psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
            if (psi_F2_new < 0)
                psi_F2_new += PI;
            if (psi_F2_new > PI)
                psi_F2_new -= PI;

            float v2p_temp = cos(nHar * PhiAsso_new[trk] - nHar * psi_F2_new) * 100;
            float v2pe_temp = cos(nHar * PhiAsso_new[trk] - nHar * TPC_EP_for_new) * 100;
            float v2pw_temp = cos(nHar * PhiAsso_new[trk] - nHar * TPC_EP_bac_new) * 100;
            float v2pe_EPD_temp = cos(nHar * PhiAsso_new[trk] - nHar * EPD_EP_east_new) * 100;
            float v2pw_EPD_temp = cos(nHar * PhiAsso_new[trk] - nHar * EPD_EP_west_new) * 100;
            float v2p_EPD1_temp = cos(nHar * PhiAsso_new[trk] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;

            float v2_final[3] = {0.};
            if (EtaAsso_tmpv2 > 0)
                v2_final[0] = v2pw_temp;
            else if (EtaAsso_tmpv2 < 0)
                v2_final[0] = v2pe_temp;
            // v2_final[0] = (v2pe_temp + v2pw_temp)/2.0;
            v2_final[1] = (v2pe_EPD_temp + v2pw_EPD_temp) / 2.0;
            v2_final[2] = v2p_EPD1_temp;

            for (int xk = 0; xk < 3; xk++)
            {
                Hist_v2_pt_obs1_alltrks[xk]->Fill(PtAsso_tmpv2, v2_final[xk]);
                Hist_v2_pt_obs2_alltrks[xk]->Fill(PtAsso_tmpv2, v2_final[xk], Eweight);
            }

            pTemp_v2_chargedhadrons->Fill(3, v2pe_temp);
            pTemp_v2_chargedhadrons->Fill(4, v2pw_temp);
            pTemp_v2_chargedhadrons->Fill(5, v2pe_EPD_temp);
            pTemp_v2_chargedhadrons->Fill(6, v2pw_EPD_temp);
            pTemp_v2_chargedhadrons->Fill(7, v2p_EPD1_temp);

            for (int trk2 = trk + 1; trk2 < NPTracks; trk2++)
            {
                if (debug_4)
                {
                    std::cout << "trk2 = " << trk2 << endl;
                    std::cout << "NPTracks = " << NPTracks << endl;
                }

                TVector3 trk_mom_tmp2(all_px->at(trk2), all_py->at(trk2), all_pz->at(trk2));
                if (!IsGoodAsso(trk_mom_tmp2.Pt(), trk_mom_tmp2.Eta(), all_dcaglobal->at(trk2)))
                    continue;

                TVector3 trk_mom_tmp_shifted(gv_gamma::trk_mom_temp->Pt() * TMath::Cos(PhiAsso_new[trk]), gv_gamma::trk_mom_temp->Pt() * TMath::Sin(PhiAsso_new[trk]), gv_gamma::trk_mom_temp->Pz());
                TVector3 trk_mom_tmp_shifted2(trk_mom_tmp2.Pt() * TMath::Cos(PhiAsso_new[trk2]), trk_mom_tmp2.Pt() * TMath::Sin(PhiAsso_new[trk2]), trk_mom_tmp2.Pz());
                TVector3 trk_mom_parent = trk_mom_tmp_shifted + trk_mom_tmp_shifted2;

                float mQx1_tmp = mQx1;
                float mQy1_tmp = mQy1;
                float mQx2_tmp = mQx2;
                float mQy2_tmp = mQy2;
                bool checked1 = false, checked2 = false;
                for (int xyz = 0; xyz < mQ1_list.size(); xyz++)
                {
                    if ((!checked1) && (trk == mQ1_list.at(xyz)))
                    {
                        mQx1_tmp -= PtAsso_tmpv2 * cos(PhiAsso_new[trk] * nHar);
                        mQy1_tmp -= PtAsso_tmpv2 * sin(PhiAsso_new[trk] * nHar);

                        checked1 = true;
                    }

                    if ((!checked2) && (trk2 == mQ1_list.at(xyz)))
                    {
                        mQx1_tmp -= trk_mom_tmp2.Pt() * cos(PhiAsso_new[trk2] * nHar);
                        mQy1_tmp -= trk_mom_tmp2.Pt() * sin(PhiAsso_new[trk2] * nHar);

                        checked2 = true;
                    }

                    if (checked1 && checked2)
                        break;
                }
                checked1 = false, checked2 = false;
                for (int xyz2 = 0; xyz2 < mQ2_list.size(); xyz2++)
                {
                    if ((!checked1) && (trk == mQ2_list.at(xyz2)))
                    {
                        mQx2_tmp -= PtAsso_tmpv2 * cos(PhiAsso_new[trk] * nHar);
                        mQy2_tmp -= PtAsso_tmpv2 * sin(PhiAsso_new[trk] * nHar);

                        checked1 = true;
                    }

                    if ((!checked2) && (trk2 == mQ2_list.at(xyz2)))
                    {
                        mQx2_tmp -= trk_mom_tmp2.Pt() * cos(PhiAsso_new[trk2] * nHar);
                        mQy2_tmp -= trk_mom_tmp2.Pt() * sin(PhiAsso_new[trk2] * nHar);

                        checked2 = true;
                    }

                    if (checked1 && checked2)
                        break;
                }

                TVector2 mQ1_tmp(mQx1_tmp, mQy1_tmp);
                float TPC_EP_east_temp = mQ1_tmp.Phi() / nHar;
                float TPC_EP_east_temp_new = TPC_EP_east_temp;
                for (int jj = 0; jj < order; jj++)
                    TPC_EP_east_temp_new += -2 * PsiMean_E[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_east_temp) / nHar / (jj + 1) + 2 * PsiMean_E[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_east_temp) / nHar / (jj + 1);
                if (TPC_EP_east_temp_new < 0)
                    TPC_EP_east_temp_new += PI;
                if (TPC_EP_east_temp_new > PI)
                    TPC_EP_east_temp_new -= PI;

                TVector2 mQ2_tmp(mQx2_tmp, mQy2_tmp);
                float TPC_EP_west_temp = mQ2_tmp.Phi() / nHar;
                float TPC_EP_west_temp_new = TPC_EP_west_temp;
                for (int jj = 0; jj < order; jj++)
                    TPC_EP_west_temp_new += -2 * PsiMean_W[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_west_temp) / nHar / (jj + 1) + 2 * PsiMean_W[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_west_temp) / nHar / (jj + 1);
                if (TPC_EP_west_temp_new < 0)
                    TPC_EP_west_temp_new += PI;
                if (TPC_EP_west_temp_new > PI)
                    TPC_EP_west_temp_new -= PI;

                float v2p_parent_temp = cos(nHar * trk_mom_parent.Phi() - nHar * psi_F2_new) * 100;
                float v2pe_parent_temp = cos(nHar * trk_mom_parent.Phi() - nHar * TPC_EP_east_temp_new) * 100;
                float v2pw_parent_temp = cos(nHar * trk_mom_parent.Phi() - nHar * TPC_EP_west_temp_new) * 100;
                float v2pe_parent_EPD_temp = cos(nHar * trk_mom_parent.Phi() - nHar * EPD_EP_east_new) * 100;
                float v2pw_parent_EPD_temp = cos(nHar * trk_mom_parent.Phi() - nHar * EPD_EP_west_new) * 100;
                float v2p_parent_EPD1_temp = cos(nHar * trk_mom_parent.Phi() - EPD_EP1_east_new - EPD_EP1_west_new) * 100;

                pTemp_v2_chargedhadrons->Fill(1, v2pe_parent_temp);
                pTemp_v2_chargedhadrons->Fill(2, v2pw_parent_temp);
                pTemp_v2_chargedhadrons->Fill(11, v2p_parent_EPD1_temp);
                pTemp_v2_chargedhadrons->Fill(12, v2pe_parent_EPD_temp);
                pTemp_v2_chargedhadrons->Fill(13, v2pw_parent_EPD_temp);

                if (debug_4)
                    std::cout << "Ending this new loop?" << endl;
            }
        }*/

        /*for (int idx_pion = 0; idx_pion < p_pion_list->size(); idx_pion++)
        {
            TVector3 p_pion_list_tmp1 = p_pion_list->at(idx_pion);

            for (int idx_pion2 = idx_pion + 1; idx_pion2 < p_pion_list->size(); idx_pion2++)
            {
                TVector3 p_pion_list_tmp2 = p_pion_list->at(idx_pion2);

                TVector3 pair_pion = p_pion_list_tmp1 + p_pion_list_tmp2;

                float mQx_j = mQx;
                float mQy_j = mQy;
                TVector2 mQ_j(mQx_j, mQy_j);
                float psi_F2 = mQ_j.Phi() / nHar;
                float psi_F2_new = psi_F2;
                for (int jj = 0; jj < order; jj++)
                    psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
                if (psi_F2_new < 0)
                    psi_F2_new += PI;
                if (psi_F2_new > PI)
                    psi_F2_new -= PI;

                float v2p_pion_temp = cos(nHar * pair_pion.Phi() - nHar * psi_F2_new) * 100;
                float v2pe_pion_temp = cos(nHar * pair_pion.Phi() - nHar * TPC_EP_east_new) * 100;
                float v2pw_pion_temp = cos(nHar * pair_pion.Phi() - nHar * TPC_EP_west_new) * 100;
                float v2pe_pion_EPD_temp = cos(nHar * pair_pion.Phi() - nHar * EPD_EP_east_new) * 100;
                float v2pw_pion_EPD_temp = cos(nHar * pair_pion.Phi() - nHar * EPD_EP_west_new) * 100;
                float v2p_pion_EPD1_temp = cos(nHar * pair_pion.Phi() - EPD_EP1_east_new - EPD_EP1_west_new) * 100;

                pTemp_v2_chargedhadrons->Fill(14, v2pe_pion_temp);
                pTemp_v2_chargedhadrons->Fill(15, v2pw_pion_temp);
                pTemp_v2_chargedhadrons->Fill(16, v2p_pion_EPD1_temp);
                pTemp_v2_chargedhadrons->Fill(17, v2pe_pion_EPD_temp);
                pTemp_v2_chargedhadrons->Fill(18, v2pw_pion_EPD_temp);

                // float v2_final_pion[3] = {0.};
                // v2_final_pion[0] = v2p_pion_temp;
                // v2_final_pion[1] = (v2pe_pion_EPD_temp + v2pw_pion_EPD_temp) * 0.5;
                // v2_final_pion[2] = v2p_pion_EPD1_temp;
                float v2_final_pion = 0;
                v2_final_pion = cos(nHar * (p_pion_list_tmp1.Phi() - p_pion_list_tmp2.Phi())) * 100.0;

                for (int mk = 0; mk < 3; mk++)
                {
                    // Hist_v2pion_pt_obs5[mk]->Fill(pair_pion.Pt(), v2_final_pion[mk], Eweight);
                    Hist_v2pion_pt_obs5[mk]->Fill(pair_pion.Pt(), v2_final_pion, Eweight);
                }

                v2_2_pion->Fill(0.5, v2_final_pion, Eweight);
            }
        }*/
    }

    WrapUpESE();

    ///////////////////////////////// Round 2: BKG Tracks /////////////////////////////////

    if (debug_1)
        std::cout << "Bkg Track Loop" << endl;
    int n_lam_rot_used = 0, n_proton_rot_used = 0, n_gamma_rot = 0, n_p_rot_used = 0, n_ap_rot_used = 0;

    int Qcount_parent_rot = 0;
    mQQx_parent_rot = 0., mQQy_parent_rot = 0., Q2_parent_rot = 0.;
    int Qcount_parent_QQcut_rot = 0;
    mQQx_parent_QQcut_rot = 0., mQQy_parent_QQcut_rot = 0., Q2_parent_QQcut_rot = 0.;

    // loop for the real analysis
    for (int trki = 0; trki < n_particle1_rot; trki++)
    {
        if (debug_1)
            std::cout << "Begin Bkg Track Loop" << endl;

        gv_gamma::trk_mom->SetXYZ(particle1_rot->at(trki).px, particle1_rot->at(trki).py, particle1_rot->at(trki).pz);
        Eta_rot = gv_gamma::trk_mom->Eta();
        Pt_rot = gv_gamma::trk_mom->Pt();
        Charge_rot = particle1_rot->at(trki).Charge; // 1 for lambda, -1 for antilambdas
        Phi_rot = gv_gamma::trk_mom->Phi();
        Theta_rot = 2. * atan(exp(-Eta));
        eff_rot = 1; // Efficiency information later
        float eff_tof = 1;
        eff *= eff_tof;

        ////// for checks systematic error cuts //////
        nHitsFitQA->Fill(particle1_dau_rot->at(2 * trki).nHitsFit);
        nHitsFitQA->Fill(particle1_dau_rot->at(2 * trki + 1).nHitsFit);
        nSigma_dauQA->Fill(particle1_dau_rot->at(2 * trki).nSigma);
        nSigma_dauQA->Fill(particle1_dau_rot->at(2 * trki + 1).nSigma);
        dca_protonQA->Fill(particle1_dau_rot->at(2 * trki).dcaglobal);
        dca_pionQA->Fill(particle1_dau_rot->at(2 * trki + 1).dcaglobal);
        dca_LambdaQA->Fill(particle1_rot->at(trki).dcaglobal);
        //////////////////////////////////////////////////////////////////

        if ((particle_option1 == 0) && (Pt_rot < 2.2))
        {
            if (debug_1)
                std::cout << "Reading eff bkg" << endl;
            if (Charge_rot == 1)
            {
                temp_eff_rot = efficiency[Centrality][(int)floor((Pt_rot) / 0.1)];

                if (debug_1)
                    std::cout << "efficiency[" << Centrality << "][" << (int)floor((Pt_rot) / 0.1) << "] = " << efficiency[Centrality][(int)floor((Pt_rot) / 0.1)] << endl;
                if (debug_1)
                    std::cout << "temp eff bkg = " << temp_eff_rot << endl;
            }
            else if (Charge_rot == -1)
            {
                temp_eff_rot = efficiency_anti[Centrality][(int)floor((Pt_rot) / 0.1)];

                if (debug_1)
                    std::cout << "efficiency_anti[" << Centrality << "][" << (int)floor((Pt_rot) / 0.1) << "] = " << efficiency_anti[Centrality][(int)floor((Pt_rot) / 0.1)] << endl;
                if (debug_1)
                    std::cout << "temp eff bkg = " << temp_eff_rot << endl;
            }

            if (debug_1)
                std::cout << "Charge_rot = " << Charge_rot << endl;
        }

        if (debug_1)
            std::cout << "temp eff bkg = " << temp_eff_rot << endl;

        if (temp_eff_rot < 0)
            continue;
        //////////////////////////////////////////////////////////////////

        int lam_index_tmp = 0;
        if (Pt_rot < 2.0)
            lam_index_tmp = (int)floor((Pt_rot - 0.5) / 0.1);
        else if (Pt_rot == 2.0)
            lam_index_tmp = 14;

        double mass_temp = 0;
        if ((particle_option1 == 0))
            mass_temp = particle1_rot->at(trki).mass;
        else if (particle_option1 == 1)
            mass_temp = 0.93827;
        TLorentzVector trk_lorentz(particle1_rot->at(trki).px, particle1_rot->at(trki).py, particle1_rot->at(trki).pz, sqrt(pow(gv_gamma::trk_mom->Mag(), 2) + pow(mass_temp, 2)));
        TLorentzVector lA;
        lA.SetPxPyPzE(Pt_rot * cos(Phi_new_rot[trki]), Pt_rot * sin(Phi_new_rot[trki]), Pt_rot * sinh(Eta_rot), sqrt(mass_temp * mass_temp + Pt_rot * Pt_rot * cosh(Eta_rot) * cosh(Eta_rot)));

        V0Mass_gamma112_rot->Fill(particle1_rot->at(trki).mass);
        hEtaPt_rot_Dist->Fill(Eta_rot, Pt_rot, Eweight);
        Hist_Pt_rot->Fill(Pt_rot, Eweight);

        float mQx_i = mQx, mQy_i = mQy;
        TVector2 mQ_i(mQx_i, mQy_i);
        float psi_F = mQ_i.Phi() / nHar;
        Hist_TPC_EP_full_m1_rot->Fill(psi_F);
        float psi_F_new = psi_F;
        for (int jj = 0; jj < order; jj++)
            psi_F_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F) / nHar / (jj + 1);
        Hist_TPC_EP_full_m1_flat_rot->Fill(psi_F_new);

        ////// obtaining v2 for Lambdas //////
        v2_sub_rot = (Eta_rot > 0) ? cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_bac_new) * 100 : cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_for_new) * 100;
        v2_rot = cos(nHar * Phi_new_rot[trki] - nHar * psi_F_new) * 100;
        v2e_rot = cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_for_new) * 100;
        v2w_rot = cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_bac_new) * 100;
        v2epos_rot = cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_for_pos_new) * 100;
        v2eneg_rot = cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_for_neg_new) * 100;
        v2wpos_rot = cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_bac_pos_new) * 100;
        v2wneg_rot = cos(nHar * Phi_new_rot[trki] - nHar * TPC_EP_bac_neg_new) * 100;
        v2_EPD1_rot = cos(2. * Phi_new_rot[trki] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
        v2_EPDe_rot = cos(nHar * Phi_new_rot[trki] - nHar * EPD_EP_east_new) * 100;
        v2_EPDw_rot = cos(nHar * Phi_new_rot[trki] - nHar * EPD_EP_west_new) * 100;

        Fillv2_rot();

        float v2_final = 0;
        if (Eta_rot > 0)
        {
            v2_final = v2w_rot;

            if (Charge_rot > 0)
                v2_sub_ss_rot = v2wpos_rot;
            else if (Charge_rot < 0)
                v2_sub_os_rot = v2wneg_rot;
        }
        else if (Eta_rot < 0)
        {
            v2_final = v2e_rot;

            if (Charge_rot > 0)
                v2_sub_ss_rot = v2epos_rot;
            else if (Charge_rot < 0)
                v2_sub_os_rot = v2eneg_rot;
        }

        Hist_v2_pt_obs1_rot->Fill(Pt_rot, v2_final);
        Hist_v2_pt_obs2_rot->Fill(Pt_rot, v2_final, Eweight);
        Hist_v2_eta_obs1_rot->Fill(Eta_rot, v2_final, 1.);
        Hist_v2_eta_obs2_rot->Fill(Eta_rot, v2_final, Eweight);
        Hist_v2_pt_EPD1_obs_rot->Fill(Pt_rot, v2_EPD1_rot, Eweight);
        Hist_v2_pt_EPD_obs_rot->Fill(Pt_rot, 0.5 * (v2_EPDe_rot + v2_EPDw_rot), Eweight);
        Hist_v2_eta_EPD1_obs_rot->Fill(Eta_rot, v2_EPD1_rot, Eweight / eff);
        Hist_v2_eta_EPD_obs_rot->Fill(Eta_rot, 0.5 * (v2_EPDe_rot + v2_EPDw_rot), Eweight / eff);
        Hist_v2_eta_EPDe_obs_rot->Fill(Eta_rot, v2_EPDe_rot, Eweight / eff);
        Hist_v2_eta_EPDw_obs_rot->Fill(Eta_rot, v2_EPDw_rot, Eweight / eff);
        //////////////////////////////////////////////////////////////////

        n_lam_rot_used++;

        if (opt_weight == 1)
            continue;
        for (int trkj = 0; trkj < n_particle2; trkj++)
        {
            Charge2 = particle2->at(trkj).Charge; // 1 for protons, -1 for antiprotons
            gv_gamma::trk_mom2->SetXYZ(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz);
            Eta2 = gv_gamma::trk_mom2->Eta();
            Pt2 = gv_gamma::trk_mom2->Pt();
            DCAGlobal2 = particle2->at(trkj).dcaglobal;
            Phi2 = gv_gamma::trk_mom2->Phi();
            eff2 = 1; // efficiency information acquired later
            float eff_tof2 = 1;
            eff2 *= eff_tof2;

            if (debug_2)
                std::cout << "Starting Track 2" << endl;

            if ((particle_option1 == 0) && ((particle1_dau_rot->at(2 * trki).trk_id == particle2->at(trkj).trk_id) || (particle1_dau_rot->at(2 * trki + 1).trk_id == particle2->at(trkj).trk_id)))
                continue;

            nSigma_prim_protonQA->Fill(particle2->at(trkj).nsigma);

            if ((particle_option1 == 1) && (particle_option2 == 0) && (particle1->at(trki).trk_id == particle2->at(trkj).trk_id))
                continue;

            TLorentzVector trk_lorentz2(particle2->at(trkj).px, particle2->at(trkj).py, particle2->at(trkj).pz, sqrt(pow(gv_gamma::trk_mom2->Mag(), 2) + pow(particle2->at(trkj).mass, 2)));
            rapidity2 = trk_lorentz2.Rapidity();
            TLorentzVector lB;
            lB.SetPxPyPzE(Pt2 * cos(Phi_pnew[trkj]), Pt2 * sin(Phi_pnew[trkj]), Pt2 * sinh(Eta2), sqrt(particle2->at(trkj).mass * particle2->at(trkj).mass + Pt2 * Pt2 * cosh(Eta2) * cosh(Eta2)));
            if (debug_2)
                std::cout << "rapidity2 = " << rapidity2 << endl;

            ///////////// Lambda Proton Invariant Mass Calculations!! /////////////
            TLorentzVector trk_lorentz3;
            trk_lorentz3 = lA + lB;
            //////////////////////////////////////////////////////////////////////////////

            if (debug_2)
                std::cout << "Pt2 = " << Pt2 << endl;
            hDpt_rot->Fill(fabs(Pt_rot - Pt2), Eweight);

            float mQx_j = mQx_i, mQy_j = mQy_i;
            TVector2 mQ_j(mQx_j, mQy_j);
            float psi_F2 = mQ_j.Phi() / nHar;
            Hist_TPC_EP_full_m2_rot->Fill(psi_F2);
            float psi_F2_new = psi_F2;
            for (int jj = 0; jj < order; jj++)
                psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
            if (psi_F2_new < 0)
                psi_F2_new += PI;
            if (psi_F2_new > PI)
                psi_F2_new -= PI;
            Hist_TPC_EP_full_m2_flat_rot->Fill(psi_F2_new);

            if (particle_option2 == 0)
            {
                temp_proton_eff_rot = fit_eff_pT_final[Centrality]->Eval(Pt2);
                temp_proton_eff_rot *= fit_tof_eff_pT_final[Centrality]->Eval(Pt2);
            }
            else
            {
                temp_proton_eff_rot = 1;
            }

            // correlations
            correlator3_rot = cos(Phi_new_rot[trki] - Phi_pnew[trkj]) * 100;

            // TPC
            correlator4e_rot[0] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - 2 * TPC_EP_east_new) * 100;
            correlator4w_rot[0] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - 2 * TPC_EP_west_new) * 100;
            correlator4_rot[0] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - 2 * psi_F2_new) * 100;
            correlator0_rot[0] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + 2 * psi_F2_new) * 100;
            correlator0e_rot[0] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + 2 * TPC_EP_east_new) * 100;
            correlator0w_rot[0] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + 2 * TPC_EP_west_new) * 100;
            correlator0_alt_rot[0] = cos(Phi_pnew[trkj] - 3 * Phi_new_rot[trki] + 2 * psi_F2_new) * 100;
            correlator0e_alt_rot[0] = cos(Phi_pnew[trkj] - 3 * Phi_new_rot[trki] + 2 * TPC_EP_east_new) * 100;
            correlator0w_alt_rot[0] = cos(Phi_pnew[trkj] - 3 * Phi_new_rot[trki] + 2 * TPC_EP_west_new) * 100;

            // EPD
            correlator4e_rot[1] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - 2 * EPD_EP_east_new) * 100;
            correlator4w_rot[1] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - 2 * EPD_EP_west_new) * 100;
            correlator4_rot[1] = 0.5 * (correlator4e_rot[1] + correlator4w_rot[1]);
            correlator0e_rot[1] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + 2 * EPD_EP_east_new) * 100;
            correlator0w_rot[1] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + 2 * EPD_EP_west_new) * 100;
            correlator0_rot[1] = 0.5 * (correlator0e_rot[1] + correlator0w_rot[1]);
            correlator0e_alt_rot[1] = cos(Phi_pnew[trkj] - 3 * Phi_new_rot[trki] + 2 * EPD_EP_east_new) * 100;
            correlator0w_alt_rot[1] = cos(Phi_pnew[trkj] - 3 * Phi_new_rot[trki] + 2 * EPD_EP_west_new) * 100;
            correlator0_alt_rot[1] = 0.5 * (correlator0e_alt_rot[1] + correlator0w_alt_rot[1]);

            // EPD 1st Order
            correlator4_rot[2] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
            correlator4e_rot[2] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
            correlator4w_rot[2] = cos(Phi_new_rot[trki] + Phi_pnew[trkj] - EPD_EP1_east_new - EPD_EP1_west_new) * 100;
            correlator0_rot[2] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0e_rot[2] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0w_rot[2] = cos(Phi_new_rot[trki] - 3 * Phi_pnew[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0_alt_rot[2] = cos(Phi_pnew[trki] - 3 * Phi_new_rot[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0e_alt_rot[2] = cos(Phi_pnew[trki] - 3 * Phi_new_rot[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;
            correlator0w_alt_rot[2] = cos(Phi_pnew[trki] - 3 * Phi_new_rot[trkj] + EPD_EP1_east_new + EPD_EP1_west_new) * 100;

            float coscos = cos(Phi_new_rot[trki] - psi_F2_new) * cos(Phi_pnew[trkj] - psi_F2_new);
            float sinsin = sin(Phi_new_rot[trki] - psi_F2_new) * sin(Phi_pnew[trkj] - psi_F2_new);

            if (debug_2)
                cout << "Gamma112 = (correlator4_rot) " << correlator4_rot << endl;
            if (debug_2)
                cout << "Gamma112 = (correlator4e_rot) " << correlator4e_rot << endl;
            if (debug_2)
                cout << "Gamma112 = (correlator4w_rot) " << correlator4w_rot << endl;

            if (trki == 0)
                Hist_Pt2_rot->Fill(Pt2, Eweight);

            ///////////// Lambda-Proton Parent v2 calculations for ESE /////////////
            double v2e_parent = cos(nHar * trk_lorentz3.Phi() - nHar * TPC_EP_for_new) * 100;
            double v2w_parent = cos(nHar * trk_lorentz3.Phi() - nHar * TPC_EP_bac_new) * 100;
            double v2e_parent_EPD = cos(nHar * trk_lorentz3.Phi() - nHar * EPD_EP_east_new) * 100;
            double v2w_parent_EPD = cos(nHar * trk_lorentz3.Phi() - nHar * EPD_EP_west_new) * 100;
            double v2_parent_EPD1 = cos(nHar * trk_lorentz3.Phi() - EPD_EP1_west_new - EPD_EP1_east_new) * 100;

            double v2_final_parent[3] = {0.};
            double phi_final_parent = 0;

            if (trk_lorentz3.Eta() > 0)
            {
                v2_final_parent[0] = v2w_parent;
                phi_final_parent = trk_lorentz3.Phi() - TPC_EP_bac_new;
            }
            else
            {
                v2_final_parent[0] = v2e_parent;
                phi_final_parent = trk_lorentz3.Phi() - TPC_EP_for_new;
            }
            v2_final_parent[1] = 0.5 * (v2e_parent_EPD + v2w_parent_EPD);
            v2_final_parent[2] = v2_parent_EPD1;

            if (phi_final_parent < -PI)
                phi_final_parent = phi_final_parent + 2 * PI;
            ///////////////////////////////////////////////////////////////////////////////////////////

            TVector3 kstar((trk_lorentz2.Px() - trk_lorentz.Px()) / 2.0, (trk_lorentz2.Py() - trk_lorentz.Py()) / 2.0, (trk_lorentz2.Pz() - trk_lorentz.Pz()) / 2.0);
            QQ_rot = pow(trk_lorentz2.E() - trk_lorentz.E(), 2) - pow(trk_lorentz2.Px() - trk_lorentz.Px(), 2) - pow(trk_lorentz2.Py() - trk_lorentz.Py(), 2) - pow(trk_lorentz2.Pz() - trk_lorentz.Pz(), 2);
            if (QQ_rot <= 0)
                QQ_rot = sqrt(-QQ_rot);
            else
                QQ_rot = -99;

            ////////////////////////// Parent Q2 Analysis //////////////////////////
            mQQx_parent_rot += cos(trk_lorentz3.Phi() * nHar);
            mQQy_parent_rot += sin(trk_lorentz3.Phi() * nHar);
            Qcount_parent_rot++;

            if (QQ >= 0.8)
            {
                mQQx_parent_QQcut_rot += cos(trk_lorentz3.Phi() * nHar);
                mQQy_parent_QQcut_rot += sin(trk_lorentz3.Phi() * nHar);
                Qcount_parent_QQcut_rot++;
            }

            pTemp_v2_parent_rot->Fill(1, v2e_parent, 1. / temp_eff_rot / temp_proton_eff_rot);
            pTemp_v2_parent_rot->Fill(2, v2w_parent, 1. / temp_eff_rot / temp_proton_eff_rot);

            if (fabs(Pt_rot - Pt2) > 0.15 && fabs(Eta_rot - Eta2) > 0.15)
            {
                pTemp_v2_parent_rot->Fill(3, v2e_parent, 1. / temp_eff_rot / temp_proton_eff_rot);
                pTemp_v2_parent_rot->Fill(4, v2w_parent, 1. / temp_eff_rot / temp_proton_eff_rot);
            }

            pTemp_v2_parent_rot->Fill(5, v2e_parent_EPD, 1. / temp_eff_rot / temp_proton_eff_rot);
            pTemp_v2_parent_rot->Fill(6, v2w_parent_EPD, 1. / temp_eff_rot / temp_proton_eff_rot);
            pTemp_v2_parent_rot->Fill(7, v2_parent_EPD1, 1. / temp_eff_rot / temp_proton_eff_rot);

            Hist_v2parent_eta_obs5_rot->Fill(trk_lorentz3.Eta(), v2_final_parent[0], Eweight / temp_eff_rot / temp_proton_eff_rot);
            for (int kl = 0; kl < 3; kl++)
                Hist_v2parent_pt_obs5_rot[kl]->Fill(trk_lorentz3.Pt(), v2_final_parent[kl], Eweight / temp_eff_rot / temp_proton_eff_rot);

            n_gamma_rot++;
            ////////////////////////// End of Parent Q2 Analysis //////////////////////////

            double delta_phi = Phi_new_rot[trki] - Phi_pnew[trkj];
            if (delta_phi < -PI)
                delta_phi = delta_phi + 2 * PI;
            if (delta_phi > PI)
                delta_phi = delta_phi - 2 * PI;

            double v2pe_tmp = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_for_new) * 100;
            double v2pw_tmp = cos(nHar * Phi_pnew[trkj] - nHar * TPC_EP_bac_new) * 100;

            float v2_final2 = 0;

            if (Eta2 > 0)
                v2_final2 = v2pw_tmp;
            else if (Eta2 < 0)
                v2_final2 = v2pe_tmp;

            if (Charge2 == 1)
                n_p_rot_used++;
            else if (Charge2 == -1)
                n_ap_rot_used++;

            ////////////////////////// Filling Gamma for Background //////////////////////////
            if (debug_2)
                std::cout << "Before FillGamma_rot " << endl;

            if (Charge_rot > 0 && Charge2 > 0)
                FillGamma_rot(1);
            if (Charge_rot < 0 && Charge2 < 0)
                FillGamma_rot(2);
            if (Charge_rot * Charge2 > 0)
                FillGamma_rot(3);
            if (Charge_rot * Charge2 < 0)
                FillGamma_rot(4);

            if (debug_2)
                std::cout << "After FillGamma_rot " << endl;

            trk_lorentz2.Clear();

            if (debug_2)
                std::cout << "Ending Track 2" << endl;
        } // 2nd track

        trk_lorentz.Clear();
    } // 1st Track

    Q2_parent_rot = double(mQQx_parent_rot * mQQx_parent_rot + mQQy_parent_rot * mQQy_parent_rot) / 100.0 / (double(Qcount_parent_rot) / 100.0 + pow((double(Qcount_parent_rot) / 10.0 * v2_parent_averaged_rot[Centrality]), 2));
    Hist_Q2_parent_rot->Fill(Q2_parent_rot);
    Q2_parent_QQcut_rot = double(mQQx_parent_QQcut_rot * mQQx_parent_QQcut_rot + mQQy_parent_QQcut_rot * mQQy_parent_QQcut_rot) / 100.0 / (double(Qcount_parent_QQcut_rot) / 100.0 + pow((double(Qcount_parent_QQcut_rot) / 10.0 * v2_parent_averaged_rot[Centrality]), 2));
    Hist_Q2_parent_QQcut_rot->Fill(Q2_parent_QQcut_rot);
    Hist_RefMult_Q2_parent_rot->Fill(Q2_parent_rot, Fcount_parent);

    WrapUpESE_rot();

    num_gamma_ss_pre->Fill(n_ss_pairs_pre);
    num_gamma_os_pre->Fill(n_os_pairs_pre);
    num_gamma_ss_final->Fill(n_ss_pairs);
    num_gamma_os_final->Fill(n_os_pairs);

    num_lam_final->Fill(n_lam_used);
    num_proton_final->Fill(n_tp_used);
    num_gamma_final->Fill(n_gamma);

    // num_proton_used->Fill(n_p_used);
    num_proton_used->Fill(n_proton_used);
    num_antiproton_used->Fill(n_ap_used);

    PhiAsso_new.clear();
    Phi_new.clear();
    Phi_pnew.clear();
    Phi_new_rot.clear();
    pTemp_v2->Reset();
    pTemp_v2_noHBT->Reset();
    pTemp_parity_e->Reset();
    pTemp_parity_e_noHBT->Reset();
    pTemp_parity_w->Reset();
    pTemp_parity_w_noHBT->Reset();
    pTemp_delta->Reset();
    pTemp_delta_noHBT->Reset();
    pTemp_v2_chargedhadrons->Reset();

    pTemp_v2_parent->Reset();

    pTemp_v2_rot->Reset();
    pTemp_parity_e_rot->Reset();
    pTemp_parity_w_rot->Reset();
    pTemp_delta_rot->Reset();

    pTemp_v2_parent_rot->Reset();

    return n_gamma;
}

void initialize_event()
{
    n_os_pairs = 0;
    n_ss_pairs = 0;
    n_os_pairs_pre = 0;
    n_ss_pairs_pre = 0;

    if (particle1->size() != 0)
        particle1->clear();
    if (particle1_dau->size() != 0)
        particle1_dau->clear();
    if (particle2->size() != 0)
        particle2->clear();
    if (particle1_rot->size() != 0)
        particle1_rot->clear();
    if (particle1_dau_rot->size() != 0)
        particle1_dau_rot->clear();

    mQ1_list.clear();
    mQ2_list.clear();
}

void determine_particlelists(const int sys_err_opt)
{
    if (debug_1)
        std::cout << "determine_particlelists" << endl;

    //////////////////////////////////////////////////////////////// Setting First Track ////////////////////////////////////////////////////////////////
    if (particle_option1 == 0) // first track is LAMBDA
    {
        if (debug_1)
            std::cout << "(details->nLambda / 2) = " << n_lam << endl;

        for (int part1 = 0; part1 < n_lam; part1++)
        {
            particle lam_tmp(lambda_px->at(part1), lambda_py->at(part1), lambda_pz->at(part1), lambda_Charge->at(part1), lambda_dcaglobal->at(part1), lambda_mass->at(part1), lambda_trk_id->at(part1), lambda_nsigma->at(part1), lambda_hits_ratio->at(part1), lambda_nhitsfit->at(part1), lambda_nhitsmax->at(part1));
            particle_dau dau1(dau_px->at(2 * part1), dau_py->at(2 * part1), dau_pz->at(2 * part1), dau_trk_id->at(2 * part1), dau_nHitsFit->at(2 * part1), dau_nSigma->at(2 * part1), dau_dcaglobal->at(2 * part1), dau_nHitsMax->at(2 * part1));
            particle_dau dau2(dau_px->at(2 * part1 + 1), dau_py->at(2 * part1 + 1), dau_pz->at(2 * part1 + 1), dau_trk_id->at(2 * part1 + 1), dau_nHitsFit->at(2 * part1 + 1), dau_nSigma->at(2 * part1 + 1), dau_dcaglobal->at(2 * part1 + 1), dau_nHitsMax->at(2 * part1 + 1));

            if ((lam_tmp.mass >= 1.113) && (lam_tmp.mass <= 1.119))
            {
                if ((sys_err_opt == 2) && ((dau1.nHitsFit <= 20) || (dau2.nHitsFit <= 20)))
                    continue;
                if ((sys_err_opt == 3) && ((dau1.nSigma > 3) || (dau2.nSigma > 3)))
                    continue;
                if ((sys_err_opt == 4) && ((dau1.dcaglobal <= 1.0) || (dau2.dcaglobal <= 2.0)))
                    continue;
                if ((sys_err_opt == 5) && ((lam_tmp.dcaglobal >= 2.0)))
                    continue;

                TVector3 trk_mom_tmp(lam_tmp.px, lam_tmp.py, lam_tmp.pz);
                TLorentzVector trk_lorentz(lam_tmp.px, lam_tmp.py, lam_tmp.pz, sqrt(pow(trk_mom_tmp.Mag(), 2) + pow(Lambda_Mass, 2)));
                float rapi = trk_lorentz.Rapidity();
                float pti = trk_mom_tmp.Perp();

                if (!IsGoodLambda(pti, rapi))
                    continue;

                particle1->push_back(lam_tmp);
                particle1_dau->push_back(dau1);
                particle1_dau->push_back(dau2);

                trk_mom_tmp.Clear();
                trk_lorentz.Clear();
            }

            if (!((lam_tmp.mass < 1.09) || (lam_tmp.mass > 1.14)))
            {
                if (!((lam_tmp.mass < 1.125) && (lam_tmp.mass > 1.105)))
                {
                    if ((sys_err_opt == 2) && ((dau1.nHitsFit <= 20) || (dau2.nHitsFit <= 20)))
                        continue;
                    if ((sys_err_opt == 3) && ((dau1.nSigma > 3) || (dau2.nSigma > 3)))
                        continue;
                    if ((sys_err_opt == 4) && ((dau1.dcaglobal <= 1.0) || (dau2.dcaglobal <= 2.0)))
                        continue;
                    if ((sys_err_opt == 5) && ((lam_tmp.dcaglobal >= 2.0)))
                        continue;

                    double lammass_temp = lam_tmp.mass;
                    TVector3 trk_mom_tmp(lam_tmp.px, lam_tmp.py, lam_tmp.pz);
                    TLorentzVector trk_lorentz(lam_tmp.px, lam_tmp.py, lam_tmp.pz, sqrt(pow(trk_mom_tmp.Mag(), 2) + pow(lammass_temp, 2)));
                    float rapi = trk_lorentz.Rapidity();
                    float pti = trk_mom_tmp.Perp();

                    if (!IsGoodLambda(pti, rapi))
                        continue;

                    particle1_rot->push_back(lam_tmp);
                    particle1_dau_rot->push_back(dau1);
                    particle1_dau_rot->push_back(dau2);

                    trk_mom_tmp.Clear();
                    trk_lorentz.Clear();
                }
            }
        }
    }
    else if (particle_option1 == 1)
    {
        for (int part1 = 0; part1 < NPTracks; part1++)
        {
            if (!all_is_proton->at(part1)) continue;

            particle proton_tmp(all_px->at(part1), all_py->at(part1), all_pz->at(part1), all_Charge->at(part1), all_dcaglobal->at(part1), 0.93827, all_trk_id->at(part1), all_nSigmaProton->at(part1), (float)all_nhitsfit->at(part1) / (float)all_nhitsmax->at(part1), all_nhitsfit->at(part1), all_nhitsmax->at(part1));
            
            proton_tmp.mass = 0.93827;
            TVector3 trk_mom2_tmp(proton_tmp.px, proton_tmp.py, proton_tmp.pz);
            TLorentzVector trk_lorentz2(proton_tmp.px, proton_tmp.py, proton_tmp.pz, sqrt(pow(trk_mom2_tmp.Mag(), 2) + pow(proton_tmp.mass, 2)));
            float rap2 = trk_lorentz2.Rapidity();
            if (((trk_mom2_tmp.Perp() < 0.4) || (trk_mom2_tmp.Mag() > 2.0)))
                continue;
            if (!IsGoodPrimary(rap2, proton_tmp.dcaglobal))
                continue;
            trk_lorentz2.Clear();
            trk_mom2_tmp.Clear();

            particle1->push_back(proton_tmp);
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////// Setting Second Track ////////////////////////////////////////////////////////////////
    if (particle_option2 == 0) // second track is PROTON
    {
        for (int part2 = 0; part2 < details->n_Proton; part2++)
        {
            if (!all_is_proton->at(part2)) continue;

            particle proton_tmp(all_px->at(part2), all_py->at(part2), all_pz->at(part2), all_Charge->at(part2), all_dcaglobal->at(part2), 0.93827, all_trk_id->at(part2), all_nSigmaProton->at(part2), (float)all_nhitsfit->at(part2) / (float)all_nhitsmax->at(part2), all_nhitsfit->at(part2), all_nhitsmax->at(part2));
            
            TVector3 trk_mom2_tmp(proton_tmp.px, proton_tmp.py, proton_tmp.pz);
            if ((sys_err_opt == 6) && (trk_mom2_tmp.Pt() <= 0.5))
                continue;
            if ((sys_err_opt == 7) && (proton_tmp.nsigma >= 1.5))
                continue;

            TLorentzVector trk_lorentz2(proton_tmp.px, proton_tmp.py, proton_tmp.pz, sqrt(pow(trk_mom2_tmp.Mag(), 2) + pow(proton_tmp.mass, 2)));
            float rap2 = trk_lorentz2.Rapidity();
            if (((trk_mom2_tmp.Perp() < 0.4) || (trk_mom2_tmp.Mag() > 2.0)))
                continue;
            if (!IsGoodPrimary(rap2, proton_tmp.dcaglobal))
                continue;
            trk_lorentz2.Clear();
            trk_mom2_tmp.Clear();

            proton_tmp.mass = 0.93827;
            particle2->push_back(proton_tmp);
        }
    }
    else if (particle_option2 == 1) // second track is PION
    {
        for (int part2 = 0; part2 < details->n_Pion; part2++)
        {
            if (!all_is_proton->at(part2)) continue;

            particle pion_tmp(all_px->at(part2), all_py->at(part2), all_pz->at(part2), all_Charge->at(part2), all_dcaglobal->at(part2), 0.13957, all_trk_id->at(part2), all_nSigmaPion->at(part2), (float)all_nhitsfit->at(part2) / (float)all_nhitsmax->at(part2), all_nhitsfit->at(part2), all_nhitsmax->at(part2));
            
            TVector3 trk_mom2_tmp(pion_tmp.px, pion_tmp.py, pion_tmp.pz);
            TLorentzVector trk_lorentz2(pion_tmp.px, pion_tmp.py, pion_tmp.pz, sqrt(pow(trk_mom2_tmp.Mag(), 2) + pow(pion_tmp.mass, 2)));
            float rap2 = trk_lorentz2.Rapidity();

            if (((trk_mom2_tmp.Perp() < 0.4) || (trk_mom2_tmp.Mag() > 2.0)))
                continue;

            if (!IsGoodPrimary(rap2, pion_tmp.dcaglobal))
                continue;
            trk_lorentz2.Clear();
            trk_mom2_tmp.Clear();

            pion_tmp.mass = 0.13957;
            particle2->push_back(pion_tmp);
        }
    }
}

bool IsGoodLambda(float p, float rap)
{
    if ((p < 0.5) || (p > 2.0))
        return false;
    if (fabs(rap) > 1)
        return false;

    return true;
}

bool IsGoodPrimary(float r, float d)
{
    if (d > DcaCut)
        return false;
    if (fabs(r) > 1)
        return false;

    return true;
}

bool IsGoodEvent(int c)
{
    hVzDiff->Fill(pVz - VPDvz);
    hVertexZ->Fill(pVz);
    hMult_Vz_new->Fill(RefMult, pVz, Eweight);
    Ref_TOF->Fill(RefMult, TOFMult);
    Ref_Day3->Fill(Day3, RefMult);
    TOF_Day3->Fill(Day3, TOFMult);
    NPT_Day3->Fill(Day3, NPTracks);
    hCentrality->Fill(Centrality); //, Eweight);
    hMult_Vz->Fill(RefMult, TOFMult);
    if ((Centrality != c))
        return false;

    return true;
}

bool CountCharge()
{
    Ntof = 0, Npos = 0, Nneg = 0, Fcount = 0;

    for (int trk = 0; trk < NPTracks; trk++)
    {
        gv_gamma::trk_mom_temp->SetXYZ(all_px->at(trk), all_py->at(trk), all_pz->at(trk));
        EtaAsso = gv_gamma::trk_mom_temp->Eta();
        PtAsso = gv_gamma::trk_mom_temp->Pt();
        ChargeAsso = all_Charge->at(trk);
        DCAGlobalAsso = all_dcaglobal->at(trk);
        nSigma_p = all_nSigmaProton->at(trk);

        if (ChargeAsso > 0 && fabs(EtaAsso) < 1 && PtAsso > 0.15 && DCAGlobalAsso < 1 && !(fabs(nSigma_p) < 3 && PtAsso < 0.4))
            Npos++;
        if (ChargeAsso < 0 && fabs(EtaAsso) < 1 && PtAsso > 0.15 && DCAGlobalAsso < 1 && !(fabs(nSigma_p) < 3 && PtAsso < 0.4))
            Nneg++;
        if (all_isTofTrack->at(trk))
            Ntof++;
        Hist_DCA->Fill(DCAGlobalAsso);
        if (!IsGoodAsso(PtAsso, EtaAsso, DCAGlobalAsso))
            continue;
        Fcount++;
    }
    if (Ntof < 2)
        return false; // at least 2 tracks match TOF
    NPA_Day3->Fill(Day3, Fcount);
    hTally->Fill(10);
    int net_Nch = Npos - Nneg;
    net_Nch_Asym = -99;
    if ((Npos + Nneg) > 0)
        net_Nch_Asym = (Npos - Nneg) / float(Npos + Nneg);
    Hist_positive->Fill(Npos);
    Hist_negative->Fill(Nneg);
    Hist_Ch->Fill(Npos + Nneg);
    Hist_netCh->Fill(net_Nch);
    if (net_Nch_Asym > -99)
    {
        Hist_netChAsym->Fill(net_Nch_Asym);
        p_netChAsym_RefMult->Fill(net_Nch_Asym, RefMult);
        if (net_Nch_Asym < (MeanNetChargeAsym - 0.8 * StdDevNetChargeAsym))
        {
            Hist_netChAsym_bin->Fill(1, net_Nch_Asym);
            net_charge_asym_bin = 0;
        }
        else if (net_Nch_Asym < (MeanNetChargeAsym - (0.3 * StdDevNetChargeAsym)))
        {
            Hist_netChAsym_bin->Fill(2, net_Nch_Asym);
            net_charge_asym_bin = 1;
        }
        else if (net_Nch_Asym < (MeanNetChargeAsym + (0.2 * StdDevNetChargeAsym)))
        {
            Hist_netChAsym_bin->Fill(3, net_Nch_Asym);
            net_charge_asym_bin = 2;
        }
        else if (net_Nch_Asym < (MeanNetChargeAsym + 0.8 * StdDevNetChargeAsym))
        {
            Hist_netChAsym_bin->Fill(4, net_Nch_Asym);
            net_charge_asym_bin = 3;
        }
        else
        {
            Hist_netChAsym_bin->Fill(5, net_Nch_Asym);
            net_charge_asym_bin = 4;
        }
    }
    return true;
}

void MakeTPC_EP()
{
    TVector2 mQ, mQ1, mQ2, mQ3, mQ4, mQ5, mQ6, mQ7, mQ8;
    mQx = 0., mQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0.;
    float mQQx = 0., mQQy = 0., mQx3 = 0., mQy3 = 0., mQx4 = 0., mQy4 = 0., mQx5 = 0., mQy5 = 0., mQx6 = 0., mQy6 = 0., mQx7 = 0., mQy7 = 0., mQx8 = 0., mQy8 = 0., mQQx_Q = 0., mQQy_Q = 0., mQQx_Q_pion = 0., mQQy_Q_pion = 0.;
    Fcount = 0;
    int Qcount = 0, proton_use_count = 0, Qcount_Full = 0, Qcount_pion = 0;
    p_pion_list->clear();
    PhiAsso_new.resize(NPTracks);
    if (debug_1)
        cout << "MakeTPC_EP" << endl;
    for (int trk = 0; trk < NPTracks; trk++)
    {
        gv_gamma::trk_mom_temp->SetXYZ(all_px->at(trk), all_py->at(trk), all_pz->at(trk));
        EtaAsso = gv_gamma::trk_mom_temp->Eta();
        PtAsso = gv_gamma::trk_mom_temp->Pt();
        ChargeAsso = all_Charge->at(trk);
        DCAGlobalAsso = all_dcaglobal->at(trk);
        PhiAsso = gv_gamma::trk_mom_temp->Phi();

        if (!IsGoodAsso(PtAsso, EtaAsso, DCAGlobalAsso))
            continue;

        FillPhiAsso(); // ChargeAsso is needed here
        ShiftPhiAsso(trk);

        // if ((p_all->at(trk).nSigmaProton < -2) || ((p_all->at(trk).nSigmaProton > 2)))
        // {
        // cout << "p_all->at(trk).px = " << p_all->at(trk).px << endl;

        mQQx_Q += cos(PhiAsso_new[trk] * nHar);
        mQQy_Q += sin(PhiAsso_new[trk] * nHar);
        Qcount_Full++;

        if (gv_gamma::iTrack[Fcount] < (Scount))
        {
            mQQx += cos(PhiAsso_new[trk] * nHar);
            mQQy += sin(PhiAsso_new[trk] * nHar);
            Qcount++;
        }

        bool proton_use = false;

        for (int ind_p = 0; ind_p < n_particle2; ind_p++)
        {
            if (all_trk_id->at(trk) == particle2->at(ind_p).trk_id)
                proton_use = true;
        }

        if (particle_option1 == 0)
        {
            for (int ind_dau = 0; ind_dau < n_particle1; ind_dau++)
                if ((all_trk_id->at(trk) == particle1_dau->at(2 * ind_dau).trk_id) || (all_trk_id->at(trk) == particle1_dau->at(2 * ind_dau + 1).trk_id))
                    proton_use = true;
            for (int ind_dau = 0; ind_dau < n_particle1_rot; ind_dau++)
                if ((all_trk_id->at(trk) == particle1_dau_rot->at(2 * ind_dau).trk_id) || (all_trk_id->at(trk) == particle1_dau_rot->at(2 * ind_dau + 1).trk_id))
                    proton_use = true;
        }
        else if (particle_option1 == 1)
        {
            for (int ind_p = 0; ind_p < n_particle1; ind_p++)
                if (all_trk_id->at(trk) == particle1->at(ind_p).trk_id)
                    proton_use = true;
        }

        if (proton_use)
        {
            proton_use_count++;
            continue;
        }

        if (all_is_pion->at(trk))
        {
            TVector3 p_pion_list_tmp(PtAsso * TMath::Cos(PhiAsso_new[trk]), PtAsso * TMath::Sin(PhiAsso_new[trk]), all_pz->at(trk));
            p_pion_list->push_back(p_pion_list_tmp);

            // cout << "is Pion!!" << endl;
        }

        if ((Scount) <= gv_gamma::iTrack[Fcount] < (2 * Scount))
        {
            mQx1 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy1 += PtAsso * sin(PhiAsso_new[trk] * nHar);
            mQ1_list.push_back(trk);
        }
        else
        {
            mQx2 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy2 += PtAsso * sin(PhiAsso_new[trk] * nHar);
            mQ2_list.push_back(trk);
        }
        // }

        mQx += PtAsso * cos(PhiAsso_new[trk] * nHar);
        mQy += PtAsso * sin(PhiAsso_new[trk] * nHar);

        if (EtaAsso > Eta_EP_Cut)
        {
            mQx3 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy3 += PtAsso * sin(PhiAsso_new[trk] * nHar);

            if (ChargeAsso > 0)
            {
                mQx5 += PtAsso * cos(PhiAsso_new[trk] * nHar);
                mQy5 += PtAsso * sin(PhiAsso_new[trk] * nHar);
            }
            else if (ChargeAsso < 0)
            {
                mQx6 += PtAsso * cos(PhiAsso_new[trk] * nHar);
                mQy6 += PtAsso * sin(PhiAsso_new[trk] * nHar);
            }
        }
        if (EtaAsso < -Eta_EP_Cut)
        {
            mQx4 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy4 += PtAsso * sin(PhiAsso_new[trk] * nHar);

            if (ChargeAsso > 0)
            {
                mQx7 += PtAsso * cos(PhiAsso_new[trk] * nHar);
                mQy7 += PtAsso * sin(PhiAsso_new[trk] * nHar);
            }
            else if (ChargeAsso < 0)
            {
                mQx8 += PtAsso * cos(PhiAsso_new[trk] * nHar);
                mQy8 += PtAsso * sin(PhiAsso_new[trk] * nHar);
            }
        }

        Fcount++;
    }

    // for (int idx_pion = 0; idx_pion < p_pion_list->size(); idx_pion++)
    // {
    //     TVector3 p_pion_list_tmp1 = p_pion_list->at(idx_pion);

    //     for (int idx_pion2 = idx_pion + 1; idx_pion2 < p_pion_list->size(); idx_pion2++)
    //     {
    //         TVector3 p_pion_list_tmp2 = p_pion_list->at(idx_pion2);

    //         TVector3 pair_pion = p_pion_list_tmp1 + p_pion_list_tmp2;

    //         mQQx_Q_pion += cos(pair_pion.Phi() * nHar);
    //         mQQy_Q_pion += sin(pair_pion.Phi() * nHar);
    //         Qcount_pion++;

    //         float v2_final_pion = 0;
    //         v2_final_pion = cos(nHar * (p_pion_list_tmp1.Phi() - p_pion_list_tmp2.Phi())) * 100.0;

    //         v2_2_pion->Fill(0.5, v2_final_pion, Eweight);
    //     }
    // }

    if (debug_1)
        cout << "After Loop" << endl;
    mQ.Set(mQx, mQy);
    mQ1.Set(mQx1, mQy1);
    mQ2.Set(mQx2, mQy2);
    mQ3.Set(mQx3, mQy3);
    mQ4.Set(mQx4, mQy4);
    mQ5.Set(mQx5, mQy5);
    mQ6.Set(mQx6, mQy6);
    mQ7.Set(mQx7, mQy7);
    mQ8.Set(mQx8, mQy8);
    TPC_EP_full = mQ.Phi() / nHar;
    TPC_EP_east = mQ1.Phi() / nHar;
    TPC_EP_west = mQ2.Phi() / nHar;
    TPC_EP_for = mQ3.Phi() / nHar;
    TPC_EP_bac = mQ4.Phi() / nHar;
    TPC_EP_for_pos = mQ5.Phi() / nHar;
    TPC_EP_for_neg = mQ6.Phi() / nHar;
    TPC_EP_bac_pos = mQ5.Phi() / nHar;
    TPC_EP_bac_neg = mQ6.Phi() / nHar;
    Q2_proper_EPD = double(mQQx_Q * mQQx_Q + mQQy_Q * mQQy_Q) / 100.0 / (double(Qcount_Full) / 100.0 + pow((double(Qcount_Full) / 10.0 * v2_averaged_TPC[Centrality]), 2));
    // Q2_proper_EPD1 = double(mQQx_Q * mQQx_Q + mQQy_Q * mQQy_Q) / 100.0 / (double(Qcount_Full) / 100.0 + pow((double(Qcount_Full) / 10.0 * v2_averaged_EPD1[Centrality]), 2));
    // Q2_proper = double(mQQx * mQQx + mQQy * mQQy) / 100.0 / (double(Qcount) / 100.0 + pow((double(Qcount) / 10.0 * v2_averaged_TPC[Centrality]), 2));
    Q2_proper_EPD1 = Q2_proper_EPD;
    Q2_proper = Q2_proper_EPD;

    // Q2_pion_EPD = double(mQQx_Q_pion * mQQx_Q_pion + mQQy_Q_pion * mQQy_Q_pion) / 100.0 / (double(Qcount_pion) / 100.0 + pow((double(Qcount_pion) / 10.0 * v2_averaged_pairpion_EPD[Centrality]), 2));
    // Q2_pion_EPD1 = double(mQQx_Q_pion * mQQx_Q_pion + mQQy_Q_pion * mQQy_Q_pion) / 100.0 / (double(Qcount_pion) / 100.0 + pow((double(Qcount_pion) / 10.0 * v2_averaged_pairpion_EPD1[Centrality]), 2));
    // Q2_pion_TPC = double(mQQx_Q_pion * mQQx_Q_pion + mQQy_Q_pion * mQQy_Q_pion) / 100.0 / (double(Qcount_pion) / 100.0 + pow((double(Qcount_pion) / 10.0 * v2_averaged_pairpion_TPC[Centrality]), 2));
    // Q2_pion_EPD1 = Q2_pion_TPC;
    // Q2_pion_EPD = Q2_pion_TPC;

    // cout << "MakeTPC Q2_pion_TPC = " << Q2_pion_TPC << endl;

    TPC_Day3_cos2->Fill(Day3, cos(2 * TPC_EP_full));
    TPC_Day3_sin2->Fill(Day3, sin(2 * TPC_EP_full));

    // Hist_Q2_EPD->Fill(Q2_proper_EPD);
    // Hist_Q2_EPD1->Fill(Q2_proper_EPD1);
    // Hist_Q2_TPC->Fill(Q2_proper);

    float tmp_Q2[3] = {Q2_proper, Q2_proper_EPD, Q2_proper_EPD1};
    float tmp_Q2_pion[3] = {Q2_pion_TPC, Q2_pion_TPC, Q2_pion_TPC};

    for (int part1 = 0; part1 < n_lam; part1++)
    {
        particle lam_tmp(lambda_px->at(part1), lambda_py->at(part1), lambda_pz->at(part1), lambda_Charge->at(part1), lambda_dcaglobal->at(part1), lambda_mass->at(part1), lambda_trk_id->at(part1), lambda_nsigma->at(part1), lambda_hits_ratio->at(part1), lambda_nhitsfit->at(part1), lambda_nhitsmax->at(part1));

        for (int lm = 0; lm < 3; lm++)
        {
            if (lam_tmp.Charge == 1)
            {
                V0Mass_Q2[lm]->Fill(tmp_Q2[lm], lam_tmp.mass);
                V0Mass_Q2_pion[lm]->Fill(tmp_Q2_pion[lm], lam_tmp.mass);
            }
            else if (lam_tmp.Charge == -1)
            {
                V0Mass_Q2_anti[lm]->Fill(tmp_Q2[lm], lam_tmp.mass);
                V0Mass_Q2_pion_anti[lm]->Fill(tmp_Q2_pion[lm], lam_tmp.mass);
            }
        }
    }

    proton_overlap_ratio->Fill((double)proton_use_count / (double)NPTracks);
}

void FillPhiAsso() // shift parameters for Particles of EP
{
    int index = 0;
    for (int kk = 0; kk < order; kk++)
    {
        if (pVz > 0)
            index = (ChargeAsso > 0) ? 1 + 8 * kk : 3 + 8 * kk;
        else
            index = (ChargeAsso > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if (EtaAsso > 0)
        {
            if (PtAsso < 0.5)
            {
                pTPCmeanPhiAsso_FF_1->Fill(index, Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_FF_1->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else if (PtAsso < 1)
            {
                pTPCmeanPhiAsso_FF_2->Fill(index, Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_FF_2->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else
            {
                pTPCmeanPhiAsso_FF_3->Fill(index, Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_FF_3->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
        }
        if (EtaAsso < 0)
        {
            if (PtAsso < 0.5)
            {
                pTPCmeanPhiAsso_RF_1->Fill(index, Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_RF_1->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else if (PtAsso < 1)
            {
                pTPCmeanPhiAsso_RF_2->Fill(index, Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_RF_2->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
            else
            {
                pTPCmeanPhiAsso_RF_3->Fill(index, Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
                pTPCmeanPhiAsso_RF_3->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
            }
        }
    }
}

void ShiftPhiAsso(int tr)
{
    int index = 0;
    if (Weight_Read && (TPCmeanAsso_FF_1->GetEntries() || TPCmeanAsso_RF_1->GetEntries()))
    {
        for (int kk = 0; kk < order; kk++)
        {
            if (pVz > 0)
                index = (ChargeAsso > 0) ? 1 + 8 * kk : 3 + 8 * kk;
            else
                index = (ChargeAsso > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if (EtaAsso > 0)
            {
                if (PtAsso < 0.5)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_FF_1->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_FF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (PtAsso < 1)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_FF_2->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_FF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_FF_3->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_FF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if (PtAsso < 0.5)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_RF_1->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_RF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (PtAsso < 1)
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_RF_2->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_RF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMeanAsso_cos[kk] = TPCmeanAsso_RF_3->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMeanAsso_sin[kk] = TPCmeanAsso_RF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    PhiAsso_new[tr] = PhiAsso;
    for (int jj = 0; jj < order; jj++)
        PhiAsso_new[tr] += -2 * PhiMeanAsso_sin[jj] * cos(jj * PhiAsso + PhiAsso) / (jj + 1) + 2 * PhiMeanAsso_cos[jj] * sin(jj * PhiAsso + PhiAsso) / (jj + 1);
}

void ShiftPsi()
{
    Hist_TPC_EP_full->Fill(TPC_EP_full, Day);
    Hist_TPC_EP_east->Fill(TPC_EP_east, Day);
    Hist_TPC_EP_west->Fill(TPC_EP_west, Day);
    Hist_TPC_EP_for->Fill(TPC_EP_for, Day);
    Hist_TPC_EP_bac->Fill(TPC_EP_bac, Day);
    Hist_TPC_EP_for_pos->Fill(TPC_EP_for_pos, Day);
    Hist_TPC_EP_for_neg->Fill(TPC_EP_for_neg, Day);
    Hist_TPC_EP_bac_pos->Fill(TPC_EP_bac_pos, Day);
    Hist_TPC_EP_bac_neg->Fill(TPC_EP_bac_neg, Day);
    Hist_EPD_EP1_east->Fill(EPD_EP1_east, Day);
    Hist_EPD_EP1_west->Fill(EPD_EP1_west, Day);
    Hist_EPD_EP_east->Fill(EPD_EP_east, Day);
    Hist_EPD_EP_west->Fill(EPD_EP_west, Day);

    for (int kk = 0; kk < order; kk++)
    {
        pTPC_EP_full->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_full), Eweight);
        pTPC_EP_full->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_full), Eweight);
        pTPC_EP_east->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_east), Eweight);
        pTPC_EP_east->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_east), Eweight);
        pTPC_EP_west->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_west), Eweight);
        pTPC_EP_west->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_west), Eweight);
        pTPC_EP_for->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_for), Eweight);
        pTPC_EP_for->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_for), Eweight);
        pTPC_EP_bac->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_bac), Eweight);
        pTPC_EP_bac->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_bac), Eweight);
        pTPC_EP_for_pos->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_for_pos), Eweight);
        pTPC_EP_for_pos->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_for_pos), Eweight);
        pTPC_EP_for_neg->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_for_neg), Eweight);
        pTPC_EP_for_neg->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_for_neg), Eweight);
        pTPC_EP_bac_pos->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_bac_pos), Eweight);
        pTPC_EP_bac_pos->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_bac_pos), Eweight);
        pTPC_EP_bac_neg->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * TPC_EP_bac_neg), Eweight);
        pTPC_EP_bac_neg->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * TPC_EP_bac_neg), Eweight);

        pEPD_EP1_east->Fill(1 + 2 * kk, Day2, cos((kk + 1) * EPD_EP1_east), Eweight);
        pEPD_EP1_east->Fill(2 + 2 * kk, Day2, sin((kk + 1) * EPD_EP1_east), Eweight);
        pEPD_EP1_west->Fill(1 + 2 * kk, Day2, cos((kk + 1) * EPD_EP1_west), Eweight);
        pEPD_EP1_west->Fill(2 + 2 * kk, Day2, sin((kk + 1) * EPD_EP1_west), Eweight);
        pEPD_EP_east->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * EPD_EP_east), Eweight);
        pEPD_EP_east->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * EPD_EP_east), Eweight);
        pEPD_EP_west->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1) * EPD_EP_west), Eweight);
        pEPD_EP_west->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1) * EPD_EP_west), Eweight);
    }

    if (Weight_Read && Read_TPC_EP_full->GetEntries())
    {
        for (int k = 0; k < 2 * order; k++)
        {
            PsiMean_F[k] = Read_TPC_EP_full->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_E[k] = Read_TPC_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_W[k] = Read_TPC_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_f[k] = Read_TPC_EP_for->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_b[k] = Read_TPC_EP_bac->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_fp[k] = Read_TPC_EP_for_pos->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_fn[k] = Read_TPC_EP_for_neg->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_bp[k] = Read_TPC_EP_bac_pos->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_bn[k] = Read_TPC_EP_bac_neg->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);

            Psi1_EPD_E[k] = Read_EPD_EP1_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi1_EPD_W[k] = Read_EPD_EP1_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_EPD_E[k] = Read_EPD_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_EPD_W[k] = Read_EPD_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
        }
    }

    TPC_EP_full_new = TPC_EP_full, TPC_EP_east_new = TPC_EP_east, TPC_EP_west_new = TPC_EP_west;
    TPC_EP_for_new = TPC_EP_for, TPC_EP_bac_new = TPC_EP_bac;
    TPC_EP_for_pos_new = TPC_EP_for_pos, TPC_EP_for_neg_new = TPC_EP_for_neg;
    TPC_EP_bac_pos_new = TPC_EP_bac_pos, TPC_EP_bac_neg_new = TPC_EP_bac_neg;

    EPD_EP_east_new = EPD_EP_east, EPD_EP_west_new = EPD_EP_west;
    EPD_EP1_east_new = EPD_EP1_east, EPD_EP1_west_new = EPD_EP1_west;

    for (int jj = 0; jj < order; jj++)
    {
        TPC_EP_full_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_full) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_full) / nHar / (jj + 1);
        TPC_EP_east_new += -2 * PsiMean_E[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_east) / nHar / (jj + 1) + 2 * PsiMean_E[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_east) / nHar / (jj + 1);
        TPC_EP_west_new += -2 * PsiMean_W[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_west) / nHar / (jj + 1) + 2 * PsiMean_W[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_west) / nHar / (jj + 1);
        TPC_EP_for_new += -2 * PsiMean_f[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_for) / nHar / (jj + 1) + 2 * PsiMean_f[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_for) / nHar / (jj + 1);
        TPC_EP_bac_new += -2 * PsiMean_b[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_bac) / nHar / (jj + 1) + 2 * PsiMean_b[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_bac) / nHar / (jj + 1);
        TPC_EP_for_pos_new += -2 * PsiMean_fp[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_for_pos) / nHar / (jj + 1) + 2 * PsiMean_fp[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_for_pos) / nHar / (jj + 1);
        TPC_EP_for_neg_new += -2 * PsiMean_fn[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_for_neg) / nHar / (jj + 1) + 2 * PsiMean_fn[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_for_neg) / nHar / (jj + 1);
        TPC_EP_bac_pos_new += -2 * PsiMean_bp[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_bac_pos) / nHar / (jj + 1) + 2 * PsiMean_bp[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_bac_pos) / nHar / (jj + 1);
        TPC_EP_bac_neg_new += -2 * PsiMean_bn[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_bac_neg) / nHar / (jj + 1) + 2 * PsiMean_bn[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_bac_neg) / nHar / (jj + 1);

        EPD_EP1_east_new += -2 * Psi1_EPD_E[1 + 2 * jj] * cos((jj + 1) * EPD_EP1_east) / (jj + 1) + 2 * Psi1_EPD_E[0 + 2 * jj] * sin((jj + 1) * EPD_EP1_east) / (jj + 1);
        EPD_EP1_west_new += -2 * Psi1_EPD_W[1 + 2 * jj] * cos((jj + 1) * EPD_EP1_west) / (jj + 1) + 2 * Psi1_EPD_W[0 + 2 * jj] * sin((jj + 1) * EPD_EP1_west) / (jj + 1);
        EPD_EP_east_new += -2 * Psi_EPD_E[1 + 2 * jj] * cos(nHar * (jj + 1) * EPD_EP_east) / nHar / (jj + 1) + 2 * Psi_EPD_E[0 + 2 * jj] * sin(nHar * (jj + 1) * EPD_EP_east) / nHar / (jj + 1);
        EPD_EP_west_new += -2 * Psi_EPD_W[1 + 2 * jj] * cos(nHar * (jj + 1) * EPD_EP_west) / nHar / (jj + 1) + 2 * Psi_EPD_W[0 + 2 * jj] * sin(nHar * (jj + 1) * EPD_EP_west) / nHar / (jj + 1);
    }

    if (TPC_EP_full_new > PI)
        TPC_EP_full_new -= PI;
    if (TPC_EP_full_new < 0)
        TPC_EP_full_new += PI;
    if (TPC_EP_east_new > PI)
        TPC_EP_east_new -= PI;
    if (TPC_EP_east_new < 0)
        TPC_EP_east_new += PI;
    if (TPC_EP_west_new > PI)
        TPC_EP_west_new -= PI;
    if (TPC_EP_west_new < 0)
        TPC_EP_west_new += PI;
    if (TPC_EP_for_new > PI)
        TPC_EP_for_new -= PI;
    if (TPC_EP_for_new < 0)
        TPC_EP_for_new += PI;
    if (TPC_EP_bac_new > PI)
        TPC_EP_bac_new -= PI;
    if (TPC_EP_bac_new < 0)
        TPC_EP_bac_new += PI;
    if (TPC_EP_for_pos_new > PI)
        TPC_EP_for_pos_new -= PI;
    if (TPC_EP_for_pos_new < 0)
        TPC_EP_for_pos_new += PI;
    if (TPC_EP_for_neg_new > PI)
        TPC_EP_for_neg_new -= PI;
    if (TPC_EP_for_neg_new < 0)
        TPC_EP_for_neg_new += PI;
    if (TPC_EP_bac_pos_new > PI)
        TPC_EP_bac_pos_new -= PI;
    if (TPC_EP_bac_pos_new < 0)
        TPC_EP_bac_pos_new += PI;
    if (TPC_EP_bac_neg_new > PI)
        TPC_EP_bac_neg_new -= PI;
    if (TPC_EP_bac_neg_new < 0)
        TPC_EP_bac_neg_new += PI;

    if (EPD_EP1_east_new > PI)
        EPD_EP1_east_new -= 2 * PI;
    if (EPD_EP1_east_new < -PI)
        EPD_EP1_east_new += 2 * PI;
    if (EPD_EP1_west_new > PI)
        EPD_EP1_west_new -= 2 * PI;
    if (EPD_EP1_west_new < -PI)
        EPD_EP1_west_new += 2 * PI;
    if (EPD_EP_east_new > 2. * PI / nHar)
        EPD_EP_east_new -= 2. * PI / nHar;
    if (EPD_EP_east_new < 0)
        EPD_EP_east_new += 2. * PI / nHar;
    if (EPD_EP_west_new > 2. * PI / nHar)
        EPD_EP_west_new -= 2. * PI / nHar;
    if (EPD_EP_west_new < 0)
        EPD_EP_west_new += 2. * PI / nHar;

    Hist_TPC_EP_east_flat->Fill(TPC_EP_east_new, Day);
    Hist_TPC_EP_west_flat->Fill(TPC_EP_west_new, Day);
    Hist_TPC_EP_for_flat->Fill(TPC_EP_for_new, Day);
    Hist_TPC_EP_bac_flat->Fill(TPC_EP_bac_new, Day);
    Hist_TPC_EP_for_pos_flat->Fill(TPC_EP_for_pos_new, Day);
    Hist_TPC_EP_for_neg_flat->Fill(TPC_EP_for_neg_new, Day);
    Hist_TPC_EP_bac_pos_flat->Fill(TPC_EP_bac_pos_new, Day);
    Hist_TPC_EP_bac_neg_flat->Fill(TPC_EP_bac_neg_new, Day);
    Hist_TPC_EP_full_flat->Fill(TPC_EP_full_new, Day);

    Hist_EPD_EP1_east_flat->Fill(EPD_EP1_east_new, Day);
    Hist_EPD_EP1_west_flat->Fill(EPD_EP1_west_new, Day);
    Hist_EPD_EP_east_flat->Fill(EPD_EP_east_new, Day);
    Hist_EPD_EP_west_flat->Fill(EPD_EP_west_new, Day);

    float EPD_EP_full = atan2(sin(nHar * EPD_EP_east_new) + sin(nHar * EPD_EP_west_new), cos(nHar * EPD_EP_east_new) + cos(nHar * EPD_EP_west_new)) / nHar;
    if (EPD_EP_full < 0)
        EPD_EP_full += 2. * PI / nHar;

    Hist_EPD_vs_TPC->Fill(TPC_EP_full_new, EPD_EP_full);
}

void FillEP_resolution()
{
    float cos_ew = cos(nHar * TPC_EP_east_new - nHar * TPC_EP_west_new);
    float cos_fb = cos(nHar * TPC_EP_for_new - nHar * TPC_EP_bac_new);
    Hist_cos->Fill(1, cos_fb, Eweight);
    Hist_cos->Fill(2, cos_ew, Eweight);
    Hist_cos->Fill(3, cos(nHar * TPC_EP_for_new - nHar * TPC_EP_full_new));
    Hist_cos->Fill(4, cos(nHar * TPC_EP_bac_new - nHar * TPC_EP_full_new));

    Hist_cos_EPD->Fill(1, cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new), Eweight);
    Hist_cos_EPD->Fill(2, cos(nHar * EPD_EP_east_new - nHar * TPC_EP_full_new), Eweight);
    Hist_cos_EPD->Fill(3, cos(nHar * EPD_EP_west_new - nHar * TPC_EP_full_new), Eweight);
    Hist_cos_EPD->Fill(4, cos(EPD_EP1_east_new - EPD_EP1_west_new), Eweight);
    Hist_cos_EPD->Fill(5, cos(nHar * TPC_EP_full_new - EPD_EP1_east_new - EPD_EP1_west_new), Eweight);
    Hist_cos_EPD->Fill(6, cos(nHar * EPD_EP_east_new - EPD_EP1_east_new - EPD_EP1_west_new), Eweight);
    Hist_cos_EPD->Fill(7, cos(nHar * EPD_EP_west_new - EPD_EP1_east_new - EPD_EP1_west_new), Eweight);

    EPD_cor_Day3->Fill(Day3, cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new));
    EPD1_cor_Day3->Fill(Day3, cos(EPD_EP1_east_new - EPD_EP1_west_new));

    float mQQx_Q_pion = 0., mQQy_Q_pion = 0.;
    int Qcount_pion = 0;
    
    for (int idx_pion = 0; idx_pion < p_pion_list->size(); idx_pion++)
    {
        TVector3 p_pion_list_tmp1 = p_pion_list->at(idx_pion);

        for (int idx_pion2 = idx_pion + 1; idx_pion2 < p_pion_list->size(); idx_pion2++)
        {
            TVector3 p_pion_list_tmp2 = p_pion_list->at(idx_pion2);

            TVector3 pair_pion = p_pion_list_tmp1 + p_pion_list_tmp2;

            mQQx_Q_pion += cos(pair_pion.Phi() * nHar);
            mQQy_Q_pion += sin(pair_pion.Phi() * nHar);
            Qcount_pion++;

            float v2_final_pion = 0;
            v2_final_pion = cos(nHar * (p_pion_list_tmp1.Phi() - p_pion_list_tmp2.Phi())) * 100.0;

            v2_2_pion->Fill(0.5, v2_final_pion, Eweight);

            float v2pe_pion_temp = cos(nHar * pair_pion.Phi() - nHar * TPC_EP_east_new) * 100;
            float v2pw_pion_temp = cos(nHar * pair_pion.Phi() - nHar * TPC_EP_west_new) * 100;
            float v2pe_pion_EPD_temp = cos(nHar * pair_pion.Phi() - nHar * EPD_EP_east_new) * 100;
            float v2pw_pion_EPD_temp = cos(nHar * pair_pion.Phi() - nHar * EPD_EP_west_new) * 100;
            float v2p_pion_EPD1_temp = cos(nHar * pair_pion.Phi() - EPD_EP1_east_new - EPD_EP1_west_new) * 100;

            pTemp_v2_chargedhadrons->Fill(14, v2pe_pion_temp);
            pTemp_v2_chargedhadrons->Fill(15, v2pw_pion_temp);
            pTemp_v2_chargedhadrons->Fill(16, v2p_pion_EPD1_temp);
            pTemp_v2_chargedhadrons->Fill(17, v2pe_pion_EPD_temp);
            pTemp_v2_chargedhadrons->Fill(18, v2pw_pion_EPD_temp);

        }
    }

    Q2_pion_TPC = double(mQQx_Q_pion * mQQx_Q_pion + mQQy_Q_pion * mQQy_Q_pion) / 100.0 / (double(Qcount_pion) / 100.0 + pow((double(Qcount_pion) / 10.0 * v2_averaged_pairpion_TPC[Centrality]), 2));
    Q2_pion_EPD1 = Q2_pion_TPC;
    Q2_pion_EPD = Q2_pion_TPC;

    // Hist_Q2->Fill(Q2_proper_EPD);
    // Hist_Q2_proper->Fill(Q2_proper);
    Hist_Q2_EPD->Fill(Q2_proper_EPD);
    Hist_Q2_EPD1->Fill(Q2_proper_EPD1);
    Hist_Q2_TPC->Fill(Q2_proper);
    p_cos_Q2->Fill(Q2_proper, cos_ew, Eweight);
    p_RefMult_Q2->Fill(Q2_proper, Fcount);
    Hist_RefMult_Q2->Fill(Q2_proper, Fcount);

    // cout << "Filling Q2_pion_TPC = " << Q2_pion_TPC << endl;

    Hist_Q2_EPD_pion->Fill(Q2_pion_EPD);
    Hist_Q2_EPD1_pion->Fill(Q2_pion_EPD1);
    Hist_Q2_TPC_pion->Fill(Q2_pion_TPC);

    p_cos_Q2_EPD->Fill(Q2_proper_EPD, cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new), Eweight);
    p_cos_Q2_EPD1->Fill(Q2_proper_EPD1, cos(EPD_EP1_east_new - EPD_EP1_west_new), Eweight);
    p_RefMult_Q2_EPD->Fill(Q2_proper_EPD, Fcount);
    Hist_RefMult_Q2_EPD->Fill(Q2_proper_EPD, Fcount);

    p_cos_Q2_TPC_pion->Fill(Q2_pion_TPC, cos_ew, Eweight);
    p_cos_Q2_EPD_pion->Fill(Q2_pion_EPD, cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new), Eweight);
    p_cos_Q2_EPD1_pion->Fill(Q2_pion_EPD1, cos(EPD_EP1_east_new - EPD_EP1_west_new), Eweight);

    Fcount_parent = Fcount;

    int net_Nch_Asym_index = -99;
    if (net_Nch_Asym > -99)
    {
        p_netChAsym_cos->Fill(net_Nch_Asym, cos_fb);
        if (net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym))
            net_Nch_Asym_index = 1;
        else if (net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym)))
            net_Nch_Asym_index = 2;
        else if (net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym)))
            net_Nch_Asym_index = 3;
        else if (net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym))
            net_Nch_Asym_index = 4;
        else
            net_Nch_Asym_index = 5;

        Hist_cos_Ach->Fill(net_Nch_Asym_index, cos_fb, Eweight);

        Hist_cos->Fill((5 + 2 * (net_Nch_Asym_index - 1)), cos(nHar * TPC_EP_for_pos_new - nHar * TPC_EP_bac_pos_new));
        Hist_cos->Fill((6 + 2 * (net_Nch_Asym_index - 1)), cos(nHar * TPC_EP_for_neg_new - nHar * TPC_EP_bac_neg_new));
    }
}

void FillPhiPOI() // shift parameters for Particles of Interest
{
    int index = 0, index2 = 0;
    if (pVz > 0)
        index2 = (Charge > 0) ? 1 : 2;
    else
        index2 = (Charge > 0) ? 3 : 4;
    for (int kk = 0; kk < order; kk++)
    {
        if (pVz > 0)
            index = (Charge > 0) ? 1 + 8 * kk : 3 + 8 * kk;
        else
            index = (Charge > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if (Eta > 0)
        {
            if (Pt < 0.5)
            {
                pTPCmeanPhi_FF_1->Fill(index, Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_FF_1->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_1->Fill(Phi, index2, Eweight);
            }
            else if (Pt < 1)
            {
                pTPCmeanPhi_FF_2->Fill(index, Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_FF_2->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_2->Fill(Phi, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_FF_3->Fill(index, Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_FF_3->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_3->Fill(Phi, index2, Eweight);
            }
        }
        if (Eta < 0)
        {
            if (Pt < 0.5)
            {
                pTPCmeanPhi_RF_1->Fill(index, Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_RF_1->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_1->Fill(Phi, index2, Eweight);
            }
            else if (Pt < 1)
            {
                pTPCmeanPhi_RF_2->Fill(index, Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_RF_2->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_2->Fill(Phi, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_RF_3->Fill(index, Day2, cos(kk * Phi + Phi), Eweight);
                pTPCmeanPhi_RF_3->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_3->Fill(Phi, index2, Eweight);
            }
        }
    }
}

void ShiftPhiPOI(int tr)
{
    int index = 0, index2 = 0;
    if (pVz > 0)
        index2 = (Charge > 0) ? 1 : 2;
    else
        index2 = (Charge > 0) ? 3 : 4;
    if (Weight_Read && (TPCmean_FF_1->GetEntries() || TPCmean_RF_1->GetEntries()))
    {
        for (int kk = 0; kk < order; kk++)
        {
            if (pVz > 0)
                index = (Charge > 0) ? 1 + 8 * kk : 3 + 8 * kk;
            else
                index = (Charge > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if (Eta > 0)
            {
                if (Pt < 0.5)
                {
                    PhiMean_cos[kk] = TPCmean_FF_1->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_FF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (Pt < 1)
                {
                    PhiMean_cos[kk] = TPCmean_FF_2->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_FF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos[kk] = TPCmean_FF_3->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_FF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if (Pt < 0.5)
                {
                    PhiMean_cos[kk] = TPCmean_RF_1->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_RF_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (Pt < 1)
                {
                    PhiMean_cos[kk] = TPCmean_RF_2->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_RF_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos[kk] = TPCmean_RF_3->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin[kk] = TPCmean_RF_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    Phi_new[tr] = Phi; // store the shifted angles
    for (int jj = 0; jj < order; jj++)
        Phi_new[tr] += -2 * PhiMean_sin[jj] * cos(jj * Phi + Phi) / (jj + 1) + 2 * PhiMean_cos[jj] * sin(jj * Phi + Phi) / (jj + 1);
    if (Phi_new[tr] > PI)
        Phi_new[tr] -= 2 * PI;
    if (Phi_new[tr] < -PI)
        Phi_new[tr] += 2 * PI;
    if (Eta > 0)
    {
        if (Pt < 0.5)
            Hist_Phi_FF_new_1->Fill(Phi_new[tr], index2, Eweight);
        else if (Pt < 1)
            Hist_Phi_FF_new_2->Fill(Phi_new[tr], index2, Eweight);
        else
            Hist_Phi_FF_new_3->Fill(Phi_new[tr], index2, Eweight);
    }
    else
    {
        if (Pt < 0.5)
            Hist_Phi_RF_new_1->Fill(Phi_new[tr], index2, Eweight);
        else if (Pt < 1)
            Hist_Phi_RF_new_2->Fill(Phi_new[tr], index2, Eweight);
        else
            Hist_Phi_RF_new_3->Fill(Phi_new[tr], index2, Eweight);
    }
}

void FillPhiPOI_rot() // shift parameters for Particles of Interest
{
    int index = 0, index2 = 0;
    if (pVz > 0)
        index2 = (Charge_rot > 0) ? 1 : 2;
    else
        index2 = (Charge_rot > 0) ? 3 : 4;
    for (int kk = 0; kk < order; kk++)
    {
        if (pVz > 0)
            index = (Charge_rot > 0) ? 1 + 8 * kk : 3 + 8 * kk;
        else
            index = (Charge_rot > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if (Eta_rot > 0)
        {
            if (Pt_rot < 0.5)
            {
                pTPCmeanPhi_FF_1_rot->Fill(index, Day2, cos(kk * Phi_rot + Phi_rot), Eweight);
                pTPCmeanPhi_FF_1_rot->Fill(index + 1, Day2, sin(kk * Phi_rot + Phi_rot), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_1_rot->Fill(Phi_rot, index2, Eweight);
            }
            else if (Pt_rot < 1)
            {
                pTPCmeanPhi_FF_2_rot->Fill(index, Day2, cos(kk * Phi_rot + Phi_rot), Eweight);
                pTPCmeanPhi_FF_2_rot->Fill(index + 1, Day2, sin(kk * Phi_rot + Phi_rot), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_2_rot->Fill(Phi_rot, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_FF_3_rot->Fill(index, Day2, cos(kk * Phi_rot + Phi_rot), Eweight);
                pTPCmeanPhi_FF_3_rot->Fill(index + 1, Day2, sin(kk * Phi_rot + Phi_rot), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_3_rot->Fill(Phi_rot, index2, Eweight);
            }
        }
        if (Eta_rot < 0)
        {
            if (Pt_rot < 0.5)
            {
                pTPCmeanPhi_RF_1_rot->Fill(index, Day2, cos(kk * Phi_rot + Phi_rot), Eweight);
                pTPCmeanPhi_RF_1_rot->Fill(index + 1, Day2, sin(kk * Phi_rot + Phi_rot), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_1_rot->Fill(Phi_rot, index2, Eweight);
            }
            else if (Pt_rot < 1)
            {
                pTPCmeanPhi_RF_2_rot->Fill(index, Day2, cos(kk * Phi_rot + Phi_rot), Eweight);
                pTPCmeanPhi_RF_2_rot->Fill(index + 1, Day2, sin(kk * Phi_rot + Phi_rot), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_2_rot->Fill(Phi_rot, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_RF_3_rot->Fill(index, Day2, cos(kk * Phi_rot + Phi_rot), Eweight);
                pTPCmeanPhi_RF_3_rot->Fill(index + 1, Day2, sin(kk * Phi_rot + Phi_rot), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_3_rot->Fill(Phi_rot, index2, Eweight);
            }
        }
    }
}

void ShiftPhiPOI_rot(int tr)
{
    int index = 0, index2 = 0;
    if (pVz > 0)
        index2 = (Charge_rot > 0) ? 1 : 2;
    else
        index2 = (Charge_rot > 0) ? 3 : 4;
    if (Weight_Read && (TPCmean_FF_1_rot->GetEntries() || TPCmean_RF_1_rot->GetEntries()))
    {
        for (int kk = 0; kk < order; kk++)
        {
            if (pVz > 0)
                index = (Charge_rot > 0) ? 1 + 8 * kk : 3 + 8 * kk;
            else
                index = (Charge_rot > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if (Eta_rot > 0)
            {
                if (Pt_rot < 0.5)
                {
                    PhiMean_cos_rot[kk] = TPCmean_FF_1_rot->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_rot[kk] = TPCmean_FF_1_rot->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (Pt_rot < 1)
                {
                    PhiMean_cos_rot[kk] = TPCmean_FF_2_rot->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_rot[kk] = TPCmean_FF_2_rot->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos_rot[kk] = TPCmean_FF_3_rot->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_rot[kk] = TPCmean_FF_3_rot->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if (Pt_rot < 0.5)
                {
                    PhiMean_cos_rot[kk] = TPCmean_RF_1_rot->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_rot[kk] = TPCmean_RF_1_rot->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (Pt_rot < 1)
                {
                    PhiMean_cos_rot[kk] = TPCmean_RF_2_rot->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_rot[kk] = TPCmean_RF_2_rot->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos_rot[kk] = TPCmean_RF_3_rot->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_rot[kk] = TPCmean_RF_3_rot->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    Phi_new_rot[tr] = Phi_rot; // store the shifted angles
    for (int jj = 0; jj < order; jj++)
        Phi_new_rot[tr] += -2 * PhiMean_sin_rot[jj] * cos(jj * Phi_rot + Phi_rot) / (jj + 1) + 2 * PhiMean_cos_p[jj] * sin(jj * Phi_rot + Phi_rot) / (jj + 1);
    if (Phi_new_rot[tr] > PI)
        Phi_new_rot[tr] -= 2 * PI;
    if (Phi_new_rot[tr] < -PI)
        Phi_new_rot[tr] += 2 * PI;
    if (Eta_rot > 0)
    {
        if (Pt_rot < 0.5)
            Hist_Phi_FF_new_1_rot->Fill(Phi_new_rot[tr], index2, Eweight);
        else if (Pt_rot < 1)
            Hist_Phi_FF_new_2_rot->Fill(Phi_new_rot[tr], index2, Eweight);
        else
            Hist_Phi_FF_new_3_rot->Fill(Phi_new_rot[tr], index2, Eweight);
    }
    else
    {
        if (Pt_rot < 0.5)
            Hist_Phi_RF_new_1_rot->Fill(Phi_new_rot[tr], index2, Eweight);
        else if (Pt_rot < 1)
            Hist_Phi_RF_new_2_rot->Fill(Phi_new_rot[tr], index2, Eweight);
        else
            Hist_Phi_RF_new_3_rot->Fill(Phi_new_rot[tr], index2, Eweight);
    }
}

void ShiftPhiPOI_p(int tr)
{
    int index = 0, index2 = 0;
    if (pVz > 0)
        index2 = (Charge2 > 0) ? 1 : 2;
    else
        index2 = (Charge2 > 0) ? 3 : 4;
    if (Weight_Read && (TPCmean_FF_1_p->GetEntries() || TPCmean_RF_1_p->GetEntries()))
    {
        for (int kk = 0; kk < order; kk++)
        {
            if (pVz > 0)
                index = (Charge2 > 0) ? 1 + 8 * kk : 3 + 8 * kk;
            else
                index = (Charge2 > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
            if (Eta2 > 0)
            {
                if (Pt2 < 0.5)
                {
                    PhiMean_cos_p[kk] = TPCmean_FF_1_p->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_FF_1_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (Pt2 < 1)
                {
                    PhiMean_cos_p[kk] = TPCmean_FF_2_p->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_FF_2_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos_p[kk] = TPCmean_FF_3_p->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_FF_3_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
            else
            {
                if (Pt2 < 0.5)
                {
                    PhiMean_cos_p[kk] = TPCmean_RF_1_p->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_RF_1_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else if (Pt2 < 1)
                {
                    PhiMean_cos_p[kk] = TPCmean_RF_2_p->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_RF_2_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
                else
                {
                    PhiMean_cos_p[kk] = TPCmean_RF_3_p->GetBinContent(index, Day2 - run_sta / 10 + 1);
                    PhiMean_sin_p[kk] = TPCmean_RF_3_p->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
                }
            }
        }
    }

    Phi_pnew[tr] = Phi2; // store the shifted angles
    for (int jj = 0; jj < order; jj++)
        Phi_pnew[tr] += -2 * PhiMean_sin_p[jj] * cos(jj * Phi2 + Phi2) / (jj + 1) + 2 * PhiMean_cos_p[jj] * sin(jj * Phi2 + Phi2) / (jj + 1);
    if (Phi_pnew[tr] > PI)
        Phi_pnew[tr] -= 2 * PI;
    if (Phi_pnew[tr] < -PI)
        Phi_pnew[tr] += 2 * PI;
    if (Eta2 > 0)
    {
        if (Pt2 < 0.5)
            Hist_Phi_FF_new_1_p->Fill(Phi_pnew[tr], index2, Eweight);
        else if (Pt2 < 1)
            Hist_Phi_FF_new_2_p->Fill(Phi_pnew[tr], index2, Eweight);
        else
            Hist_Phi_FF_new_3_p->Fill(Phi_pnew[tr], index2, Eweight);
    }
    else
    {
        if (Pt2 < 0.5)
            Hist_Phi_RF_new_1_p->Fill(Phi_pnew[tr], index2, Eweight);
        else if (Pt2 < 1)
            Hist_Phi_RF_new_2_p->Fill(Phi_pnew[tr], index2, Eweight);
        else
            Hist_Phi_RF_new_3_p->Fill(Phi_pnew[tr], index2, Eweight);
    }
}

void FillPhiPOI_p() // shift parameters for Particles of Interest (Protons)
{
    int index = 0, index2 = 0;
    if (pVz > 0)
        index2 = (Charge2 > 0) ? 1 : 2;
    else
        index2 = (Charge2 > 0) ? 3 : 4;
    for (int kk = 0; kk < order; kk++)
    {
        if (pVz > 0)
            index = (Charge2 > 0) ? 1 + 8 * kk : 3 + 8 * kk;
        else
            index = (Charge2 > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if (Eta2 > 0)
        {
            if (Pt2 < 0.5)
            {
                pTPCmeanPhi_FF_1_p->Fill(index, Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_FF_1_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_1_p->Fill(Phi2, index2, Eweight);
            }
            else if (Pt2 < 1)
            {
                pTPCmeanPhi_FF_2_p->Fill(index, Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_FF_2_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_2_p->Fill(Phi2, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_FF_3_p->Fill(index, Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_FF_3_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if (kk == 0)
                    Hist_Phi_FF_3_p->Fill(Phi2, index2, Eweight);
            }
        }
        if (Eta2 < 0)
        {
            if (Pt2 < 0.5)
            {
                pTPCmeanPhi_RF_1_p->Fill(index, Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_RF_1_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_1_p->Fill(Phi2, index2, Eweight);
            }
            else if (Pt2 < 1)
            {
                pTPCmeanPhi_RF_2_p->Fill(index, Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_RF_2_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_2_p->Fill(Phi2, index2, Eweight);
            }
            else
            {
                pTPCmeanPhi_RF_3_p->Fill(index, Day2, cos(kk * Phi2 + Phi2), Eweight);
                pTPCmeanPhi_RF_3_p->Fill(index + 1, Day2, sin(kk * Phi2 + Phi2), Eweight);
                if (kk == 0)
                    Hist_Phi_RF_3_p->Fill(Phi2, index2, Eweight);
            }
        }
    }
}

void FillCMW()
{
    p_v2_Ach->Fill(net_Nch_Asym, v2_sub);
    if (net_Nch_Asym <= -99)
        return;
    if (net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym))
    {
        if (Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(1, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(1, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(1, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(1, Pt);
        }

        Hist_pt_Ach[0]->Fill(Pt);
    }
    else if (net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym)))
    {
        if (Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(2, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(2, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(2, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(2, Pt);
        }

        Hist_pt_Ach[1]->Fill(Pt);
    }
    else if (net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym)))
    {
        if (Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(3, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(3, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(3, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(3, Pt);
        }

        Hist_pt_Ach[2]->Fill(Pt);
    }
    else if (net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym))
    {
        if (Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(4, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(4, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(4, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(4, Pt);
        }

        Hist_pt_Ach[3]->Fill(Pt);
    }
    else
    {
        if (Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(5, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(5, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(5, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(5, Pt);
        }

        Hist_pt_Ach[4]->Fill(Pt);
    }
}

void Fillv2()
{
    pTemp_v2->Fill(1, v2e, 1. / temp_eff);
    pTemp_v2->Fill(2, v2w, 1. / temp_eff);
    if (Charge > 0)
    {
        pTemp_v2->Fill(3, v2e, 1. / temp_eff);
        pTemp_v2->Fill(4, v2w, 1. / temp_eff);
    }
    if (Charge < 0)
    {
        pTemp_v2->Fill(5, v2e, 1. / temp_eff);
        pTemp_v2->Fill(6, v2w, 1. / temp_eff);
    }
    pTemp_v2->Fill(7, v2, 1. / temp_eff);

    pTemp_v2->Fill(11, v2_EPD1, 1. / temp_eff);
    pTemp_v2->Fill(12, v2_EPDe, 1. / temp_eff);
    pTemp_v2->Fill(13, v2_EPDw, 1. / temp_eff);
}

void finish_Gamma112(int cen, int opt_weight, TString JobIDName)
{
    rc = (TH1D *)Hist_Pt_TOF->Clone();
    rc->SetName("rc");
    rc->Divide(Hist_Pt);
    WriteHistogram(cen, opt_weight, JobIDName);
    if (opt_weight == 1)
    {
        WriteWeight(fname_new);
    }

    // delete mEpFinder;
    // mEpFinder = NULL;
    // delete mEpdHits;
    // mEpdHits = NULL;

    //    debugfile.close();
}

//////////////////////////////////

/////////////////////////////////////////////////
bool IsGoodAsso(float p, float e, float d)
{
    if (p > pt_asso_up || p < pt_asso_lo)
        return false;
    if (d > DcaCut)
        return false;
    if (e > EtaCut || e < -EtaCut)
        return false;
    return true;
}
/////////////////////////////////////////////////


/////////////////////////////////////////////////////
/*bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt)
{
    // hTall->Fill(1);
    if (p->gDCA(pV).Mag() > 1)
        return false;
    // hTall->Fill(2);
    float eta_i = p->pMom().Eta();
    if (fabs(eta_i) > 0.9)
        return false;
    // hTall->Fill(3);
    float pt_i = p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if (pt_i < 0.2 || p_i > 1.6)
        return false;
    // hTall->Fill(4);
    float nSig_i = p->nSigmaPion();
    float ndEdx_i = p->nHitsDedx();
    if (ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2)
        return false;
    // hTall->Fill(5);
    if (opt == 0)
        return true;
    if (!(p->isTofTrack()))
        return false;
    // hTall->Fill(6);
    StPicoBTofPidTraits *trait = d->btofPidTraits(p->bTofPidTraitsIndex());
    if (!trait)
        return false;
    // hTall->Fill(7);
    if (trait->btof() <= 0)
        return false;
    // hTall->Fill(8);
    if (fabs(trait->btofYLocal()) > 1.8)
        return false;
    // hTall->Fill(9);
    float beta_i = trait->btofBeta();
    if (beta_i == 0)
        return false;
    // hTall->Fill(10);
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if (mass2_i < -0.01 || mass2_i > 0.1)
        return false;
    // hTall->Fill(11);
    return true;
}
///////////////////////////////////////////////////
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt)
{
    if (p->gDCA(pV).Mag() > 1)
        return false;
    float eta_i = p->pMom().Eta();
    if (fabs(eta_i) > 0.9)
        return false;
    float pt_i = p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if (pt_i < 0.2 || p_i > 1.6)
        return false;

    float nSig_i = p->nSigmaKaon();
    float ndEdx_i = p->nHitsDedx();
    if (ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2)
        return false;

    if (opt == 0)
        return true;
    if (!(p->isTofTrack()))
        return false;
    StPicoBTofPidTraits *trait = d->btofPidTraits(p->bTofPidTraitsIndex());
    if (!trait)
        return false;
    if (trait->btof() <= 0)
        return false;
    if (fabs(trait->btofYLocal()) > 1.8)
        return false;
    float beta_i = trait->btofBeta();
    if (beta_i == 0)
        return false;
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if (mass2_i < 0.2 || mass2_i > 0.35)
        return false;
    return true;
}
////////////////////////////////////////////////////
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt)
{
    if (p->gDCA(pV).Mag() > 1)
        return false;
    float eta_i = p->pMom().Eta();
    if (fabs(eta_i) > 0.9)
        return false;
    float pt_i = p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if (pt_i < 0.4 || p_i > 2)
        return false;

    float nSig_i = p->nSigmaProton();
    float ndEdx_i = p->nHitsDedx();
    if (ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2)
        return false;

    if (opt == 0)
        return true;
    if (!(p->isTofTrack()))
        return false;
    StPicoBTofPidTraits *trait = d->btofPidTraits(p->bTofPidTraitsIndex());
    if (!trait)
        return false;
    if (trait->btof() <= 0)
        return false;
    if (fabs(trait->btofYLocal()) > 1.8)
        return false;
    float beta_i = trait->btofBeta();
    if (beta_i == 0)
        return false;
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if (mass2_i < 0.8 || mass2_i > 1)
        return false;
    return true;
}*/
////////////////////////////////////

/////////////////////////////////
void WriteHistogram(int c, int o, TString JobIDName)
{
    char fname_out[200];
    TString Name2 = "sched";
    Name2.Append(JobIDName);
    if (o != 1)
        sprintf(fname_out, "cen%d.gamma112_fullEP_eff_pT02_module.root", c);
    if (o == 1)
        sprintf(fname_out, "cen%d.v2_fullEP_eff_pT02_module.root", c);
    Name2.Append(fname_out);
    TFile *fout = new TFile(Name2, "RECREATE");

    //        gROOT->GetList()->ls();
    TList *list = gROOT->GetList(); // GetListOfKeys();
    TIter next(list);
    TKey *key;
    TObject *obj;
    while ((key = (TKey *)next()))
    {
        TString tempStr(key->GetName());
        if (tempStr.Contains("Temp"))
            continue;
        if (o == 1 && (tempStr.Contains("Parity") || tempStr.Contains("Delta")))
            continue;
        obj = gROOT->Get(key->GetName());
        if (!obj)
            continue;
        if (obj->IsA() == TDirectory::Class())
        {
            delete obj;
            obj = NULL;
            continue;
        }
        obj->Write();
    }

    if (o != 1)
    {
        for (int i = 0; i < 5; i++)
        {
            sprintf(fname_out, "Hist_v2_pt_obs2_caysm_%d", i);
            Hist_v2_pt_obs2_caysm[i]->Write(fname_out, TObject::kOverwrite);
            sprintf(fname_out, "Hist_v2_pt_obs2_caysm_os_%d", i);
            Hist_v2_pt_obs2_caysm_os[i]->Write(fname_out, TObject::kOverwrite);
            sprintf(fname_out, "Hist_v2_pt_obs2_p_caysm_%d", i);
            Hist_v2_pt_obs2_p_caysm[i]->Write(fname_out, TObject::kOverwrite);
            sprintf(fname_out, "Hist_pt_Ach_%d", i);
            Hist_pt_Ach[i]->Write(fname_out, TObject::kOverwrite);
        }

        num_lam_tree->Write();
        num_antilam_tree->Write();
        num_lam_rot_tree->Write();
        num_antilam_rot_tree->Write();
        num_proton_tree->Write();
        num_antiproton_tree->Write();
        num_lam_final->Write();
        num_proton_final->Write();
        num_gamma_final->Write();

        num_proton_used->Write();
        num_antiproton_used->Write();

        // for (int indl = 0; indl < 6; indl++)
        // {
        //     for (int f = 0; f < 8; f++)
        //     {
        //         for (int l = 0; l < 15; l++)
        //         {
        //             pParity_int_obs3_splitpt[indl][l][f]->Write();
        //         }
        //     }
        // }

        Hist_112ss_v2_parent_obs2->Write("Hist_112ss_v2_parent_obs2", TObject::kOverwrite);
        Hist_112os_v2_parent_obs2->Write("Hist_112os_v2_parent_obs2", TObject::kOverwrite);
        Hist_132ss_v2_parent_obs2->Write("Hist_132ss_v2_parent_obs2", TObject::kOverwrite);
        Hist_132os_v2_parent_obs2->Write("Hist_132os_v2_parent_obs2", TObject::kOverwrite);
        Hist_112ss_pt_parent_obs2->Write("Hist_112ss_pt_parent_obs2", TObject::kOverwrite);
        Hist_112os_pt_parent_obs2->Write("Hist_112os_pt_parent_obs2", TObject::kOverwrite);
        Hist_132ss_pt_parent_obs2->Write("Hist_132ss_pt_parent_obs2", TObject::kOverwrite);
        Hist_132os_pt_parent_obs2->Write("Hist_132os_pt_parent_obs2", TObject::kOverwrite);
        Hist_112ss_v2_parent_obs2_low->Write("Hist_112ss_v2_parent_obs2_low", TObject::kOverwrite);
        Hist_112os_v2_parent_obs2_low->Write("Hist_112os_v2_parent_obs2_low", TObject::kOverwrite);
        Hist_132ss_v2_parent_obs2_low->Write("Hist_132ss_v2_parent_obs2_low", TObject::kOverwrite);
        Hist_132os_v2_parent_obs2_low->Write("Hist_132os_v2_parent_obs2_low", TObject::kOverwrite);
        Hist_112ss_v2_parent_obs2_high->Write("Hist_112ss_v2_parent_obs2_high", TObject::kOverwrite);
        Hist_112os_v2_parent_obs2_high->Write("Hist_112os_v2_parent_obs2_high", TObject::kOverwrite);
        Hist_132ss_v2_parent_obs2_high->Write("Hist_132ss_v2_parent_obs2_high", TObject::kOverwrite);
        Hist_132os_v2_parent_obs2_high->Write("Hist_132os_v2_parent_obs2_high", TObject::kOverwrite);

        Hist_112ss_m_parent_obs2->Write("Hist_112ss_m_parent_obs2", TObject::kOverwrite);
        Hist_112os_m_parent_obs2->Write("Hist_112os_m_parent_obs2", TObject::kOverwrite);
        Hist_v2_parent_lam_pT->Write("Hist_v2_parent_lam_pT", TObject::kOverwrite);
        Hist_v2_parent_parent_pT->Write("Hist_v2_parent_parent_pT", TObject::kOverwrite);
        Hist_v2_parent_parent_pT_ss->Write("Hist_v2_parent_parent_pT_ss", TObject::kOverwrite);
        Hist_v2_parent_parent_pT_os->Write("Hist_v2_parent_parent_pT_os", TObject::kOverwrite);
        Hist_v2_lam_parent_pT->Write("Hist_v2_lam_parent_pT", TObject::kOverwrite);
        Hist_v2_p_parent_pT->Write("Hist_v2_p_parent_pT", TObject::kOverwrite);

        Hist_parent_phi_low_pT->Write("Hist_parent_phi_low_pT", TObject::kOverwrite);
        Hist_parent_phi_high_pT->Write("Hist_parent_phi_high_pT", TObject::kOverwrite);

        Hist_delta_phi_parent_pT_1->Write("Hist_delta_phi_parent_pT_1", TObject::kOverwrite);
        Hist_delta_phi_parent_pT_2->Write("Hist_delta_phi_parent_pT_2", TObject::kOverwrite);
        Hist_delta_phi_parent_pT_3->Write("Hist_delta_phi_parent_pT_3", TObject::kOverwrite);
        Hist_delta_phi_parent_pT_4->Write("Hist_delta_phi_parent_pT_4", TObject::kOverwrite);
        Hist_delta_phi_parent_pT_5->Write("Hist_delta_phi_parent_pT_5", TObject::kOverwrite);
        Hist_delta_phi_parent_pT_6->Write("Hist_delta_phi_parent_pT_6", TObject::kOverwrite);
        Hist_delta_phi_parent_pT_7->Write("Hist_delta_phi_parent_pT_7", TObject::kOverwrite);

        Hist_QQ_dist->Write("Hist_QQ_dist", TObject::kOverwrite);
        Hist_QQ_dist_coupling->Write("Hist_QQ_dist_coupling", TObject::kOverwrite);

        Hist_gamma112ss_QQ->Write();
        Hist_gamma112os_QQ->Write();
        Hist_gamma112ss_k->Write();
        Hist_gamma112os_k->Write();
        Hist_gamma132ss_QQ->Write();
        Hist_gamma132os_QQ->Write();
        Hist_gamma112ss_lpt->Write();
        Hist_gamma112os_lpt->Write();
        Hist_gamma112ss_ppt->Write();
        Hist_gamma112os_ppt->Write();
        Hist_gamma112ss_pi_QQ->Write();
        Hist_gamma112os_pi_QQ->Write();
        Hist_gamma112ss_pibar_QQ->Write();
        Hist_gamma112os_pibar_QQ->Write();
        Hist_lam_v2_low_parent_v2->Write();
        Hist_p_v2_low_parent_v2->Write();
        Hist_lam_pT_low_parent_v2->Write();
        Hist_p_pT_low_parent_v2->Write();
        V0Mass_gamma112->Write();
        V0Mass_gamma112_rot->Write();
        Hist_Pt->Write();
        Hist_Pt_rot->Write();
        Hist_Pt2->Write();
        Hist_Pt2_rot->Write();
        Hist_Pt2_anti->Write();
        Hist_Pt2_part->Write();

        Hist_charge1_dist->Write();
        Hist_charge2_dist->Write();
        Hist_charge1_ss_dist->Write();
        Hist_charge2_ss_dist->Write();
        Hist_charge1_os_dist->Write();
        Hist_charge2_os_dist->Write();

        num_gamma_ss_pre->Write();
        num_gamma_os_pre->Write();
        num_gamma_ss_final->Write();
        num_gamma_os_final->Write();

        nHitsFitQA->Write();
        nSigma_dauQA->Write();
        dca_protonQA->Write();
        dca_pionQA->Write();
        dca_LambdaQA->Write();
        nSigma_prim_protonQA->Write();

        proton_overlap_ratio->Write();

        hEtaPtDist->Write();
        hEtaPtDist_anti->Write();

        Hist_check_trksplitting_mag_ss->Write();
        Hist_check_trksplitting_mag_os->Write();
        Hist_check_trksplitting_pmom_ss->Write();
        Hist_check_trksplitting_pmom_os->Write();
        Hist_check_trksplitting_ptrkid_ss->Write();
        Hist_check_trksplitting_ptrkid_os->Write();
        Hist_check_trksplitting_ptrkid2_ss->Write();
        Hist_check_trksplitting_ptrkid2_os->Write();

        Hist_Q2_parent_tp->Write();
        Hist_Q2_parent_p->Write();
        Hist_Q2_parent_ap->Write();
        Hist_Q2_parent_lam->Write();
        Hist_ngamma_wrange->Write();
        Hist_nlam_wrange->Write();
        Hist_np_wrange->Write();
        Hist_nap_wrange->Write();
        Hist_ntp_wrange->Write();
        Hist_check_trksplitting_flowvector_ss->Write();
        Hist_check_trksplitting_flowvector_os->Write();
        Hist_Q2_parent_diff->Write();
        Hist_Q2_parent_pair_ss->Write();
        Hist_Q2_parent_pair_os->Write();

        Hist_check_trksplitting_phi_ss->Write();
        Hist_check_trksplitting_phi_os->Write();

        Hist_check_trksplitting_phi_all_ss->Write();
        Hist_check_trksplitting_phi_all_os->Write();

        Hist_check_trksplitting_phi_p->Write();
        Hist_check_trksplitting_phi_ap->Write();

        Hist_phi_below1_ss->Write();
        Hist_phi_after1_ss->Write();
        Hist_phi_both1_ss->Write();
        Hist_phi_below1_os->Write();
        Hist_phi_after1_os->Write();
        Hist_phi_both1_os->Write();

        Hist_Q2_parent_vs_Qcount_parent->Write();
        nHits_ratio_under_peak->Write();

        // Hist_v2_pt_obs1_alltrks->Write();
        // Hist_v2_pt_obs2_alltrks->Write();
        // Hist_v2_pt_obs1_alltrks_1stOrder->Write();
        // Hist_v2_pt_obs2_alltrks_1stOrder->Write();

        Hist_Q2_TPC->Write();
        Hist_Q2_EPD->Write();
        Hist_Q2_EPD1->Write();
        Hist_cos_EPD->Write();
        Hist_cos->Write();

        for (int lm = 0; lm < 3; lm++)
        {
            Hist_v2_pt_obs1_alltrks[lm]->Write();
            Hist_v2_pt_obs2_alltrks[lm]->Write();
            Hist_v2parent_pt_obs5[lm]->Write();
            Hist_v2parent_pt_obs5_rot[lm]->Write();

            Hist_v2pion_pt_obs5[lm]->Write();
        }

        Hist_v2_pt_obs1->Write();
        Hist_v2_pt_obs2->Write();
        Hist_v2_eta_obs1->Write();
        Hist_v2_eta_obs2->Write();
        Hist_v2_pt_EPD1_obs->Write();
        Hist_v2_pt_EPD_obs->Write();
        Hist_v2_eta_EPD1_obs->Write();
        Hist_v2_eta_EPD_obs->Write();
        Hist_v2_eta_EPDe_obs->Write();
        Hist_v2_eta_EPDw_obs->Write();
        Hist_v2_pt_obs1_rot->Write();
        Hist_v2_pt_obs2_rot->Write();
        Hist_v2_eta_obs1_rot->Write();
        Hist_v2_eta_obs2_rot->Write();
        Hist_v2_pt_EPD1_obs_rot->Write();
        Hist_v2_pt_EPD_obs_rot->Write();
        Hist_v2_eta_EPD1_obs_rot->Write();
        Hist_v2_eta_EPD_obs_rot->Write();
        Hist_v2_eta_EPDe_obs_rot->Write();
        Hist_v2_eta_EPDw_obs_rot->Write();

        for (int b = 0; b < 6; b++)
        {
            p_v2_Q2_obs1[b]->Write();
            p_v2_Q2_obs2[b]->Write();
            p_v2_Q2_obs1_rot[b]->Write();
            p_v2_Q2_obs2_rot[b]->Write();
        }

        pParity_e_Q2_obs1->Write();
        pParity_e_Q2_obs2->Write();
        pParity_e_Q2_obs1_rot->Write();
        pParity_e_Q2_obs2_rot->Write();
        pParity_w_Q2_obs1->Write();
        pParity_w_Q2_obs2->Write();
        pParity_w_Q2_obs1_rot->Write();
        pParity_w_Q2_obs2_rot->Write();

        Hist_v2_v2parent_3->Write();
        Hist_v2_v2parent_3_rot->Write();

        p_v2parente_Q2parent_obs1->Write();
        p_v2parente_Q2parent_obs2->Write();
        p_v2parentw_Q2parent_obs1->Write();
        p_v2parentw_Q2parent_obs2->Write();
        p_v2parente_Q2parent_obs1_rot->Write();
        p_v2parente_Q2parent_obs2_rot->Write();
        p_v2parentw_Q2parent_obs1_rot->Write();
        p_v2parentw_Q2parent_obs2_rot->Write();

        p_v2parente_Q2parent_EPD_obs1->Write();
        p_v2parente_Q2parent_EPD_obs2->Write();
        p_v2parentw_Q2parent_EPD_obs1->Write();
        p_v2parentw_Q2parent_EPD_obs2->Write();
        p_v2parente_Q2parent_EPD1_obs1->Write();
        p_v2parente_Q2parent_EPD1_obs2->Write();
        p_v2parentw_Q2parent_EPD1_obs1->Write();
        p_v2parentw_Q2parent_EPD1_obs2->Write();
        p_v2parente_Q2parent_EPD_obs1_rot->Write();
        p_v2parente_Q2parent_EPD_obs2_rot->Write();
        p_v2parentw_Q2parent_EPD_obs1_rot->Write();
        p_v2parentw_Q2parent_EPD_obs2_rot->Write();
        p_v2parente_Q2parent_EPD1_obs1_rot->Write();
        p_v2parente_Q2parent_EPD1_obs2_rot->Write();
        p_v2parentw_Q2parent_EPD1_obs1_rot->Write();
        p_v2parentw_Q2parent_EPD1_obs2_rot->Write();

        p_v2e_Q2parent_EPD_obs1->Write();
        p_v2e_Q2parent_EPD_obs2->Write();
        p_v2w_Q2parent_EPD_obs1->Write();
        p_v2w_Q2parent_EPD_obs2->Write();
        p_v2e_Q2parent_EPD1_obs1->Write();
        p_v2e_Q2parent_EPD1_obs2->Write();
        p_v2w_Q2parent_EPD1_obs1->Write();
        p_v2w_Q2parent_EPD1_obs2->Write();

        p_Parity_v2e_parent_obs1->Write();
        p_Parity_v2e_parent_obs2->Write();
        p_Parity_v2e_parent_obs3->Write();
        p_Parity_v2e_parent_obs4->Write();
        p_Parity_v2w_parent_obs1->Write();
        p_Parity_v2w_parent_obs2->Write();
        p_Parity_v2w_parent_obs3->Write();
        p_Parity_v2w_parent_obs4->Write();
        p_Parity_v2e_parent_obs1_rot->Write();
        p_Parity_v2e_parent_obs2_rot->Write();
        p_Parity_v2e_parent_obs3_rot->Write();
        p_Parity_v2e_parent_obs4_rot->Write();
        p_Parity_v2w_parent_obs1_rot->Write();
        p_Parity_v2w_parent_obs2_rot->Write();
        p_Parity_v2w_parent_obs3_rot->Write();
        p_Parity_v2w_parent_obs4_rot->Write();

        Hist_parentv2_event->Write();
        Hist_parentv2_event_rot->Write();

        pParity_e_Q2parent_obs1->Write();
        pParity_e_Q2parent_obs2->Write();
        pParity_w_Q2parent_obs1->Write();
        pParity_w_Q2parent_obs2->Write();
        pParity_e_Q2parent_obs1_rot->Write();
        pParity_e_Q2parent_obs2_rot->Write();
        pParity_w_Q2parent_obs1_rot->Write();
        pParity_w_Q2parent_obs2_rot->Write();

        for (int k = 0; k < 3; k++)
        {
            V0Mass_Q2[k]->Write();
            V0Mass_Q2_anti[k]->Write();
            V0Mass_Q2_pion[k]->Write();
            V0Mass_Q2_pion_anti[k]->Write();
        }

        p_v2e_Q2parent_obs1->Write();
        p_v2e_Q2parent_obs2->Write();
        p_v2w_Q2parent_obs1->Write();
        p_v2w_Q2parent_obs2->Write();
        p_v2e_Q2parent_EPD_obs1->Write();
        p_v2e_Q2parent_EPD_obs2->Write();
        p_v2w_Q2parent_EPD_obs1->Write();
        p_v2w_Q2parent_EPD_obs2->Write();
        p_v2e_Q2parent_EPD1_obs1->Write();
        p_v2e_Q2parent_EPD1_obs2->Write();

        p_v2e_Q2_obs1->Write();
        p_v2e_Q2_obs2->Write();
        p_v2w_Q2_obs1->Write();
        p_v2w_Q2_obs2->Write();
        p_v2e_Q2_EPD_obs1->Write();
        p_v2e_Q2_EPD_obs2->Write();
        p_v2w_Q2_EPD_obs1->Write();
        p_v2w_Q2_EPD_obs2->Write();
        p_v2e_Q2_EPD1_obs1->Write();
        p_v2e_Q2_EPD1_obs2->Write();

        p_pionv2e_Q2_obs1->Write();
        p_pionv2e_Q2_obs2->Write();
        p_pionv2w_Q2_obs1->Write();
        p_pionv2w_Q2_obs2->Write();
        p_pionv2e_Q2_EPD_obs1->Write();
        p_pionv2e_Q2_EPD_obs2->Write();
        p_pionv2w_Q2_EPD_obs1->Write();
        p_pionv2w_Q2_EPD_obs2->Write();
        p_pionv2e_Q2_EPD1_obs1->Write();
        p_pionv2e_Q2_EPD1_obs2->Write();

        p_cos_Q2->Write();
        p_cos_Q2_EPD->Write();
        p_cos_Q2_EPD1->Write();
        p_cos_Q2_TPC_pion->Write();
        p_cos_Q2_EPD_pion->Write();
        p_cos_Q2_EPD1_pion->Write();

        Hist_Q2_EPD_pion->Write();
        Hist_Q2_EPD1_pion->Write();
        Hist_Q2_TPC_pion->Write();

        v2_single_lamda->Write();
        v2_single_proton->Write();
        v2_single_lamp->Write();
        v2_single_lamda_err->Write();
        v2_single_proton_err->Write();
        v2_single_lamp_err->Write();

        v2_2_pion->Write();

        pDelta_int_ss_obs1->Write();
        pDelta_int_ss_obs1_rot->Write();
        pDelta_int_ss_obs3->Write();
        pDelta_int_ss_obs3_rot->Write();
        pDelta_int_ss_obs1_anti->Write();
        pDelta_int_ss_obs1_anti_rot->Write();
        pDelta_int_ss_obs3_anti->Write();
        pDelta_int_ss_obs3_anti_rot->Write();
        pDelta_int_ss_obs1_QQcut->Write();
        pDelta_int_ss_obs1_rot_QQcut->Write();
        pDelta_int_ss_obs3_QQcut->Write();
        pDelta_int_ss_obs3_rot_QQcut->Write();
        pDelta_int_ss_obs1_QQcut_anti->Write();
        pDelta_int_ss_obs1_rot_QQcut_anti->Write();
        pDelta_int_ss_obs3_QQcut_anti->Write();
        pDelta_int_ss_obs3_rot_QQcut_anti->Write();
    }

    fout->Write();
    fout->Close();
}
//////////////////////////////////////////////////
void WriteWeight(TString OutFileName)
{
    TFile *fWgtNew = new TFile(OutFileName, "UPDATE");
    Hist_netChAsym->Write();
    pTPCmeanPhi_FF_1->Write();
    pTPCmeanPhi_RF_1->Write();
    pTPCmeanPhi_FF_1_p->Write();
    pTPCmeanPhi_RF_1_p->Write();
    pTPCmeanPhi_FF_1_rot->Write();
    pTPCmeanPhi_RF_1_rot->Write();
    pTPCmeanPhiAsso_FF_1->Write();
    pTPCmeanPhiAsso_RF_1->Write();
    pTPCmeanPhi_FF_2->Write();
    pTPCmeanPhi_RF_2->Write();
    pTPCmeanPhi_FF_2_p->Write();
    pTPCmeanPhi_RF_2_p->Write();
    pTPCmeanPhi_FF_2_rot->Write();
    pTPCmeanPhi_RF_2_rot->Write();
    pTPCmeanPhiAsso_FF_2->Write();
    pTPCmeanPhiAsso_RF_2->Write();
    pTPCmeanPhi_FF_3->Write();
    pTPCmeanPhi_RF_3->Write();
    pTPCmeanPhi_FF_3_p->Write();
    pTPCmeanPhi_RF_3_p->Write();
    pTPCmeanPhi_FF_3_rot->Write();
    pTPCmeanPhi_RF_3_rot->Write();
    pTPCmeanPhiAsso_FF_3->Write();
    pTPCmeanPhiAsso_RF_3->Write();
    pTPC_EP_east->Write();
    pTPC_EP_west->Write();
    pTPC_EP_for->Write();
    pTPC_EP_bac->Write();
    pTPC_EP_for_pos->Write();
    pTPC_EP_bac_pos->Write();
    pTPC_EP_for_neg->Write();
    pTPC_EP_bac_neg->Write();
    pTPC_EP_full->Write();
    pEPD_EP1_east->Write();
    pEPD_EP1_west->Write();
    pEPD_EP_east->Write();
    pEPD_EP_west->Write();

    //  rc->Write();
    fWgtNew->Close();
}
/////////////////////////////////////////////////////

//////////////////////////////////////////////////
void FillGamma(int ord)
{
    if (ord == 4)
    {
        pParity_int_ss_oppo_run->Fill(Day2, correlator4[0], Eweight / eff / eff2);
        pDelta_int_ss_oppo_run->Fill(Day2, correlator3, Eweight / eff / eff2);
    }
    if (ord == 3)
    {
        pParity_int_ss_same_run->Fill(Day2, correlator4[0], Eweight / eff / eff2);
        pDelta_int_ss_same_run->Fill(Day2, correlator3, Eweight / eff / eff2);
    }

    int qq_index = 0, lambda_index = 0, charge_index = -1;

    if (Pt < 2.0)
        lambda_index = (int)floor((Pt - 0.5) / 0.1);
    else if (Pt == 2.0)
        lambda_index = 14;

    if (Charge == 1)
        charge_index = 0;
    else if (Charge == -1)
        charge_index = 4;

    pParity_int_obs3_splitpt[2][lambda_index][charge_index]->Fill(ord, correlator3, Eweight / temp_eff / temp_proton_eff);
    pParity_int_obs3_splitpt[5][lambda_index][charge_index]->Fill(ord, correlator3);

    for (int tp = 0; tp < 3; tp++)
    {
        pParity_int_obs3_splitpt[0][lambda_index][charge_index]->Fill(ord + 4 + (8 * tp), correlator0_alt[tp], Eweight / temp_eff / temp_proton_eff);
        pParity_int_obs3_splitpt[0][lambda_index][charge_index]->Fill(ord + (8 * tp), correlator0[tp], Eweight / temp_eff / temp_proton_eff);
        pParity_int_obs3_splitpt[1][lambda_index][charge_index]->Fill(ord + (4 * tp), correlator4[tp], Eweight / temp_eff / temp_proton_eff);
        pParity_int_obs3_splitpt[3][lambda_index][charge_index]->Fill(ord + 4 + (8 * tp), correlator0_alt[tp]);
        pParity_int_obs3_splitpt[3][lambda_index][charge_index]->Fill(ord + (8 * tp), correlator0[tp]);
        pParity_int_obs3_splitpt[4][lambda_index][charge_index]->Fill(ord + (4 * tp), correlator4[tp]);
    }

    if (Charge == 1)
    {
        pDelta_int_ss_obs1->Fill(ord, correlator3);
        pDelta_int_ss_obs3->Fill(ord, correlator3, Eweight / temp_eff / temp_proton_eff);

        for (int tp1 = 0; tp1 < 3; tp1++)
        {
            pParity_int_obs1->Fill(ord + (8 * tp1), correlator0[tp1]);
            pParity_int_obs1->Fill(ord + 4 + (8 * tp1), correlator0_alt[tp1]);
            pParity_int_obs3->Fill(ord + (8 * tp1), correlator0[tp1], Eweight / temp_eff / temp_proton_eff);
            pParity_int_obs3->Fill(ord + 4 + (8 * tp1), correlator0_alt[tp1], Eweight / temp_eff / temp_proton_eff);
            pParity_int_ss_obs1->Fill(ord + (4 * tp1), correlator4[tp1]);
            pParity_int_ss_obs3->Fill(ord + (4 * tp1), correlator4[tp1], Eweight / temp_eff / temp_proton_eff);
        }
    }
    else if (Charge == -1)
    {
        pDelta_int_ss_obs1_anti->Fill(ord, correlator3);
        pDelta_int_ss_obs3_anti->Fill(ord, correlator3, Eweight / temp_eff / temp_proton_eff);

        for (int tp2 = 0; tp2 < 3; tp2++)
        {
            pParity_int_obs1_anti->Fill(ord + (8 * tp2), correlator0[tp2]);
            pParity_int_obs1_anti->Fill(ord + 4 + (8 * tp2), correlator0_alt[tp2]);
            pParity_int_obs3_anti->Fill(ord + (8 * tp2), correlator0[tp2], Eweight / temp_eff / temp_proton_eff);
            pParity_int_obs3_anti->Fill(ord + 4 + (8 * tp2), correlator0_alt[tp2], Eweight / temp_eff / temp_proton_eff);
            pParity_int_ss_obs1_anti->Fill(ord + (4 * tp2), correlator4[tp2]);
            pParity_int_ss_obs3_anti->Fill(ord + (4 * tp2), correlator4[tp2], Eweight / temp_eff / temp_proton_eff);
        }
    }

    pParity_eta_ss_obs1->Fill(ord, 0.5 * (Eta + Eta2), correlator4[0]);
    pParity_eta_ss_obs3->Fill(ord, 0.5 * (Eta + Eta2), correlator4[0], Eweight / eff / eff2);
    pParity_Deta_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4[0]);
    pParity_Deta_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4[0], Eweight / eff / eff2);

    // new graphs
    pParity_pt_ss_obs1->Fill(ord, Pt, correlator0[0], Eweight / eff / eff2); // Gamma132
    pParity_pt_ss_obs3->Fill(ord, Pt, correlator4[0], Eweight / eff / eff2); // Gamma112
    pDelta_pt_ss_obs3->Fill(ord, Pt, correlator3, Eweight / eff / eff2);

    pParity_Dpt_ss_obs1->Fill(ord, fabs(Pt - Pt2), correlator4[0]);
    pParity_Dpt_ss_obs3->Fill(ord, fabs(Pt - Pt2), correlator4[0], Eweight / eff / eff2);
    pDelta_eta_ss_obs1->Fill(ord, 0.5 * (Eta + Eta2), correlator3);
    pDelta_eta_ss_obs3->Fill(ord, 0.5 * (Eta + Eta2), correlator3, Eweight / eff / eff2);
    pDelta_Deta_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator3);
    pDelta_Deta_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator3, Eweight / eff / eff2);
    pDelta_pt_ss_obs1->Fill(ord, 0.5 * (Pt + Pt2), correlator3);
    pDelta_Dpt_ss_obs1->Fill(ord, fabs(Pt - Pt2), correlator3);
    pDelta_Dpt_ss_obs3->Fill(ord, fabs(Pt - Pt2), correlator3, Eweight / eff / eff2);

    for (int tp3 = 0; tp3 < 3; tp3++)
    {
        pTemp_parity_e->Fill(ord + (12 * tp3), correlator4e[tp3], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_w->Fill(ord + (12 * tp3), correlator4w[tp3], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_e->Fill(ord + 4 + (12 * tp3), correlator0e[tp3], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_w->Fill(ord + 4 + (12 * tp3), correlator0w[tp3], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_e->Fill(ord + 8 + (12 * tp3), correlator0e_alt[tp3], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_w->Fill(ord + 8 + (12 * tp3), correlator0w_alt[tp3], 1. / temp_eff / temp_proton_eff);
    }

    pTemp_delta->Fill(ord, correlator3, 1. / temp_eff / temp_proton_eff);

    if (fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
    {
        pTemp_parity_e_noHBT->Fill(ord, correlator4e[0], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_w_noHBT->Fill(ord, correlator4w[0], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_e_noHBT->Fill(ord + 4, correlator4e[0], 1. / temp_eff / temp_proton_eff);
        pTemp_parity_w_noHBT->Fill(ord + 4, correlator4w[0], 1. / temp_eff / temp_proton_eff);
        pTemp_delta_noHBT->Fill(ord, correlator3, 1. / temp_eff / temp_proton_eff);
    }

    if (fabs(Pt - Pt2) > 0.15)
    {
        pParity_Deta_highDpt_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4[0]);
        pParity_Deta_highDpt_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4[0], Eweight / eff / eff2);
        pDelta_Deta_highDpt_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4[0]);
        pDelta_Deta_highDpt_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4[0], Eweight / eff / eff2);
    }
    if (fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
    {
        pParity_noHBT_ss_obs1->Fill(ord, correlator4[0]);
        pParity_noHBT_ss_obs3->Fill(ord, correlator4[0], Eweight / eff / eff2);
        pDelta_noHBT_ss_obs1->Fill(ord, correlator3);
        pDelta_noHBT_ss_obs3->Fill(ord, correlator3, Eweight / eff / eff2);
    }

    if (ord == 3)
    {
        pParity_int_pt_ss_obs3->Fill(Pt, correlator4[0], Eweight / eff / eff2);
        pParity_int_pt_ssb_obs3->Fill(Pt, correlator0[0], Eweight / eff / eff2);
    }
    if (ord == 4)
    {
        pParity_int_pt_os_obs3->Fill(Pt, correlator4[0], Eweight / eff / eff2);
        pParity_int_pt_osb_obs3->Fill(Pt, correlator0[0], Eweight / eff / eff2);
    }

    if ((QQ <= 0.8))
        return;

    int lambda_index_QQcut = 0;

    if (Pt < 2.0)
        lambda_index_QQcut = (int)floor((Pt - 0.5) / 0.1);
    else if (Pt == 2.0)
        lambda_index_QQcut = 14;

    pParity_int_obs3_splitpt[2][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator3, Eweight / temp_eff / temp_proton_eff);
    pParity_int_obs3_splitpt[5][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator3);

    pParity_int_obs3_splitpt[0][lambda_index_QQcut][charge_index + 2]->Fill(ord + 4, correlator0_alt[0], Eweight / temp_eff / temp_proton_eff);
    pParity_int_obs3_splitpt[0][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator0[0], Eweight / temp_eff / temp_proton_eff);
    pParity_int_obs3_splitpt[1][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator4[0], Eweight / temp_eff / temp_proton_eff);
    pParity_int_obs3_splitpt[3][lambda_index_QQcut][charge_index + 2]->Fill(ord + 4, correlator0_alt[0]);
    pParity_int_obs3_splitpt[3][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator0[0]);
    pParity_int_obs3_splitpt[4][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator4[0]);

    if (Charge == 1)
    {
        pDelta_int_ss_obs1_QQcut->Fill(ord, correlator3);
        pDelta_int_ss_obs3_QQcut->Fill(ord, correlator3, Eweight / temp_eff / temp_proton_eff);

        for (int tp1 = 0; tp1 < 3; tp1++)
        {
            pParity_int_obs1_QQcut->Fill(ord + (8 * tp1), correlator0[tp1]);
            pParity_int_obs1_QQcut->Fill(ord + 4 + (8 * tp1), correlator0_alt[tp1]);
            pParity_int_obs3_QQcut->Fill(ord + (8 * tp1), correlator0[tp1], Eweight / temp_eff / temp_proton_eff);
            pParity_int_obs3_QQcut->Fill(ord + 4 + (8 * tp1), correlator0_alt[tp1], Eweight / temp_eff / temp_proton_eff);
            pParity_int_ss_obs1_QQcut->Fill(ord + (4 * tp1), correlator4[tp1]);
            pParity_int_ss_obs3_QQcut->Fill(ord + (4 * tp1), correlator4[tp1], Eweight / temp_eff / temp_proton_eff);
        }
    }
    else if (Charge == -1)
    {
        pDelta_int_ss_obs1_QQcut_anti->Fill(ord, correlator3);
        pDelta_int_ss_obs3_QQcut_anti->Fill(ord, correlator3, Eweight / temp_eff / temp_proton_eff);

        for (int tp2 = 0; tp2 < 3; tp2++)
        {
            pParity_int_obs1_QQcut_anti->Fill(ord + (8 * tp2), correlator0[tp2]);
            pParity_int_obs1_QQcut_anti->Fill(ord + 4 + (8 * tp2), correlator0_alt[tp2]);
            pParity_int_obs3_QQcut_anti->Fill(ord + (8 * tp2), correlator0[tp2], Eweight / temp_eff / temp_proton_eff);
            pParity_int_obs3_QQcut_anti->Fill(ord + 4 + (8 * tp2), correlator0_alt[tp2], Eweight / temp_eff / temp_proton_eff);
            pParity_int_ss_obs1_QQcut_anti->Fill(ord + (4 * tp2), correlator4[tp2]);
            pParity_int_ss_obs3_QQcut_anti->Fill(ord + (4 * tp2), correlator4[tp2], Eweight / temp_eff / temp_proton_eff);
        }
    }

    pTemp_parity_e_QQcut->Fill(ord, correlator4e[0], 1. / temp_eff / temp_proton_eff);
    pTemp_parity_w_QQcut->Fill(ord, correlator4w[0], 1. / temp_eff / temp_proton_eff);
    pTemp_delta_QQcut->Fill(ord, correlator3, 1. / temp_eff / temp_proton_eff);
}

/////////////////////////////////////////////////////
void Fillv2_rot()
{
    pTemp_v2_rot->Fill(1, v2e_rot, 1. / temp_eff_rot);
    pTemp_v2_rot->Fill(2, v2w_rot, 1. / temp_eff_rot);
    if (Charge_rot > 0)
    {
        pTemp_v2_rot->Fill(3, v2e_rot, 1. / temp_eff_rot);
        pTemp_v2_rot->Fill(4, v2w_rot, 1. / temp_eff_rot);
    }
    if (Charge_rot < 0)
    {
        pTemp_v2_rot->Fill(5, v2e_rot, 1. / temp_eff_rot);
        pTemp_v2_rot->Fill(6, v2w_rot, 1. / temp_eff_rot);
    }
    pTemp_v2_rot->Fill(7, v2_rot, 1. / temp_eff_rot);

    pTemp_v2_rot->Fill(11, v2_EPD1_rot, 1. / temp_eff_rot);
    pTemp_v2_rot->Fill(12, v2_EPDe_rot, 1. / temp_eff_rot);
    pTemp_v2_rot->Fill(13, v2_EPDw_rot, 1. / temp_eff_rot);
}

//////////////////////////////////////////////////

//////////////////////////////////////////////////
void FillGamma_rot(int ord)
{
    if (debug_2)
        std::cout << "temp_eff_rot = " << temp_eff_rot << endl;
    if (debug_2)
        std::cout << "temp_proton_eff_rot = " << temp_proton_eff_rot << endl;

    int qq_index = 0, lambda_index = 0, charge_index = -1;

    if (Pt_rot < 2.0)
        lambda_index = (int)floor((Pt_rot - 0.5) / 0.1);
    else if (Pt_rot == 2.0)
        lambda_index = 14;

    if (Charge_rot == 1)
        charge_index = 1;
    else if (Charge_rot == -1)
        charge_index = 5;

    pParity_int_obs3_splitpt[2][lambda_index][charge_index]->Fill(ord, correlator3_rot, Eweight / temp_eff_rot / temp_proton_eff_rot);
    pParity_int_obs3_splitpt[5][lambda_index][charge_index]->Fill(ord, correlator3_rot);

    for (int tp = 0; tp < 3; tp++)
    {
        pParity_int_obs3_splitpt[0][lambda_index][charge_index]->Fill(ord + 4 + (8 * tp), correlator0_alt_rot[tp], Eweight / temp_eff_rot / temp_proton_eff_rot);
        pParity_int_obs3_splitpt[0][lambda_index][charge_index]->Fill(ord + (8 * tp), correlator0_rot[tp], Eweight / temp_eff_rot / temp_proton_eff_rot);
        pParity_int_obs3_splitpt[1][lambda_index][charge_index]->Fill(ord + (4 * tp), correlator4_rot[tp], Eweight / temp_eff_rot / temp_proton_eff_rot);
        pParity_int_obs3_splitpt[3][lambda_index][charge_index]->Fill(ord + 4 + (8 * tp), correlator0_alt_rot[tp]);
        pParity_int_obs3_splitpt[3][lambda_index][charge_index]->Fill(ord + (8 * tp), correlator0_rot[tp]);
        pParity_int_obs3_splitpt[4][lambda_index][charge_index]->Fill(ord + (4 * tp), correlator4_rot[tp]);
    }

    if (Charge_rot == 1)
    {
        for (int tp1 = 0; tp1 < 3; tp1++)
        {
            pParity_int_obs1_rot->Fill(ord + (8 * tp1), correlator0_rot[tp1]);
            pParity_int_obs1_rot->Fill(ord + 4 + (8 * tp1), correlator0_alt_rot[tp1]);
            pParity_int_obs3_rot->Fill(ord + (8 * tp1), correlator0_rot[tp1], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_obs3_rot->Fill(ord + 4 + (8 * tp1), correlator0_alt_rot[tp1], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_ss_obs1_rot->Fill(ord + (4 * tp1), correlator4_rot[tp1]);
            pParity_int_ss_obs3_rot->Fill(ord + (4 * tp1), correlator4_rot[tp1], Eweight / temp_eff_rot / temp_proton_eff_rot);
        }

        pDelta_int_ss_obs1_rot->Fill(ord, correlator3_rot);
        pDelta_int_ss_obs3_rot->Fill(ord, correlator3_rot, Eweight / temp_eff_rot / temp_proton_eff_rot);
    }
    else if (Charge_rot == -1)
    {
        for (int tp2 = 0; tp2 < 3; tp2++)
        {
            pParity_int_obs1_anti_rot->Fill(ord + (8 * tp2), correlator0_rot[tp2]);
            pParity_int_obs1_anti_rot->Fill(ord + 4 + (8 * tp2), correlator0_alt_rot[tp2]);
            pParity_int_obs3_anti_rot->Fill(ord + (8 * tp2), correlator0_rot[tp2], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_obs3_anti_rot->Fill(ord + 4 + (8 * tp2), correlator0_alt_rot[tp2], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_ss_obs1_anti_rot->Fill(ord + (4 * tp2), correlator4_rot[tp2]);
            pParity_int_ss_obs3_anti_rot->Fill(ord + (4 * tp2), correlator4_rot[tp2], Eweight / temp_eff_rot / temp_proton_eff_rot);
        }

        pDelta_int_ss_obs1_anti_rot->Fill(ord, correlator3_rot);
        pDelta_int_ss_obs3_anti_rot->Fill(ord, correlator3_rot, Eweight / temp_eff_rot / temp_proton_eff_rot);
    }

    for (int tp3 = 0; tp3 < 3; tp3++)
    {
        pTemp_parity_e_rot->Fill(ord + (12 * tp3), correlator4e_rot[tp3], 1. / temp_eff_rot / temp_proton_eff_rot);
        pTemp_parity_w_rot->Fill(ord + (12 * tp3), correlator4w_rot[tp3], 1. / temp_eff_rot / temp_proton_eff_rot);
        pTemp_parity_e_rot->Fill(ord + 4 + (12 * tp3), correlator0e_rot[tp3], 1. / temp_eff_rot / temp_proton_eff_rot);
        pTemp_parity_w_rot->Fill(ord + 4 + (12 * tp3), correlator0w_rot[tp3], 1. / temp_eff_rot / temp_proton_eff_rot);
        pTemp_parity_e_rot->Fill(ord + 8 + (12 * tp3), correlator0e_alt_rot[tp3], 1. / temp_eff_rot / temp_proton_eff_rot);
        pTemp_parity_w_rot->Fill(ord + 8 + (12 * tp3), correlator0w_alt_rot[tp3], 1. / temp_eff_rot / temp_proton_eff_rot);
    }

    pTemp_delta_rot->Fill(ord, correlator3_rot, 1. / temp_eff_rot / temp_proton_eff_rot);

    if ((QQ_rot <= 0.8))
        return;

    int lambda_index_QQcut = 0;

    if (Pt_rot < 2.0)
        lambda_index_QQcut = (int)floor((Pt_rot - 0.5) / 0.1);
    else if (Pt_rot == 2.0)
        lambda_index_QQcut = 14;

    pParity_int_obs3_splitpt[2][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator3_rot, Eweight / temp_eff_rot / temp_proton_eff_rot);
    pParity_int_obs3_splitpt[5][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator3_rot);

    pParity_int_obs3_splitpt[0][lambda_index_QQcut][charge_index + 2]->Fill(ord + 4, correlator0_alt_rot[0], Eweight / temp_eff_rot / temp_proton_eff_rot);
    pParity_int_obs3_splitpt[0][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator0_rot[0], Eweight / temp_eff_rot / temp_proton_eff_rot);
    pParity_int_obs3_splitpt[1][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator4_rot[0], Eweight / temp_eff_rot / temp_proton_eff_rot);
    pParity_int_obs3_splitpt[3][lambda_index_QQcut][charge_index + 2]->Fill(ord + 4, correlator0_alt_rot[0]);
    pParity_int_obs3_splitpt[3][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator0_rot[0]);
    pParity_int_obs3_splitpt[4][lambda_index_QQcut][charge_index + 2]->Fill(ord, correlator4_rot[0]);

    if (Charge_rot == 1)
    {
        for (int tp1 = 0; tp1 < 3; tp1++)
        {
            pParity_int_obs1_rot_QQcut->Fill(ord + (8 * tp1), correlator0_rot[tp1]);
            pParity_int_obs1_rot_QQcut->Fill(ord + 4 + (8 * tp1), correlator0_alt_rot[tp1]);
            pParity_int_obs3_rot_QQcut->Fill(ord + (8 * tp1), correlator0_rot[tp1], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_obs3_rot_QQcut->Fill(ord + 4 + (8 * tp1), correlator0_alt_rot[tp1], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_ss_obs1_rot_QQcut->Fill(ord + (4 * tp1), correlator4_rot[tp1]);
            pParity_int_ss_obs3_rot_QQcut->Fill(ord + (4 * tp1), correlator4_rot[tp1], Eweight / temp_eff_rot / temp_proton_eff_rot);
        }

        pDelta_int_ss_obs1_rot_QQcut->Fill(ord, correlator3_rot);
        pDelta_int_ss_obs3_rot_QQcut->Fill(ord, correlator3_rot, Eweight / temp_eff_rot / temp_proton_eff_rot);
    }
    else if (Charge_rot == -1)
    {
        for (int tp2 = 0; tp2 < 3; tp2++)
        {
            pParity_int_obs1_rot_QQcut_anti->Fill(ord + (8 * tp2), correlator0_rot[tp2]);
            pParity_int_obs1_rot_QQcut_anti->Fill(ord + 4 + (8 * tp2), correlator0_alt_rot[tp2]);
            pParity_int_obs3_rot_QQcut_anti->Fill(ord + (8 * tp2), correlator0_rot[tp2], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_obs3_rot_QQcut_anti->Fill(ord + 4 + (8 * tp2), correlator0_alt_rot[tp2], Eweight / temp_eff_rot / temp_proton_eff_rot);
            pParity_int_ss_obs1_rot_QQcut_anti->Fill(ord + (4 * tp2), correlator4_rot[tp2]);
            pParity_int_ss_obs3_rot_QQcut_anti->Fill(ord + (4 * tp2), correlator4_rot[tp2], Eweight / temp_eff_rot / temp_proton_eff_rot);
        }

        pDelta_int_ss_obs1_rot_QQcut_anti->Fill(ord, correlator3_rot);
        pDelta_int_ss_obs3_rot_QQcut_anti->Fill(ord, correlator3_rot, Eweight / temp_eff_rot / temp_proton_eff_rot);
    }

    pTemp_parity_e_rot_QQcut->Fill(ord, correlator4e_rot[0], 1. / temp_eff_rot / temp_proton_eff_rot);
    pTemp_parity_w_rot_QQcut->Fill(ord, correlator4w_rot[0], 1. / temp_eff_rot / temp_proton_eff_rot);
    pTemp_delta_rot_QQcut->Fill(ord, correlator3_rot, 1. / temp_eff_rot / temp_proton_eff_rot);
}

//////////////////////////////////////////////////////
void WrapUpESE()
{
    //////////////////////// TPC Results ////////////////////////
    // v2 vs. Q2
    for (int tmp1 = 1; tmp1 <= 6; tmp1++)
    {
        if (pTemp_v2->GetBinContent(tmp1) != 0)
        {
            p_v2_Q2_obs1[tmp1 - 1]->Fill(Q2_proper, pTemp_v2->GetBinContent(tmp1));
            p_v2_Q2_obs2[tmp1 - 1]->Fill(Q2_proper, pTemp_v2->GetBinContent(tmp1), Eweight);
        }
    }

    // Gamma vs. Q2
    for (int tmp2 = 1; tmp2 <= 8; tmp2++)
    {
        if (pTemp_parity_e->GetBinContent(tmp2) != 0)
        {
            pParity_e_Q2_obs1->Fill(tmp2, Q2_proper, pTemp_parity_e->GetBinContent(tmp2));
            pParity_e_Q2_obs2->Fill(tmp2, Q2_proper, pTemp_parity_e->GetBinContent(tmp2), Eweight);
        }
        if (pTemp_parity_w->GetBinContent(tmp2) != 0)
        {
            pParity_w_Q2_obs1->Fill(tmp2, Q2_proper, pTemp_parity_w->GetBinContent(tmp2));
            pParity_w_Q2_obs2->Fill(tmp2, Q2_proper, pTemp_parity_w->GetBinContent(tmp2), Eweight);
        }
    }

    // Delta vs. Q2
    for (int tmp3 = 1; tmp3 <= 4; tmp3++)
    {
        if (pTemp_delta->GetBinContent(tmp3) != 0)
        {
            pDelta_Q2_obs1->Fill(tmp3, Q2_proper, pTemp_delta->GetBinContent(tmp3));
            pDelta_Q2_obs2->Fill(tmp3, Q2_proper, pTemp_delta->GetBinContent(tmp3), Eweight);
        }
    }

    ////////////////////////// Parent v2 vs. All Charged Hadrons Q2 //////////////////////////
    float Temp_v2e = pTemp_v2_parent->GetBinContent(1);
    float Temp_v2w = pTemp_v2_parent->GetBinContent(2);
    float Temp_v2d_e = pTemp_v2->GetBinContent(1);
    float Temp_v2d_w = pTemp_v2->GetBinContent(2);

    if ((Temp_v2e != 0) && (Temp_v2w != 0) && (Temp_v2d_e != 0) && (Temp_v2d_w != 0))
    {
        Hist_v2_v2parent_3->Fill(0.5 * (Temp_v2e + Temp_v2w), 0.5 * (Temp_v2d_e + Temp_v2d_w));
    }

    // TPC
    if (Temp_v2e != 0)
    {
        p_v2parente_Q2parent_obs1->Fill(Q2_proper, Temp_v2e);
        p_v2parente_Q2parent_obs2->Fill(Q2_proper, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2parentw_Q2parent_obs1->Fill(Q2_proper, Temp_v2w);
        p_v2parentw_Q2parent_obs2->Fill(Q2_proper, Temp_v2w, Eweight);
    }

    // EPD
    Temp_v2e = pTemp_v2_parent->GetBinContent(5);
    Temp_v2w = pTemp_v2_parent->GetBinContent(6);

    if (Temp_v2e != 0)
    {
        p_v2parente_Q2parent_EPD_obs1->Fill(Q2_proper_EPD, Temp_v2e);
        p_v2parente_Q2parent_EPD_obs2->Fill(Q2_proper_EPD, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2parentw_Q2parent_EPD_obs1->Fill(Q2_proper_EPD, Temp_v2w);
        p_v2parentw_Q2parent_EPD_obs2->Fill(Q2_proper_EPD, Temp_v2w, Eweight);
    }

    // EPD1
    Temp_v2e = pTemp_v2_parent->GetBinContent(7);

    if (Temp_v2e != 0)
    {
        p_v2parente_Q2parent_EPD1_obs1->Fill(Q2_proper_EPD1, Temp_v2e);
        p_v2parente_Q2parent_EPD1_obs2->Fill(Q2_proper_EPD1, Temp_v2e, Eweight);
    }
    //////////////////////////////////////////////////////////////////////////////

    ////////////////////////// Pair Hadrons v2 vs. All Charged Hadrons Q2 //////////////////////////
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(1);
    Temp_v2w = pTemp_v2_chargedhadrons->GetBinContent(2);

    // TPC
    if (Temp_v2e != 0)
    {
        p_v2e_Q2parent_obs1->Fill(Q2_proper, Temp_v2e);
        p_v2e_Q2parent_obs2->Fill(Q2_proper, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2w_Q2parent_obs1->Fill(Q2_proper, Temp_v2w);
        p_v2w_Q2parent_obs2->Fill(Q2_proper, Temp_v2w, Eweight);
    }

    // EPD
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(12);
    Temp_v2w = pTemp_v2_chargedhadrons->GetBinContent(13);

    if (Temp_v2e != 0)
    {
        p_v2e_Q2parent_EPD_obs1->Fill(Q2_proper_EPD, Temp_v2e);
        p_v2e_Q2parent_EPD_obs2->Fill(Q2_proper_EPD, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2w_Q2parent_EPD_obs1->Fill(Q2_proper_EPD, Temp_v2w);
        p_v2w_Q2parent_EPD_obs2->Fill(Q2_proper_EPD, Temp_v2w, Eweight);
    }

    // EPD1
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(11);

    if (Temp_v2e != 0)
    {
        p_v2e_Q2parent_EPD1_obs1->Fill(Q2_proper_EPD1, Temp_v2e);
        p_v2e_Q2parent_EPD1_obs2->Fill(Q2_proper_EPD1, Temp_v2e, Eweight);
    }
    //////////////////////////////////////////////////////////////////////////////

    ////////////////////////// Single Hadrons v2 vs. All Charged Hadrons Q2 //////////////////////////
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(3);
    Temp_v2w = pTemp_v2_chargedhadrons->GetBinContent(4);

    // TPC
    if (Temp_v2e != 0)
    {
        p_v2e_Q2_obs1->Fill(Q2_proper, Temp_v2e);
        p_v2e_Q2_obs2->Fill(Q2_proper, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2w_Q2_obs1->Fill(Q2_proper, Temp_v2w);
        p_v2w_Q2_obs2->Fill(Q2_proper, Temp_v2w, Eweight);
    }

    // EPD
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(5);
    Temp_v2w = pTemp_v2_chargedhadrons->GetBinContent(6);

    if (Temp_v2e != 0)
    {
        p_v2e_Q2_EPD_obs1->Fill(Q2_proper_EPD, Temp_v2e);
        p_v2e_Q2_EPD_obs2->Fill(Q2_proper_EPD, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2w_Q2_EPD_obs1->Fill(Q2_proper_EPD, Temp_v2w);
        p_v2w_Q2_EPD_obs2->Fill(Q2_proper_EPD, Temp_v2w, Eweight);
    }

    // EPD1
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(7);

    if (Temp_v2e != 0)
    {
        p_v2e_Q2_EPD1_obs1->Fill(Q2_proper_EPD1, Temp_v2e);
        p_v2e_Q2_EPD1_obs2->Fill(Q2_proper_EPD1, Temp_v2e, Eweight);
    }
    //////////////////////////////////////////////////////////////////////////////

    ////////////////////////// Pair Pions v2 vs. Pair Pions Q2 //////////////////////////
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(14);
    Temp_v2w = pTemp_v2_chargedhadrons->GetBinContent(15);

    // TPC
    if (Temp_v2e != 0)
    {
        p_pionv2e_Q2_obs1->Fill(Q2_pion_TPC, Temp_v2e);
        p_pionv2e_Q2_obs2->Fill(Q2_pion_TPC, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_pionv2w_Q2_obs1->Fill(Q2_pion_TPC, Temp_v2w);
        p_pionv2w_Q2_obs2->Fill(Q2_pion_TPC, Temp_v2w, Eweight);
    }

    // EPD
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(16);
    Temp_v2w = pTemp_v2_chargedhadrons->GetBinContent(17);

    if (Temp_v2e != 0)
    {
        p_pionv2e_Q2_EPD_obs1->Fill(Q2_pion_EPD, Temp_v2e);
        p_pionv2e_Q2_EPD_obs2->Fill(Q2_pion_EPD, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_pionv2w_Q2_EPD_obs1->Fill(Q2_pion_EPD, Temp_v2w);
        p_pionv2w_Q2_EPD_obs2->Fill(Q2_pion_EPD, Temp_v2w, Eweight);
    }

    // EPD1
    Temp_v2e = pTemp_v2_chargedhadrons->GetBinContent(18);

    if (Temp_v2e != 0)
    {
        p_pionv2e_Q2_EPD1_obs1->Fill(Q2_pion_EPD1, Temp_v2e);
        p_pionv2e_Q2_EPD1_obs2->Fill(Q2_pion_EPD1, Temp_v2e, Eweight);
    }
    //////////////////////////////////////////////////////////////////////////////

    ////////////////////////// Parent Gamma112/Gamma132 vs. Parent v2 //////////////////////////
    // TPC
    Temp_v2e = pTemp_v2_parent->GetBinContent(1);
    Temp_v2w = pTemp_v2_parent->GetBinContent(2);
    float Temp_parity3e = pTemp_parity_e->GetBinContent(3);
    float Temp_parity4e = pTemp_parity_e->GetBinContent(4);
    float Temp_parity3w = pTemp_parity_w->GetBinContent(3);
    float Temp_parity4w = pTemp_parity_w->GetBinContent(4);
    float Temp_parity3e_132 = (float)(pTemp_parity_e->GetBinContent(7) + pTemp_parity_e->GetBinContent(11)) / 2.0;
    float Temp_parity4e_132 = (float)(pTemp_parity_e->GetBinContent(8) + pTemp_parity_e->GetBinContent(12)) / 2.0;
    float Temp_parity3w_132 = (float)(pTemp_parity_w->GetBinContent(7) + pTemp_parity_w->GetBinContent(11)) / 2.0;
    float Temp_parity4w_132 = (float)(pTemp_parity_w->GetBinContent(8) + pTemp_parity_w->GetBinContent(12)) / 2.0;

    if ((Temp_v2e != 0) && (Temp_parity3e != 0))
    {
        p_Parity_v2e_parent_obs1->Fill(Temp_v2e, Temp_parity3e);
        p_Parity_v2e_parent_obs2->Fill(Temp_v2e, Temp_parity3e, Eweight);
    }
    if ((Temp_v2e != 0) && (Temp_parity4e != 0))
    {
        p_Parity_v2e_parent_obs3->Fill(Temp_v2e, Temp_parity4e);
        p_Parity_v2e_parent_obs4->Fill(Temp_v2e, Temp_parity4e, Eweight);
    }
    if ((Temp_v2w != 0) && (Temp_parity3w != 0))
    {
        p_Parity_v2w_parent_obs1->Fill(Temp_v2w, Temp_parity3w);
        p_Parity_v2w_parent_obs2->Fill(Temp_v2w, Temp_parity3w, Eweight);
    }
    if ((Temp_v2w != 0) && (Temp_parity4w != 0))
    {
        p_Parity_v2w_parent_obs3->Fill(Temp_v2w, Temp_parity4w);
        p_Parity_v2w_parent_obs4->Fill(Temp_v2w, Temp_parity4w, Eweight);
    }
    if ((Temp_v2e != 0) && (Temp_v2w != 0))
    {
        Hist_parentv2_event->Fill(0.5 * (Temp_v2e + Temp_v2w));
        pv2parent_CosFW->Fill(0.5 * (Temp_v2e + Temp_v2w), cos(nHar * TPC_EP_for_new - nHar * TPC_EP_bac_new), Eweight / eff);
    }

    ////////////////////////// Parent Gamma112/Gamma132 vs. All Charged Hadrons Q2 //////////////////////////
    // tmp4: 0-1 TPC, 2-3 EPD, 4-5 EPD1, pion: 6-7 TPC, 8-9 EPD, 10-11 EPD1
    float tmp_Q2[6] = {Q2_proper, Q2_proper_EPD, Q2_proper_EPD1, Q2_pion_TPC, Q2_pion_EPD, Q2_pion_EPD1};
    // int index_for_data[18] = {3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 36};
    int index_for_data[24] = {3, 4, 7, 8, 15, 16, 19, 20, 27, 28, 31, 32,
                              39, 40, 43, 44, 51, 52, 55, 56, 63, 64, 67, 68};

    for (int tmp4 = 0; tmp4 < 12; tmp4++)
    {
        float for_Temp_parity3e = 0, for_Temp_parity4e = 0, for_Temp_parity3w = 0, for_Temp_parity4w = 0;
        if (tmp4 % 2 == 0)
        {
            for_Temp_parity3e = pTemp_parity_e->GetBinContent(index_for_data[(tmp4 % 6) * 2]);
            for_Temp_parity4e = pTemp_parity_e->GetBinContent(index_for_data[(tmp4 % 6) * 2 + 1]);
            for_Temp_parity3w = pTemp_parity_w->GetBinContent(index_for_data[(tmp4 % 6) * 2]);
            for_Temp_parity4w = pTemp_parity_w->GetBinContent(index_for_data[(tmp4 % 6) * 2 + 1]);
        }
        else
        {
            for_Temp_parity3e = (float)(pTemp_parity_e->GetBinContent(index_for_data[(tmp4 % 6) * 2]) + pTemp_parity_e->GetBinContent(index_for_data[(tmp4 % 6) * 2] + 4)) / 2.0;
            for_Temp_parity4e = (float)(pTemp_parity_e->GetBinContent(index_for_data[(tmp4 % 6) * 2 + 1]) + pTemp_parity_e->GetBinContent(index_for_data[(tmp4 % 6) * 2 + 1] + 4)) / 2.0;
            for_Temp_parity3w = (float)(pTemp_parity_w->GetBinContent(index_for_data[(tmp4 % 6) * 2]) + pTemp_parity_w->GetBinContent(index_for_data[(tmp4 % 6) * 2] + 4)) / 2.0;
            for_Temp_parity4w = (float)(pTemp_parity_w->GetBinContent(index_for_data[(tmp4 % 6) * 2 + 1]) + pTemp_parity_w->GetBinContent(index_for_data[(tmp4 % 6) * 2 + 1] + 4)) / 2.0;
        }

        if (for_Temp_parity3e != 0)
        {
            pParity_e_Q2parent_obs1->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3e);
            pParity_e_Q2parent_obs2->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3e, Eweight);
        }
        if (for_Temp_parity4e != 0)
        {
            pParity_e_Q2parent_obs1->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4e);
            pParity_e_Q2parent_obs2->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4e, Eweight);
        }
        if (for_Temp_parity3w != 0)
        {
            pParity_w_Q2parent_obs1->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3w);
            pParity_w_Q2parent_obs2->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3w, Eweight);
        }
        if (for_Temp_parity4w != 0)
        {
            pParity_w_Q2parent_obs1->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4w);
            pParity_w_Q2parent_obs2->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4w, Eweight);
        }
    }
    //////////////////////////////////////////////////////////////////////////////
    // QQ Cut
    //  Temp_v2e = pTemp_v2_parent->GetBinContent(1);
    //  Temp_v2w = pTemp_v2_parent->GetBinContent(2);
    //  Temp_parity3e = pTemp_parity_e_QQcut->GetBinContent(3);
    //  Temp_parity4e = pTemp_parity_e_QQcut->GetBinContent(4);
    //  Temp_parity3w = pTemp_parity_w_QQcut->GetBinContent(3);
    //  Temp_parity4w = pTemp_parity_w_QQcut->GetBinContent(4);
    //  if((Temp_v2e != 0) && (Temp_parity3e != 0))
    //  {
    //      p_Parity_v2e_parent_QQcut_obs1->Fill(Temp_v2e, Temp_parity3e);
    //      p_Parity_v2e_parent_QQcut_obs2->Fill(Temp_v2e, Temp_parity3e, Eweight);
    //  }
    //  if((Temp_v2e != 0) && (Temp_parity4e != 0))
    //  {
    //      p_Parity_v2e_parent_QQcut_obs3->Fill(Temp_v2e, Temp_parity4e);
    //      p_Parity_v2e_parent_QQcut_obs4->Fill(Temp_v2e, Temp_parity4e, Eweight);
    //  }
    //  if((Temp_v2w != 0) && (Temp_parity3w != 0))
    //  {
    //      p_Parity_v2w_parent_QQcut_obs1->Fill(Temp_v2w, Temp_parity3w);
    //      p_Parity_v2w_parent_QQcut_obs2->Fill(Temp_v2w, Temp_parity3w, Eweight);
    //  }
    //  if((Temp_v2w != 0) && (Temp_parity4w != 0))
    //  {
    //      p_Parity_v2w_parent_QQcut_obs3->Fill(Temp_v2w, Temp_parity4w);
    //      p_Parity_v2w_parent_QQcut_obs4->Fill(Temp_v2w, Temp_parity4w, Eweight);
    //  }
    //  if(Temp_v2e != 0)
    //  {
    //      p_v2parente_Q2parent_QQcut_obs1->Fill(Q2_parent_QQcut, Temp_v2e);
    //      p_v2parente_Q2parent_QQcut_obs2->Fill(Q2_parent_QQcut, Temp_v2e, Eweight);
    //  }
    //  if(Temp_v2w != 0)
    //  {
    //      p_v2parentw_Q2parent_QQcut_obs1->Fill(Q2_parent_QQcut, Temp_v2w);
    //      p_v2parentw_Q2parent_QQcut_obs2->Fill(Q2_parent_QQcut, Temp_v2w, Eweight);
    //  }
    //  if(Temp_parity3e != 0)
    //  {
    //      pParity_e_Q2parent_QQcut_obs1->Fill(3, Q2_parent_QQcut, Temp_parity3e);
    //      pParity_e_Q2parent_QQcut_obs2->Fill(3, Q2_parent_QQcut, Temp_parity3e, Eweight);
    //  }
    //  if(Temp_parity4e != 0)
    //  {
    //      pParity_e_Q2parent_QQcut_obs1->Fill(4, Q2_parent_QQcut, Temp_parity4e);
    //      pParity_e_Q2parent_QQcut_obs2->Fill(4, Q2_parent_QQcut, Temp_parity4e, Eweight);
    //  }
    //  if(Temp_parity3w != 0)
    //  {
    //      pParity_w_Q2parent_QQcut_obs1->Fill(3, Q2_parent_QQcut, Temp_parity3w);
    //      pParity_w_Q2parent_QQcut_obs2->Fill(3, Q2_parent_QQcut, Temp_parity3w, Eweight);
    //  }
    //  if(Temp_parity4w != 0)
    //  {
    //      pParity_w_Q2parent_QQcut_obs1->Fill(4, Q2_parent_QQcut, Temp_parity4w);
    //      pParity_w_Q2parent_QQcut_obs2->Fill(4, Q2_parent_QQcut, Temp_parity4w, Eweight);
    //  }
}

//////////////////////////////////////////////////////
void WrapUpESE_rot()
{
    //////////////////////// TPC Results ////////////////////////
    // v2 vs. Q2
    for (int tmp1 = 1; tmp1 <= 6; tmp1++)
    {
        if (pTemp_v2->GetBinContent(tmp1) != 0)
        {
            p_v2_Q2_obs1_rot[tmp1 - 1]->Fill(Q2_proper, pTemp_v2_rot->GetBinContent(tmp1));
            p_v2_Q2_obs2_rot[tmp1 - 1]->Fill(Q2_proper, pTemp_v2_rot->GetBinContent(tmp1), Eweight);
        }
    }

    // Gamma vs. Q2
    for (int tmp2 = 1; tmp2 <= 8; tmp2++)
    {
        if (pTemp_parity_e_rot->GetBinContent(tmp2) != 0)
        {
            pParity_e_Q2_obs1_rot->Fill(tmp2, Q2_proper, pTemp_parity_e_rot->GetBinContent(tmp2));
            pParity_e_Q2_obs2_rot->Fill(tmp2, Q2_proper, pTemp_parity_e_rot->GetBinContent(tmp2), Eweight);
        }
        if (pTemp_parity_w_rot->GetBinContent(tmp2) != 0)
        {
            pParity_w_Q2_obs1_rot->Fill(tmp2, Q2_proper, pTemp_parity_w_rot->GetBinContent(tmp2));
            pParity_w_Q2_obs2_rot->Fill(tmp2, Q2_proper, pTemp_parity_w_rot->GetBinContent(tmp2), Eweight);
        }
    }

    // Delta vs. Q2
    for (int tmp3 = 1; tmp3 <= 4; tmp3++)
    {
        if (pTemp_delta_rot->GetBinContent(tmp3) != 0)
        {
            pDelta_Q2_obs1_rot->Fill(tmp3, Q2_proper, pTemp_delta_rot->GetBinContent(tmp3));
            pDelta_Q2_obs2_rot->Fill(tmp3, Q2_proper, pTemp_delta_rot->GetBinContent(tmp3), Eweight);
        }
    }

    ////////////////////////// Parent v2 vs. All Charged Hadrons Q2 //////////////////////////
    float Temp_v2e = pTemp_v2_parent_rot->GetBinContent(1);
    float Temp_v2w = pTemp_v2_parent_rot->GetBinContent(2);
    float Temp_v2d_e = pTemp_v2_rot->GetBinContent(1);
    float Temp_v2d_w = pTemp_v2_rot->GetBinContent(2);

    if ((Temp_v2e != 0) && (Temp_v2w != 0) && (Temp_v2d_e != 0) && (Temp_v2d_w != 0))
    {
        Hist_v2_v2parent_3_rot->Fill(0.5 * (Temp_v2e + Temp_v2w), 0.5 * (Temp_v2d_e + Temp_v2d_w));
    }

    // TPC
    if (Temp_v2e != 0)
    {
        p_v2parente_Q2parent_obs1_rot->Fill(Q2_proper, Temp_v2e);
        p_v2parente_Q2parent_obs2_rot->Fill(Q2_proper, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2parentw_Q2parent_obs1_rot->Fill(Q2_proper, Temp_v2w);
        p_v2parentw_Q2parent_obs2_rot->Fill(Q2_proper, Temp_v2w, Eweight);
    }

    // EPD
    Temp_v2e = pTemp_v2_parent_rot->GetBinContent(5);
    Temp_v2w = pTemp_v2_parent_rot->GetBinContent(6);

    if (Temp_v2e != 0)
    {
        p_v2parente_Q2parent_EPD_obs1_rot->Fill(Q2_proper_EPD, Temp_v2e);
        p_v2parente_Q2parent_EPD_obs2_rot->Fill(Q2_proper_EPD, Temp_v2e, Eweight);
    }
    if (Temp_v2w != 0)
    {
        p_v2parentw_Q2parent_EPD_obs1_rot->Fill(Q2_proper_EPD, Temp_v2w);
        p_v2parentw_Q2parent_EPD_obs2_rot->Fill(Q2_proper_EPD, Temp_v2w, Eweight);
    }

    // EPD1
    Temp_v2e = pTemp_v2_parent_rot->GetBinContent(7);

    if (Temp_v2e != 0)
    {
        p_v2parente_Q2parent_EPD1_obs1_rot->Fill(Q2_proper_EPD1, Temp_v2e);
        p_v2parente_Q2parent_EPD1_obs2_rot->Fill(Q2_proper_EPD1, Temp_v2e, Eweight);
    }
    //////////////////////////////////////////////////////////////////////////////

    ////////////////////////// Parent Gamma112/Gamma132 vs. Parent v2 //////////////////////////
    // TPC
    Temp_v2e = pTemp_v2_parent_rot->GetBinContent(1);
    Temp_v2w = pTemp_v2_parent_rot->GetBinContent(2);
    float Temp_parity3e = pTemp_parity_e_rot->GetBinContent(3);
    float Temp_parity4e = pTemp_parity_e_rot->GetBinContent(4);
    float Temp_parity3w = pTemp_parity_w_rot->GetBinContent(3);
    float Temp_parity4w = pTemp_parity_w_rot->GetBinContent(4);
    float Temp_parity3e_132 = (float)(pTemp_parity_e_rot->GetBinContent(7) + pTemp_parity_e_rot->GetBinContent(11)) / 2.0;
    float Temp_parity4e_132 = (float)(pTemp_parity_e_rot->GetBinContent(8) + pTemp_parity_e_rot->GetBinContent(12)) / 2.0;
    float Temp_parity3w_132 = (float)(pTemp_parity_w_rot->GetBinContent(7) + pTemp_parity_w_rot->GetBinContent(11)) / 2.0;
    float Temp_parity4w_132 = (float)(pTemp_parity_w_rot->GetBinContent(8) + pTemp_parity_w_rot->GetBinContent(12)) / 2.0;

    if ((Temp_v2e != 0) && (Temp_parity3e != 0))
    {
        p_Parity_v2e_parent_obs1_rot->Fill(Temp_v2e, Temp_parity3e);
        p_Parity_v2e_parent_obs2_rot->Fill(Temp_v2e, Temp_parity3e, Eweight);
    }
    if ((Temp_v2e != 0) && (Temp_parity4e != 0))
    {
        p_Parity_v2e_parent_obs3_rot->Fill(Temp_v2e, Temp_parity4e);
        p_Parity_v2e_parent_obs4_rot->Fill(Temp_v2e, Temp_parity4e, Eweight);
    }
    if ((Temp_v2w != 0) && (Temp_parity3w != 0))
    {
        p_Parity_v2w_parent_obs1_rot->Fill(Temp_v2w, Temp_parity3w);
        p_Parity_v2w_parent_obs2_rot->Fill(Temp_v2w, Temp_parity3w, Eweight);
    }
    if ((Temp_v2w != 0) && (Temp_parity4w != 0))
    {
        p_Parity_v2w_parent_obs3_rot->Fill(Temp_v2w, Temp_parity4w);
        p_Parity_v2w_parent_obs4_rot->Fill(Temp_v2w, Temp_parity4w, Eweight);
    }
    if ((Temp_v2e != 0) && (Temp_v2w != 0))
    {
        Hist_parentv2_event_rot->Fill(0.5 * (Temp_v2e + Temp_v2w));
        pv2parent_CosFW->Fill(0.5 * (Temp_v2e + Temp_v2w), cos(nHar * TPC_EP_for_new - nHar * TPC_EP_bac_new), Eweight / eff);
    }

    ////////////////////////// Parent Gamma112/Gamma132 vs. All Charged Hadrons Q2 //////////////////////////
    // tmp4: 0-1 TPC, 2-3 EPD, 4-5 EPD1, pion: 6-7 TPC, 8-9 EPD, 10-11 EPD1
    float tmp_Q2[6] = {Q2_proper, Q2_proper_EPD, Q2_proper_EPD1, Q2_pion_TPC, Q2_pion_EPD, Q2_pion_EPD1};
    // int index_for_data[18] = {3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 36};
    int index_for_data[24] = {3, 4, 7, 8, 15, 16, 19, 20, 27, 28, 31, 32,
                              39, 40, 43, 44, 51, 52, 55, 56, 63, 64, 67, 68};

    for (int tmp4 = 0; tmp4 < 12; tmp4++)
    {
        float for_Temp_parity3e = 0, for_Temp_parity4e = 0, for_Temp_parity3w = 0, for_Temp_parity4w = 0;
        if (tmp4 % 2 == 0)
        {
            for_Temp_parity3e = pTemp_parity_e_rot->GetBinContent(index_for_data[tmp4 * 2]);
            for_Temp_parity4e = pTemp_parity_e_rot->GetBinContent(index_for_data[tmp4 * 2 + 1]);
            for_Temp_parity3w = pTemp_parity_w_rot->GetBinContent(index_for_data[tmp4 * 2]);
            for_Temp_parity4w = pTemp_parity_w_rot->GetBinContent(index_for_data[tmp4 * 2 + 1]);
        }
        else
        {
            for_Temp_parity3e = (float)(pTemp_parity_e_rot->GetBinContent(index_for_data[tmp4 * 2]) + pTemp_parity_e_rot->GetBinContent(index_for_data[tmp4 * 2] + 4)) / 2.0;
            for_Temp_parity4e = (float)(pTemp_parity_e_rot->GetBinContent(index_for_data[tmp4 * 2 + 1]) + pTemp_parity_e_rot->GetBinContent(index_for_data[tmp4 * 2 + 1] + 4)) / 2.0;
            for_Temp_parity3w = (float)(pTemp_parity_w_rot->GetBinContent(index_for_data[tmp4 * 2]) + pTemp_parity_w_rot->GetBinContent(index_for_data[tmp4 * 2] + 4)) / 2.0;
            for_Temp_parity4w = (float)(pTemp_parity_w_rot->GetBinContent(index_for_data[tmp4 * 2 + 1]) + pTemp_parity_w_rot->GetBinContent(index_for_data[tmp4 * 2 + 1] + 4)) / 2.0;
        }

        if (for_Temp_parity3e != 0)
        {
            pParity_e_Q2parent_obs1_rot->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3e);
            pParity_e_Q2parent_obs2_rot->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3e, Eweight);
        }
        if (for_Temp_parity4e != 0)
        {
            pParity_e_Q2parent_obs1_rot->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4e);
            pParity_e_Q2parent_obs2_rot->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4e, Eweight);
        }
        if (for_Temp_parity3w != 0)
        {
            pParity_w_Q2parent_obs1_rot->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3w);
            pParity_w_Q2parent_obs2_rot->Fill(index_for_data[tmp4 * 2], tmp_Q2[tmp4 / 2], for_Temp_parity3w, Eweight);
        }
        if (for_Temp_parity4w != 0)
        {
            pParity_w_Q2parent_obs1_rot->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4w);
            pParity_w_Q2parent_obs2_rot->Fill(index_for_data[tmp4 * 2 + 1], tmp_Q2[tmp4 / 2], for_Temp_parity4w, Eweight);
        }
    }
    //////////////////////////////////////////////////////////////////////////////
}
