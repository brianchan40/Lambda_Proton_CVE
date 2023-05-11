#include "event.h"

event::event(int in_cent, int in_num_trk, int in_nLambda, int in_nLambdaRot, int in_Run, int in_TOFMult, int in_RefMult, int in_n_Proton, int in_n_Pion, float in_VPDvz, float in_PVtxz, float in_PVtxx, float in_PVtxy, float in_Eweight, float in_Magn, int in_EventID, float in_EPD_EP1_east, float in_EPD_EP1_west, float in_EPD_EP_east, float in_EPD_EP_west)
{
    cent = in_cent;
    num_trk = in_num_trk;
    nLambda = in_nLambda;
    nLambdaRot = in_nLambdaRot;
    Run = in_Run;
    TOFMult = in_TOFMult;
    RefMult = in_RefMult;
    n_Proton = in_n_Proton;
    n_Pion = in_n_Pion;
    VPDvz = in_VPDvz;
    PVtxz = in_PVtxz;
    PVtxx = in_PVtxx;
    PVtxy = in_PVtxy;
    Eweight = in_Eweight;
    
    Magn = in_Magn;
    
    EventID = in_EventID;

    EPD_EP1_east = in_EPD_EP1_east;
    EPD_EP1_west = in_EPD_EP1_west;
    EPD_EP_east = in_EPD_EP_east;
    EPD_EP_west = in_EPD_EP_west;
}
