#ifndef event_hh
#define event_hh

class event
{
public:
    event()
    {
        cent = -999;
        num_trk = -999;
        nLambda = -999;
        nLambdaRot = -999;
        Run = -999;
        TOFMult = -999;
        RefMult = -999;
        n_Proton = -999;
        n_Pion = -999;
        VPDvz = -3.14999;
        PVtxz = -3.14999;
        PVtxx = -3.14999;
        PVtxy = -3.14999;
        Eweight = -3.14999;
        EventID = -999;

        Magn = -3.14999;

        EPD_EP1_east = -3.14999;
        EPD_EP1_west = -3.14999;
        EPD_EP_east = -3.14999;
        EPD_EP_west = -3.14999;
    }
    event(int in_cent, int in_num_trk, int in_nLambda, int in_nLambdaRot, int in_Run, int in_TOFMult, int in_RefMult, int in_n_Proton, int in_n_Pion, float in_VPDvz, float in_PVtxz, float in_PVtxx, float in_PVtxy, float in_Eweight, float in_Magn, int in_EventID, float in_EPD_EP1_east, float in_EPD_EP1_west, float in_EPD_EP_east, float in_EPD_EP_west);

    virtual       ~event() { }

    int cent, num_trk, nLambda, nLambdaRot, Run, TOFMult, RefMult, n_Proton, n_Pion, EventID;
    float VPDvz;
    float PVtxz, PVtxx, PVtxy, Eweight, EPD_EP1_east, EPD_EP1_west, EPD_EP_east, EPD_EP_west;
    
    float Magn;
};



#endif
