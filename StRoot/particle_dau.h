#ifndef particle_dau_hh
#define particle_dau_hh

class particle_dau
{
public:
    particle_dau() {
        px = -3.14999;
        py = -3.14999;
        pz = -3.14999;
//        x = -3.14999;
//        y = -3.14999;
//        z = -3.14999;
//        Charge = -3.14999;
//        TOFflag = -999;
        dcaglobal = -3.14999;
//        prim = -999;
        nSigma = -3.14999;
//        isTofTrack = -999;
//        mass = -3.14999;
        trk_id = -999;
//        ToF = -3.14999;
//        BTOFYLocal = -3.14999;
//        Mass2 = -3.14999;
        nHitsFit = -999;
        nHitsMax = -999;
    }
    particle_dau(double in_px, double in_py, double in_pz, int in_trk_id, float in_nHitsFit, float in_nSigma, float in_dcaglobal, float in_nHitsMax);

    virtual       ~particle_dau() { }

    float px, py, pz;
//    float x, y, z;
//    float Charge;
//    Int_t TOFflag;
    float dcaglobal;
//    Int_t prim;
    float nSigma;
//    Int_t isTofTrack;
//    double mass;
    int trk_id;
//    double ToF;
//    double BTOFYLocal;
//    float Mass2;
    float nHitsFit;
    float nHitsMax;
};



#endif
