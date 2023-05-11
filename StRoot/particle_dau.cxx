#include "particle_dau.h"

particle_dau::particle_dau(double in_px, double in_py, double in_pz, int in_trk_id, float in_nHitsFit, float in_nSigma, float in_dcaglobal, float in_nHitsMax)
{
    px = in_px;
    py = in_py;
    pz = in_pz;
//    x = in_x;
//    y = in_y;
//    z = in_z;
//    Charge = in_charge;
//    TOFflag = in_TOFflag;
    dcaglobal = in_dcaglobal;
//    prim = in_prim;
    nSigma = in_nSigma;
//    isTofTrack = in_isTofTrack;
//    mass = in_mass;
    trk_id = in_trk_id;
//    ToF = in_ToF;
//    BTOFYLocal = in_BTOFYLocal;
//    Mass2 = in_Mass2;
    nHitsFit = in_nHitsFit;
    nHitsMax = in_nHitsMax;
}
