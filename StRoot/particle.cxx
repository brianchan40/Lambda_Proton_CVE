#include "particle.h"

particle::particle(float in_px, float in_py, float in_pz, float in_charge, float in_dcaglobal, float in_mass, int in_trk_id, float in_nsigma, float in_hits_ratio, float in_nhitsfit, float in_nhitsmax)
{
    px = in_px;
    py = in_py;
    pz = in_pz;
//    x = in_x;
//    y = in_y;
//    z = in_z;
    Charge = in_charge;
//    TOFflag = in_TOFflag;
    dcaglobal = in_dcaglobal;
//    prim = in_prim;
//    isTofTrack = in_isTofTrack;
    mass = in_mass;
    trk_id = in_trk_id;
//    ToF = in_ToF;
//    BTOFYLocal = in_BTOFYLocal;
//    Mass2 = in_Mass2;
    nsigma = in_nsigma;
    hits_ratio = in_hits_ratio;
    nhitsfit = in_nhitsfit;
    nhitsmax = in_nhitsmax;
}
