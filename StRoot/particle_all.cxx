#include "particle_all.h"

particle_all::particle_all(float in_px, float in_py, float in_pz, float in_charge, float in_dcaglobal, float in_nSigmaProton, bool in_isTofTrack, int in_trk_id, bool in_is_pion)
{
    px = in_px;
    py = in_py;
    pz = in_pz;
    Charge = in_charge;
    dcaglobal = in_dcaglobal;
    nSigmaProton = in_nSigmaProton;
    isTofTrack = in_isTofTrack;
    trk_id = in_trk_id;
    is_pion = in_is_pion;
}
