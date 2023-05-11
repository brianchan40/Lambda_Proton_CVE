#ifndef particle_all_hh
#define particle_all_hh

class particle_all
{
public:
    particle_all() {
        px = -3.14999;
        py = -3.14999;
        pz = -3.14999;
        Charge = -3.14999;
        dcaglobal = -3.14999;
        nSigmaProton = -3.14999;
        isTofTrack = false;
        trk_id = -999;
        is_pion = false;
    }
    particle_all(float in_px, float in_py, float in_pz, float in_charge, float in_dcaglobal, float in_nSigmaProton, bool in_isTofTrack, int in_trk_id, bool in_is_pion);

    virtual       ~particle_all() { }

    float px, py, pz;
    float Charge;
    float dcaglobal;
    float nSigmaProton;
    bool isTofTrack;
    int trk_id;
    bool is_pion;
};



#endif
