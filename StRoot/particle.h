#ifndef particle_hh
#define particle_hh

class particle
{
public:
    particle() {
        px = -3.14999;
        py = -3.14999;
        pz = -3.14999;
//        x = -3.14999;
//        y = -3.14999;
//        z = -3.14999;
        Charge = -3.14999;
//        TOFflag = -999;
        dcaglobal = -3.14999;
        nsigma = -3.14999;
//        prim = -999;
//        isTofTrack = -999;
        mass = -3.14999;
        trk_id = -999;
//        ToF = -3.14999;
//        BTOFYLocal = -3.14999;
//        Mass2 = -3.14999;
        hits_ratio = -3.149999;
        nhitsfit = -3.1499999;
        nhitsmax = -3.1499999;

    }
    particle(float in_px, float in_py, float in_pz, float in_charge, float in_dcaglobal, float in_mass, int in_trk_id, float in_nsigma, float in_hits_ratio, float in_nhitsfit, float in_nhitsmax);

    virtual       ~particle() { }

    float px, py, pz;
    float Charge;
    float dcaglobal;
    float nsigma;
    float mass;
    int trk_id;
    float hits_ratio;
    float nhitsfit;
    float nhitsmax;
};



#endif
