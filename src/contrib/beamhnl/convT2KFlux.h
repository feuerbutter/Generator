#ifndef convT2kFlux_h
#define convT2kFlux_h

#include <map>
#include <cstdio>
#include <iostream>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

namespace convT2K{

    // Utility functions
    double boostNu(double px, double py, double pz, double E, double px_p, double py_p, double pz_p, double E_p);

    double getMassPar(int pid);

    void calCMKin(double m_parent, double m_dg1, double m_dg2);

    std::vector<std::pair< Double_t, Double_t> > calOriGeomWeight(const TVector3& P3d_parent, const TVector3& pos_decay, const TVector3& pos_det);


    // const variables
    const int maxArray = 30;
    const int maxC     = 100;
    const float nd_x   = -322.2292;
    const float nd_y   = -814.557;
    const float nd_z   = 28010.0;

    const double m_mu  = 0.105658;     //muon mass in GeV
    float m_dg[2] = {0.105658, 0.0};            // the mass of the two decay products 
    // m_dg[0] = m_mu;
    // m_dg[1] = 0.0;

    // < nd280HNLSim.DetConfig.xPosND = -322.2292 >  Deault X-coordinate of ND center in cm
    // < nd280HNLSim.DetConfig.yPosND = -814.557  >  Deault Y-coordinate of ND center in cm
    // < nd280HNLSim.DetConfig.zPosND = 28010.0   >  Deault Z-coordinate of ND center in cm

    // t2k-flux input variables
    float       xpi[3];
    float       npi[3];          // direction of parent
    float       ppi;             // momentum of parent
    int         mode;            // decay mode of the parent hadron
    int         ppid;            // PID of parent based on Geant3
    int         ppid_pdg;        // Parent PDG code
    float       Enu;             // neutrino energy in lab frame
    float       nnu[3];          // direction of neutrino 
    float       xnu;             // x-position of nu hitting the dector
    float       ynu;             // y-position of nu hitting the dector
    float       t2k_norm;        // the combined weight input from the t2k flux

    // calculated variables
    float       m_parent;        // mass of parent particle
    float       beta_dg_cm[2];      // the speed for the two decay daugthers in the cm frame
    float       E_dg_cm[2];         // the energy for the two decay daugthers in the cm frame
    float       p_dg_cm;            // momentum of the decay products in the cm

    std::map< int, int > g32pdg = { { 1, 22 }, { 2, -11 }, { 3, 11 }, { 4, 12 }, { 5, -13 }, { 6, 13 }, { 7, 111 }, { 8, 211 }, { 9, -211 }, { 10, 130 }, { 11, 321 }, { 12, -321 }, { 13, 2112 }, { 14, 2212 }, { 15, -2212 }, { 16, 130 } }; // there are more but w/e

    // int         dArSize           = 0;         // Size of location arrays
    // int         dAnArSize         = 0;         // Size of ancestor arrays
    // int         dTrArSize         = 0;         // Size of traj     arrays
    // // - - -
    // int         dJob              = 0;         // Simulation job number
    // double      dPotnum           = 0.0;       // N POT for this v (0 job-beginning, total POT job-end)
    // double      dPpvx             = 0.0;       // x component of v parent production vertex in cm
    // double      dPpvy             = 0.0;       // y component of v parent production vertex in cm
    // double      dPpvz             = 0.0;       // z component of v parent production vertex in cm
    // // - - -
    // int         dDecayDotNorig    = 0;         // v origin (1,2,3 = tgt/baffle, muon decay, other)
    // int         dDecayDotNdecay   = 0;         // Decay code of decay that produced v
    // int         dDecayDotNtype    = 0;         // GEANT particle code of v
    // double      dDecayDotVx       = 0.0;       // x component of v vertex position in cm
    // double      dDecayDotVy       = 0.0;       // y component of v vertex position in cm
    // double      dDecayDotVz       = 0.0;       // z component of v vertex position in cm
    // double      dDecayDotPdpx     = 0.0;       // x component of final parent momentum in GeV
    // double      dDecayDotPdpy     = 0.0;       // y component of final parent momentum in GeV
    // double      dDecayDotPdpz     = 0.0;       // z component of final parent momentum in GeV
    // double      dDecayDotPpdxdz   = 0.0;       // px/pz of parent at parent production point
    // double      dDecayDotPpdydz   = 0.0;       // py/pz of parent at parent production point
    // double      dDecayDotPppz     = 0.0;       //    pz of parent at parent production point in GeV
    // double      dDecayDotPpenergy = 0.0;       //     E of parent at parent production point in GeV
    // int         dDecayDotPpmedium = 0;         // empty branch
    // int         dDecayDotPtype    = 0;         // GEANT particle code of parent
    // double      dDecayDotMuparpx  = 0.0;       // (parent == mu) ? grandparent px in GeV : -99999
    // double      dDecayDotMuparpy  = 0.0;       // (parent == mu) ? grandparent py in GeV : -99999
    // double      dDecayDotMuparpz  = 0.0;       // (parent == mu) ? grandparent pz in GeV : -99999
    // double      dDecayDotMupare   = 0.0;       // (parent == mu) ? grandparent  E in GeV : -99999
    // double      dDecayDotNecm     = 0.0;       // v E in parent rest frame in GeV
    // double      dDecayDotNimpwt   = 0.0;       // Importance weight
    // // - - -
    // double      dNurayDotPx[maxArray];         // v px in GeV for each location in meta
    // double      dNurayDotPy[maxArray];         // v py in GeV for each location in meta
    // double      dNurayDotPz[maxArray];         // v pz in GeV for each location in meta
    // double      dNurayDotE[maxArray];          // v  E in GeV for each location in meta
    // double      dNurayDotWgt[maxArray];        // weights to make v flux spectra for each location in meta
    // // - - -
    // int         dAncestorDotPdg[maxArray];     // PDG code of ancestor
    // double      dAncestorDotStartx[maxArray];  // x component of ancestor start position in cm
    // double      dAncestorDotStarty[maxArray];  // y component of ancestor start position in cm
    // double      dAncestorDotStartz[maxArray];  // z component of ancestor start position in cm
    // double      dAncestorDotStartt[maxArray];  // t component of ancestor start position in (ns ?)
    // double      dAncestorDotStartpx[maxArray]; // ancestor initial px in GeV
    // double      dAncestorDotStartpy[maxArray]; // ancestor initial py in GeV
    // double      dAncestorDotStartpz[maxArray]; // ancestor initial pz in GeV
    // double      dAncestorDotStoppx[maxArray];  // ancestor final   px in GeV 
    // double      dAncestorDotStoppy[maxArray];  // ancestor final   py in GeV
    // double      dAncestorDotStoppz[maxArray];  // ancestor final   pz in GeV
    // double      dAncestorDotPolx[maxArray];    // empty branch
    // double      dAncestorDotPoly[maxArray];    // empty branch
    // double      dAncestorDotPolz[maxArray];    // empty branch
    // double      dAncestorDotPprodpx[maxArray]; // parent px prior to secondaries production (meaning?)
    // double      dAncestorDotPprodpy[maxArray]; // parent py prior to secondaries production (meaning?)
    // double      dAncestorDotPprodpz[maxArray]; // parent pz prior to secondaries production (meaning?)
    // int         dAncestorDotNucleus[maxArray]; // PDG code of nucleus where the interaction happened
    // char   dAncestorDotProc[maxArray*maxC];    // Describes processes that created each ancestor
    // char   dAncestorDotIvol[maxArray*maxC];    // Describes volume   where each ancestor was created
    // char   dAncestorDotImat[maxArray*maxC];    // Describes material where each ancestor was created
    // // - - -
    // double      dTgtexitDotTvx    = 0.0;       // x position of parent target exit in cm
    // double      dTgtexitDotTvy    = 0.0;       // y position of parent target exit in cm
    // double      dTgtexitDotTvz    = 0.0;       // z position of parent target exit in cm
    // double      dTgtexitDotTpx    = 0.0;       // parent px at target exit in GeV
    // double      dTgtexitDotTpy    = 0.0;       // parent py at target exit in GeV
    // double      dTgtexitDotTpz    = 0.0;       // parent pz at target exit in GeV
    // int         dTgtexitDotTptype = 0;         // GEANT particle code of ancestor that exited target
    // int         dTgtexitDotTgen   = 0;         // Generation number of v
    // // - - -
    // double      dTrajDotTrkx[maxArray];        // ?
    // double      dTrajDotTrky[maxArray];        // ?
    // double      dTrajDotTrkz[maxArray];        // ?
    // double      dTrajDotTrkpx[maxArray];       // ?
    // double      dTrajDotTrkpy[maxArray];       // ?
    // double      dTrajDotTrkpz[maxArray];       // ?

    // define a map to store the mapping of PID to mass
    std::map<int, double> pid_mass_map = {
        {1, 0.0},
        {2, 0.000511},
        {3, 0.000511},
        {4, 0.0},
        {5, 0.1057},
        {6, 0.1057},
        {7, 0.135},
        {8, 0.1396},
        {9, 0.1396},
        {10, 0.4977},
        {11, 0.4937},
        {12, 0.4937},
        {13, 0.9396},
        {14, 0.9383},
        {15, 0.9383},
        {16, 0.4977},
        {17, 0.5475},
        {18, 1.116},
        {19, 1.189},
        {20, 1.193},
        {21, 1.197},
        {22, 1.315},
        {23, 1.321},
        {24, 1.672},
        {25, 0.9396},
        {26, 1.116},
        {27, 1.189},
        {28, 1.193},
        {29, 1.197},
        {30, 1.315},
        {31, 1.321},
        {32, 1.672},
        {45, 1.876},
        {46, 2.809},
        {47, 3.727},
        {48, 0.0},
        {49, 2.809},
        {50, 0.0}
    };

    double getMassPar(int pid) {
        // check if the pid is in the map
        if (pid_mass_map.count(pid) == 0) {
            // if pid is not in the map, return -1
            return -1.0;
        } else {
            // if pid is in the map, return the corresponding mass
            return pid_mass_map[pid];
        }
    }

    
    void calCMKin(double m_parent, double m_dg1, double m_dg2){
        E_dg_cm[0] =(m_parent*m_parent - m_dg2*m_dg2 + m_dg1*m_dg1)/(2.0*m_parent);
        E_dg_cm[1] = m_parent - E_dg_cm[0]; 
        p_dg_cm = sqrt(E_dg_cm[1]*E_dg_cm[1] - m_dg2*m_dg2);
	
	    beta_dg_cm[0] = beta_dg_cm[1] = 0.;
	
	    if(E_dg_cm[0]!=0)
	      beta_dg_cm[0] = p_dg_cm/E_dg_cm[0];
	
	    if(E_dg_cm[1]!=0)
	      beta_dg_cm[1] = p_dg_cm/E_dg_cm[1];
    }
        
}


#endif
