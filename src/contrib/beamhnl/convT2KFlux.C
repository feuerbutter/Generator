#include "convT2KFlux.h"

using namespace convT2K;

void convT2KFlux(std::string T2Kflux_filename){
    // pass as argument "/cvmfs/t2k.egi.eu/flux/nu.13a_nom_ND6_250ka_flukain.all.root"
    TFile * fin = new TFile(T2Kflux_filename.c_str(),"READ");
    TTree* t_t2k = (TTree*) fin->Get("h3002");

    t_t2k->SetBranchAddress("xpi",      &xpi);
    t_t2k->SetBranchAddress("npi",      &npi);
    t_t2k->SetBranchAddress("ppi",      &ppi);
    t_t2k->SetBranchAddress("ppid",     &ppid);
    t_t2k->SetBranchAddress("Enu",      &Enu);
    t_t2k->SetBranchAddress("nnu",      &nnu);
    
    int n_evt = t_t2k->GetEntries();
    std::cout << "This file has " << n_evt << " entries." << std::endl;

    //------------variables-------------

    int n_conv = 100; // # of entries to convert
    double pots = (double) n_conv;

    int job    = 1;
    // global coordinates of parent meson decay vertex (i.e. v emission) [cm]
    double decay_vx,decay_vy,decay_vz;

    // momentum [GeV] of parent at its decay vertex in global coordinates
    double decay_pdpx,decay_pdpy,decay_pdpz;

    int decay_ptype; //decay_pdpx

    // SM v energy in CM frame. Needed to calculate acceptance correction
    double decay_necm; 

    double decay_nimpwt = 1.0; // stat multiplier as parents with similar kinematics not simulated in dk2nu -- can leave 1 for now

    //-- open the output ROOT file
    // TString fout_name = Form("genie_lt160momNmu1pi/genie_lt160momNmu1pi_%d.root",file_ind);
    TString fout_name = "dk2nu_t2k.root";
    TFile* fout = new TFile(fout_name,"RECREATE");

    //-- create the output tree
    TTree* dkTree = new TTree("dkTree", "dkTree");
    TTree* dkMeta = new TTree("dkMeta", "dkMeta");
    dkMeta->Branch( "job",              &job,                  "job/I"                             );
    dkMeta->Branch( "pots",              &pots,                  "pots/D"                             );

    dkTree->Branch( "job",              &job,                  "job/I"                             );
    dkTree->Branch( "potnum",           &pots,               "potnum/D"                          );
    dkTree->Branch( "decay_ptype",      &ppid,        "decay_ptype/I"                     );
    dkTree->Branch( "decay_vx",         &decay_vx,           "decay_vx/D"                        );
    dkTree->Branch( "decay_vy",         &decay_vy,           "decay_vy/D"                        );
    dkTree->Branch( "decay_vz",         &decay_vz,           "decay_vz/D"                        );
    dkTree->Branch( "decay_pdpx",       &decay_pdpx,         "decay_pdpx/D"                      );
    dkTree->Branch( "decay_pdpy",       &decay_pdpy,         "decay_pdpy/D"                      );
    dkTree->Branch( "decay_pdpz",       &decay_pdpz,         "decay_pdpz/D"                      );
    dkTree->Branch( "decay_necm",       &decay_necm,         "decay_necm/D"                      );
    dkTree->Branch( "decay_nimpwt",     &decay_nimpwt,       "decay_nimpwt/D"                    );



    for (int edx=0;edx<n_conv;edx++){
        t_t2k->GetEntry(edx);
        std::cout << "=======Checking event " << edx << "====================" << std::endl;
        decay_vx = (double) xpi[0];
        decay_vy = (double) xpi[1];
        decay_vz = (double) xpi[2];
        printf("Decay vertex at : (%6.3f,%6.3f,%6.3f)\n",decay_vx,decay_vy,decay_vz);

        decay_pdpx = (double) ppi * npi[0];
        decay_pdpy = (double) ppi * npi[1];
        decay_pdpz = (double) ppi * npi[2];

        double m_par = convT2K::getMassPar(ppid);
        double E_par = std::sqrt((decay_pdpx*decay_pdpx) + (decay_pdpy*decay_pdpy)+(decay_pdpz*decay_pdpz) +m_par*m_par );

        decay_necm = boostNu(Enu*nnu[0],Enu*nnu[1],Enu*nnu[2],Enu,decay_pdpx,decay_pdpy,decay_pdpz,E_par);

        printf("Parent id is %d, with mass %6.3f GeV, and energy %6.3f GeV \n",ppid,m_par,E_par);
        printf("parent momentum is : (%6.3f,%6.3f,%6.3f)\n",decay_pdpx,decay_pdpy,decay_pdpz);
        printf("neutrino momentum is : (%6.3f,%6.3f,%6.3f,%6.3f)\n",Enu*nnu[0],Enu*nnu[1],Enu*nnu[2],Enu);

        dkTree->Fill();
        dkMeta->Fill();

    }

    fin->Close();
    fout->Write();
    fout->Close();
}

double convT2K::boostNu(double px, double py, double pz, double E, double px_p, double py_p, double pz_p, double E_p) {
    // Create the four-momentum vectors for the neutrino and parent meson
    TLorentzVector neutrino(px, py, pz, E);
    TLorentzVector parent(px_p, py_p, pz_p, E_p);

    // Boost the neutrino back to the rest frame of the parent meson
    TLorentzVector boosted_neutrino = neutrino;
    boosted_neutrino.Boost(-parent.BoostVector());
    
    // you can check the values of the four-momentum in the rest frame of the parent meson
    cout << "Boosted neutrino four-momentum: " << boosted_neutrino.X() << " " << boosted_neutrino.Y() << " " << boosted_neutrino.Z() << " " << boosted_neutrino.T() << endl;

    return boosted_neutrino.T();
}

