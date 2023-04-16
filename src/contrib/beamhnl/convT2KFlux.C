#include "convT2KFlux.h"

using namespace convT2K;

void convT2KFlux(std::string T2Kflux_filename){
    // pass as argument "/cvmfs/t2k.egi.eu/flux/nu.13a_nom_ND6_250ka_flukain.all.root"
    TFile * fin = new TFile(T2Kflux_filename.c_str(),"READ");
    TTree* t_t2k = (TTree*) fin->Get("h3002");

    t_t2k->SetBranchAddress("xpi",      &xpi);
    t_t2k->SetBranchAddress("npi",      &npi);
    t_t2k->SetBranchAddress("ppi",      &ppi);
    t_t2k->SetBranchAddress("mode",     &mode);
    t_t2k->SetBranchAddress("ppid",     &ppid);
    t_t2k->SetBranchAddress("Enu",      &Enu);
    t_t2k->SetBranchAddress("nnu",      &nnu);
    t_t2k->SetBranchAddress("xnu",      &xnu);
    t_t2k->SetBranchAddress("ynu",      &ynu);
    t_t2k->SetBranchAddress("norm",     &t2k_norm);
    
    int n_evt = t_t2k->GetEntries();
    std::cout << "This file has " << n_evt << " entries." << std::endl;

    //------------variables-------------

    int n_conv = 1000; // # of entries to convert
    // int n_conv = n_evt; // # of entries to convert
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

    double t2k_geom_weight;
    double t2k_ungeom_norm;

    //-- open the output ROOT file
    // TString fout_name = Form("genie_lt160momNmu1pi/genie_lt160momNmu1pi_%d.root",file_ind);
    TString fout_name = "/home/weijunli/Apps/GENIE/flux/conv_cm/dk2nu_t2k_konly.root";
    // TString fout_name = "test.root";
    TFile* fout = new TFile(fout_name,"RECREATE");

    //-- create the output tree
    TTree* dkTree = new TTree("dkTree", "dkTree");
    TTree* dkMeta = new TTree("dkMeta", "dkMeta");
    dkMeta->Branch( "job",              &job,                  "job/I"                             );
    dkMeta->Branch( "pots",              &pots,                  "pots/D"                             );

    dkTree->Branch( "potnum",           &pots,               "potnum/D"                          );
    dkTree->Branch( "decay_ptype",      &ppid_pdg,        "decay_ptype/I"                     );
    dkTree->Branch( "decay_vx",         &decay_vx,           "decay_vx/D"                        );
    dkTree->Branch( "decay_vy",         &decay_vy,           "decay_vy/D"                        );
    dkTree->Branch( "decay_vz",         &decay_vz,           "decay_vz/D"                        );
    dkTree->Branch( "decay_pdpx",       &decay_pdpx,         "decay_pdpx/D"                      );
    dkTree->Branch( "decay_pdpy",       &decay_pdpy,         "decay_pdpy/D"                      );
    dkTree->Branch( "decay_pdpz",       &decay_pdpz,         "decay_pdpz/D"                      );
    dkTree->Branch( "decay_necm",       &decay_necm,         "decay_necm/D"                      );
    dkTree->Branch( "decay_nimpwt",     &decay_nimpwt,       "decay_nimpwt/D"                    );


    dkMeta->Branch( "t2k_geom_weight",  &t2k_geom_weight,    "t2k_geom_weight/D"                 );
    dkMeta->Branch( "t2k_ungeom_norm",  &t2k_ungeom_norm,    "t2k_ungeom_norm/D"                 );

    for (int edx=0;edx<n_conv;edx++){
      t_t2k->GetEntry(edx);

      if (! (mode==12)) continue;
	    if( edx % 100 == 0 ) std::cout << "=======Checking event " << edx << "====================" << std::endl;
      decay_vx = (double) xpi[0]/10.;
      decay_vy = (double) xpi[1]/10.;
      decay_vz = (double) xpi[2]/10.;
      printf("Decay vertex at : (%6.3f,%6.3f,%6.3f) [cm]\n",decay_vx,decay_vy,decay_vz);

      decay_pdpx = (double) ppi * npi[0];
      decay_pdpy = (double) ppi * npi[1];
      decay_pdpz = (double) ppi * npi[2];

  	  ppid_pdg = (*convT2K::g32pdg.find(ppid)).second;

      double m_par = convT2K::getMassPar(ppid);
      double E_par = std::sqrt((decay_pdpx*decay_pdpx) + (decay_pdpy*decay_pdpy)+(decay_pdpz*decay_pdpz) +m_par*m_par );

      decay_necm = boostNu(Enu*nnu[0],Enu*nnu[1],Enu*nnu[2],Enu,decay_pdpx,decay_pdpy,decay_pdpz,E_par);

      //printf("Parent id is %d, with mass %6.3f GeV, and energy %6.3f GeV \n",ppid,m_par,E_par);
      //printf("parent momentum is : (%6.3f,%6.3f,%6.3f)\n",decay_pdpx,decay_pdpy,decay_pdpz);
      //printf("neutrino momentum is : (%6.3f,%6.3f,%6.3f,%6.3f)\n",Enu*nnu[0],Enu*nnu[1],Enu*nnu[2],Enu);

      //calculation of the t2k geom weight

      //original position of the neutrino point in global coordinate 


      m_parent = m_par;
      TVector3 p_par(decay_pdpx,decay_pdpy,decay_pdpz);
      TVector3 pos_det(nd_x+xnu, nd_y+ynu, nd_z);
      TVector3 pos_decay(decay_vx,decay_vy,decay_vz);
      // posDet.SetXYZ(nd_x+xnu, nd_y+ynu, nd_z);
      convT2K::calCMKin(m_parent,m_dg[0],m_dg[1]);

      // std::cout << "p_par from t2kinput " << ppi << " p_par from vector " << p_par.Mag() << std::endl; // this shows npi is a unit vector

      std::vector<std::pair<Double_t, Double_t> > weight_energy1;
      weight_energy1 = convT2K::calOriGeomWeight(p_par, pos_decay, pos_det);

      Double_t weight1  = weight_energy1[0].first; //exactly one solution for mass-less neutrino
      Double_t wcal_Enu = weight_energy1[0].second; 

      // printf("direct Enu is %6.3f, boosted_Enu is %6.3f, Enu from weight cal %6.3f\n ", Enu, decay_necm, wcal_Enu);

      t2k_geom_weight = weight1;
      t2k_ungeom_norm = t2k_norm / t2k_geom_weight;

      // 5.5653498e+20 obtained manually by getting a first converted file
      decay_nimpwt = t2k_ungeom_norm / 5.5653498e+20;
      
      dkTree->Fill();
      dkMeta->Fill();

    }

    fin->Close();
    fout->Write();
    fout->Close();


    // TFile* fout_up = new TFile(fout_name,"update");
    // TTree* dk = (TTree*) fout_up->Get("dkTree");
    // TTree* dm = (TTree*) fout_up->Get("dkMeta");
    // dk->SetBranchStatus("*", 0);
    // // dk->SetBranchAddress("decay_nimpwt",      &decay_nimpwt);
    // TBranch* br_new = dk->Branch( "decay_nimpwt",     &decay_nimpwt,       "decay_nimpwt/D"                    );
    // // dk->SetBranchStatus("decay_nimpwt",1);
    // dm->SetBranchAddress("t2k_ungeom_norm",      &t2k_ungeom_norm);

    // Double_t min_norm = dm->GetMinimum("t2k_ungeom_norm");

    // for (Long64_t iEntry = 0; iEntry < dk->GetEntries(); ++iEntry) {
    //   // Get the entry
    //   dk->GetEntry(iEntry);
    //   dm->GetEntry(iEntry);

    //   // cout << "weight bf update" << decay_nimpwt << endl;
    //   decay_nimpwt = t2k_ungeom_norm / min_norm;

    //   // cout << "weight af update" << decay_nimpwt << endl;
    //   br_new->Fill();
    // }
    // dk->Write("",TObject::kOverwrite);
    // fout_up->Close();



}

double convT2K::boostNu(double px, double py, double pz, double E, double px_p, double py_p, double pz_p, double E_p) {
    // Create the four-momentum vectors for the neutrino and parent meson
    TLorentzVector neutrino(px, py, pz, E);
    TLorentzVector parent(px_p, py_p, pz_p, E_p);

    // Boost the neutrino back to the rest frame of the parent meson
    TLorentzVector boosted_neutrino = neutrino;
    boosted_neutrino.Boost(-parent.BoostVector());
    
    // you can check the values of the four-momentum in the rest frame of the parent meson
    // cout << "Boosted neutrino four-momentum: " << boosted_neutrino.X() << " " << boosted_neutrino.Y() << " " << boosted_neutrino.Z() << " " << boosted_neutrino.T() << endl;

    return boosted_neutrino.T();
}

 std::vector<std::pair< Double_t, Double_t> > convT2K::calOriGeomWeight(const TVector3& P3d_parent, const TVector3& pos_decay, const TVector3& pos_det){
	
	std::vector<std::pair< Double_t, Double_t> >  weight_energy;
	
	Double_t E_prod2 = 0.;
	
	//calculate angle and distance to the detector point
	TVector3 ddir(pos_det - pos_decay);
	
	Double_t L  = ddir.Mag();
	ddir        = ddir.Unit();
	
	//parent energy and momentum
	Double_t Plab_parent        = P3d_parent.Mag();
	Double_t Elab_parent        = sqrt(Plab_parent*Plab_parent + m_parent * m_parent);
	Double_t beta_parent        = Plab_parent/Elab_parent;
	
	
	Double_t gamma_parent = Elab_parent/m_parent;
	
	TVector3 dir(P3d_parent.Unit());
	
	// theta is the angle between parent direction and  to detector (forced here for second product) 
	Double_t theta_lab = dir.Angle(ddir); 
	
	Double_t cos_theta_lab = cos(theta_lab); //cos of theta in lab
	
	// calculate energy in lab frame
	Double_t tmp    = gamma_parent*gamma_parent*beta_parent*beta_parent*cos_theta_lab*cos_theta_lab;
	Double_t A      = gamma_parent*gamma_parent - tmp;
	Double_t B      = gamma_parent*E_dg_cm[1];
	Double_t C      = E_dg_cm[1]*E_dg_cm[1] + tmp*m_dg[1]*m_dg[1];
	
	
	// check that the solution is possible 
	if(B*B-A*C < 0.){
	  weight_energy.push_back(std::make_pair( 0., 0.));
	  return weight_energy;
	}
	
	// two solutions hold so have to be careful here 
	// need to consider two main cases:
	// v_cm > v_parent_lab --> take the largest solution since the smallest corresponds to negative cos_theta_lab 
	// v_cm < v_parent_lab --> have to deal with both, each gives its own probability value
	
	std::vector<int> bin_tmp;
	bin_tmp.push_back(1);  
	
	
	if(beta_dg_cm[1]<beta_parent)
	  bin_tmp.push_back(-1);
	
	for (int i = 0; i<bin_tmp.size(); i++){
	  Double_t E2       = (B + bin_tmp[i]*sqrt(B*B - A*C))/A;
	  Double_t P2       = sqrt(E2*E2 - m_dg[1]*m_dg[1]);
	  Double_t beta2    = P2/E2;
	  Double_t weight   = 1.0; 
	
	  /*!
	   * calculate w(E_lab, Omega_lab) = (P_lab*E_cm)*w(P_cm, Omega_cm)/(P_cm*P_cm)
	   * w(P_cm, Omega_cm) = 1/(4*TMath::Pi())*delta(P-P_cm)
	   * in order to get the distribution in lab frame need to change the variable under the delta function
	   * which basically involves retrieving d(P_cm)/d(P_lab)
	   * having two solutions leads to having two final delta-functions in lab frame
	   */
	  
	  // g_cm 
	  Double_t g_cm   = beta_parent/beta_dg_cm[1];
	  
	  // D function
	  Double_t D      = 1 + gamma_parent*gamma_parent*(1-g_cm*g_cm)*tan(theta_lab)*tan(theta_lab); 
	
	  // 1/|d(P_cm)/d(P_lab)|
	  Double_t tmp_P  = cos_theta_lab*(g_cm + bin_tmp[i]*sqrt(D))/(gamma_parent*(1 - beta_parent*beta_parent*cos_theta_lab*cos_theta_lab));
	  Double_t my = 1/g_cm;
	
	  // final distribution in the lab frame
	  Double_t w_E_Omega = fabs(P2*E_dg_cm[1]/(p_dg_cm*p_dg_cm)*tmp_P);   
	
	  //here assume the surface iz Z oriented!!!!!
	  weight *= w_E_Omega*ddir.z()/(4*TMath::Pi()*L*L);
	  weight_energy.push_back(std::make_pair(weight, E2));
	
	}
	
	return weight_energy;
}