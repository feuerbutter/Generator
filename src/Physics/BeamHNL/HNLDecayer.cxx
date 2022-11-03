//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
 John Plows <komninos-john.plows \at physics.ox.ac.uk>
 University of Oxford
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GMode.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/BeamHNL/HNLDecayer.h"
#include "Physics/BeamHNL/HNLDecayUtils.h"
#include "Physics/BeamHNL/HNLDecayMode.h"

using namespace genie;
using namespace genie::HNL;

//____________________________________________________________________________
HNLDecayer::HNLDecayer() :
EventRecordVisitorI("genie::HNL::HNLDecayer")
{

}
//____________________________________________________________________________
HNLDecayer::HNLDecayer(string config) :
EventRecordVisitorI("genie::HNL::HNLDecayer",config)
{

}
//____________________________________________________________________________
HNLDecayer::~HNLDecayer()
{

}
//____________________________________________________________________________
void HNLDecayer::ProcessEventRecord(GHepRecord * event) const
{

  LOG("HNL", pDEBUG)
    << "Entering ProcessEventRecord...";

  Interaction * interaction = event->Summary();

  fCurrInitStatePdg = interaction->InitState().ProbePdg();
  fCurrDecayMode = (HNLDecayMode_t) interaction->ExclTag().DecayMode();

  LOG("HNL", pNOTICE)
    << "Simulating HNL decay " << utils::hnl::AsString(fCurrDecayMode)
    << " for an initial state with PDG code " << fCurrInitStatePdg;

  this->AddInitialState(event);
  //this->GenerateDecayPosition(event);
  this->GenerateDecayProducts(event);
  this->UpdateEventRecord(event);

}
//____________________________________________________________________________
void HNLDecayer::AddInitialState(GHepRecord * event) const
{
  std::vector< double > * prodVtx = 0;

  Interaction * interaction = event->Summary();
  InitialState * init_state = interaction->InitStatePtr();

  TLorentzVector p4;
  if( event->Particle(0) ){
    // p4 was already set using HNLFluxCreator. No action needed.
    // Read event vertex == HNL production vertex. We will find the decay vertex later.
    p4 = *( init_state->GetProbeP4() );
    LOG( "HNL", pDEBUG ) << "\nHere's the p4 seen at InitialState(): " << utils::print::P4AsString( &p4 )
			 << "\nand the v4 seen at InitialState(): " << utils::print::X4AsString( event->Vertex() ) << " [cm, ns]";

    prodVtx = new std::vector< double >();
    prodVtx->emplace_back( event->Vertex()->X() );
    prodVtx->emplace_back( event->Vertex()->Y() );
    prodVtx->emplace_back( event->Vertex()->Z() );
    prodVtx->emplace_back( event->Vertex()->T() );
  } else {
    std::vector< double > * p3HNL = this->GenerateMomentum( event );

    double px = p3HNL->at(0);
    double py = p3HNL->at(1);
    double pz = p3HNL->at(2);
    double E = interaction->InitState().ProbeE(kRfLab);

    p4 = TLorentzVector( px, py, pz, E );

    if( !event->Vertex() || (event->Vertex()->Vect()).Mag() == 0.0 )
      prodVtx = this->GenerateDecayPosition( event );
    else{
      prodVtx = new std::vector< double >();
      prodVtx->emplace_back( event->Vertex()->X() );
      prodVtx->emplace_back( event->Vertex()->Y() );
      prodVtx->emplace_back( event->Vertex()->Z() );
      prodVtx->emplace_back( event->Vertex()->T() );
    }
  }

  // RETHERE don't sample production vtx if user isn't asking for geom! It's pointless.
  TLorentzVector v4( prodVtx->at(0), prodVtx->at(1), prodVtx->at(2), prodVtx->at(3) );

  init_state->SetProbeP4( p4 );

  LOG( "HNL", pDEBUG )
    << "\nProbe p4 = " << utils::print::P4AsString( &p4 )
    << "\nProd vtx = " << utils::print::X4AsString( &v4 );

  int hpdg = interaction->InitState().ProbePdg();
  if( !event->Particle(0) )
    event->AddParticle(hpdg, kIStInitialState, 0,-1,-1,-1, p4, v4);
}
//____________________________________________________________________________
void HNLDecayer::GenerateDecayProducts(GHepRecord * event) const
{
  LOG("HNL", pINFO) << "Generating decay...";
  fDecLepPdg = 0; fDecHadPdg = 0;

  // do we have nubar?
  PDGCodeList pdgv0 = utils::hnl::DecayProductList(fCurrDecayMode);
  int typeMod = ( fCurrInitStatePdg >= 0 ) ? 1 : -1; 
  if( fCurrDecayMode == kHNLDcyPi0Nu ){
    fDecHadPdg = typeMod * kPdgPi0;
  } else if( fCurrDecayMode != kHNLDcyNuNuNu ){
    fDecHadPdg = typeMod * kPdgPiP;
  }

  if( event->Particle(0) && !fIsMajorana ){
    typeMod = ( event->Particle(0)->Pdg() >= 0 ) ? 1 : -1;
  } else if( event->Particle(0) && fIsMajorana ){
    // equal probability to decay to 1 of 2 charge-conjugated states
    RandomGen * Rng = RandomGen::Instance();
    double rnd = Rng->RndGen().Uniform(0.0, 1.0);
    typeMod = ( rnd >= 0.5 ) ? 1.0 : -1.0;
  }
  PDGCodeList pdgv(true);
  for( std::vector<int>::iterator it = pdgv0.begin(); it != pdgv0.end(); ++it ){
    int pdgc = *it; 
    int newpdgc = ( pdgc == genie::kPdgPi0 ) ? pdgc : typeMod * pdgc; // pi-0 is its own antiparticle
    LOG("HNL", pDEBUG) << "Adding " << pdgc << " --> " << newpdgc;
    pdgv.push_back( newpdgc );
  }

  LOG("HNL", pINFO) << "Decay product IDs: " << pdgv;
  assert ( pdgv.size() > 1);

  // if user wants to include polarisation effects, start prep now
  double fPolDirMag = 0.0;
  LOG( "HNL", pDEBUG ) << "fPolDir.size() = " << fPolDir.size();
  if( fPolDir.size() == 3 ){
    fPolDirMag = std::sqrt( ( fPolDir.at(0) * fPolDir.at(0) ) +
			    ( fPolDir.at(1) * fPolDir.at(1) ) + 
			    ( fPolDir.at(2) * fPolDir.at(2) ) );
    LOG( "HNL", pDEBUG ) << "fPolDir = ( " << fPolDir.at(0) << ", " << fPolDir.at(1) << ", " << fPolDir.at(2) << " )";
  }
  bool doPol = fDoPol;

  std::ostringstream asts;
  if( !doPol ) asts << "Performing a phase space decay...";
  else asts << "Performing a polarised decay with polarisation vector:"
	    << " ( " << fPolDir.at(0) << ", " << fPolDir.at(1) << ", " << fPolDir.at(2) << " )";
  LOG("HNL", pINFO) << asts.str();

  // Get the decay product masses

  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[i++] = m;
    sum += m;
  }

  LOG("HNL", pINFO)
    << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;

  int hnl_id = 0;
  GHepParticle * hnl = event->Particle(hnl_id);
  assert(hnl);
  TLorentzVector * p4d = hnl->GetP4(); TVector3 HNLBVec = p4d->BoostVector();
  TLorentzVector * p4d_rest = (TLorentzVector *) p4d->Clone(); p4d_rest->Boost( -HNLBVec );
  TLorentzVector * v4d = hnl->GetX4();

  LOG("HNL", pINFO)
    << "\nDecaying system p4 = " << utils::print::P4AsString(p4d) 
    << "\nin rest frame, p4 = " << utils::print::P4AsString(p4d_rest);

  // Set the decay
  TGenPhaseSpace fPhaseSpaceGenerator;
  //bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d_rest, pdgv.size(), mass);
  if(!permitted) {
     LOG("HNL", pERROR)
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(p4d);
     // clean-up
     delete [] mass;
     delete p4d, p4d_rest;
     delete v4d;
     // throw exception
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Decay not permitted kinematically");
     exception.SwitchOnFastForward();
     throw exception;
  }

  // Get the maximum weight
  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int idec=0; idec<200; idec++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("HNL", pNOTICE)
     << "Max phase space gen. weight @ current hadronic system: " << wmax;

  // Generate a decay
  bool decayFailed = false;
  if( doPol && fCurrDecayMode != kHNLDcyNuNuNu ){
    // if polarisation is on we must calculate the polarisation modulus for comparison
    // for now, assume 2-body production and 2-body decay describes the pol modulus
    // see arXiv:1805.06419[hep-ph]
    TVector3 vecPolDir( fPolDir.at(0), fPolDir.at(1), fPolDir.at(2) );
    LOG( "HNL", pDEBUG ) << "Doing a polarised decay...";
    this->PolarisedDecay( fPhaseSpaceGenerator, pdgv, wmax, vecPolDir, decayFailed );
  } else {
    if( doPol && fCurrDecayMode == kHNLDcyNuNuNu ){
      // no charged lepton here... warn the user about it, though
      LOG( "HNL", pWARN )
	<< "Polarisation for uncharged FS not implemented yet, defaulting to phase-space decay...";
    }

    LOG( "HNL", pDEBUG ) << "Doing a phase-space decay...";
    this->UnpolarisedDecay( fPhaseSpaceGenerator, pdgv, wmax, decayFailed );
  }
  if( decayFailed ){
    // clean up
    delete [] mass;
    delete p4d;
    delete v4d;
    // throw exception
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Couldn't select decay after N attempts");
    exception.SwitchOnFastForward();
    throw exception;
  }

  // Insert final state products into a TClonesArray of GHepParticle's
  TLorentzVector v4(*v4d);
  int idp = 0; int npip = 0, npi0 = 0, npim = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     if( !fGetCMFrameInstead )
       p4fin->Boost( HNLBVec ); // from rest to lab again
     GHepStatus_t ist = kIStStableFinalState;
     event->AddParticle(pdgc, ist, hnl_id,-1,-1,-1, *p4fin, v4);
     switch( pdgc ){
     case kPdgPiP: npip++; break;
     case kPdgPi0: npi0++; break;
     case kPdgPiM: npim++; break;
     }
     idp++;
  }
  Interaction * interaction = event->Summary();
  interaction->ExclTagPtr()->SetNPions( npip, npi0, npim );

  // Manually set up some mother-daughter links
  // find primary (=leading) *charged!* lepton
  double Elead = -1.0; int ilead = -1;
  std::vector< int > lpdgv = { 11, 13, 15 };
  for( std::vector<int>::iterator lit = lpdgv.begin(); lit != lpdgv.end(); ++lit ) {
    GHepParticle * tmpPart = event->FindParticle( (*lit), kIStStableFinalState, 1 );
    if( !tmpPart ) tmpPart = event->FindParticle( -1 * (*lit), kIStStableFinalState, 1 ); //antiparticle?
    if( tmpPart ){
      double tmpE = tmpPart->E();
      if( tmpE > Elead ){ Elead = tmpE; ilead = event->ParticlePosition( tmpPart ); }
    }
  }
  event->Particle( 0 )->SetFirstDaughter( ilead );
  event->Particle( 0 )->SetFirstMother(-1); // why do I need to do this explicitly?
  event->Particle( 0 )->SetLastMother(-1);

  assert( event->Probe() );
  if( !event->FinalStatePrimaryLepton() ){ // no charged lepton means invisible or pi0 nu or test
    LOG( "HNL", pWARN )
      << "No final state primary lepton for this event.";
    assert( fCurrDecayMode == kHNLDcyPi0Nu || fCurrDecayMode == kHNLDcyNuNuNu
	    || fCurrDecayMode == kHNLDcyPi0Pi0Nu || fCurrDecayMode == kHNLDcyTEST );
  }
  //assert( event->FinalStatePrimaryLepton() );
  
  // loop over all FS particles and set their mother to HNL
  int itmp = 1, ilast = 1;
  while( event->Particle( itmp ) ){
    if( event->Particle(itmp)->Status() != kIStStableFinalState ){ itmp++; continue; }
    event->Particle(itmp)->SetFirstMother(0);
    event->Particle(itmp)->SetLastMother(-1);
    if( itmp != ilead ) ilast = itmp;
    itmp++;
  }
  event->Particle( 0 )->SetLastDaughter( ilast );
  // "last daughter" of HNL means last non-primary-FS-lepton, so can be less than "first" daughter

  LOG("HNL", pNOTICE)
    << "Finished with decay products. Clean up and exit!";

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//____________________________________________________________________________
std::vector< double > * HNLDecayer::GenerateDecayPosition( GHepRecord * /* event */ ) const
{
  // let's query *where* the HNL decayed from.
  if( std::strcmp( std::getenv( "PRODVTXDIR" ), "NODIR" ) != 0 ){
    if( ( !fProdVtxHist || fProdVtxHist == 0 ) ){
      std::string pvPath = std::getenv( "PRODVTXDIR" );
      LOG( "HNL", pDEBUG ) << "pvPath = " << pvPath.c_str();
      std::string pdName = "";
      std::string pvName = "hHNLVtxPos";
      fProdVtxHist = HNLFluxReader::getFluxHist3D( pvPath, pdName, pvName );
    }
    LOG( "HNL", pDEBUG )
      << "Found production vertex histo with " << fProdVtxHist->GetEntries() << " entries. Good!";
  }
  else{
    if( !fProdVtxHist ) fProdVtxHist = new TH3D( "dummy", "dummy", 100, 0, 1, 100, 0, 1, 100, 0, 1 );
  }
  assert( fProdVtxHist );
  
  std::vector< double > * prodVtx = HNLFluxReader::generateVtx3X( fProdVtxHist );
  LOG( "HNL", pDEBUG )
    << "Production vertex at: ( " << prodVtx->at(0) << ", " << prodVtx->at(1) << ", " << prodVtx->at(2) << ") [cm]";
  
  TLorentzVector v4( prodVtx->at(0), prodVtx->at(1), prodVtx->at(2), 0.0 );

  SetProdVtxPosition( v4 );

  return prodVtx;
}
//____________________________________________________________________________
std::vector< double > * HNLDecayer::GenerateMomentum( GHepRecord * event ) const
{
  Interaction * interaction = event->Summary();
  double E = interaction->InitState().ProbeE(kRfLab);
  double M = PDGLibrary::Instance()->Find(kPdgHNL)->Mass();
  double p = TMath::Sqrt(E*E-M*M);

  // set some initial deviation from beam axis due to collimation effect
  double thetaDev = fAngularDeviation; // deg
  thetaDev *= genie::constants::kPi / 180.0; // rad
  RandomGen * Rng = RandomGen::Instance();
  double theta = Rng->RndGen().Gaus(0.0, thetaDev);
  if( theta < 0.0 ) theta *= -1.0;
  double phi = Rng->RndGen().Uniform(0.0, 2.0 * genie::constants::kPi);

  double px = p * std::sin(theta) * std::cos(phi);
  double py = p * std::sin(theta) * std::sin(phi);
  double pz = p * std::cos(theta);

  std::vector< double > * p3HNL = new std::vector< double >();
  p3HNL->emplace_back(px);
  p3HNL->emplace_back(py);
  p3HNL->emplace_back(pz);

  LOG( "HNL", pDEBUG )
    << "Generated momentum: ( " << px << ", " << py << ", " << pz << " )"; 

  return p3HNL;
}
//____________________________________________________________________________
void HNLDecayer::UpdateEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();

  interaction->KinePtr()->Sett( 0.0, true );
  interaction->KinePtr()->SetW( interaction->InitStatePtr()->Probe()->Mass(), true );
  TLorentzVector * p4HNL = event->Particle(0)->GetP4(); assert( p4HNL );
  // primary lepton is FirstDaughter() of Probe()
  // need Probe() as a GHepParticle(), not a TParticlePDG()!
  // get from event record position 0
  LOG( "HNL", pDEBUG ) << "Particle(0) has PDG code " << event->Particle(0)->Pdg();
  TLorentzVector * p4FSL = 0;
  if( event->FinalStatePrimaryLepton() ){
    int iFSL = event->Particle(0)->FirstDaughter();
    LOG( "HNL", pDEBUG ) << "First daughter = " << iFSL << " with status " 
			 << (int) (event->Particle( iFSL ))->Status();
    assert( event->Particle( iFSL ) );
    p4FSL = ( event->Particle( iFSL ) )->GetP4(); 
    assert( p4FSL );
    TLorentzVector p4DIF( p4HNL->Px() - p4FSL->Px(),
			  p4HNL->Py() - p4FSL->Py(),
			  p4HNL->Pz() - p4FSL->Pz(),
			  p4HNL->E() - p4FSL->E() );
    interaction->KinePtr()->SetQ2( p4DIF.M2(), true );
    
    LOG( "HNL", pDEBUG )
      << "\nHNL p4 = ( " << p4HNL->E() << ", " << p4HNL->Px() << ", " << p4HNL->Py() << ", " << p4HNL->Pz() << " )"
      << "\nFSL p4 = ( " << p4FSL->E() << ", " << p4FSL->Px() << ", " << p4FSL->Py() << ", " << p4FSL->Pz() << " )"
      << "\nDIF p4 = ( " << p4DIF.E() << ", " << p4DIF.Px() << ", " << p4DIF.Py() << ", " << p4DIF.Pz() << " )";

  }
    
  // Set probe
  //interaction->InitStatePtr()->SetProbePdg( event->Particle(0)->Pdg() );

  //interaction->InitStatePtr()->SetProbeP4( *(event->Particle(0)->P4()) );
  
  /*
  // Set target: always Particle(1)
  // This is the charged pion in channels that have it, the pi0 in N --> pi0 pi0 v,
  // and the SM neutrino in 3-lepton channels (for N --> v v v it is the one marked "nu_e")
  interaction->InitStatePtr()->SetTgtPdg( event->Particle(1)->Pdg() );
  interaction->InitStatePtr()->SetTgtP4( *(event->Particle(1)->P4()) );
  
  LOG( "HNL", pDEBUG )
    << "Target info: " << Target();
  */
  
  // clean up
  delete p4HNL;
  if(p4FSL) delete p4FSL;
}
//____________________________________________________________________________
void HNLDecayer::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void HNLDecayer::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void HNLDecayer::LoadConfig(void)
{
  LOG("HNL", pDEBUG)
    << "Loading configuration from file...";

  this->GetParam( "HNL-Mass", fMass );
  std::vector< double > U4l2s;
  this->GetParamVect( "HNL-LeptonMixing", U4l2s );
  SetHNLCouplings( U4l2s.at(0), U4l2s.at(1), U4l2s.at(2) );
  this->GetParam( "HNL-Majorana", fIsMajorana );
  this->GetParam( "HNL-Type", fType );

  this->GetParam( "HNL-angular_deviation", fAngularDeviation );

  this->GetParam( "IncludePolarisation", fDoPol );

  this->GetParamVect( "Near2User_T", fB2UTranslation );
  this->GetParamVect( "Near2Beam_R", fB2URotation );
  SetBeam2User( fB2UTranslation, fB2URotation );

  fIntChannels = {}; bool itChan = false;
  // RETHERE for now I parse the channels manually... need to add automatic recognition?
  this->GetParam( "HNL-2B_mu_pi",  itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyPiMu );
  this->GetParam( "HNL-2B_e_pi",   itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyPiE );
  this->GetParam( "HNL-2B_nu_pi0", itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyPi0Nu );
  this->GetParam( "HNL-3B_nu_nu_nu",   itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyNuNuNu );
  this->GetParam( "HNL-3B_nu_mu_mu",   itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyNuMuMu );
  this->GetParam( "HNL-3B_nu_e_e",     itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyNuEE );
  this->GetParam( "HNL-3B_nu_mu_e",    itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyNuMuE );
  this->GetParam( "HNL-3B_e_pi_pi0",   itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyPiPi0E );
  this->GetParam( "HNL-3B_mu_pi_pi0",  itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyPiPi0Mu );
  this->GetParam( "HNL-3B_nu_pi0_pi0", itChan ); if( itChan ) fIntChannels.push_back( kHNLDcyPi0Pi0Nu );

  // RETHERE remove
  this->GetParam( "GetCMFrameInstead", fGetCMFrameInstead );

  fIsConfigLoaded = true;
}
//___________________________________________________________________________
void HNLDecayer::SetHNLCouplings( double Ue42, double Um42, double Ut42 ) const
{
  fUe42 = Ue42;
  fUm42 = Um42;
  fUt42 = Ut42;
}
//___________________________________________________________________________
void HNLDecayer::SetBeam2User( std::vector< double > translation, std::vector< double > rotation ) const
{
  fTx = -1.0 * translation.at(0);
  fTy = -1.0 * translation.at(1);
  fTz = -1.0 * translation.at(2);

  fR1 = rotation.at(0);
  fR2 = rotation.at(1);
  fR3 = rotation.at(2);

  LOG( "HNL", pDEBUG )
    << "Set BEAM origin = (0,0,0) to UNROTATED USER coordinates = ( " << fTx << ", " << fTy << ", " << fTz << " ) [m]"
    << "\nSet Euler (extrinsic x-z-x) angles to ( " << fR1 << ", " << fR2 << ", " << fR3 << " ) [rad]";  
}
//___________________________________________________________________________
SimpleHNL HNLDecayer::GetHNLInstance(string config) const
{
  SimpleHNL sh = SimpleHNL( "HNLInstance", 0, genie::kPdgHNL, genie::kPdgKP,
			    fMass, fUe42, fUm42, fUt42, fIsMajorana );
  sh.SetType( fType );
  sh.SetInterestingChannelsVec( fIntChannels );
  sh.SetAngularDeviation( fAngularDeviation );
  sh.SetBeam2UserTranslation( fTx, fTy, fTz );
  sh.SetBeam2UserRotation( fR1, fR2, fR3 );
  return sh;
}
//____________________________________________________________________________
void HNLDecayer::SetProdVtxPosition(const TLorentzVector & v4) const
{
  TLorentzVector * pv4 = new TLorentzVector();
  pv4->SetXYZT( v4.X(), v4.Y(), v4.Z(), v4.T() );
  fProdVtx = pv4;
}
//____________________________________________________________________________
void HNLDecayer::ReadCreationInfo( flux::GNuMIFluxPassThroughInfo gnmf ) const
{
  LOG( "HNL", pDEBUG )
    << "Reading creation info...";

  if( fPolDir.size() > 0 ) fPolDir.clear();
  fPolDir.emplace_back( gnmf.ppvx );
  fPolDir.emplace_back( gnmf.ppvy );
  fPolDir.emplace_back( gnmf.ppvz );

  fParentPdg = gnmf.ptype;
  fProdLepPdg = gnmf.ppmedium;
}
//____________________________________________________________________________
void HNLDecayer::UnpolarisedDecay( TGenPhaseSpace & fPSG, PDGCodeList pdgv, double wm, bool failed = false ) const
{

  RandomGen * rnd = RandomGen::Instance();
  
  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay) {
    itry++;
    
    if(itry > controls::kMaxUnweightDecayIterations) {
      // report and return
      LOG("HNL", pWARN)
	<< "Couldn't generate an unweighted phase space decay after "
	<< itry << " attempts";
      failed = true;
      return;
    }
    double w  = fPSG.Generate();
    if(w > wm) {
      LOG("HNL", pWARN)
	<< "Decay weight = " << w << " > max decay weight = " << wm;
    }
    double gw = wm * rnd->RndHadro().Rndm();
    accept_decay = (gw<=w);
    
    /*
      LOG("HNL", pDEBUG)
      << "Decay weight = " << w << " / R = " << gw
      << " - accepted: " << accept_decay;
    */
    
  } //!accept_decay
  
}
//____________________________________________________________________________
void HNLDecayer::PolarisedDecay( TGenPhaseSpace & fPSG, PDGCodeList pdgv, double wm, TVector3 vPolDir, bool failed = false ) const
{ 
  // calculate polarisation modulus
  PDGLibrary * pdgl = PDGLibrary::Instance();
  double MHNL = pdgl->Find( kPdgHNL )->Mass();
  double polMag = this->CalcPolMag( fParentPdg, fProdLepPdg, MHNL );
  double polMod = -999.99;
  
  // do decays until weight \in Uniform[0.0, 2.0] < 1 \mp polMod * cos\theta
  unsigned int iUPD = 0;
  double polWgt = -999.9;
  RandomGen * rnd = RandomGen::Instance();
  double rwgt = 0.0;
  bool isAccepted = false;

  // first, check to see if we have pi0 + v. Then let the neutrino be a QLep.
  bool isPi0Nu = ( pdgv.size() == 3 && 
		   std::abs( pdgv.at(1) ) == kPdgPi0 &&
		   std::abs( pdgv.at(2) ) == kPdgNuMu );
  
  while( !isAccepted && iUPD < controls::kMaxUnweightDecayIterations ){
    this->UnpolarisedDecay( fPSG, pdgv, wm, failed );
    
    // find charged lepton of FS. If two, take the leading one.
    // For now, this method doesn't handle vvv invisible decay mode.
    
    TVector3 lepDir;
    vector<int>::const_iterator pdg_iter;
    int idc = 0; double Elead = -1.0;
    for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
      int pdgc = *pdg_iter; Elead = -1.0;
      bool isQLep = ( std::abs( pdgc ) == kPdgElectron ||
		      std::abs( pdgc ) == kPdgMuon ||
		      std::abs( pdgc ) == kPdgTau ||
		      ( isPi0Nu && std::abs( pdgc ) == kPdgNuMu ) );
      if( isQLep ){
	TLorentzVector * p4lep = fPSG.GetDecay(idc);
	LOG( "HNL", pDEBUG )
	  << "Found QLep " << pdgc << " at idc = " << idc << " with E = " << p4lep->E();
	if( p4lep->E() > Elead ){
	  Elead = p4lep->E();
	  lepDir = p4lep->Vect();
	  fDecLepPdg = pdgc;
	  if( std::abs(fDecLepPdg) != std::abs(pdgc) || polMod < -1.0 ){
	    // update polarisation modulus for new leading lepton
	    polMod = this->CalcPolMod( polMag, fDecLepPdg, fDecHadPdg, MHNL );
	    LOG( "HNL", pDEBUG ) << "polMod = " << polMod;
	  } // std::abs(fDecLepPdg) != std::abs(pdgc) || polMod == -999.9
	} // p4lep->E() > Elead
      } // isQLep
      
      idc++;
    }

    LOG( "HNL", pDEBUG )
    << "\nfParent, ProdLep, DecLep, DecHad Pdg = " << fParentPdg
    << ", " << fProdLepPdg << ", " << fDecLepPdg << ", " << fDecHadPdg
    << "\npolMag, polMod = " << polMag << ", " << polMod
    << "\nvPolDir = " << utils::print::Vec3AsString( &vPolDir );

    // find angle \theta of leading FS QLep with vPolDir
    // assume differential decay rate \propto ( 1 \mp pMod * cos\theta )

    double theta = vPolDir.Angle( lepDir ); // rad
    double ctheta = TMath::Cos( theta );

    rwgt = rnd->RndGen().Uniform(0.0, 2.0);
    int typeMod = ( *(pdgv.begin()) > 0 ) ? 1 : -1;
    polWgt = 1 - typeMod * polMod * ctheta;

    isAccepted = ( rwgt >= polWgt );

    LOG( "HNL", pDEBUG )
      << "*** For polarised decay attempt " << iUPD << ":"
      << "\nLeading lepton has direction " << utils::print::Vec3AsString( &lepDir )
      << "\npolDir = " << utils::print::Vec3AsString( &vPolDir ) << ", angle = "
      << theta * TMath::RadToDeg() << " [deg]"
      << "\nrwgt, polWgt = " << rwgt << ", " << polWgt << ", isAccepted = " << isAccepted;

    iUPD++;
  } // while( rwgt >= polWgt && iUPD < controls::kMaxUnweightDecayIterations )

  if( iUPD == controls::kMaxUnweightDecayIterations ){
    // report and return
    LOG("HNL", pWARN)
      << "Couldn't generate a polarised decay after "
      << iUPD << " attempts";
    failed = true;
    return;
  }

}
//____________________________________________________________________________
double HNLDecayer::CalcPolMag( int parPdg, int lepPdg, double M ) const
{
  LOG( "HNL", pDEBUG )
    << "\nCalcPolMag( parPdg = " << parPdg << ", lepPdg = " << lepPdg << ", M = " << M << " )" ;

  PDGLibrary * pdgl = PDGLibrary::Instance();
  double mPar = pdgl->Find( std::abs( parPdg ) )->Mass();
  double mLep = pdgl->Find( std::abs( lepPdg ) )->Mass();

  double num1 = mLep * mLep - M * M;
  double num2 = TMath::Sqrt(utils::hnl::Kallen( mPar*mPar, M*M, mLep*mLep ));
  
  double den1 = mPar*mPar*( mLep*mLep + M*M );
  double den2 = mLep*mLep - M*M;

  // pMag is a modulus, not a magnitude... not positive semi-definite. See Fig.4 in 1805.06419[hep-ph]
  double pMag = -1.0 * num1*num2 / ( den1 - den2*den2 );

  LOG( "HNL", pDEBUG )
    << "\nmPar, mLep, M = " << mPar << ", " << mLep << ", " << M << " GeV"
    << "\nnum1, num2, den1, den2 = " << num1 << ", " << num2 << ", " << den1 << ", " << den2;
  
  return pMag;
}
//____________________________________________________________________________
double HNLDecayer::CalcPolMod( double polMag, int lepPdg, int hadPdg, double M ) const
{
  LOG( "HNL", pDEBUG )
    << "\nCalcPolMod( polMag = " << polMag << ", lepPdg = " << lepPdg << ", hadPdg = " << hadPdg << ", M = " << M << " )" ;

  PDGLibrary * pdgl = PDGLibrary::Instance();
  double mLep = pdgl->Find( std::abs( lepPdg ) )->Mass();
  double mHad = pdgl->Find( std::abs( hadPdg ) )->Mass();
  
  double num1 = M*M - mLep*mLep;
  double num2 = TMath::Sqrt(utils::hnl::Kallen( M*M, mLep*mLep, mHad*mHad ));
  double num3 = polMag;

  double den1 = M*M - mLep*mLep;
  double den2 = mHad*mHad*( M*M + mLep*mLep );

  double pMod = num1*num2*num3 / ( den1*den1 - den2 );

  LOG( "HNL", pDEBUG )
    << "\nM, mLep, mHad = " << M << ", " << mLep << ", " << mHad << " GeV"
    << "\nnum1, num2, num3, den1, den2 = " << num1 << ", " << num2 << ", "
    << num3 << ", " << den1 << ", " << den2;
  
  return pMod;
}