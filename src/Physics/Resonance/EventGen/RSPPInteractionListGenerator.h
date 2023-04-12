//____________________________________________________________________________
/*!

\class    genie::RSPPInteractionListGenerator

\brief    Creates a list of all the interactions that can be generated by the
          SPP thread (generates exclusive inelastic 1 pion reactions proceeding
          through resonance neutrinoproduction).
          Concrete implementations of the InteractionListGeneratorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 13, 2005

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _RSPP_INTERACTION_LIST_GENERATOR_H_
#define _RSPP_INTERACTION_LIST_GENERATOR_H_

#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/Interaction/SppChannel.h"

namespace genie {

class RSPPInteractionListGenerator : public InteractionListGeneratorI {

public :
  RSPPInteractionListGenerator();
  RSPPInteractionListGenerator(string config);
 ~RSPPInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void AddFinalStateInfo (Interaction * i, SppChannel_t chan) const;
  void LoadConfigData(void);

  bool          fIsCC;
  bool          fIsNC;
};

}      // genie namespace
#endif // _RSPP_INTERACTION_LIST_GENERATOR_H_
