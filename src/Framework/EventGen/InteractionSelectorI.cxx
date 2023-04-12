//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionSelectorI.h"
#include "Framework/Interaction/Interaction.h"

using namespace genie;

//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI() :
Algorithm()
{

}
//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
InteractionSelectorI::~InteractionSelectorI()
{

}
//___________________________________________________________________________
