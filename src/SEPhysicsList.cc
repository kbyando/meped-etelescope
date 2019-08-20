//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Physics List for Space Electronics Applications                //
//  Date:   20 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#include "SEPhysicsList.hh"

#include "SEDecayPhysics.hh"
#include "SEBosonPhysics.hh"
#include "SELeptonPhysics.hh"
#include "SENeutronPhysics.hh"
#include "SEHadronPhysics.hh"
#include "SEIonPhysics.hh"

SEPhysicsList::SEPhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 0.01*mm; // 10 micron (yikes!)
//  defaultCutValue = 1.0*mm;
  // SetVerboseLevel(1);

  // Particle decays
  RegisterPhysics( new SEDecayPhysics("decay"));

  // Bosons (gamma + geantinos)
  RegisterPhysics( new SEBosonPhysics("boson"));

  // Leptons
  RegisterPhysics( new SELeptonPhysics("lepton"));

  // Neutron Physics
  RegisterPhysics( new SENeutronPhysics("neutron"));

  // Hadron Physics
  RegisterPhysics( new SEHadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new SEIonPhysics("ion"));
}

SEPhysicsList::~SEPhysicsList()
{;}

void SEPhysicsList::SetCuts()
{
  // Use default cut values gamma and e processes
  SetCutsWithDefault();   
}
