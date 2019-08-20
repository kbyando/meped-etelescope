// NOAA POES Monte Carlo Simulation v1.0, 14/04/2008e
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_PrimaryGeneratorAction.hh" - (generic)

#ifndef MEPED_PrimaryGeneratorAction_h
#define MEPED_PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;
class MEPED_DetectorConstruction;
class MEPED_EventAction; 
class MEPED_PrimaryGeneratorMessenger;
class MEPED_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MEPED_PrimaryGeneratorAction(MEPED_DetectorConstruction*, MEPED_EventAction*); 
   ~MEPED_PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4GeneralParticleSource* particleGun;
    MEPED_DetectorConstruction* MEPED_Detector;
    MEPED_EventAction* eventaction;

    MEPED_PrimaryGeneratorMessenger* gunMessenger;
};

#endif
