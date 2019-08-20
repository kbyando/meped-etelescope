// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
// 
// "MEPED_PrimaryGeneratorMessenger.hh" - (generic)

#ifndef MEPED_PrimaryGeneratorMessenger_h
#define MEPED_PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class MEPED_PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit; // (0.4.6)

class MEPED_PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    MEPED_PrimaryGeneratorMessenger(MEPED_PrimaryGeneratorAction*);
   ~MEPED_PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);

  private:
    MEPED_PrimaryGeneratorAction* MEPED_Action;
    G4UIdirectory*               gunDir; 
};

#endif

