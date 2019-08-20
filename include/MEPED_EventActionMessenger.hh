// NOAA POES Monte Carlo Simulation v1.0, 14/04/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_EventActionMessenger.hh"

#ifndef MEPED_EventActionMessenger_h
#define MEPED_EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MEPED_EventAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;

class MEPED_EventActionMessenger: public G4UImessenger
{
  public:
    MEPED_EventActionMessenger(MEPED_EventAction*);
   ~MEPED_EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MEPED_EventAction*     eventAction;
    G4UIdirectory*        eventDir;   
    G4UIcmdWithAnInteger* PrintCmd;    
};

#endif

