// NOAA POES Monte Carlo Simulation v1.0, 14/04/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_SteppingAction.hh"

#ifndef MEPED_SteppingAction_h
#define MEPED_SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class MEPED_DetectorConstruction;
class MEPED_EventAction;

class MEPED_SteppingAction : public G4UserSteppingAction
{
  public:
    MEPED_SteppingAction(MEPED_DetectorConstruction*, MEPED_EventAction*);
   ~MEPED_SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    MEPED_DetectorConstruction* detector;
    MEPED_EventAction*		eventaction;
};

#endif
