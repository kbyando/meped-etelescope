// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_DetectorMessenger.hh" - (generic)

#ifndef MEPED_DetectorMessenger_h
#define MEPED_DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MEPED_DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class MEPED_DetectorMessenger: public G4UImessenger
{
  public:
    MEPED_DetectorMessenger(MEPED_DetectorConstruction* );
   ~MEPED_DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MEPED_DetectorConstruction* MEPED_Detector;
    
    G4UIdirectory*             MEPEDDir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        DeteMaterCmd;
    G4UIcmdWithAString*        ShieldMaterCmd;
    G4UIcmdWithAString*        HouseMaterCmd;
    G4UIcmdWithAString*        ContactMaterCmd;
    G4UIcmdWithAString*        BelvilleMaterCmd;
    G4UIcmdWithAString*        InsulatorMaterCmd;
    G4UIcmdWithAString*	       FoilMaterCmd;
    G4UIcmdWithADoubleAndUnit* DetectorThickCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

