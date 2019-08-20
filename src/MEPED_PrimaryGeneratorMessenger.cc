// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_PrimaryGeneratorMessenger.cc" - (generic)

#include "MEPED_PrimaryGeneratorMessenger.hh"

#include "MEPED_PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh" // (0.4.6)

MEPED_PrimaryGeneratorMessenger::MEPED_PrimaryGeneratorMessenger(
        MEPED_PrimaryGeneratorAction* MEPEDGun)
:MEPED_Action(MEPEDGun)
{
  gunDir = new G4UIdirectory("/MEPED/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
// NOTE (09/16/2008): v1.2 and earlier contained old code that provided
//  functionality since superceded by the use of the GPS module.

}

MEPED_PrimaryGeneratorMessenger::~MEPED_PrimaryGeneratorMessenger()
{
}

void MEPED_PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
}

