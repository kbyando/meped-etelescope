// NOAA POES Monte Carlo Simulation v0.2.1, 22/01/2007
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2007
//
// "MEPED_EventActionMessenger.cc"

#include "MEPED_EventActionMessenger.hh"

#include "MEPED_EventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

MEPED_EventActionMessenger::MEPED_EventActionMessenger(MEPED_EventAction* EvAct)
:eventAction(EvAct)
{
  eventDir = new G4UIdirectory("/MEPED/event/");
  eventDir->SetGuidance("event control");
   
  PrintCmd = new G4UIcmdWithAnInteger("/MEPED/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
}

MEPED_EventActionMessenger::~MEPED_EventActionMessenger()
{
  delete PrintCmd;
  delete eventDir;   
}

void MEPED_EventActionMessenger::SetNewValue(
                                        G4UIcommand* command,G4String newValue)
{ 
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}
}

