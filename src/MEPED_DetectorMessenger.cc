// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_DetectorMessenger.cc" - (generic)

#include "MEPED_DetectorMessenger.hh"

#include "MEPED_DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

MEPED_DetectorMessenger::MEPED_DetectorMessenger(
                                           MEPED_DetectorConstruction* MEPEDDet)

:MEPED_Detector(MEPEDDet)
{ 
  MEPEDDir = new G4UIdirectory("/MEPED/");
  MEPEDDir->SetGuidance("UI commands of this example");
  
  detDir = new G4UIdirectory("/MEPED/det/");
  detDir->SetGuidance("detector control");
       
  DeteMaterCmd = new G4UIcmdWithAString("/MEPED/det/setDetMat",this);
  DeteMaterCmd->SetGuidance("Select Material of the Detector.");
  DeteMaterCmd->SetParameterName("choice",false);
  DeteMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ShieldMaterCmd = new G4UIcmdWithAString("/MEPED/det/setShieldMat",this);
  ShieldMaterCmd->SetGuidance("Select Material of the Shielding.");
  ShieldMaterCmd->SetParameterName("choice",false);
  ShieldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HouseMaterCmd = new G4UIcmdWithAString("/MEPED/det/setHouseMat",this);
  HouseMaterCmd->SetGuidance("Select Material of the Housing.");
  HouseMaterCmd->SetParameterName("choice",false);
  HouseMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ContactMaterCmd = new G4UIcmdWithAString("/MEPED/det/setContactMat",this);
  ContactMaterCmd->SetGuidance("Select Material of the Contacts.");
  ContactMaterCmd->SetParameterName("choice",false);
  ContactMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BelvilleMaterCmd = new G4UIcmdWithAString("/MEPED/det/setBelvMat",this);
  BelvilleMaterCmd->SetGuidance("Select Material of the Belville Washers.");
  BelvilleMaterCmd->SetParameterName("choice",false);
  BelvilleMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  InsulatorMaterCmd = new G4UIcmdWithAString("/MEPED/det/setInsulMat",this);
  InsulatorMaterCmd->SetGuidance("Select Material of the Insulator.");
  InsulatorMaterCmd->SetParameterName("choice",false);
  InsulatorMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FoilMaterCmd = new G4UIcmdWithAString("/MEPED/det/setFoilMat", this);
  FoilMaterCmd->SetGuidance("Select Material of the Cross-Aperture Foil.");
  FoilMaterCmd->SetParameterName("choice",false);
  FoilMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  DetectorThickCmd = new G4UIcmdWithADoubleAndUnit("/MEPED/det/setDeteThick",this);
  DetectorThickCmd->SetGuidance("Set Thickness of the Detector.");
  DetectorThickCmd->SetParameterName("Size",false);
  DetectorThickCmd->SetRange("Size>=0.");
  DetectorThickCmd->SetUnitCategory("Length");
  DetectorThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/MEPED/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);      
}

MEPED_DetectorMessenger::~MEPED_DetectorMessenger()
{
  delete DeteMaterCmd;
  delete ShieldMaterCmd;
  delete HouseMaterCmd;
  delete ContactMaterCmd;
  delete BelvilleMaterCmd;
  delete InsulatorMaterCmd;
  delete FoilMaterCmd;
  delete DetectorThickCmd;
  delete UpdateCmd;
  delete detDir;
  delete MEPEDDir;  
}

void MEPED_DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == DeteMaterCmd )
   { MEPED_Detector->SetDetectorMaterial(newValue);}
   
  if( command == ShieldMaterCmd )
   { MEPED_Detector->SetShieldingMaterial(newValue);}
  
  if( command == HouseMaterCmd )
   { MEPED_Detector->SetHousingMaterial(newValue);}
  
  if( command == ContactMaterCmd )
   { MEPED_Detector->SetContactMaterial(newValue);}
  
  if( command == BelvilleMaterCmd )
   { MEPED_Detector->SetBelvilleMaterial(newValue);}
  
  if( command == InsulatorMaterCmd )
   { MEPED_Detector->SetInsulatorMaterial(newValue);}
 
  if( command == FoilMaterCmd )
   { MEPED_Detector->SetFoilMaterial(newValue);}
 
  if( command == DetectorThickCmd )
   { MEPED_Detector->SetDetectorThickness(DetectorThickCmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { MEPED_Detector->UpdateGeometry(); }
}

