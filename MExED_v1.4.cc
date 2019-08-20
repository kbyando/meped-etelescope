// NOAA POES Monte Carlo Simulation v1.4, 03/07/2010e
// MEPED - Electron Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2010
//
// Modeled on the Geant4 package example N03 (G. Cosmo... John Allison)

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4VUIWin32.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "MEPED_DetectorConstruction.hh"
#include "SEPhysicsList.hh" // Courtesy SLAC
#include "MEPED_PrimaryGeneratorAction.hh"
#include "MEPED_RunAction.hh"
#include "MEPED_EventAction.hh"
#include "MEPED_SteppingAction.hh"
#include "MEPED_SteppingVerbose.hh"

int main(int argc,char** argv)
{
  // Choose Random Engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // User Verbose Output Class
  G4VSteppingVerbose* verbosity = new MEPED_SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);

  // construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory user initialization classes
  MEPED_DetectorConstruction* detector = new MEPED_DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new SEPhysicsList;
  runManager->SetUserInitialization(physics);
 
  G4UIsession* session=0;

  if (argc==1)  // Define UI session for interactive mode
    {
       // G4UIterminal is a dumb terminal
#if defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
      session = new G4UIterminal();
#endif
    }

// Visualization Options
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Set User Action Classes
  // Mandatory Classes
  //
  MEPED_RunAction* run_action = new MEPED_RunAction;
  runManager->SetUserAction(run_action);
  //
  MEPED_EventAction* event_action = new MEPED_EventAction(run_action);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action = new MEPED_SteppingAction(detector, event_action);
  runManager->SetUserAction(stepping_action);    

  // Moved down from the top of this list (03/27/2008)
  G4VUserPrimaryGeneratorAction* gen_action = new MEPED_PrimaryGeneratorAction(detector,event_action);
  runManager->SetUserAction(gen_action);

  // Initialize G4 Kernel
  runManager->Initialize();

  // Get the pointer to the UI manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (session)    // Define (G)UI for interactive mode
  {
     // G4UIterminal is a dumb terminal
     UI->ApplyCommand("/control/execute vis.mac");
     session->SessionStart();
     delete session;
   }
 else         // Batch Mode
   {
#ifdef G4VIS_USE
      visManager->SetVerboseLevel("quiet");
#endif
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  // job validation
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete verbosity;

  return 0;
}




