// NOAA POES Monte Carlo Simulation v1.4, 03/07/2010e
// NOAA POES Monte Carlo Simulation v1.3, 21/10/2008e
// MEPED :: Electron Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2010
//
// "MEPED_EventAction.cc"

#include "MEPED_EventAction.hh"

#include "MEPED_RunAction.hh"
#include "MEPED_EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

MEPED_EventAction::MEPED_EventAction(MEPED_RunAction* run)
:runAct(run),printModulo(1),eventMessenger(0)
{
  eventMessenger = new MEPED_EventActionMessenger(this);
}

MEPED_EventAction::~MEPED_EventAction()
{
  delete eventMessenger;
}

void MEPED_EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) {
//   CLHEP::HepRandom::showEngineStatus();
 }

 // Event Initialization
  EnergyDete1 = 0.;		// energy deposited in Dete1 (Dete3)
  TrackLDete1 = 0.;		// length of particle track in ""

  eTotal = 0.;			// [incident] particle energy

  // Note: these values are updated from SteppingAction using the methods
  //   StoreE(), AddDete1() defined in EventAction.hh
}

void MEPED_EventAction::EndOfEventAction(const G4Event* evt)
{
  // Accumulate Statistics
  runAct->fillPerEvent(EnergyDete1, TrackLDete1);

  // Number Per Event (modulo n)
  G4int evtNb = evt->GetEventID();

  // Output simulation data
  if (evtNb%printModulo == 0) {

    // this is where one inserts code to save the seed

    if (EnergyDete1 > 0.) {

      // hits were scored: report event information
      runAct->incrHitCount();  

      G4int prec = G4cout.precision(8);

      G4cout
	<< evtNb << ", "	// eventID
   	<< x / mm << ", " 	// x-coordinate of position, in mm
   	<< y / mm << ", " 	// y-coordinate of position, in mm
   	<< z / mm << ", "	// z-coordinate of position, in mm
   	<< px << ", " 		// x-direction of momentum (normalized)
   	<< py << ", "		// y-direction of momentum (normalized)
   	<< pz << ", "		// z-direction of momentum (normalized)
   	<< eTotal / keV << ", "	// Incident Energy, in keV
   	<< EnergyDete1 / keV << ", " // Energy Deposited in D1 (D3), in keV
   	<< -999 		// Dummy Variable (retain compatibility)
   	<< G4endl;
    }
    else {
      // no hits were scored
      runAct->incrNullCount();
    }
  }
} 

