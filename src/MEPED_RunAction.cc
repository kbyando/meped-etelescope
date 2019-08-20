// NOAA POES Monte Carlo Simulation v1.4, 03/07/2010e
// NOAA POES Monte Carlo Simulation v1.3, 15/09/2008e
// MEPED :: Proton Telescope 
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
//"MEPED_RunAction.cc"

#include "MEPED_RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

MEPED_RunAction::MEPED_RunAction()
{}

MEPED_RunAction::~MEPED_RunAction()
{}

void MEPED_RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "% ### Run " << aRun->GetRunID() << " start." << G4endl;

  // Print Header Informations (previously in EventAction.cc)
  G4cout << "\n%DATA BEGIN: eventID, (x,y,z)[mm], (px,py,pz), E_(incident,depD3,-999)[keV]" << G4endl;

  // Store random numbers (per event; written to CurrentEvent.rndm)
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // Initialize Cumulative Quantities
  sumEActiveD1 = sum2EActiveD1 = 0.;
  sumLActiveD1 = sum2LActiveD1 = 0.;

  // Initialize Counters
  hitCount = nullCount = 0;
}

void MEPED_RunAction::fillPerEvent(G4double EActiveD1, G4double LActiveD1)
{
  // Accumulate Statistics
  sumEActiveD1 += EActiveD1; sum2EActiveD1 += EActiveD1*EActiveD1;

  sumLActiveD1 += LActiveD1; sum2LActiveD1 += LActiveD1*LActiveD1;
}

void MEPED_RunAction::incrHitCount()
{ // increments int "hitCount"
  hitCount++;
}

void MEPED_RunAction::incrNullCount()
{ // increments int "nullCount"
  nullCount++;
}

void MEPED_RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

//  // report final seed status
//  G4RunManager::GetRunManager()->rndmSaveThisEvent();

  //compute statistics: mean and rms
  //

  sumEActiveD1 /= NbOfEvents; sum2EActiveD1 /= NbOfEvents;
  G4double rmsEActiveD1 = sum2EActiveD1 - sumEActiveD1*sumEActiveD1;
  if (rmsEActiveD1 >0.) rmsEActiveD1 = std::sqrt(rmsEActiveD1); else rmsEActiveD1 = 0.;

  sumLActiveD1 /= NbOfEvents; sum2LActiveD1 /= NbOfEvents;
  G4double rmsLActiveD1 = sum2LActiveD1 - sumLActiveD1*sumLActiveD1;
  if (rmsLActiveD1 >0.) rmsLActiveD1 = std::sqrt(rmsLActiveD1); else rmsLActiveD1 = 0.;

  //print
  //
  G4cout
     << "%DATA END------------End of Run------------------------------\n"
     << "\n% :: Energy deposited by " << hitCount << " hits, arising from " 
     << NbOfEvents << " events (" << nullCount << " null primaries).  "
     << "Not all hits will register after logic."
     << "\n% mean Energy in Detector 1 : " << G4BestUnit(sumEActiveD1,"Energy")
     << " +- "                          << G4BestUnit(rmsEActiveD1,"Energy")
     << G4endl;

  G4cout
     << "\n% mean trackLength in Detector 1 : " << G4BestUnit(sumLActiveD1,"Length")
     << " +- "                               << G4BestUnit(rmsLActiveD1,"Length")
     << G4endl;
 }


