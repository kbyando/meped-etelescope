// NOAA POES Monte Carlo Simulation v1.4, 03/07/2010e
// NOAA POES Monte Carlo Simulation v1.3, 21/10/2008e
// MEPED :: Electron Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_EventAction.hh"

#ifndef MEPED_EventAction_h
#define MEPED_EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class MEPED_RunAction;
class MEPED_EventActionMessenger;

class MEPED_EventAction : public G4UserEventAction
{
  public:
    MEPED_EventAction(MEPED_RunAction*);
   ~MEPED_EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

    // method "StoreE()": works in conjunction with MEPED_SteppingAction.cc and 
    // MEPED_EventAction.cc to store initial energy of particle (taken at step 0)
    void StoreE(G4double etot) {eTotal += etot;};

    // methods to receive intial position and momentum data (v1.2+)
    void StoreQ(G4double dx, G4double dy, G4double dz) 
	{x = dx; y = dy; z = dz;};
    void StoreP(G4double dPx, G4double dPy, G4double dPz)
	{px = dPx; py = dPy; pz = dPz;};

    // Tally deposited energies
    void AddDete1(G4double de, G4double dl) {EnergyDete1 += de; TrackLDete1 +=dl;};

    void SetPrintModulo(G4int	val) {printModulo = val;};

  private:
    MEPED_RunAction* runAct;

    G4double eTotal;
    G4double x, y, z, px, py, pz;

   // Energy and Tracklength Collection for D1 (D3) (active areas only)
    G4double EnergyDete1, TrackLDete1;

    G4int	printModulo;

    MEPED_EventActionMessenger* eventMessenger;
};

#endif
