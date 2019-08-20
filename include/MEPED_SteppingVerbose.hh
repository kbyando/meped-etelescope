// NOAA POES Monte Carlo Simulation v1.0, 14/04/2008e
// MEPED :: Proton Telescope
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2008
//
// "MEPED_SteppingVerbose.hh"


class MEPED_SteppingVerbose;

#ifndef MEPED_SteppingVerbose_h
#define MEPED_SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class MEPED_SteppingVerbose : public G4SteppingVerbose
{
 public:   

   MEPED_SteppingVerbose();
  ~MEPED_SteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};

#endif
