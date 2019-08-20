// NOAA POES Monte Carlo Simulation v1.3, 17/09/2008e
// MEPED (Medium Energy Proton and Electron Detector)
// Karl Yando, Professor Robyn Millan
// Dartmouth College, 2007
//
// "MEPED_PrimaryGeneratorAction.cc" - (generic)

#include "MEPED_PrimaryGeneratorAction.hh"

#include "MEPED_DetectorConstruction.hh"
#include "MEPED_EventAction.hh"
#include "MEPED_PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

MEPED_PrimaryGeneratorAction::MEPED_PrimaryGeneratorAction(MEPED_DetectorConstruction* MEPED_DC, MEPED_EventAction* evt)
:MEPED_Detector(MEPED_DC),eventaction(evt)
{
  // Instantiate a GPS-Class Object{
  particleGun = new G4GeneralParticleSource();

  // Create a Messenger of this Class
  gunMessenger = new MEPED_PrimaryGeneratorMessenger(this);

// Particle Lookup Routine
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);

// Specify Number of Particles
  G4int n_particle = 1;
  particleGun->SetNumberOfParticles(n_particle);
}
  
MEPED_PrimaryGeneratorAction::~MEPED_PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

void MEPED_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}

