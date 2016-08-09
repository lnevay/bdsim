#include "BDSCutsAndLimits.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"
#include "G4ProcessManager.hh"

BDSCutsAndLimits::BDSCutsAndLimits():G4VPhysicsConstructor("BDSCutsAndLimits"),_wasActivated(false)
{
  stepLimiter = new G4StepLimiter;
  specialCuts = new G4UserSpecialCuts;
}

BDSCutsAndLimits::~BDSCutsAndLimits()
{
  delete stepLimiter;
  delete specialCuts;
}

void BDSCutsAndLimits::ConstructParticle(){
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();    
  G4Proton::Proton();
  return;
}

void BDSCutsAndLimits::ConstructProcess(){
  if(_wasActivated) return;
  _wasActivated=true;

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();

    if((particle->GetParticleName()=="gamma")||
       (particle->GetParticleName()=="e-")||
       (particle->GetParticleName()=="e+")||
       (particle->GetParticleName()=="proton")){
      particle->SetApplyCutsFlag(true);
    }
    ph->RegisterProcess(stepLimiter,particle);
#ifndef NOUSERSPECIALCUTS
    ph->RegisterProcess(specialCuts,particle);
#endif
  }
  return;
}

