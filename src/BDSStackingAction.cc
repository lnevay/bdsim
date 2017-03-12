#include "BDSGlobalConstants.hh"
#include "BDSRunManager.hh"
#include "BDSStackingAction.hh"

#include "globals.hh" // geant4 globals / types
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

BDSStackingAction::BDSStackingAction():
  maxTracksPerEvent(BDSGlobalConstants::Instance()->MaximumTracksPerEvent())
{
  auto globals = BDSGlobalConstants::Instance();
  killNeutrinos     = globals->KillNeutrinos();
  stopSecondaries   = globals->StopSecondaries();
  stopTracks        = globals->StopTracks();
}

BDSStackingAction::~BDSStackingAction()
{;}

G4ClassificationOfNewTrack BDSStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;

#ifdef BDSDEBUG 
  G4cout<<"StackingAction: ClassifyNewtrack "<<aTrack->GetTrackID()<<
    " "<<aTrack->GetDefinition()->GetParticleName()<<G4endl;
  G4StackManager* SM = G4EventManager::GetEventManager()->GetStackManager();
  G4cout<<"N total tracks :     " << SM->GetNTotalTrack()     << G4endl;
  G4cout<<"N waiting tracks :   " << SM->GetNWaitingTrack()   << G4endl;
  G4cout<<"N urgent tracks :    " << SM->GetNUrgentTrack()    << G4endl;
  G4cout<<"N postponed tracks : " << SM->GetNPostponedTrack() << G4endl;
  G4cout<<"Events to process : "<<
    BDSRunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed()<<G4endl;
  G4cout<<"Number of event : "<<
    BDSRunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent()<<G4endl;
#endif

  // '0' is default -> no action
  if (maxTracksPerEvent > 0 && (aTrack->GetTrackID() > maxTracksPerEvent))
    {return fKill;}

  // Kill all neutrinos
  if(killNeutrinos)
    {
      G4int pdgNr = std::abs(aTrack->GetParticleDefinition()->GetPDGEncoding());
      if( pdgNr == 12 || pdgNr == 14 || pdgNr == 16)
	{return fKill;}
    }

  // kill secondaries
  if(stopSecondaries && (aTrack->GetParentID() > 0) )
    {return fKill;}
  
  if(stopTracks) // if tracks killed after interaction
    {
      G4int parentID = aTrack->GetParentID();
      // note this definition getter will be deprecated in G4 soon.
      G4ParticleDefinition* definition = aTrack->GetDefinition();
      // kill secondary electrons      
      if( (parentID > 0) && (definition == G4Electron::ElectronDefinition() ) )
	{return fKill;}
      
      // kill secondary photons      
      if( (parentID > 0) && (definition == G4Gamma::GammaDefinition()))
	{return fKill;}
      
      // kill secondary positrons
      if((parentID > 0) && (definition == G4Positron::PositronDefinition()))
	{return fKill;}

      // kill secondary protons/antiprotons
      if( (parentID > 0) && 
	  ( (definition == G4Proton::ProtonDefinition() ) ||
	    (definition == G4AntiProton::AntiProtonDefinition()) ) )
	{return fKill;}
    }
  return classification;
}

void BDSStackingAction::NewStage()
{
  // urgent stack empty, looking into the waiting stack
 
#ifdef BDSDEBUG
  G4cout<<"StackingAction: New stage"<<G4endl;
#endif

  return;
}
    
void BDSStackingAction::PrepareNewEvent()
{;}


