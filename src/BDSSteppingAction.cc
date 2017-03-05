#include "BDSSteppingAction.hh"

#include "globals.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"


BDSSteppingAction::BDSSteppingAction():
  verboseStep(false),
  verboseEventNumber(-1)	     
{;}

BDSSteppingAction::BDSSteppingAction(G4bool verboseStepIn,
				     G4int  verboseEventNumberIn):
  verboseStep(verboseStepIn),
  verboseEventNumber(verboseEventNumberIn)
{;}

BDSSteppingAction::~BDSSteppingAction()
{;}

void BDSSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int event_number = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  if(verboseStep || event_number == verboseEventNumber)
    {VerboseSteppingAction(step);}
}

void BDSSteppingAction::VerboseSteppingAction(const G4Step* step)
{ 
  //output in case of verbose step
  G4Track* track        = step->GetTrack();
  int ID                = track->GetTrackID();
  G4VPhysicalVolume* pv = track->GetVolume();
  G4LogicalVolume* lv   = pv->GetLogicalVolume();
  G4ThreeVector pos     = track->GetPosition();
  G4ThreeVector mom     = track->GetMomentum() / CLHEP::GeV;
  G4String materialName = track->GetMaterial()->GetName();
  
  int G4precision = G4cout.precision();
  G4cout.precision(10);
  G4cout << "Physical volume = " << pv->GetName() << G4endl;
  G4cout << "Logical volume  = " << lv->GetName() << G4endl;
  G4cout << "ID="        << ID
	 << " part="     << track->GetDefinition()->GetParticleName()
	 << " Energy="   << track->GetTotalEnergy()/CLHEP::GeV
	 << " position=" << pos << "mm "
	 << " momentum=" << mom
	 << " material=" << materialName
	 << G4endl;
	    
  auto preProcess  = step->GetPreStepPoint()->GetProcessDefinedStep();
  auto postProcess = step->GetPostStepPoint()->GetProcessDefinedStep();
  
  if (preProcess)
    {G4cout << "Pre-step process= " << preProcess->GetProcessName() << G4endl;}
  if (postProcess)
    {G4cout << "Post-step process= " << postProcess->GetProcessName() << G4endl;}

  // set precision back
  G4cout.precision(G4precision);
}
