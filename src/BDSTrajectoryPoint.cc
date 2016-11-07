#include "G4Allocator.hh"
#include "G4ProcessType.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "globals.hh" // geant4 types / globals
#include "G4ThreeVector.hh"
#include "G4NavigationHistory.hh"

#include "BDSAuxiliaryNavigator.hh"
#include "BDSDebug.hh"
#include "BDSGlobalConstants.hh"
#include "BDSPhysicalVolumeInfoRegistry.hh"
#include "BDSPhysicalVolumeInfo.hh"
#include "BDSProcessMap.hh"
#include "BDSTrajectoryPoint.hh"

#include <ostream>

G4Allocator<BDSTrajectoryPoint> bdsTrajectoryPointAllocator;

// Don't use transform caching in the aux navigator as it's used for all over the geometry here.
BDSAuxiliaryNavigator* BDSTrajectoryPoint::auxNavigator = new BDSAuxiliaryNavigator();
G4int BDSTrajectoryPoint::numberOfPoints = 0;

BDSTrajectoryPoint::BDSTrajectoryPoint():
  G4TrajectoryPoint(G4ThreeVector())
{
  InitialiseVariables();
}

BDSTrajectoryPoint::BDSTrajectoryPoint(const G4Track* track):
  G4TrajectoryPoint(track->GetPosition())
{
  InitialiseVariables();

  // Need to store the creator process
  const G4VProcess* creatorProcess = track->GetCreatorProcess();
  if(creatorProcess) {
    preProcessType = track->GetCreatorProcess()->GetProcessType();
    preProcessSubType = track->GetCreatorProcess()->GetProcessSubType();
    postProcessType = preProcessType;
    postProcessSubType = preProcessSubType;
  }

  preWeight  = track->GetWeight();
  postWeight = preWeight;

}

BDSTrajectoryPoint::BDSTrajectoryPoint(const G4Step* step):
  G4TrajectoryPoint(step->GetPostStepPoint()->GetPosition())
{
  InitialiseVariables();
  
  const G4StepPoint* prePoint    = step->GetPreStepPoint();
  const G4StepPoint* postPoint   = step->GetPostStepPoint();
  const G4VProcess*  preProcess  = prePoint->GetProcessDefinedStep();
  const G4VProcess*  postProcess = postPoint->GetProcessDefinedStep();

  if(preProcess)
    {
      preProcessType    = preProcess->GetProcessType();
      preProcessSubType = preProcess->GetProcessSubType();
    }
  
  if(postProcess)
    {
      postProcessType    = postProcess->GetProcessType();
      postProcessSubType = postProcess->GetProcessSubType();
    }

  preWeight  = prePoint->GetWeight();
  postWeight = postPoint->GetWeight();
  energy     = step->GetTotalEnergyDeposit();
  preEnergy  = prePoint->GetKineticEnergy();
  postEnergy = postPoint->GetKineticEnergy();

  // s position for pre and post step point
  G4VPhysicalVolume* curvilinearVol = auxNavigator->LocateGlobalPointAndSetup(step);
  BDSPhysicalVolumeInfo* info = BDSPhysicalVolumeInfoRegistry::Instance()->GetInfo(curvilinearVol);

#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << BDSProcessMap::Instance()->GetProcessName(postProcessType, postProcessSubType) << G4endl;
#endif
  if (info)
    {
      prePosLocal  = auxNavigator->ConvertToLocal(prePoint->GetPosition());
      postPosLocal = auxNavigator->ConvertToLocal(postPoint->GetPosition());
      
      G4double sCentre = info->GetSPos();
      preS             = sCentre + prePosLocal.z();
      postS            = sCentre + postPosLocal.z();
      beamlineIndex    = info->GetBeamlineIndex();
      turnstaken       = BDSGlobalConstants::Instance()->TurnsTaken();
    }
}

BDSTrajectoryPoint::~BDSTrajectoryPoint()
{
  if (numberOfPoints == 0)
    {delete auxNavigator; auxNavigator = nullptr;}
  numberOfPoints--;
}

void BDSTrajectoryPoint::InitialiseVariables()
{
  numberOfPoints++;
  preProcessType     = -1;
  preProcessSubType  = -1;
  postProcessType    = -1;
  postProcessSubType = -1;
  preWeight          = -1.;
  postWeight         = -1.;
  preEnergy          = -1.;
  postEnergy         = -1.;
  energy             = -1.;
  preS               = -1000;
  postS              = -1000;
  beamlineIndex      = -1;
  turnstaken         = 0;
  prePosLocal        = G4ThreeVector();
  postPosLocal       = G4ThreeVector();
}

std::ostream& operator<< (std::ostream& out, BDSTrajectoryPoint const &p)
{
  out << p.GetPosition();
  return out;
}
