#include "BDSTrajectory.hh"

#include "BDSDebug.hh"
#include "BDSGlobalConstants.hh"
#include "BDSTrajectoryPoint.hh"

#include "globals.hh" // geant4 globals / types
#include "G4Allocator.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4TrajectoryContainer.hh"  // also provides TrajectoryVector type(def)

#include <map>
#include <ostream>

G4Allocator<BDSTrajectory> bdsTrajectoryAllocator;

BDSTrajectory* BDS::GetPrimaryTrajectory(G4TrajectoryContainer* trajCont)
{
  TrajectoryVector* trajVec = trajCont->GetVector();
  BDSTrajectory*    primary = nullptr;
  for (const auto iT1 : *trajVec)
    {
      BDSTrajectory* traj = static_cast<BDSTrajectory*>(iT1);
      if(traj->GetParentID() == 0)
	{primary = traj; break;}
    }
  return primary;
}

BDSTrajectory::BDSTrajectory(const G4Track* aTrack, G4bool interactiveIn):
  G4Trajectory(aTrack), interactive(interactiveIn)
{
  const G4VProcess *proc = aTrack->GetCreatorProcess();
  if(proc)
  {
    creatorProcessType    = aTrack->GetCreatorProcess()->GetProcessType();
    creatorProcessSubType = aTrack->GetCreatorProcess()->GetProcessSubType();
  }
  else
  {
    creatorProcessType    = -1;
    creatorProcessSubType = -1;
  }
  weight = aTrack->GetWeight();

  fParentIndex = -1;
  fpBDSPointsContainer = new BDSTrajectoryPointsContainer();
  // this is for the first point of the track
  (*fpBDSPointsContainer).push_back(new BDSTrajectoryPoint(aTrack));
}

BDSTrajectory::~BDSTrajectory()
{
  // clean points container
  for (auto i : *fpBDSPointsContainer)
    {delete i;}
  fpBDSPointsContainer->clear();
  delete fpBDSPointsContainer;
}

void BDSTrajectory::AppendStep(const G4Step* aStep)
{
  // we do not use G4Trajectory::AppendStep here as that would
  // duplicate position information in its own vector of positions
  // which we prevent access to be overrideing GetPoint

  // TODO filter transportation steps if storing trajectory and batch

  if(BDSGlobalConstants::Instance()->TrajNoTransportation() && !interactive )
    {
      // decode aStep and if on storage.
      auto preStepPoint  = aStep->GetPreStepPoint();
      auto postStepPoint = aStep->GetPostStepPoint();
      
      // add step
      const G4VProcess *preProcess = preStepPoint->GetProcessDefinedStep();
      const G4VProcess *postProcess = postStepPoint->GetProcessDefinedStep();
      
      if (preProcess && postProcess)
	{
	  G4int preProcessType = preProcess->GetProcessType();
	  G4int postProcessType = postProcess->GetProcessType();
	  if((preProcessType  != 1   /* transportation */ &&
	      preProcessType  != 10 /* parallel world */) ||
	     (postProcessType != 1   /* transportation */ &&
	      postProcessType != 10 /* parallel world */) )
	    {fpBDSPointsContainer->push_back(new BDSTrajectoryPoint(aStep));}
	}
      /*
      else
	{
	  fpBDSPointsContainer->push_back(new BDSTrajectoryPoint(aStep));
	  }*/
    }
  else
    {
      fpBDSPointsContainer->push_back(new BDSTrajectoryPoint(aStep));
    }
}

void BDSTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif

  if(!secondTrajectory)
    {return;}
  
  BDSTrajectory* second = (BDSTrajectory*)secondTrajectory;
  G4int ent = second->GetPointEntries();
  // initial point of the second trajectory should not be merged
  for(G4int i = 1; i < ent; ++i)
    {fpBDSPointsContainer->push_back((*(second->fpBDSPointsContainer))[i]);}
  delete (*second->fpBDSPointsContainer)[0];
  second->fpBDSPointsContainer->clear();
}

BDSTrajectoryPoint* BDSTrajectory::FirstInteraction()const
{
  // loop over trajectory to find non transportation step
  for (G4int i=0; i < GetPointEntries(); ++i)
    {
      BDSTrajectoryPoint* point = static_cast<BDSTrajectoryPoint*>(GetPoint(i));
      if (point->IsScatteringPoint())
	{
	  return point;
	}
    }
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "no interaction" << G4endl;
#endif
  return static_cast<BDSTrajectoryPoint*>(GetPoint(0));
}

BDSTrajectoryPoint* BDSTrajectory::LastInteraction()const
{
  // loop over trajectory backwards to find non transportation step
  for (G4int i = GetPointEntries()-1; i >= 0; --i)
  {
    BDSTrajectoryPoint* point = static_cast<BDSTrajectoryPoint*>(GetPoint(i));
    if (point->IsScatteringPoint())
      {
	return point;
      }
  }
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "no interaction" << G4endl;
#endif
  return static_cast<BDSTrajectoryPoint*>(GetPoint(GetPointEntries()-1));
}

std::ostream& operator<< (std::ostream& out, BDSTrajectory const& t)
{
  for(G4int i = 0; i < t.GetPointEntries(); i++)
    {out << *(static_cast<BDSTrajectoryPoint*>(t.GetPoint(i))) << G4endl;}
  return out;
}

