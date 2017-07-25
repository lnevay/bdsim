#include "BDSAuxiliaryNavigator.hh"
#include "BDSEnergyCounterHit.hh"
#include "BDSEnergyCounterSD.hh"
#include "BDSDebug.hh"
#include "BDSGlobalConstants.hh"
#include "BDSPhysicalVolumeInfo.hh"
#include "BDSPhysicalVolumeInfoRegistry.hh"
#include "BDSStep.hh"
#include "BDSUtilities.hh"

#include "globals.hh" // geant4 types / globals
#include "G4AffineTransform.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

BDSEnergyCounterSD::BDSEnergyCounterSD(G4String name):
  G4VSensitiveDetector("energy_counter/"+name),
  colName(name),
  energyCounterCollection(nullptr),
  HCIDe(-1),
  enrg(0.0),
  weight(0.0),
  X(0.0),
  Y(0.0),
  Z(0.0),
  sBefore(0.0),
  sAfter(0.0),
  x(0.0),
  y(0.0),
  z(0.0),
  stepLength(0.0),
  ptype(0),
  trackID(-1),
  parentID(-1),
  volName(""),
  turnstaken(0),
  eventnumber(0),
  auxNavigator(new BDSAuxiliaryNavigator())
{
  verbose = BDSGlobalConstants::Instance()->Verbose();
  collectionName.insert(colName);
}

BDSEnergyCounterSD::~BDSEnergyCounterSD()
{
  delete auxNavigator;
}

void BDSEnergyCounterSD::Initialize(G4HCofThisEvent* HCE)
{
  energyCounterCollection = new BDSEnergyCounterHitsCollection(GetName(),colName);
  if (HCIDe < 0)
    {HCIDe = G4SDManager::GetSDMpointer()->GetCollectionID(energyCounterCollection);}
  HCE->AddHitsCollection(HCIDe,energyCounterCollection);
  
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "Energy Counter SD Hits Collection ID: " << HCIDe << G4endl;
#endif
}

G4bool BDSEnergyCounterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  /// TBC - if the vacuum is sensitive this will stop all tracks including primaries!
  /// TBC - we can cache this stop tracks bool in the class for speed.
  G4bool stopTracks = BDSGlobalConstants::Instance()->StopTracks();
  if(stopTracks)
    {enrg = (aStep->GetTrack()->GetTotalEnergy() - aStep->GetTotalEnergyDeposit());} // Why subtract the energy deposit of the step? Why not add?
  //this looks like accounting for conservation of energy when you're killing a particle
  //which may normally break energy conservation for the whole event
  //see developer guide 6.2.2...
  else
    {enrg = aStep->GetTotalEnergyDeposit();}
#ifdef BDSDEBUG
  G4cout << "BDSEnergyCounterSD> energy = " << enrg << G4endl;
#endif
  //if the energy is 0, don't do anything
  if (!BDS::IsFinite(enrg))
    {return false;}

  // avoid double getting pv
  auto hitMassWorldPV = aStep->GetPreStepPoint()->GetPhysicalVolume();
  volName             = hitMassWorldPV->GetName();
  G4int nCopy         = hitMassWorldPV->GetCopyNo();
  
  // attribute the energy deposition to a uniformly random position along the step - correct!
  // random distance - store to use twice to ensure global and local represent the same point
  G4double randDist = G4UniformRand();
  
  // global coordinate positions of the step
  G4ThreeVector posbefore = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector posafter  = aStep->GetPostStepPoint()->GetPosition();
  G4ThreeVector eDepPos   = posbefore + randDist*(posafter - posbefore);

  // calculate local coordinates
  BDSStep stepLocal = auxNavigator->ConvertToLocal(aStep);
  const G4ThreeVector& posbeforelocal = stepLocal.PreStepPoint();
  const G4ThreeVector& posafterlocal  = stepLocal.PostStepPoint();
  G4ThreeVector eDepPosLocal = posbeforelocal + randDist*(posafterlocal - posbeforelocal);
  stepLength = (posafterlocal - posbeforelocal).mag();
  
  // global
  X = eDepPos.x();
  Y = eDepPos.y();
  Z = eDepPos.z();
  // local
  x = eDepPosLocal.x();
  y = eDepPosLocal.y();
  z = eDepPosLocal.z();
  
  // get the s coordinate (central s + local z)
  // volume is from curvilinear coordinate parallel geometry
  BDSPhysicalVolumeInfo* theInfo = BDSPhysicalVolumeInfoRegistry::Instance()->GetInfo(stepLocal.VolumeForTransform());
  G4int beamlineIndex = -1;
  
  // declare lambda for updating parameters if info found (avoid duplication of code)
  auto UpdateParams = [&](BDSPhysicalVolumeInfo* info)
    {
      G4double sCentre = info->GetSPos();
      sAfter           = sCentre + posafterlocal.z();
      sBefore          = sCentre + posbeforelocal.z();
      beamlineIndex    = info->GetBeamlineIndex();
    };
  
  if (theInfo)
    {UpdateParams(theInfo);}
  else
    {
      // Try again but with the pre step point only
      G4ThreeVector unitDirection = (posafter - posbefore).unit();
      BDSStep stepLocal2 = auxNavigator->ConvertToLocal(posbefore, unitDirection);
      theInfo = BDSPhysicalVolumeInfoRegistry::Instance()->GetInfo(stepLocal2.VolumeForTransform());
      if (theInfo)
	{UpdateParams(theInfo);}
      else
	{
	  // Try yet again with just a slight shift (100um is bigger than any padding space).
	  G4ThreeVector shiftedPos = posbefore + 0.1*CLHEP::mm*unitDirection;
	  stepLocal2 = auxNavigator->ConvertToLocal(shiftedPos, unitDirection);
	  theInfo = BDSPhysicalVolumeInfoRegistry::Instance()->GetInfo(stepLocal2.VolumeForTransform());
	  if (theInfo)
	    {UpdateParams(theInfo);}
	  else
	    {
	      G4cerr << "No volume info for ";
	      auto vol = stepLocal.VolumeForTransform();
	      if (vol)
		{G4cerr << vol->GetName() << G4endl;}
	      else
		{G4cerr << "Unknown" << G4endl;}
	      // unphysical default value to allow easy identification in output
	      sAfter        = -1000;
	      sBefore       = -1000;
	      beamlineIndex = -2;
	    }
	}
    }
  
  G4double sHit = sBefore + randDist*(sAfter - sBefore);
  
  eventnumber = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  
  if(verbose)
    {
      G4cout << "BDSEnergyCounterSD: Current Volume: " 
	     << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() 
	     << "\tEvent:  " << eventnumber 
	     << "\tEnergy: " << enrg/CLHEP::GeV 
	     << "GeV\tPosition: " << sAfter/CLHEP::m <<" m"<< G4endl;
    }
  
  weight     = aStep->GetTrack()->GetWeight();
  ptype      = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  trackID    = aStep->GetTrack()->GetTrackID();
  parentID   = aStep->GetTrack()->GetParentID();
  turnstaken = BDSGlobalConstants::Instance()->TurnsTaken();
  
  //create hits and put in hits collection of the event
  BDSEnergyCounterHit* ECHit = new BDSEnergyCounterHit(nCopy,
                                                       enrg,
                                                       X, Y, Z,
                                                       sBefore,
                                                       sAfter,
                                                       sHit,
                                                       x, y, z,
                                                       volName,
                                                       ptype,
                                                       trackID,
                                                       parentID,
                                                       weight,
                                                       turnstaken,
                                                       eventnumber,
                                                       stepLength,
                                                       beamlineIndex);
  
  // don't worry, won't add 0 energy tracks as filtered at top by if statement
  energyCounterCollection->insert(ECHit);

  // TBC - this will kill all particles - both primaries and secondaries, but if it's being
  // recorded in an SD that means it's hit something, so ok
  // BUT, we can make the vacuum sensitive too for ionisation energy loss
  if(stopTracks)
    {aStep->GetTrack()->SetTrackStatus(fStopAndKill);}
   
  return true;
}
