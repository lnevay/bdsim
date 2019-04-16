/*
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway,
University of London 2001 - 2019.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "BDSLaserPhotoDetachment.hh"
#include "BDSLogicalVolumeLaser.hh"
#include "BDSLaser.hh"
#include "BDSPhotoDetachmentEngine.hh"
#include "BDSStep.hh"
#include "BDSGlobalConstants.hh"

#include "globals.hh"
#include "G4AffineTransform.hh"
#include "G4Electron.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessType.hh"
#include "G4Proton.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

BDSLaserPhotoDetachment::BDSLaserPhotoDetachment(const G4String& processName):
  G4VDiscreteProcess(processName, G4ProcessType::fUserDefined),
  auxNavigator(new BDSAuxiliaryNavigator()),
  photoDetachmentEngine(new BDSPhotoDetachmentEngine())
{;}

BDSLaserPhotoDetachment::~BDSLaserPhotoDetachment()
{
  delete auxNavigator;
  delete photoDetachmentEngine;
}

/*
G4double BDSLaserPhotoDetachment::GetMeanFreePath(const G4Track& track,
                                                  G4double previousStepSize,
						  G4ForceCondition* forceCondition)
{
  G4LogicalVolume* lv = track.GetVolume()->GetLogicalVolume();
  if (!lv->IsExtended())
    {// not extended so can't be a laser logical volume
      return DBL_MAX;
    }
  BDSLogicalVolumeLaser* lvv = dynamic_cast<BDSLogicalVolumeLaser*>(lv);
  if (!lvv)
    {// it's an extended volume but not ours (could be a crystal)
      return DBL_MAX;
    }

  // else proceed
  const BDSLaser* laser = lvv->Laser();

  aParticleChange.Initialize(track);
  const G4DynamicParticle* ion = track.GetDynamicParticle();

  if (ion->GetCharge()==-1)
  {
    //get particle coordinates for global and local
    G4ThreeVector particlePositionGlobal = track.GetPosition();
    G4ThreeVector particleDirectionMomentumGlobal = track.GetMomentumDirection();
    const G4RotationMatrix* rot = track.GetTouchable()->GetRotation();
    const G4AffineTransform transform = track.GetTouchable()->GetHistory()->GetTopTransform();
    G4ThreeVector particlePositionLocal = transform.TransformPoint(particlePositionGlobal);
    G4ThreeVector particleDirectionMomentumLocal = transform.TransformPoint(particleDirectionMomentumGlobal).unit();

    // create photon
    G4ThreeVector photonUnit(0,0,1);
    photonUnit.transform(*rot);
    G4double photonE = (CLHEP::h_Planck*CLHEP::c_light)/laser->Wavelength();
    G4ThreeVector photonVector = photonUnit*photonE;
    G4LorentzVector photonLorentz = G4LorentzVector(photonVector,photonE);

    // ion information to boost photon energy to the right frame
    G4double ionEnergy = ion->GetTotalEnergy();
    G4ThreeVector ionMomentum = ion->GetMomentum();
    G4double ionMass = ion->GetMass();
    G4ThreeVector ionBeta = ionMomentum/ionEnergy;
    //G4double ionGamma = ionEnergy/ionMass;
    //G4ThreeVector ionVelocity = ionMomentum/(ionMass*ionGamma);
    //G4double ionVelocityMag = ionVelocity.mag();
    photonLorentz.boost(ionBeta);
    G4double photonEnergy = photonLorentz.e();
    G4double crossSection = photoDetachmentEngine->CrossSection(photonEnergy)*CLHEP::m2;

    //calculating length of particle path through laser volume
    G4double safety = BDSGlobalConstants::Instance()->LengthSafety();
    auto transportMgr = G4TransportationManager::GetTransportationManager();
    auto fLinearNavigator = transportMgr->GetNavigatorForTracking();
    G4double linearStepLength = fLinearNavigator->ComputeStep(particlePositionGlobal,
							    particleDirectionMomentumGlobal,
							    9.0e99,
							    safety);

    //things needed for loop to sum over photon density
    G4double stepSize = linearStepLength/10;// hard coded for now will later be based on max intensity and width
    G4double count = 0;
    G4double totalPhotonDensity = 0;
    G4double maxPhotonDensity = 0;
    G4double sum=0;

    for (G4double i = 0; i <= linearStepLength; i = i+stepSize)
     {
        G4ThreeVector particlePositionLocalStep=particlePositionLocal-i*particleDirectionMomentumLocal;
        G4double photonDensityStep = laser->Intensity(particlePositionLocalStep, 0)/(photonEnergy);
        totalPhotonDensity = totalPhotonDensity + photonDensityStep;
        sum += (photonDensityStep*stepSize);
        if(photonDensityStep >= maxPhotonDensity)
          {maxPhotonDensity = photonDensityStep;}
        count += 1.0;
     }
    G4double averagePhotonDensity = (totalPhotonDensity/count);
    G4double mfp = 1.0/(crossSection*sum)*CLHEP::c_light;
    return mfp;
  }
  else
    {
      return DBL_MAX;
    }
}
*/

G4double BDSLaserPhotoDetachment::GetMeanFreePath(const G4Track& track,
                                                  G4double /*previousStepSize*/,
                                                  G4ForceCondition* /*forceCondition*/)
{
  G4LogicalVolume* lv = track.GetVolume()->GetLogicalVolume();
  if (!lv->IsExtended())
  {// not extended so can't be a laser logical volume
    return DBL_MAX;
  }
  BDSLogicalVolumeLaser* lvv = dynamic_cast<BDSLogicalVolumeLaser*>(lv);
  if (!lvv)
  {// it's an extended volume but not ours (could be a crystal)
    return DBL_MAX;
  }

  // else proceed
  const BDSLaser* laser = lvv->Laser();

  aParticleChange.Initialize(track);
  const G4DynamicParticle* ion = track.GetDynamicParticle();

  if (ion->GetCharge()==-1)
  {

    //get particle coordinates for global and local
    G4ThreeVector particlePositionGlobal = track.GetPosition();
    G4ThreeVector particleDirectionMomentumGlobal = track.GetMomentumDirection();
    const G4RotationMatrix* rot = track.GetTouchable()->GetRotation();
    const G4AffineTransform transform = track.GetTouchable()->GetHistory()->GetTopTransform();
    G4ThreeVector particlePositionLocal = transform.TransformPoint(particlePositionGlobal);
    G4ThreeVector particleDirectionMomentumLocal = transform.TransformPoint(particleDirectionMomentumGlobal).unit();

    // create photon
    G4ThreeVector photonUnit(0,0,1);
    photonUnit.transform(*rot);
    G4double photonE = (CLHEP::h_Planck*CLHEP::c_light)/laser->Wavelength();
    G4ThreeVector photonVector = photonUnit*photonE;
    G4LorentzVector photonLorentz = G4LorentzVector(photonVector,photonE);

    // ion information to boost photon energy to the right frame
    G4double ionEnergy = ion->GetTotalEnergy();
    G4ThreeVector ionMomentum = ion->GetMomentum();
    G4ThreeVector ionBeta = ionMomentum/ionEnergy;
    photonLorentz.boost(ionBeta);
    G4double photonEnergy = photonLorentz.e();

    G4int copyNo = track.GetTouchable()->GetCopyNumber();

    G4double crossSection = photoDetachmentEngine->CrossSection(photonEnergy)*CLHEP::m2;
    G4double photonDensity = laser->LaserIntensities(copyNo)/photonEnergy;// units per mm^2 per ns

    //G4double ionVelocity = ionMomentum.mag()/(ion->GetMass()*(1.0/std::sqrt(1.0-ionBeta*ionBeta)));

    //G4double mft = 1.0/(crossSection*photonDensity);
    G4double mfp = 1.0/(crossSection*photonDensity)*CLHEP::c_light;
    return mfp;
  }
  else
  {
    return DBL_MAX;
  }

}

G4VParticleChange* BDSLaserPhotoDetachment::PostStepDoIt(const G4Track& track,
							 const G4Step&  step)
{
  // get coordinates for photon desity calculations
  aParticleChange.Initialize(track);
  aParticleChange.SetNumberOfSecondaries(1);

  const G4DynamicParticle* ion = track.GetDynamicParticle();
  G4double ionKe = ion->GetKineticEnergy();
  G4ThreeVector ionMomentum = ion->GetMomentum();
  G4double ionMass = ion->GetMass();

  
  G4double hydrogenMass = (CLHEP::electron_mass_c2+CLHEP::proton_mass_c2);
  aParticleChange.ProposeMass(hydrogenMass);
  G4ThreeVector hydrogenMomentum = (ionMomentum / ionMass) * hydrogenMass;
  aParticleChange.ProposeMomentumDirection(hydrogenMomentum.unit());

  G4ThreeVector electronMomentum = (ionMomentum / ionMass)*CLHEP::electron_mass_c2;  
  G4double electronKe = ionKe*(CLHEP::electron_mass_c2/ionMass);
  G4DynamicParticle* electron = new G4DynamicParticle(G4Electron::ElectronDefinition(),electronMomentum.unit(),electronKe);
  aParticleChange.AddSecondary(electron);
  aParticleChange.ProposeCharge(0);
  
  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

