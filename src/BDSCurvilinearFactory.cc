#include "BDSCurvilinearFactory.hh"
#include "BDSExtent.hh"
#include "BDSGlobalConstants.hh"
#include "BDSSimpleComponent.hh"
#include "BDSTiltOffset.hh"
#include "BDSUtilities.hh"

#include "globals.hh"
#include "G4CutTubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

BDSCurvilinearFactory::BDSCurvilinearFactory():
  lengthSafety(BDSGlobalConstants::Instance()->LengthSafety())
{;}

BDSCurvilinearFactory::~BDSCurvilinearFactory()
{;}

BDSSimpleComponent* BDSCurvilinearFactory::CreateCurvilinearVolume(const G4String name,
								   const G4double chordLength,
								   const G4double radius)
{
  G4double halfLength = chordLength * 0.5 - lengthSafety;
  G4Tubs* solid = new G4Tubs(name + "_cl_solid", // name
			     0,                  // inner radius
			     radius,             // outer radius
			     halfLength,         // z half width
			     0,                  // start angle
			     CLHEP::twopi);      // sweep angle

  G4ThreeVector inputFaceNormal  = G4ThreeVector(0, 0,-1);
  G4ThreeVector outputFaceNormal = G4ThreeVector(0, 0, 1);

  return CommonConstruction(name, chordLength, chordLength, radius,
			    solid, inputFaceNormal, outputFaceNormal, 0);
}


BDSSimpleComponent* BDSCurvilinearFactory::CreateCurvilinearVolume(const G4String       name,
								   const G4double       arcLength,
								   const G4double       chordLength,
								   const G4double       radius,
								   const G4double       angle,
								   const G4ThreeVector  inputFaceNormal,
								   const G4ThreeVector  outputFaceNormal,
								   const BDSTiltOffset* tiltOffset)
{
  // angle is finite!
  // factor of 0.8 here is arbitrary tolerance as g4 cut tubs seems to fail
  // with cutting entranace / exit planes close to limit.
  // s = r*theta -> r = s/theta
  G4double radiusFromAngleLength =  std::abs(chordLength / angle) * 0.8;
  G4double radiusLocal = std::min(radius, radiusFromAngleLength);

  // copy in case we need to modify in the case of tilt offset
  G4ThreeVector inputface  = inputFaceNormal;
  G4ThreeVector outputface = outputFaceNormal;

  if (tiltOffset)
    {// could be nullptr
      G4double tilt = tiltOffset->GetTilt();
      if (BDS::IsFinite(tilt))
	{// rotate normal faces
	  inputface = inputface.rotateZ(tilt);
	  outputface = outputface.rotateZ(tilt);
	}
    }
  
  G4double halfLength = chordLength * 0.5 - lengthSafety;
  G4CutTubs* solid = new G4CutTubs(name + "_cl_solid", // name
				   0,                  // inner radius
				   radiusLocal,        // outer radius
				   halfLength,         // half length (z)
				   0,                  // rotation start angle
				   CLHEP::twopi,       // rotation sweep angle
				   inputface,          // input face normal vector
				   outputface);        // output face normal vector

  return CommonConstruction(name, arcLength, chordLength, radiusLocal, solid, inputface, outputface, angle);
}

BDSSimpleComponent* BDSCurvilinearFactory::CreateCurvilinearVolume(const G4String       name,
								   const G4double       arcLength,
								   const G4double       chordLength,
								   const G4double       radius,
								   const G4double       angle,
								   const BDSTiltOffset* tiltOffset)
{
  std::pair<G4ThreeVector,G4ThreeVector> faces = BDS::CalculateFaces(-0.5*angle, -0.5*angle);
  G4ThreeVector inputFaceNormal  = faces.first;
  G4ThreeVector outputFaceNormal = faces.second;

  return CreateCurvilinearVolume(name, arcLength, chordLength, radius, angle, inputFaceNormal, outputFaceNormal, tiltOffset);
}

BDSSimpleComponent* BDSCurvilinearFactory::CommonConstruction(const G4String      name,
							      const G4double      arcLength,
							      const G4double      chordLength,
							      const G4double      radius,
							      G4VSolid*           solid,
							      const G4ThreeVector inputFaceNormal,
							      const G4ThreeVector outputFaceNormal,
							      const G4double      angle)
{
  // nullptr for material ONLY ok in parallel world!
  G4LogicalVolume* lv =  new G4LogicalVolume(solid,            // solid
					     nullptr,          // material
					     name + "_cl_lv"); // name

  // always debug visualisation for read out geometry - only viewed via explicit commands
  lv->SetVisAttributes(BDSGlobalConstants::Instance()->GetVisibleDebugVisAttr());

  BDSSimpleComponent* result = new BDSSimpleComponent(name + "_cl",
						      arcLength,
						      angle,
						      solid,
						      lv,
						      inputFaceNormal,
						      outputFaceNormal);

  BDSExtent extent = BDSExtent(radius, radius, chordLength*0.5);
  result->SetExtent(extent);
  return result;
}
