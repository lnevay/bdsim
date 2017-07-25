#include "BDSCavity.hh"
#include "BDSCavityFactoryPillBox.hh"
#include "BDSCavityInfo.hh"

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"

BDSCavityFactoryPillBox::BDSCavityFactoryPillBox()
{;}

BDSCavityFactoryPillBox::~BDSCavityFactoryPillBox()
{;}

G4double BDSCavityFactoryPillBox::CreateSolids(G4String             name,
					       G4double             totalChordLength,
					       const BDSCavityInfo* info)
{
  G4double chordLength  = totalChordLength;
  G4double cavityRadius = info->equatorRadius;
  G4VSolid* outerSolid = new G4Tubs(name + "_outer_solid",          // name
				    info->irisRadius,               // inner radius
				    cavityRadius + info->thickness, // outer radius
				    0.5*chordLength - lengthSafety, // half length
				    0.0,                            // start angle
				    CLHEP::twopi);                  // sweep angle
  allSolids.push_back(outerSolid);
   
  // creates a cylinder  to subtract from larger cylinder. 
  G4VSolid* innerSolid = new G4Tubs(name + "_inner_solid",            // name
				    0.0,                              // inner radius 
				    cavityRadius,                     // outer radius
				    0.5*chordLength - info->thickness,// half length
				    0.0,                              // star angle
				    CLHEP::twopi);                    // sweep angle
  allSolids.push_back(innerSolid);
  
  // subtraction
  cavitySolid = new G4SubtractionSolid(name + "_cavity_solid",    //name
				       outerSolid,                //solid1
				       innerSolid);               //minus solid2
  allSolids.push_back(cavitySolid);
  
  // vacuum: union of two solid - one cylinder (VacuumInnerCavity) to fill the centre,
  // and a longer, thinner cylinder (vacuumAperture) to fill the ends provided by the thickness.
  G4double vacuumInnerRadius = cavityRadius - info->thickness - lengthSafety;
  G4double vacuumHalfLength  = 0.5*chordLength - info->thickness - lengthSafety;
  G4VSolid* vacuumInnerCavity = new G4Tubs(name + "_vacuum_inner_cavity_solid",// name
					   0.0,                                // inner radius
					   vacuumInnerRadius,                  // outer radius
					   vacuumHalfLength,                   // half length
					   0.0,                                // start angle
					   CLHEP::twopi);                      // sweep angle

  G4VSolid* vacuumAperture = new G4Tubs(name + "_vacuum_aperture_solid",   // name
					0.0,                               // inner radius
					info->irisRadius - lengthSafety,   // outer radius
					0.5*chordLength - lengthSafety,    // half length
					0.0,                               // start angle
					CLHEP::twopi);                     // sweep angle
  
  // create the vacuum as a union of the two solides defined prior
  vacuumSolid = new G4UnionSolid(name + "_vacuum_solid",  // name
				 vacuumInnerCavity,       // solid one
				 vacuumAperture);         // added to solid two.
  allSolids.push_back(vacuumInnerCavity);
  allSolids.push_back(vacuumAperture);
  allSolids.push_back(vacuumSolid);

  G4double containerRadius = cavityRadius + info->thickness + lengthSafety;
  containerSolid = new G4Tubs(name + "_container_solid",   // name
			      0.0,                         // innerRadius
			      containerRadius,             // outerRadius
			      chordLength*0.5,             // half length
			      0.0,                         // starting angle
			      CLHEP::twopi);               // sweep angle
  return containerRadius;
}
