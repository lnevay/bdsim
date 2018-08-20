/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

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
#include "BDSMagnetOuterFactoryCylindrical.hh"

#include "BDSBeamPipe.hh"
#include "BDSColours.hh"
#include "BDSDebug.hh"
#include "BDSExtent.hh"
#include "BDSGeometryComponent.hh"
#include "BDSGlobalConstants.hh"
#include "BDSMagnetOuter.hh"
#include "BDSMagnetOuterInfo.hh"
#include "BDSMaterials.hh"
#include "BDSUtilities.hh"                 // for calculateorientation

#include "globals.hh"                      // geant4 globals / types
#include "G4CutTubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"

#include <cmath>
#include <utility>                         // for std::pair
#include <algorithm>                       // for std::max

BDSMagnetOuterFactoryCylindrical* BDSMagnetOuterFactoryCylindrical::instance = nullptr;

BDSMagnetOuterFactoryCylindrical* BDSMagnetOuterFactoryCylindrical::Instance()
{
  if (instance == nullptr)
    {instance = new BDSMagnetOuterFactoryCylindrical();}
  return instance;
}

BDSMagnetOuterFactoryCylindrical::BDSMagnetOuterFactoryCylindrical()
{;}

BDSMagnetOuterFactoryCylindrical::~BDSMagnetOuterFactoryCylindrical()
{
  instance = nullptr;
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateSectorBend(G4String     name,
								   G4double     length,
								   const BDSBeamPipe* beamPipe,
								   G4double     containerLength,
								   const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  G4double angleIn          = recipe->angleIn;
  G4double angleOut         = recipe->angleOut;
  
  // clear up variables
  CleanUp();
  
  // test input parameters - set global options as default if not specified
  TestInputParameters(beamPipe,horizontalWidth);
    
  // Simple cylinder if no poleface rotation, otherwise angled.
  if (!BDS::IsFinite(angleIn) && !BDS::IsFinite(angleOut))
    {
      CreateCylindricalSolids(name,length, beamPipe, containerLength, horizontalWidth);
      G4double magnetContainerRadius = (0.5 * horizontalWidth) + lengthSafety;
      BuildMagnetContainerSolidStraight(name, containerLength, magnetContainerRadius);
    }
  else
    {
      std::pair<G4ThreeVector,G4ThreeVector> faces = BDS::CalculateFaces(angleIn,angleOut);
      inputFaceNormal = faces.first;
      outputFaceNormal = faces.second;

      CreateCylindricalSolidsAngled(name, length, beamPipe, containerLength, horizontalWidth);
    
      // build the container for the whole magnet object - this horizontal width should be
      // larger than the magnet outer piece width which is just 'horizontalWidth' wide.
      G4double magnetContainerRadius = (0.5 * horizontalWidth) + lengthSafety;
      BuildMagnetContainerSolidAngled(name, containerLength, magnetContainerRadius);
    }
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateRectangularBend(G4String     name,
									G4double     length,
									const BDSBeamPipe* beamPipe,
									G4double     containerLength,
									const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  G4double angleIn          = recipe->angleIn;
  G4double angleOut         = recipe->angleOut;
  
  //rectangular bends currently just make a shorter straight volume, so ignore angle for now
  // clear up variables
  CleanUp();
  
  // test input parameters - set global options as default if not specified
  TestInputParameters(beamPipe, horizontalWidth);
    
  // Simple cylinder if no poleface rotation, otherwise angled.
  if ((!BDS::IsFinite(angleIn)) && !BDS::IsFinite(angleOut))
    {
      CreateCylindricalSolids(name,length, beamPipe, containerLength, horizontalWidth);
      BuildMagnetContainerSolidStraight(name,containerLength,0.5*horizontalWidth);
    }
  else
    {
      std::pair<G4ThreeVector,G4ThreeVector> faces = BDS::CalculateFaces(angleIn,angleOut);
      inputFaceNormal = faces.first;
      outputFaceNormal = faces.second;

      CreateCylindricalSolidsAngled(name, length, beamPipe, containerLength, horizontalWidth);
      BuildMagnetContainerSolidAngled(name, containerLength, 0.5*horizontalWidth);
    }

  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateQuadrupole(G4String     name,
								   G4double     length,
								   BDSBeamPipe* beamPipe,
								   G4double     containerLength,
								   const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateSextupole(G4String     name,
								  G4double     length,
								  BDSBeamPipe* beamPipe,
								  G4double     containerLength,
								  const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateOctupole(G4String     name,
								 G4double     length,
								 BDSBeamPipe* beamPipe,
								 G4double     containerLength,
								 const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateDecapole(G4String     name,
								 G4double     length,
								 BDSBeamPipe* beamPipe,
								 G4double     containerLength,
								 const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateSolenoid(G4String     name,
								 G4double     length,
								 BDSBeamPipe* beamPipe,
								 G4double     containerLength,
								 const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateMultipole(G4String     name,
								  G4double     length,
								  BDSBeamPipe* beamPipe,
								  G4double     containerLength,
								  const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateRfCavity(G4String     name,
								 G4double     length,
								 BDSBeamPipe* beamPipe,
								 G4double     containerLength,
								 const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateMuonSpoiler(G4String     name,
								    G4double     length,
								    BDSBeamPipe* beamPipe,
								    G4double     containerLength,
								    const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth = recipe->horizontalWidth;
  CleanUp();
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CreateKicker(G4String                  name,
							       G4double                  length,
							       const BDSBeamPipe*        beamPipe,
							       G4double                  containerLength,
							       const BDSMagnetOuterInfo* recipe,
							       G4bool                    /*vertical*/)
{
  G4double horizontalWidth  = recipe->horizontalWidth;
  CleanUp();
  // in this factory, h and v kickers will look the same so ignore bool vertical
  // have to retain it though for virtual base class compatability
  CreateCylindricalSolids(name, length, beamPipe, containerLength, horizontalWidth);
  return CommonFinalConstructor(name, length, recipe);}

/// functions below here are private to this particular factory

void BDSMagnetOuterFactoryCylindrical::CreateCylindricalSolids(G4String     name,
							       G4double     length,
							       const BDSBeamPipe* beamPipe,
							       G4double     magnetContainerLength,
							       G4double     horizontalWidth)
{
  // build the container for the whole magnet object - this horizontal width should be
  // larger than the magnet outer piece width which is just 'horizontalWidth' wide.
  G4double magnetContainerRadius = (0.5 * horizontalWidth) + lengthSafety;
  BuildMagnetContainerSolidStraight(name, magnetContainerLength, magnetContainerRadius);
  
  if (beamPipe->ContainerIsCircular())
    {
      //circular beampipe so we can simply use its radius
      yokeSolid = new G4Tubs(name + "_yoke_solid",        // name
			     beamPipe->GetContainerRadius() + 2*lengthSafety, // inner radius
			     horizontalWidth*0.5,           // outer radius
			     length*0.5-2*lengthSafety,   // half length
			     0,                           // rotation start angle
			     CLHEP::twopi);               // rotation finish angle

      //container is similar but slightly wider and hollow (to allow placement of beampipe)
      containerSolid = new G4Tubs(name + "_container_solid",      // name
				  beamPipe->GetContainerRadius() + lengthSafety, // inner radius
				  horizontalWidth*0.5 + lengthSafety,// outer radius
				  length*0.5,                      // half length
				  0,                               // rotation start angle
				  CLHEP::twopi);                   // rotation finish angle
    }
  else
    {
      G4VSolid* yokeSolidCylinder = new G4Tubs(name + "_yoke_solid_cylinder",// name
					       0,  // inner radius - for unambiguous subtraction
					       horizontalWidth*0.5,            // outer radius
					       length*0.5-2*lengthSafety,    // half length
					       0,                            // rotation start angle
					       CLHEP::twopi);                // rotation finish angle
      allSolids.push_back(yokeSolidCylinder);
      yokeSolid = new G4SubtractionSolid(name + "_yoke_solid",
					 yokeSolidCylinder,
					 beamPipe->GetContainerSubtractionSolid());

      //container is similar but slightly wider
      G4VSolid* containerSolidCylinder = new G4Tubs(name + "_container_solid_cylinder",// name
						    0,  // inner radius - for unambiguous subtraction
						    horizontalWidth*0.5 + lengthSafety,  // outer radius
						    length*0.5,                        // half length
						    0,                                 // rotation start angle
						    CLHEP::twopi);                     // rotation finish angle
      allSolids.push_back(containerSolidCylinder);
      containerSolid = new G4SubtractionSolid(name + "_container_solid",
					      containerSolidCylinder,
					      beamPipe->GetContainerSubtractionSolid());
    }
  allSolids.push_back(yokeSolid);
}

// Function for cylinder with angled faces - for pole face rotation in dipoles.
void BDSMagnetOuterFactoryCylindrical::CreateCylindricalSolidsAngled(G4String     name,
								     G4double     length,
								     const BDSBeamPipe* beamPipe,
								     G4double     magnetContainerLength,
								     G4double     horizontalWidth)
{ 
  // build the container for the whole magnet object - this horizontal width should be
  // larger than the magnet outer piece width which is just 'horizontalWidth' wide.
  G4double magnetContainerRadius = (0.5 * horizontalWidth) + lengthSafety;
  BuildMagnetContainerSolidStraight(name, magnetContainerLength, magnetContainerRadius);
  
  if (beamPipe->ContainerIsCircular())
    {
      //circular beampipe so we can simply use its radius
      yokeSolid = new G4CutTubs(name + "_yoke_solid",        // name
			     beamPipe->GetContainerRadius() + 2*lengthSafety, // inner radius
			     horizontalWidth*0.5,           // outer radius
			     length*0.5-2*lengthSafety,   // half length
			     0,                           // rotation start angle
			     CLHEP::twopi,                // rotation finish angle
			     inputFaceNormal,             // input face normal
			     outputFaceNormal);           // output face normal);

      //container is similar but slightly wider and hollow (to allow placement of beampipe)
      containerSolid = new G4CutTubs(name + "_container_solid",    // name
				  beamPipe->GetContainerRadius() + lengthSafety, // inner radius
				  horizontalWidth*0.5 + lengthSafety,// outer radius
				  length*0.5,         // half length
				  0,                  // rotation start angle
				  CLHEP::twopi,       // rotation finish angle
				  inputFaceNormal,    // input face normal
				  outputFaceNormal);  // output face normal);
    }
  else
    {
      G4VSolid* yokeSolidCylinder = new G4CutTubs(name + "_yoke_solid_cylinder",// name
						  0,  // inner radius - for unambiguous subtraction
						  horizontalWidth*0.5,        // outer radius
						  length*0.5-2*lengthSafety,// half length
						  0,                        // rotation start angle
						  CLHEP::twopi,             // rotation finish angle
						  inputFaceNormal,          // input face normal
						  outputFaceNormal);              // output face normal);
      allSolids.push_back(yokeSolidCylinder);
      yokeSolid = new G4SubtractionSolid(name + "_yoke_solid",
					 yokeSolidCylinder,
					 beamPipe->GetContainerSubtractionSolid());

      //container is similar but slightly wider
      G4VSolid* containerSolidCylinder = new G4CutTubs(name + "_container_solid_cylinder",// name
						       0,  // inner radius - for unambiguous subtraction
						       horizontalWidth*0.5 + lengthSafety,// outer radius
						       length*0.5,       // half length
						       0,                // rotation start angle
						       CLHEP::twopi,     // rotation finish angle
						       inputFaceNormal,  // input face normal
						       outputFaceNormal);// output face normal);
      allSolids.push_back(containerSolidCylinder);
      containerSolid = new G4SubtractionSolid(name + "_container_solid",
					      containerSolidCylinder,
					      beamPipe->GetContainerSubtractionSolid());
    }
  allSolids.push_back(yokeSolid);
}

void BDSMagnetOuterFactoryCylindrical::TestInputParameters(const BDSBeamPipe* beamPipe,
							   G4double&    horizontalWidth)
{
  // ensure box size is bigger than the beampipe
  if (beamPipe->ContainerIsCircular())
    {
      // if it's circular, just check radius
      if (horizontalWidth < 2*(beamPipe->GetContainerRadius()) )
	{horizontalWidth = 2*(beamPipe->GetContainerRadius()) + 1*CLHEP::mm;}
    }
  else
    {
      // it's not circular - have a look at extents
      // +ve - -ve
      G4double extentX = beamPipe->GetExtentX().second - beamPipe->GetExtentX().first;
      G4double extentY = beamPipe->GetExtentY().second - beamPipe->GetExtentY().first;
      if ( (horizontalWidth < extentX) || (horizontalWidth < extentY) )
	{
	  // horizontalWidth isn't sufficient for range in x or y
	  horizontalWidth = std::max(extentX,extentY) + 1*CLHEP::mm;
	}
    }
}

BDSMagnetOuter* BDSMagnetOuterFactoryCylindrical::CommonFinalConstructor(G4String    name,
									 G4double    length,
									 const BDSMagnetOuterInfo* recipe)
{
  G4double horizontalWidth = recipe->horizontalWidth;
  G4Material* outerMaterial = recipe->outerMaterial;
  if (!outerMaterial)
    {outerMaterial = BDSMaterials::Instance()->GetMaterial(BDSGlobalConstants::Instance()->EmptyMaterial());}
  G4Colour* colour = recipe->colour;
  
  CreateLogicalVolumes(name, colour, outerMaterial);
  SetUserLimits();
  CreateMagnetContainerComponent();

  // PLACEMENT
  // place the components inside the container
  // note we don't need the pointer for anything - it's registered upon construction with g4
  yokePV = new G4PVPlacement((G4RotationMatrix*)nullptr,   // no rotation
			     (G4ThreeVector)0,             // position
			     yokeLV,                       // lv to be placed
			     name + "_yoke_pv",            // name
			     containerLV,                  // mother lv to be placed in
			     false,                        // no boolean operation
			     0,                            // copy number
			     checkOverlaps); // whether to check overlaps

  // record extents
  // container radius is the same for all methods as all cylindrical
  G4double containerRadius = horizontalWidth + lengthSafety;
  BDSExtent ext = BDSExtent(containerRadius, containerRadius, length*0.5);
  
  // build the BDSMagnetOuter instance and return it
  BDSMagnetOuter* outer = new BDSMagnetOuter(containerSolid,
					     containerLV, ext,
					     magnetContainer);

  SetFaceNormals(outer);
  
  // register all objects that go with the final geometry component (from internal vectors)
  outer->RegisterSolid(allSolids);
  outer->RegisterLogicalVolume(yokeLV);
  outer->RegisterSensitiveVolume(yokeLV);
  outer->RegisterPhysicalVolume(yokePV);
  outer->RegisterVisAttributes(outerVisAttributes);
  // no rotation matrices used in this factory
  return outer;
}


