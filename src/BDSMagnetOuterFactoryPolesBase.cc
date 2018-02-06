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
#include "BDSMagnetOuterFactoryPolesBase.hh"

#include "BDSBeamPipe.hh"
#include "BDSColours.hh"
#include "BDSDebug.hh"
#include "BDSExtent.hh"
#include "BDSSimpleComponent.hh"
#include "BDSGlobalConstants.hh"
#include "BDSMagnetOuter.hh"
#include "BDSMagnetOuterFactoryCylindrical.hh" // for default geometry
#include "BDSMaterials.hh"
#include "BDSUtilities.hh"

#include "globals.hh"                      // geant4 globals / types
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4CutTubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4TwoVector.hh"
#include "G4VisAttributes.hh"
#include "G4VSolid.hh"

#include <algorithm>                       // for std::max, std::swap
#include <cmath>
#include <string>                          // for std::to_string
#include <utility>                         // for std::pair
#include <vector>

BDSMagnetOuterFactoryPolesBase::BDSMagnetOuterFactoryPolesBase():
  BDSMagnetOuterFactoryPolesBase(1.1)
{
  CleanUp();
}

BDSMagnetOuterFactoryPolesBase::BDSMagnetOuterFactoryPolesBase(G4double poleStopFactorIn):
  poleFraction(0.7),
  poleAngularFraction(0.7),
  poleTipFraction(0.2),
  poleAnnulusFraction(0.1),
  bendHeightFraction(0.7),
  poleStopFactor(poleStopFactorIn)
{
  // now the base class constructor should be called first which
  // should call clean up (in the derived class) which should initialise
  // the variables I think, but doing here just to be sure.
  CleanUp();
}

void BDSMagnetOuterFactoryPolesBase::CleanUp()
{
  BDSMagnetOuterFactoryBase::CleanUp();

  yokeStartRadius       = 0;
  yokeFinishRadius      = 0;
  magnetContainerRadius = 0;
  buildPole             = true;
  poleStartRadius       = 0;
  poleFinishRadius      = 0;
  poleSquareWidth       = 0;
  poleSquareStartRadius = 0;
  segmentAngle          = 0;
  poleAngle             = 0;
  poleTranslation       = G4ThreeVector(0,0,0);
  coilHeight            = 0;
  coilCentreRadius      = 0;
  endPieceLength        = 0;
  endPieceInnerR        = 0;
  endPieceOuterR        = 0;
  poleIntersectionSolid = nullptr;
  coilLeftSolid         = nullptr;
  coilRightSolid        = nullptr;
  endPieceContainerSolid = nullptr;
  coilLeftLV            = nullptr;
  coilRightLV           = nullptr;
  endPieceCoilLV        = nullptr;
  endPieceContainerLV   = nullptr;
  endPiece              = nullptr;
  leftPoints.clear();
  rightPoints.clear();
  endPiecePoints.clear();
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateSectorBend(G4String      name,
								 G4double      length,
								 BDSBeamPipe*  beamPipe,
								 G4double      outerDiameter,
								 G4double      containerLength,
								 G4double      angleIn,
								 G4double      angleOut,
								 G4bool        yokeOnLeft,
								 G4bool        hStyle,
								 G4Material*   outerMaterial,
								 G4bool        buildEndPiece,
								 G4double      vhRatio,
								 G4double      coilWidthFraction,
								 G4double      coilHeightFraction)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  auto colour = BDSColours::Instance()->GetColour("sectorbend");
  if (hStyle)
    {
      return CreateDipoleH(name, length, beamPipe, outerDiameter, containerLength, angleIn,
			   angleOut, outerMaterial, colour, false, buildEndPiece, vhRatio,
			   coilWidthFraction, coilHeightFraction);
    }
  else
    {
      return CreateDipoleC(name, length, beamPipe, outerDiameter, containerLength, angleIn,
			   angleOut, outerMaterial, yokeOnLeft, colour, false, buildEndPiece,
			   vhRatio, coilWidthFraction, coilHeightFraction);
    }
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateRectangularBend(G4String      name,
								      G4double      length,
								      BDSBeamPipe*  beamPipe,
								      G4double      outerDiameter,
								      G4double      containerLength,
								      G4double      angleIn,
								      G4double      angleOut,
								      G4bool        yokeOnLeft,
								      G4bool        hStyle,
								      G4Material*   outerMaterial,
								      G4bool        buildEndPiece,
								      G4double      vhRatio,
								      G4double      coilWidthFraction,
								      G4double      coilHeightFraction)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  auto colour = BDSColours::Instance()->GetColour("rectangularbend");
  if (hStyle)
    {
      return CreateDipoleH(name, length, beamPipe, outerDiameter, containerLength, angleIn,
			   angleOut, outerMaterial, colour, false, buildEndPiece, vhRatio,
			   coilWidthFraction, coilHeightFraction);
    }
  else
    {
      return CreateDipoleC(name, length, beamPipe, outerDiameter, containerLength, angleIn,
			   angleOut, outerMaterial, yokeOnLeft, colour, false, buildEndPiece,
			   vhRatio, coilWidthFraction, coilHeightFraction);
    }
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateQuadrupole(G4String      name,
								 G4double      length,
								 BDSBeamPipe*  beamPipe,
								 G4double      outerDiameter,
								 G4double      containerLength,
								 G4Material*   outerMaterial,
								 G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return CommonConstructor(name, length, beamPipe, 2, outerDiameter, outerMaterial, containerLength, buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateSextupole(G4String      name,
								G4double      length,
								BDSBeamPipe*  beamPipe,
								G4double      outerDiameter,
								G4double      containerLength,
								G4Material*   outerMaterial,
								G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return CommonConstructor(name, length, beamPipe, 3, outerDiameter, outerMaterial, containerLength, buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateOctupole(G4String      name,
							       G4double      length,
							       BDSBeamPipe*  beamPipe,
							       G4double      outerDiameter,
							       G4double      containerLength,
							       G4Material*   outerMaterial,
							       G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return CommonConstructor(name, length, beamPipe, 4, outerDiameter, outerMaterial, containerLength, buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateDecapole(G4String      name,
							       G4double      length,
							       BDSBeamPipe*  beamPipe,
							       G4double      outerDiameter,
							       G4double      containerLength,
							       G4Material*   outerMaterial,
							       G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return CommonConstructor(name, length, beamPipe, 5, outerDiameter, outerMaterial, containerLength, buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateSolenoid(G4String      name,
							       G4double      length,
							       BDSBeamPipe*  beamPipe,
							       G4double      outerDiameter,
							       G4double      containerLength,
							       G4Material*   outerMaterial,
							       G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return BDSMagnetOuterFactoryCylindrical::Instance()->CreateSolenoid(name,length,beamPipe,outerDiameter,
								      containerLength,outerMaterial,buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateMultipole(G4String      name,
								G4double      length,
								BDSBeamPipe*  beamPipe,
								G4double      outerDiameter,
								G4double      containerLength,
								G4Material*   outerMaterial,
								G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return BDSMagnetOuterFactoryCylindrical::Instance()->CreateMultipole(name,length,beamPipe,outerDiameter,
								       containerLength,outerMaterial,buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateRfCavity(G4String      name,
							       G4double      length,
							       BDSBeamPipe*  beamPipe,
							       G4double      outerDiameter,
							       G4double      containerLength,
							       G4Material*   outerMaterial,
							       G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return BDSMagnetOuterFactoryCylindrical::Instance()->CreateRfCavity(name,length,beamPipe,outerDiameter,
								      containerLength,outerMaterial,buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateMuSpoiler(G4String      name,
								G4double      length,
								BDSBeamPipe*  beamPipe,
								G4double      outerDiameter,
								G4double      containerLength,
								G4Material*   outerMaterial,
								G4bool        buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  return BDSMagnetOuterFactoryCylindrical::Instance()->CreateMuSpoiler(name,length,beamPipe,outerDiameter,
								       containerLength,outerMaterial,buildEndPiece);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateKicker(G4String      name,
							     G4double      length,
							     BDSBeamPipe*  beamPipe,
							     G4double      outerDiameter,
							     G4double      containerLength,
							     G4bool        vertical,
							     G4Material*   outerMaterial,
							     G4bool        buildEndPiece,
							     G4bool        hStyle,
							     G4double      vhRatio,
							     G4double      coilWidthFraction,
							     G4double      coilHeightFraction)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  G4String colourName = (vertical) ? "vkicker" : "hkicker";
  auto colour = BDSColours::Instance()->GetColour(colourName);
  if (hStyle)
    {
      return CreateDipoleH(name, length, beamPipe, outerDiameter, containerLength, 0, 0,
			   outerMaterial, colour, vertical, buildEndPiece, vhRatio,
			   coilWidthFraction, coilHeightFraction);
    }
  else
    {
      return CreateDipoleC(name, length, beamPipe, outerDiameter, containerLength, 0, 0,
			   outerMaterial, true, colour, vertical, buildEndPiece, vhRatio,
			   coilWidthFraction, coilHeightFraction);
    }
}

/// functions below here are private to this particular factory
BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CommonConstructor(G4String     name,
								  G4double     length,
								  BDSBeamPipe* beamPipe,
								  G4int        order,
								  G4double     outerDiameter,
								  G4Material*  outerMaterial,
								  G4double     magnetContainerLength,
								  G4bool       buildEndPiece)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  // reset all pointers and variables to protect against bugs using previous use of factory
  CleanUp();
  
  CalculatePoleAndYoke(outerDiameter, beamPipe, order);
  if (buildPole)
    {
      CreatePoleSolid(name, length, order);
      CreateCoilSolids(name, length);
    }
  CreateYokeAndContainerSolid(name, length, order, magnetContainerLength);
  if (buildPole)
    {IntersectPoleWithYoke(name, length, order);}
  G4Colour* magnetColour = BDSColours::Instance()->GetMagnetColour(order);
  CreateLogicalVolumes(name, magnetColour, outerMaterial);
  CreateMagnetContainerComponent();
  if (buildPole && buildEndPiece)
    {CreateEndPiece(name);}
  PlaceComponents(name, order); //returns vector of PVs
  if (buildPole)
    {PlaceComponentsCoils(name, order);}
  
  // record extents
  // container radius is just outerDiamter as yoke is circular
  G4double containerRadius = outerDiameter + lengthSafety;
  BDSExtent ext = BDSExtent(containerRadius, containerRadius, length*0.5);
  
  // build the BDSMagnetOuter instance and return it
  BDSMagnetOuter* outer = new BDSMagnetOuter(containerSolid,
					     containerLV, ext,
					     magnetContainer);

  outer->RegisterRotationMatrix(allRotationMatrices);
  outer->RegisterPhysicalVolume(allPhysicalVolumes);
  outer->RegisterVisAttributes(allVisAttributes);
  
  // Register logical volumes and set sensitivity
  // test if poleLV exists first as some derived classes use their own vector of
  // poles and don't use this one - therefore it'll be null
  if (poleLV)
    {
      outer->RegisterLogicalVolume(poleLV);
      outer->RegisterSensitiveVolume(poleLV);
    }
  outer->RegisterLogicalVolume(yokeLV);
  outer->RegisterSensitiveVolume(yokeLV);

  outer->SetEndPieceBefore(endPiece);
  outer->SetEndPieceAfter(endPiece);

  SetFaceNormals(outer);

  return outer;
}

void BDSMagnetOuterFactoryPolesBase::CalculatePoleAndYoke(G4double     outerDiameter,
							  BDSBeamPipe* beamPipe,
							  G4int        order)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  G4double bpRadius = beamPipe->GetContainerRadius();
  
  // check parameters are valid
  if (outerDiameter*0.5 < bpRadius)
    {
      G4cerr << __METHOD_NAME__
	     << "outerDiameter (" << outerDiameter << ") must be greater than 2*beampipe radius ("
	     << 2*bpRadius << ")" << G4endl;
      exit(1);
    }

  // layout markers for the pole and yoke - radially out from centre
  poleStartRadius       = bpRadius + lengthSafety;
  yokeFinishRadius      = outerDiameter*0.5 - lengthSafety;
  G4double totalLength  = yokeFinishRadius - poleStartRadius;
  poleFinishRadius      = poleStartRadius + poleFraction*totalLength;
  yokeStartRadius       = poleFinishRadius + lengthSafety;
  magnetContainerRadius = yokeFinishRadius + lengthSafetyLarge;
  poleSquareWidth       = (yokeFinishRadius - yokeStartRadius)*1.3;

  G4int nPoles = 2*order;
  // full circle is divided into segments for each pole
  segmentAngle = CLHEP::twopi / (G4double)nPoles;
  poleAngle = segmentAngle * poleAngularFraction;

  G4double minimumPoleFraction = 0.05;
  // If the space occupied by the yoke / pole is < 5% of outerDiameter, don't bother with
  // pole and coil - it's clearly unphysical.
  if ( (yokeFinishRadius - yokeStartRadius) < (minimumPoleFraction * outerDiameter) )
    {
        buildPole = false;
        yokeStartRadius = poleStartRadius; // make the magnet all yoke
    }
}

void BDSMagnetOuterFactoryPolesBase::CreatePoleSolid(G4String name,
						     G4double length,
						     G4int    order)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  // calculate geometrical parameters first.
  // pole is ellipse at tip, then (possibly) tapered section going outwards and then
  // a section that is straight - ie constant width going outwards. This is later
  // intersected with the yoke solid to give the matching shape. Therefore, the pole
  // is made slightly too long (ie > poleFinishRadius).
    
  // make some % of pole length the curved ellipsoidal surface at the pole tip
  G4double poleLength      = poleFinishRadius - poleStartRadius - 2*lengthSafety;
  G4double ellipsoidHeight = poleTipFraction*poleLength; // full height of an ellipse (2*a)

  // calculate the ellipsoid centre radius from the centre of the magnet.
  G4double ellipsoidCentreY = poleStartRadius + ellipsoidHeight*0.5;
  
  // calculate the width (2*b) of the ellipse given its height - at the centre of the
  // ellipsoid, go out by the pole angle and that's the half width, x2.
  G4double ellipsoidWidth  = 2 * ellipsoidCentreY * tan(poleAngle*0.5);
  // don't allow the pole to be wider than it is long - silly proportions - stop the scaling
  ellipsoidWidth = std::min(ellipsoidWidth, poleLength); 

  // calculate the distance from the centre of the magnet where the straight sides start
  G4double fraction = std::min(0.5, (G4double)order/10.0); // scale slightly with order
  poleSquareStartRadius = ellipsoidCentreY + (poleFinishRadius - ellipsoidCentreY) * fraction;

  // poleSquareWidth is initially calculated in CalculatePoleAndYoke - here we check it
  // ensure pole doesn't get narrower than ellipsoid tip part
  poleSquareWidth = std::max(poleSquareWidth, ellipsoidWidth);
  // but ensure the square section can't protrude outside the angular segment for a single pole
  G4double maxPoleSquareWidth = 2 * poleSquareStartRadius * tan(0.5*poleAngle);
  poleSquareWidth = std::min(maxPoleSquareWidth, poleSquareWidth);

  // points that define extruded polygon
  std::vector<G4TwoVector> points;

  // generate points on a ellipse in the positive quadrant for the pole tip
  // note, the QT visualiser can't visualise the union of an ellipsoid and extruded
  // polygon so can't use that method and this method produces a better transition
  // between the ellipse and the pole. We just generate an extruded solid for the whole pole.
  std::vector<G4double> xEllipse;
  std::vector<G4double> yEllipse;
  const G4int numAngPoints = 6;
  const G4double iterant = CLHEP::halfpi/(G4double)numAngPoints;
  for (G4double angle = 0; angle < CLHEP::halfpi; angle += iterant)
    {
      xEllipse.push_back(0.5*ellipsoidWidth*sin(angle));
      yEllipse.push_back(0.5*ellipsoidHeight*cos(angle));
    }

  // add bottom left quadrant - miss out first point so don't reach apex.
  for (G4int i = 0; i < (G4int)xEllipse.size()-1; i++)
    {points.emplace_back(xEllipse[i], ellipsoidCentreY-yEllipse[i]);}
  points.emplace_back(poleSquareWidth*0.5, poleSquareStartRadius);
  // top points are x poleStopFactor for unambiguous intersection later on
  points.emplace_back(poleSquareWidth*0.5, poleFinishRadius*poleStopFactor);
  points.emplace_back(-poleSquareWidth*0.5, poleFinishRadius*poleStopFactor);
  points.emplace_back(-poleSquareWidth*0.5, poleSquareStartRadius);
  // add bottom right quadrant - miss out first point so don't reach apex.
  // required start index is of course size()-1 also.
  for (G4int i = (G4int)xEllipse.size()-2; i > 0; i--)
    {points.emplace_back(-xEllipse[i], ellipsoidCentreY-yEllipse[i]);}
  
  G4TwoVector zOffsets(0,0); // the transverse offset of each plane from 0,0
  G4double zScale = 1;       // the scale at each end of the points = 1
  poleSolid = new G4ExtrudedSolid(name + "_pole_raw_solid",  // name
				  points,                    // transverse 2d coordinates
				  length*0.5 - lengthSafety, // z half length
				  zOffsets, zScale,          // dx,dy offset for each face, scaling
				  zOffsets, zScale);         // dx,dy offset for each face, scaling
  
  allSolids.push_back(poleSolid);
}

void BDSMagnetOuterFactoryPolesBase::CreateCoilSolids(G4String name,
						      G4double length)
{
  // create an extruded polygon even though a square to avoid the need for
  // individual rotated translations. This also allows the same rotation
  // to be used for the coils as the poles.
  CreateCoilPoints();

  G4TwoVector zOffsets(0,0); // the transverse offset of each plane from 0,0
  G4double zScale = 1;       // the scale at each end of the points = 1
  
  coilLeftSolid = new G4ExtrudedSolid(name + "_coil_left_solid", // name
				      leftPoints,                // transverse 2d coordinates
				      length*0.5 - lengthSafety, // z half length
				      zOffsets, zScale, // dx,dy offset for each face, scaling
				      zOffsets, zScale);// dx,dy offset for each face, scaling

  coilRightSolid = new G4ExtrudedSolid(name + "_coil_right_solid", // name
				       rightPoints,                // transverse 2d coordinates
				       length*0.5 - lengthSafety, // z half length
				       zOffsets, zScale, // dx,dy offset for each face, scaling
				       zOffsets, zScale);// dx,dy offset for each face, scaling
  
  allSolids.push_back(coilLeftSolid);
  allSolids.push_back(coilRightSolid);
}

void BDSMagnetOuterFactoryPolesBase::CreateCoilPoints()
{
  G4double innerX = 0.5*poleSquareWidth + lengthSafetyLarge;

  // Start the coil just after the pole becomes square in cross-section
  // start it just beyond this point - say + 20% square section length
  G4double lowerY = poleSquareStartRadius + 0.2*(poleFinishRadius - poleSquareStartRadius);

  // At lower Y, project out in X as far as possible within the pole segment
  G4double outerX = 0.9 * lowerY * tan (segmentAngle*0.5);

  // At outer X, we find out the maximum y before hitting the inside of the yoke and
  // back off a wee bit
  G4double upperY = 0.95 * std::abs(std::sqrt(std::pow(yokeStartRadius,2) - std::pow(outerX,2)));
  
  coilHeight  = upperY - lowerY;
  G4double dx = outerX - innerX;
  // Check proportions and make roughly square.
  if (coilHeight < dx)
    {
      G4double newXHeight = std::max(dx*0.5, dx - (dx-coilHeight));
      outerX = innerX + newXHeight;
      // recalculate upperY given the new outerX
      upperY = 0.95 * std::abs(std::sqrt(std::pow(yokeStartRadius,2) - std::pow(outerX,2)));
    }

  // update coil height as used for end piece construction later
  coilHeight = upperY - lowerY;
  
  coilCentreRadius = lowerY + 0.5*coilHeight;
  endPieceInnerR   = lowerY - lengthSafetyLarge;
  endPieceOuterR   = yokeStartRadius;
  
  leftPoints.emplace_back(innerX, lowerY);
  leftPoints.emplace_back(outerX, lowerY);
  leftPoints.emplace_back(outerX, upperY);
  leftPoints.emplace_back(innerX, upperY);

  // must be in clockwise order
  rightPoints.emplace_back(-innerX, lowerY);
  rightPoints.emplace_back(-innerX, upperY);
  rightPoints.emplace_back(-outerX, upperY);
  rightPoints.emplace_back(-outerX, lowerY);

  // this will be the eventual length along z but for now its the amplitude in y.
  // make slightly smaller version as endPieceLength used for container dimensions
  endPieceLength = (outerX - innerX) - lengthSafetyLarge;
  for (G4double angle = 0; angle <= CLHEP::halfpi; angle+= CLHEP::halfpi / 8.0)
    {
      G4double x = outerX + endPieceLength * (cos(angle) - 1.0);
      G4double y = endPieceLength * sin(angle);
      endPiecePoints.emplace_back(x,y);
    }
  for (G4double angle = 0; angle <= CLHEP::halfpi; angle+= CLHEP::halfpi / 8.0)
    {
      G4double x = -outerX - endPieceLength * (sin(angle) - 1.0);
      G4double y = endPieceLength * cos(angle);
      endPiecePoints.emplace_back(x,y);
    }
}

void BDSMagnetOuterFactoryPolesBase::CreateYokeAndContainerSolid(G4String name,
								 G4double length,
								 G4int    /*order*/,
								 G4double magnetContainerLength)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif  
  
  // circular yoke so pretty easy
  yokeSolid = new G4Tubs(name + "_yoke_solid",      // name
			 yokeStartRadius,           // start radius
			 yokeFinishRadius,          // finish radius
			 length*0.5 - lengthSafety, // z half length
			 0,                         // start angle
			 CLHEP::twopi);             // sweep angle

  poleIntersectionSolid = new G4Tubs(name + "_yoke_intersection_solid", // name
				     0,                                 // start radius
				     yokeStartRadius - lengthSafetyLarge,
				     length,  // long half length for unamibiguous intersection
				     0,                                 // start angle
				     CLHEP::twopi);                     // sweep angle
  allSolids.push_back(poleIntersectionSolid);
  

  // note container must have hole in it for the beampipe to fit in!
  // poled geometry doesn't fit tightly to beampipe so can alays use a circular aperture
  containerSolid = new G4Tubs(name + "_outer_container_solid", // name
			      poleStartRadius,                 // start radius
			      yokeFinishRadius + lengthSafety, // finish radius
			      length*0.5,                      // z half length
			      0,                               // start angle
			      CLHEP::twopi);                   // sweep angle
  allSolids.push_back(yokeSolid);

  // magnet container radius calculated when poles are calculated and assigned to
  // BDSMagnetOuterFactoryBase::magnetContainerRadius
  BuildMagnetContainerSolidStraight(name, magnetContainerLength, magnetContainerRadius); 
}

void BDSMagnetOuterFactoryPolesBase::IntersectPoleWithYoke(G4String name,
							   G4double /*length*/,
							   G4int  /*order*/)
{ 
  // cut pole here with knowledge of yoke shape.
  poleSolid = new G4IntersectionSolid(name + "_pole_solid",   // name
				      poleSolid,              // solid a
				      poleIntersectionSolid); // solid b
}

void BDSMagnetOuterFactoryPolesBase::CreateLogicalVolumes(G4String    name,
							  G4Colour*   colour,
							  G4Material* outerMaterial)
{
  BDSMagnetOuterFactoryBase::CreateLogicalVolumes(name, colour, outerMaterial);
  CreateLogicalVolumesCoil(name);
}

void BDSMagnetOuterFactoryPolesBase::CreateLogicalVolumesCoil(G4String name)
{
  if (coilLeftSolid)
    { // if one exists, all coil solids exist
      G4Material* coilMaterial = BDSMaterials::Instance()->GetMaterial("copper");
      coilLeftLV = new G4LogicalVolume(coilLeftSolid,
				       coilMaterial,
				       name + "_coil_left_lv");

      coilRightLV = new G4LogicalVolume(coilRightSolid,
					coilMaterial,
					name + "_coil_right_lv");
      
      G4Colour* coil = BDSColours::Instance()->GetColour("coil");
      G4VisAttributes* coilVisAttr = new G4VisAttributes(*coil);
      coilVisAttr->SetVisibility(true);

      coilLeftLV->SetVisAttributes(coilVisAttr);
      coilRightLV->SetVisAttributes(coilVisAttr);
      allVisAttributes.push_back(coilVisAttr);
    }
}

void BDSMagnetOuterFactoryPolesBase::TestInputParameters(BDSBeamPipe* beamPipe,
							 G4double&    outerDiameter,
							 G4Material*& outerMaterial)// reference to a pointer
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  //function arguments by reference to they can be modified in place
  //check outer material is something
  if (!outerMaterial)
    {outerMaterial = BDSMaterials::Instance()->GetMaterial(BDSGlobalConstants::Instance()->EmptyMaterial());}

  // ensure box size is bigger than the beampipe
  if (beamPipe->ContainerIsCircular())
    {
      // if it's circular, just check radius
      if (outerDiameter < 2*(beamPipe->GetContainerRadius()) )
	{outerDiameter = 2*(beamPipe->GetContainerRadius()) + 1*CLHEP::mm;}
    }
  else
    {
      // it's not circular - have a look at extents
      // +ve - -ve
      G4double extentX = beamPipe->GetExtentX().second - beamPipe->GetExtentX().first;
      G4double extentY = beamPipe->GetExtentY().second - beamPipe->GetExtentY().first;
      if ( (outerDiameter < extentX) || (outerDiameter < extentY) )
	{
	  // outerDiameter isn't sufficient for range in x or y
	  outerDiameter = std::max(extentX,extentY) + 1*CLHEP::mm;
	}
    }
}

void BDSMagnetOuterFactoryPolesBase::TestCoilFractions(G4double& coilWidthFraction,
						       G4double& coilHeightFraction)
{
  const G4double  lowerLimit   = 0.05;
  const G4double  upperLimit   = 0.98;
  std::vector<G4double*> fractions = {&coilWidthFraction, &coilHeightFraction};
  G4String names[2] = {"coilWidthFraction", "coilHeightFraction"};
  for (G4int i = 0; i < 2; i++)
    {
      if ((*(fractions[i]) < lowerLimit) && (*(fractions[i]) >= 0))
	{
	  *(fractions[i]) = lowerLimit;
	  G4cout << names[i] << " is below lower limit - limiting to " << lowerLimit << G4endl;
	}
      else if (*(fractions[i]) > upperLimit)
	{
	  *(fractions[i]) = upperLimit;
	  G4cout << names[i] << " is above upper limit - limiting to " << upperLimit << G4endl;
	}
      else
	{continue;}
    }
}

void BDSMagnetOuterFactoryPolesBase::PlaceComponents(const G4String& name,
						     G4int    order)
{
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  // place the components inside the container
  yokePV = new G4PVPlacement((G4RotationMatrix *) nullptr, // no rotation
			     G4ThreeVector(),              // position
			     yokeLV,                       // lv to be placed
			     name + "_yoke_pv",            // name
			     containerLV,                  // mother lv to be placed in
			     false,                        // no boolean operation
			     0,                            // copy number
			     checkOverlaps);               // whether to check overlaps
  allPhysicalVolumes.push_back(yokePV);
  
  // place poles
  if (!buildPole)
    {return;}
  // else continue and place poles and coils
  G4PVPlacement* aPolePV     = nullptr;  
  for (G4int n = 0; n < 2*order; ++n)
    {
      G4RotationMatrix* rm  = new G4RotationMatrix();
      allRotationMatrices.push_back(rm);
      rm->rotateZ((n+0.5)*segmentAngle + CLHEP::pi*0.5);
      G4String pvName = name + "_pole_" + std::to_string(n) + "_pv";
      // poleTranslation is by default (0,0,0)
      aPolePV = new G4PVPlacement(rm,                 // rotation
				  poleTranslation,    // position
				  poleLV,             // logical volume
				  pvName,             // name
				  containerLV,        // mother lv to be placed in
				  false,              // no boolean operation
				  n,                  // copy number
				  checkOverlaps);     // check overlaps
      allPhysicalVolumes.push_back(aPolePV);
    }
}

void BDSMagnetOuterFactoryPolesBase::PlaceComponentsCoils(G4String name,
							  G4int    order)
{
  if (!buildPole)
    {return;}
  G4PVPlacement* coilLeftPV   = nullptr;
  G4PVPlacement* coilRightPV  = nullptr;
  G4PVPlacement* endCoilPV    = nullptr;
  G4RotationMatrix* endCoilRM = new G4RotationMatrix();
  endCoilRM->rotateX(CLHEP::halfpi);

  G4ThreeVector endCoilTranslation(0,coilCentreRadius,0.5*endPieceLength);
  
  for (G4int n = 0; n < 2*order; ++n)
    {
      // prepare a new rotation matrix - must be new and can't reuse the same one
      // as the placement doesn't own it - changing the existing one will affect all
      // previously placed objects
      G4RotationMatrix* rm  = new G4RotationMatrix();
      G4RotationMatrix* ecrm = new G4RotationMatrix(*endCoilRM);
      allRotationMatrices.push_back(rm);
      G4double rotationAngle = (n+0.5)*segmentAngle + CLHEP::pi*0.5;
      rm->rotateZ((n+0.5)*segmentAngle + CLHEP::pi*0.5);

      // nominally this should be around z axis, but I guess z is now y give
      // rotation already... could do with rotate about axis in future
      ecrm->rotateY(rotationAngle);

      G4String pvName = name + "_coil_" + std::to_string(n);
      coilLeftPV = new G4PVPlacement(rm,                 // rotation
				     (G4ThreeVector)0,   // position
				     coilLeftLV,         // logical volume
				     pvName + "_left_pv",// name      
				     containerLV,        // mother lv to be placed in
				     false,              // no boolean operation
				     n,                  // copy number
				     checkOverlaps);     // check overlaps
      
      coilRightPV = new G4PVPlacement(rm,                 // rotation
				      (G4ThreeVector)0,   // position
				      coilRightLV,        // logical volume
				      pvName + "_right_pv",// name      
				      containerLV,        // mother lv to be placed in
				      false,              // no boolean operation
				      n,                  // copy number
				      checkOverlaps);     // check overlaps

      G4ThreeVector placementOffset = G4ThreeVector(endCoilTranslation);
      placementOffset.rotateZ(rotationAngle);
      endCoilPV = new G4PVPlacement(ecrm,                    // rotation
				    placementOffset,         // position
				    endPieceCoilLV,          // logical volume
				    name + "_end_piece_coil_pv", // name      
				    endPieceContainerLV,     // mother lv to be placed in
				    false,                   // no boolean operation
				    n,                       // copy number
				    checkOverlaps);          // check overlaps
      endPiece->RegisterPhysicalVolume(endCoilPV);
      endPiece->RegisterRotationMatrix(ecrm);
      
      allPhysicalVolumes.push_back(coilLeftPV);
      allPhysicalVolumes.push_back(coilRightPV);
    }
}

void BDSMagnetOuterFactoryPolesBase::CreateEndPiece(const G4String& name)
{
  // container solid
  endPieceContainerSolid = new G4Tubs(name + "_end_container_solid", // name
				      endPieceInnerR,                // inner radius
				      endPieceOuterR,                // outer radius
				      endPieceLength*0.5,            // z half length
				      0,                             // start angle
				      CLHEP::twopi);                 // sweep angle
  
  // container lv
  G4Material* emptyMaterial = BDSMaterials::Instance()->GetMaterial(BDSGlobalConstants::Instance()->EmptyMaterial());
  endPieceContainerLV = new G4LogicalVolume(endPieceContainerSolid,     // solid
					    emptyMaterial,              // material
					    name + "_end_container_lv");// name
  
  // coil end piece solid
  G4TwoVector zOffsets(0,0); // the transverse offset of each plane from 0,0
  G4double zScale = 1;       // the scale at each end of the points = 1
  G4VSolid* endPieceCoilSolid = new G4ExtrudedSolid(name + "_end_coil_solid", // name
						    endPiecePoints, // transverse 2d coordinates
						    coilHeight*0.5, // z half length
						    zOffsets, zScale,
						    zOffsets, zScale);
  
  // coil end piece lv
  G4Material* copper = BDSMaterials::Instance()->GetMaterial("copper");
  endPieceCoilLV = new G4LogicalVolume(endPieceCoilSolid,      // solid
				       copper,                 // material
				       name + "_end_coil_lv"); // name
				       
  // vis attributes - copy the coil ones that are already built
  G4VisAttributes* endPieceCoilVis = new G4VisAttributes(*(coilLeftLV->GetVisAttributes()));
  endPieceCoilLV->SetVisAttributes(endPieceCoilVis);
  // container vis attributes are held in global constants - only need copy pointer
  endPieceContainerLV->SetVisAttributes(containerLV->GetVisAttributes());
  
  // user limits - don't register as using global one
  endPieceCoilLV->SetUserLimits(BDSGlobalConstants::Instance()->DefaultUserLimits());

  // package it all up
  endPiece = new BDSSimpleComponent(name + "_end_piece",
				    endPieceLength,
				    0,
				    endPieceContainerSolid,
				    endPieceContainerLV);

  endPiece->RegisterSolid(endPieceCoilSolid);
  endPiece->RegisterLogicalVolume(endPieceCoilLV);
  endPiece->RegisterVisAttributes(endPieceCoilVis);
  endPiece->RegisterSensitiveVolume(endPieceCoilLV);
  endPiece->SetExtent(BDSExtent(endPieceOuterR, endPieceOuterR, endPieceLength*0.5));
  endPiece->SetInnerExtent(BDSExtent(endPieceInnerR, endPieceInnerR, endPieceLength*0.5));
}

void BDSMagnetOuterFactoryPolesBase::DipoleCommonPreConstruction(BDSBeamPipe*    beamPipe,
								 const G4String& name,
								 const G4double& angleIn,
								 const G4double& angleOut,
								 const G4double& length,
								 G4double&       outerDiameter,
								 G4Material*&    material,
								 G4double&       vhRatio)
{
  CleanUp();
 
  // Test input parameters
  TestInputParameters(beamPipe, outerDiameter, material);

  // Test faces
  if (BDS::WillIntersect(-angleIn, -angleOut, outerDiameter, length))
    {
      G4cout << __METHOD_NAME__ << "Error: Faces of magnet (section) named \""
	     << name << "\" will overlap!" << G4endl;
      G4cout << "Length of magnet " << length << " mm"
	     << " is too short for the angle of the pole faces: (" << angleIn << "," << angleOut << ")." << G4endl;
      exit(1);
    }
  
  // vhRatio - don't allow a ratio greater than 10:1
  if (vhRatio > 10)
    {vhRatio = 10;}
  else if (vhRatio < 0.1)
    {vhRatio = 0.1;}
}

void BDSMagnetOuterFactoryPolesBase::DipoleCalculations(const G4bool&      hStyle,
							const G4bool&      buildVertically,
							const BDSBeamPipe* beamPipe,
							const G4double&    length,
							const G4double&    outerDiameter,
							const G4double&    angleIn,
							const G4double&    angleOut,
							const G4double&    yokeThicknessFraction,
							const G4double&    vhRatio,
							const G4double&    coilWidthFraction,
							const G4double&    coilHeightFraction,
							G4double& bpHalfWidth,
							G4double& bpHalfHeight,
							G4double& poleHalfGap,
							G4double& poleWidth,
							G4double& poleHeight,
							G4double& yokeWidth,
							G4double& yokeHalfHeight,
							G4double& yokeThickness,
							G4double& yokeOverHang,
							G4double& coilWidth,
							G4double& coilHeightIn,
							G4double& coilToYokeGap,
							G4double& coilToPoleGap,
							G4double& sLength,
							G4double& containerSLength,
							G4double& intersectionRadius)
{
  // calculate any geometrical parameters
  bpHalfWidth  = beamPipe->GetExtent().MaximumX();
  bpHalfHeight = beamPipe->GetExtent().MaximumY();
  
  poleWidth = 2*(buildVertically ? bpHalfHeight : bpHalfWidth)  + 2*lengthSafetyLarge;
  poleWidth = std::max(poleWidth, outerDiameter*0.36);
  // in the case of a very wide beam pipe, we can't build a pole that matches
  if (poleWidth > 0.9*outerDiameter)
    {poleWidth = outerDiameter*0.36;}
  
  poleHalfGap = (buildVertically ? bpHalfWidth : bpHalfHeight) + lengthSafetyLarge;

  // propose outer dimensions.
  yokeWidth      = outerDiameter;
  yokeHalfHeight = 0.5 * outerDiameter * vhRatio;
  
  // ensure outer edges aren't smaller than beam pipe
  const G4double margin = 50*CLHEP::um; // minimum allowable 'yoke'
  G4double yokeWidthLowerLimit      = 2*(buildVertically ? bpHalfHeight : bpHalfWidth) + margin;
  G4double yokeHalfHeightLowerLimit = (buildVertically ? bpHalfWidth  : bpHalfHeight)+ margin;
  if (yokeWidth <= yokeWidthLowerLimit)
    {yokeWidth = yokeWidthLowerLimit + margin;}
  if (yokeHalfHeight <= yokeHalfHeightLowerLimit)
    {yokeHalfHeight = yokeHalfHeightLowerLimit + margin;}
  
  // propose yoke thickness
  yokeThickness = yokeThicknessFraction * outerDiameter;
  // if there's not enough space (given the beam pipe and outer edge)
  // for the specified fraction of yoke, don't build pole. Also coerce
  // yoke thickness.
  if (yokeThickness > yokeHalfHeight - poleHalfGap)
    {
      buildPole     = false;
      yokeThickness = yokeHalfHeight - poleHalfGap;
    }

  // don't build pole if there's not enough room - coerce yoke thickness
  G4double factor = hStyle ? 2.0 : 1.0; // 2x thickness for h style
  if (factor * yokeThickness > yokeWidth - poleWidth)
    {
      buildPole = false;
      yokeThickness = (yokeWidth - poleWidth) / factor;
    }
  // don't build pole if yoke tight around beam pipe - coerce yoke thickness
  G4double test = buildVertically ? bpHalfWidth : bpHalfHeight;
  if (factor * yokeThickness > yokeWidth - 2*test)
    {
      buildPole = false;
      yokeThickness = (yokeWidth - 2*test) / factor;
    }
  
  if (buildPole)
    {poleHeight = yokeHalfHeight - yokeThickness - poleHalfGap;}
  else
    {poleHeight = 0;} 

  // ensure the yoke isn't too thick - choose smaller of two limits
  G4double yokeThicknessLimitHorizontal = yokeWidth - yokeWidthLowerLimit - margin*0.5;
  G4double yokeThicknessLimitVertical   = yokeHalfHeight - yokeHalfHeightLowerLimit - margin;
  G4double yokeThicknessLimit           = std::min(yokeThicknessLimitHorizontal, 2*yokeThicknessLimitVertical);
  if (yokeThickness > yokeThicknessLimit)
    {yokeThickness = yokeThicknessLimit;}
  if (yokeThickness < margin)
    {yokeThickness = margin;} // ensure minimum width of yoke
  
  // prevent negative coil widths by yoke becoming too wide in the case
  // of a wide pole
  if (yokeThickness + poleWidth > outerDiameter - margin)
    {yokeThickness = outerDiameter - poleWidth - margin;}

  if (hStyle)
    {yokeOverHang = 0.5*yokeWidth - yokeThickness - 0.5*poleWidth;}
  else
    {yokeOverHang = yokeWidth - yokeThickness - poleWidth;}

  // must ensure that:
  // yoke length < outer container length < full magnet container length
  // whether straight or angled
  // yoke full length -> length - 2*lsl
  // outer container full length -> length - 2*ls
  // full magnet container full length -> container length
  // if we have angled faces, make square faced solids longer for intersection.
  if (BDS::IsFinite(angleIn) || BDS::IsFinite(angleOut))
    {
      // In the case of angled faces, calculate a length so that the straight solids
      // used in intersection are long enough to reach the edges of the angled faces.
      // Could simply do 2x length, but for short dipole sections with strongly angled
      // faces this doesn't work. Calculate extent along z for each angled face. This
      // is called the 'safe' length -> sLength
      G4double hypotenuse = std::hypot(yokeWidth, yokeHalfHeight);
      G4double dzIn  = tan(std::abs(angleIn))  * 1.2*hypotenuse; // 20% over estimation for safety
      G4double dzOut = tan(std::abs(angleOut)) * 1.2*hypotenuse;
      // take the longest of different estimations (2x and 1.5x + dZs)
      sLength = std::max(2*length, 1.5*length + dzIn + dzOut);
      containerSLength = sLength;
    }

  intersectionRadius = std::hypot(0.5*poleWidth + yokeOverHang, poleHalfGap + poleHeight);
  // if finite thickness yoke (independent of overall size)
  if (yokeThickness > 0.05*outerDiameter) // nicely round edges of outer side without cutting inner
    {intersectionRadius += 0.8 * std::hypot(yokeThickness, yokeThickness);}
  else
    {intersectionRadius *= 1.3;} // some margin

  coilWidth = yokeOverHang * coilWidthFraction;
  coilHeightIn = poleHeight * coilHeightFraction;
  coilToYokeGap = (poleHeight - coilHeightIn) * 0.5;
  coilToPoleGap = lengthSafetyLarge;
}

std::vector<G4ThreeVector> BDSMagnetOuterFactoryPolesBase::CalculateCoilDisplacements(G4double  poleHalfWidthIn,
										      G4double  poleHalfGapIn,
										      G4double  coilWidthIn,
										      G4double  coilHeightIn,
										      G4double  cDY,
										      G4double& coilDY)
{
  // T = top, B = bottom, L = left, R = right
  G4double coilDX = poleHalfWidthIn  + 0.5*coilWidthIn  + lengthSafetyLarge;
  coilDY = poleHalfGapIn + 0.5*coilHeightIn + lengthSafetyLarge + cDY;
  
  G4ThreeVector coilTLDisp = G4ThreeVector( coilDX, coilDY, 0);
  G4ThreeVector coilTRDisp = G4ThreeVector(-coilDX, coilDY, 0);
  G4ThreeVector coilBLDisp = G4ThreeVector( coilDX,-coilDY, 0);
  G4ThreeVector coilBRDisp = G4ThreeVector(-coilDX,-coilDY, 0);
  std::vector<G4ThreeVector> coilDisps;
  coilDisps.push_back(coilTLDisp);
  coilDisps.push_back(coilTRDisp);
  coilDisps.push_back(coilBLDisp);
  coilDisps.push_back(coilBRDisp);
  return coilDisps;
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateDipoleC(G4String     name,
							      G4double     length,
							      BDSBeamPipe* beamPipe,
							      G4double     outerDiameter,
							      G4double     containerLength,
							      G4double     angleIn,
							      G4double     angleOut,
							      G4Material*  material,
							      G4bool       yokeOnLeft,
							      G4Colour*    colour,
							      G4bool       buildVertically,
							      G4bool       buildEndPiece,
							      G4double     vhRatio,
							      G4double     coilWidthFraction,
							      G4double     coilHeightFraction)
{
  DipoleCommonPreConstruction(beamPipe, name, angleIn, angleOut, length, outerDiameter,
			      material, vhRatio);
  TestCoilFractions(coilWidthFraction, coilHeightFraction);

  // 1 calculations
  // 2 c shaped solid
  // 3 angled solids for intersection
  // 4 intersection solids
  // 5 logical volumes
  // 6 placement
  // general order - yoke, container, magnet container, coils

  G4double bpHalfWidth         = 0;
  G4double bpHalfHeight        = 0;
  G4double poleHalfGap         = 0;
  G4double poleWidth           = 0;
  G4double poleHeight          = 0;
  G4double yokeWidth           = 0;
  G4double yokeHalfHeight      = 0;
  G4double yokeThickness       = 0;
  G4double yokeOverHang        = 0;
  G4double coilWidth           = 0;
  // coilHeight is class member
  G4double coilToYokeGap       = 0;
  G4double coilToPoleGap       = 0;
  G4double sLength             = length; // 'safe' length fo intersection - default is normal length
  G4double containerSLength    = containerLength; // similarly for the container
  G4double intersectionRadius  = 0;
  DipoleCalculations(false, buildVertically, beamPipe, length, outerDiameter, angleIn,
		     angleOut, 0.23, vhRatio, coilWidthFraction, coilHeightFraction,
		     bpHalfWidth, bpHalfHeight, poleHalfGap, poleWidth, poleHeight,
		     yokeWidth, yokeHalfHeight, yokeThickness, yokeOverHang, coilWidth,
		     coilHeight, coilToYokeGap, coilToPoleGap,
		     sLength, containerSLength, intersectionRadius);
  
  //G4double yokeInsideX = 0.5*poleWidth - yokeWidth + yokeThickness;
  G4double yokeInsideX = -0.5*poleWidth - yokeOverHang;

  // vertical offset of coil from pole tip
  G4double cDY = (poleHeight - coilHeight)*0.5;

  // coil displacements
  std::vector<G4ThreeVector> coilDisps;
  G4double coilDY = 0;
  if (buildPole)
    {
      coilDisps= CalculateCoilDisplacements(0.5*poleWidth, poleHalfGap, coilWidth,
					    coilHeight, cDY, coilDY);
    }

  // create vectors of points for extruded solids
  // create yoke + pole (as one solid about 0,0,0)
  std::vector<G4TwoVector> yokePoints;
  std::vector<G4TwoVector> cPoints;  // container points (for magnet outer only)
  std::vector<G4TwoVector> mCPoints; // magnet container points

  // variables for extents
  G4double extXPos = 0;
  G4double extXNeg = 0;
  G4double extYPos = 0;
  G4double extYNeg = 0;
  // Typically we have a positive bend angle that (by convention) causes a
  // bend to the -ve x direction in right hadned coordinates. Also, typically,
  // a C shaped magnet has the yoke to the inside so there is an aperture for
  // any radiation to escape to the outside. Therefore, we build the yoke like this
  // and flip it if required. Points are done in clock wise order from the bottom left
  // corner of the top pole.
  const G4double lsl = lengthSafetyLarge;   // shortcut
  const G4double phw = 0.5*poleWidth;       // shortcut
  const G4double phg = poleHalfGap;         // shortcut
  const G4double yw  = yokeWidth; // shortcut
  const G4double yhh = yokeHalfHeight;   // shortcut
  yokePoints.emplace_back(phw - lsl,       phg + lsl);
  yokePoints.emplace_back(phw - lsl,       yhh - lsl);
  yokePoints.emplace_back(phw + lsl - yw,  yhh - lsl);
  yokePoints.emplace_back(phw + lsl - yw, -yhh + lsl);
  yokePoints.emplace_back(phw - lsl,      -yhh + lsl);
  yokePoints.emplace_back(phw - lsl,      -phg - lsl);
  if (buildPole)
    {
      yokePoints.emplace_back(-phw + lsl,         -phg - lsl);
      yokePoints.emplace_back(-phw + lsl,         -yhh + yokeThickness);
      yokePoints.emplace_back( yokeInsideX - lsl, -yhh + yokeThickness);
      yokePoints.emplace_back( yokeInsideX - lsl,  yhh - yokeThickness);
      yokePoints.emplace_back(-phw + lsl,          yhh - yokeThickness);
      yokePoints.emplace_back(-phw + lsl,          phg + lsl);
    }
  else
    {
      yokePoints.emplace_back(yokeInsideX - lsl, -phg - lsl);
      yokePoints.emplace_back(yokeInsideX - lsl,  phg + lsl);
    }
  
  // points for container for magnet outer only
  if (buildPole == false) // redundant point when buildPole == true
    {cPoints.emplace_back(phw + lsl, phg);}
  else
    {
      cPoints.emplace_back(phw + coilWidth + 2*lsl, phg);
      cPoints.emplace_back(phw + coilWidth + 2*lsl, phg + coilHeight + 2*lsl + cDY);
      cPoints.emplace_back(phw + lsl, phg + coilHeight + 2*lsl + cDY);
    }
  cPoints.emplace_back(phw + lsl,       yhh + lsl);
  cPoints.emplace_back(phw - yw - lsl,  yhh + lsl);
  cPoints.emplace_back(phw - yw - lsl, -yhh - lsl);
  cPoints.emplace_back(phw + lsl,      -yhh - lsl);
  if (buildPole)
    {
      cPoints.emplace_back(phw + lsl, -phg - coilHeight - 2*lsl - cDY);
      cPoints.emplace_back(phw + coilWidth + 2*lsl, -phg - coilHeight - 2*lsl - cDY);
      cPoints.emplace_back(phw + coilWidth + 2*lsl, -phg);
    }
  else // redundant point when buildPole == true
    {cPoints.emplace_back(phw + lsl, -phg);}
  cPoints.emplace_back(yokeInsideX, -phg);
  cPoints.emplace_back(yokeInsideX,  phg);

  // points for container for full magnet object including beam pipe
  // first one in y here is -lsl to cancel +lsl to phg originally
  // container can be same height as beam pipe as it's always wider
  G4double maxLeft = std::max(phw, bpHalfWidth);
  mCPoints.emplace_back(maxLeft + lsl, phg - lsl);
  if (buildPole)
    {
      mCPoints.emplace_back(phw + coilWidth + 4*lsl, phg - lsl);
      mCPoints.emplace_back(phw + coilWidth + 4*lsl, phg + coilHeight + 3*lsl + cDY);
      mCPoints.emplace_back(phw + 2*lsl, phg + coilHeight + 3*lsl + cDY);
    }
  mCPoints.emplace_back(phw + 2*lsl,       yhh + 2*lsl);
  mCPoints.emplace_back(phw - yw - 2*lsl,  yhh + 2*lsl);
  mCPoints.emplace_back(phw - yw - 2*lsl, -yhh - 2*lsl);
  mCPoints.emplace_back(phw + 2*lsl,      -yhh - 2*lsl);
  if (buildPole)
    {
      mCPoints.emplace_back(phw + 2*lsl, -phg - coilHeight - 3*lsl - cDY);
      mCPoints.emplace_back(phw + coilWidth + 4*lsl, -phg - coilHeight - 3*lsl - cDY);
      mCPoints.emplace_back(phw + coilWidth + 4*lsl, -phg + lsl);
    }
  mCPoints.emplace_back(maxLeft + lsl, -phg + lsl);
  
  // extents
  extXPos = phw + lsl;
  extXNeg = phw - yw - 2*lsl;
  extYPos = yhh + lsl;
  extYNeg = -(yhh +lsl);
  if (buildPole)
    {extXPos += coilWidth + lsl;}

  if (yokeOnLeft)
    {
      // flip x component so geometry is reflected horizontally
      for (auto& vec : yokePoints)
	{vec.setX(vec.x() * -1);}
      // reverse order so it's clockwise as geant4 requires
      std::reverse(yokePoints.begin(), yokePoints.end());

      // same for container
      for (auto& vec : cPoints)
	{vec.setX(vec.x() * -1);}
      std::reverse(cPoints.begin(), cPoints.end());

      // same for magnet container
      for (auto& vec : mCPoints)
	{vec.setX(vec.x() * -1);}
      std::reverse(mCPoints.begin(), mCPoints.end());

      // flip extents
      std::swap(extXPos, extXNeg);
    }
  
  // rotate if building vertically
  if (buildVertically)
    {
      for (auto& point : yokePoints)
	{point = BDS::Rotate(point, CLHEP::halfpi);}
      for (auto& point : cPoints)
	{point = BDS::Rotate(point, CLHEP::halfpi);}
      for (auto& point : mCPoints)
	{point = BDS::Rotate(point, CLHEP::halfpi);}
      // 'rotate' extents too
      std::swap(extXPos, extYPos);
      std::swap(extXNeg, extYNeg);
    }
  BDSExtent ext = BDSExtent(extXNeg, extXPos, extYNeg, extYPos,
			    -length*0.5, length*0.5);
  magContExtent = ext; // copy to container variable - basically the same
  
  G4TwoVector zOffsets(0,0); // the transverse offset of each plane from 0,0
  G4double zScale = 1;       // the scale at each end of the points = 1
  yokeSolid = new G4ExtrudedSolid(name + "_yoke_sq_solid",// name
				  yokePoints,             // transverse 2d coordinates
				  sLength*0.5 - lsl,      // z half length
				  zOffsets, zScale,       // dx,dy offset for each face, scaling
				  zOffsets, zScale);      // dx,dy offset for each face, scaling

  // container for magnet outer
  containerSolid = new G4ExtrudedSolid(name + "_outer_sq_container_solid", // name
				       cPoints,                    // transverse 2d coordinates
				       sLength*0.5 - lengthSafety, // z half length
				       zOffsets, zScale, // dx,dy offset for each face, scaling
				       zOffsets, zScale);// dx,dy offset for each face, scaling

  // container for full magnet
  magnetContainerSolid = new G4ExtrudedSolid(name + "_sq_container_solid", // name
					     mCPoints,                     // transverse 2d coordinates
					     containerSLength*0.5,         // z half length
					     zOffsets, zScale, // dx,dy offset for each face, scaling
					     zOffsets, zScale);// dx,dy offset for each face, scaling

  // register existing square solids as we're going to overwrite them with intersected ones
  allSolids.push_back(yokeSolid);
  allSolids.push_back(containerSolid);
  allSolids.push_back(magnetContainerSolid);
  
  return DipoleCommonConstruction(name, outerDiameter, buildEndPiece, coilWidth, length,
                                  containerLength, sLength, angleIn, angleOut,
                                  colour, material,
                                  coilDisps, buildVertically, ext, 0.5*poleWidth, poleHalfGap,
                                  cDY, coilDY, intersectionRadius);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::CreateDipoleH(G4String     name,
							      G4double     length,
							      BDSBeamPipe* beamPipe,
							      G4double     outerDiameter,
							      G4double     containerLength,
							      G4double     angleIn,
							      G4double     angleOut,
							      G4Material*  material,
							      G4Colour*    colour,
							      G4bool       buildVertically,
							      G4bool       buildEndPiece,
							      G4double     vhRatio,
							      G4double     coilWidthFraction,
							      G4double     coilHeightFraction)
{
  DipoleCommonPreConstruction(beamPipe, name, angleIn, angleOut, length, outerDiameter, material, vhRatio);
  TestCoilFractions(coilWidthFraction, coilHeightFraction);
    
  // 1 calculations
  // 2 h shaped solid
  // 3 angled solids for intersection
  // 4 intersection solids
  // 5 logical volumes
  // 6 placement
  // general order - yoke, container, magnet container, coils
  
  G4double bpHalfWidth    = 0;
  G4double bpHalfHeight   = 0;
  G4double poleHalfGap    = 0;
  G4double poleWidth      = 0;
  G4double poleHeight     = 0;
  G4double yokeWidth      = 0;
  G4double yokeHalfHeight = 0;
  G4double yokeThickness  = 0;
  G4double yokeOverHang   = 0;
  G4double coilWidth      = 0;
  // coilHeight is member variable
  G4double coilToYokeGap  = 0;
  G4double coilToPoleGap  = 0;
  G4double sLength        = length; // 'safe' length fo intersection - default is normal length
  G4double containerSLength = containerLength; // similarly for the container
  G4double intersectionRadius = 0;
  DipoleCalculations(true, buildVertically, beamPipe, length, outerDiameter, angleIn,
		     angleOut, 0.12, vhRatio, coilWidthFraction, coilHeightFraction,
		     bpHalfWidth, bpHalfHeight, poleHalfGap, poleWidth, poleHeight,
		     yokeWidth, yokeHalfHeight, yokeThickness, yokeOverHang,
		     coilWidth, coilHeight, coilToYokeGap, coilToPoleGap,
		     sLength, containerSLength, intersectionRadius);

  // distance from axis to inside of yoke horizontally
  G4double yokeInsideX = 0.5*poleWidth + yokeOverHang;
  //yokeWidth*0.5 - yokeThickness;
  //G4double yokeInsideY = yokeHalfHeight - yokeThickness;
  G4double yokeInsideY = poleHalfGap + poleHeight;
  
  // Vertical offset of coil from pole tip
  G4double cDY = (poleHeight - coilHeight)*0.5;

  const G4double lsl = lengthSafetyLarge; // shortcut
  
  // coil displacements
  std::vector<G4ThreeVector> coilDisps;
  G4double coilDY = 0;
  if (buildPole)
    {
      coilDisps = CalculateCoilDisplacements(0.5*poleWidth, poleHalfGap,
					     coilWidth, coilHeight, cDY, coilDY);
    }

  // these may be used later so need to be outside if statement below
  // create vectors of points for extruded solids
  std::vector<G4TwoVector> yokePoints;
  std::vector<G4TwoVector> cPoints;  // container points (for 'magnet outer' inner only)
  std::vector<G4TwoVector> mCPoints; // magnet container points

  // variables for extents
  G4double extXPos = 0;
  G4double extXNeg = 0;
  G4double extYPos = 0;
  G4double extYNeg = 0;

  G4VSolid* yokeInnerSolid = nullptr;
  if (buildPole)
    {// use extruded solid for inner yoke shape and box for out with subtraction
      // create yoke + pole (as one solid about 0,0,0)
      // points are done in clock wise order from the bottom right corner of the top pole.
      const G4double phw = 0.5*poleWidth;
      const G4double phg = poleHalfGap;
      yokePoints.emplace_back( phw - lsl,  phg + lsl);
      yokePoints.emplace_back( phw - lsl,  yokeInsideY    + lsl);
      yokePoints.emplace_back( yokeInsideX   + lsl,  yokeInsideY    + lsl);
      yokePoints.emplace_back( yokeInsideX   + lsl, -yokeInsideY    - lsl);
      yokePoints.emplace_back( phw - lsl, -yokeInsideY    - lsl);
      yokePoints.emplace_back( phw - lsl, -phg - lsl);
      yokePoints.emplace_back(-phw + lsl, -phg - lsl);
      yokePoints.emplace_back(-phw + lsl, -yokeInsideY    - lsl);
      yokePoints.emplace_back(-yokeInsideX   - lsl, -yokeInsideY    - lsl);
      yokePoints.emplace_back(-yokeInsideX   - lsl,  yokeInsideY    + lsl);
      yokePoints.emplace_back(-phw + lsl,  yokeInsideY    + lsl);
      yokePoints.emplace_back(-phw + lsl,  phg + lsl);

      // rotate if building vertically
      if (buildVertically)
	{
	  for (auto& point : yokePoints)
	    {point = BDS::Rotate(point, CLHEP::halfpi);}
	}

      G4TwoVector zOffsets(0,0); // the transverse offset of each plane from 0,0
      G4double zScale = 1;       // the scale at each end of the points = 1
      yokeInnerSolid = new G4ExtrudedSolid(name + "_yoke_inner_solid", // name
					   yokePoints,                 // transverse 2d coordinates
					   sLength,                    // z half length
					   zOffsets, zScale, // dx,dy offset for each face, scaling
					   zOffsets, zScale);// dx,dy offset for each face, scaling
      // note 1.0x length > 0.5 length for unambiguous subtraction
    }
  else // if build pole
    {// don't build pole - just build yoke -> box
      G4double yIX = buildVertically ? yokeInsideY : yokeInsideX;
      G4double yIY = buildVertically ? yokeInsideX : yokeInsideY;
      yokeInnerSolid = new G4Box(name + "_yoke_inner_solid", // name
				 yIX + lsl,                  // x half length
				 yIY + lsl,                  // y half length
				 sLength);                   // z half length
      // note 1.0x length > 0.5 length for unambiguous subtraction
    }

  // extents
  extXPos =  0.5*yokeWidth + lsl;
  extXNeg = -0.5*yokeWidth - lsl;
  extYPos =  0.5*yokeWidth + lsl;
  extYNeg = -0.5*yokeWidth - lsl;
  
  if (buildVertically)
    {// 'rotate' extents too
      std::swap(extXPos, extYPos);
      std::swap(extXNeg, extYNeg);
    }
  BDSExtent ext = BDSExtent(extXNeg, extXPos, extYNeg, extYPos,
			    -length*0.5, length*0.5);
  magContExtent = ext; // copy to container variable - basically the same

  G4double yokeOuterX = 0.5*yokeWidth - lsl;
  G4double yokeOuterY = yokeHalfHeight - lsl;
  G4double yOX = buildVertically ? yokeOuterY : yokeOuterX;
  G4double yOY = buildVertically ? yokeOuterX : yokeOuterY;
  G4VSolid* yokeOuterSolid = new G4Box(name + "_yoke_outer_solid", // name
				       yOX,                        // x half length
				       yOY,                        // y half length
				       sLength * 0.5 - lsl);       // z half length
  
  yokeSolid = new G4SubtractionSolid(name + "_yoke_solid", // name,
				     yokeOuterSolid,       // this
				     yokeInnerSolid);      // minus this
  
  // container for magnet outer 
  G4double containerdx = buildVertically ? bpHalfWidth : yokeInsideX;
  G4double containerdy = buildVertically ? yokeInsideX : bpHalfHeight;

  // full length for unambiguous subtraction
  G4VSolid* containerInnerSolid = new G4Box(name + "_container_inner_solid", // name
					    containerdx, containerdy, sLength);
  G4double cOX = buildVertically ? yokeHalfHeight : 0.5*yokeWidth;
  G4double cOY = buildVertically ? 0.5*yokeWidth : yokeHalfHeight;
  G4VSolid* containerOuterSolid = new G4Box(name + "_container_outer_solid", // name
					    cOX,                             // x half length
					    cOY,                             // y half length
					    sLength * 0.5);                  // z half length

  containerSolid = new G4SubtractionSolid(name + "_outer_sq_container_solid", // name
					  containerOuterSolid,  // this
					  containerInnerSolid); // minus this

  // container for full magnet
  magnetContainerSolid = new G4Box(name + "_sq_container_solid", // name
				   cOX + lsl,                    // x half length
				   cOY + lsl,                    // y half length
				   containerSLength * 0.5);      // z half length

  // register existing square solids as we're going to overwrite them with intersected ones
  allSolids.push_back(yokeInnerSolid);
  allSolids.push_back(yokeOuterSolid);
  allSolids.push_back(yokeSolid);
  allSolids.push_back(containerInnerSolid);
  allSolids.push_back(containerOuterSolid);
  allSolids.push_back(containerSolid);
  allSolids.push_back(magnetContainerSolid);

  return DipoleCommonConstruction(name, outerDiameter, buildEndPiece, coilWidth, length,
				  containerLength, sLength, angleIn, angleOut,
				  colour, material,
				  coilDisps, buildVertically, ext, poleWidth*0.5, poleHalfGap,
				  cDY, coilDY, intersectionRadius);
}

BDSMagnetOuter* BDSMagnetOuterFactoryPolesBase::DipoleCommonConstruction(G4String    name,
									 G4double    outerDiameter,
									 G4bool      buildEndPiece,
									 G4double    coilWidth,
									 G4double    length,
									 G4double    containerLength,
									 G4double    sLength,
									 G4double    angleIn,
									 G4double    angleOut,
									 G4Colour*   colour,
									 G4Material* material,
									 std::vector<G4ThreeVector>& coilDisps,
									 G4bool      buildVertically,
									 BDSExtent&  ext,
									 G4double    poleHalfWidth,
									 G4double    poleHalfGap,
									 G4double    cDY,
									 G4double    coilDY,
									 G4double    intersectionRadius)
{
  const G4double lsl = lengthSafetyLarge; // shortcut
  
  // create coil - one solid that will be placed 4 times... or
  // if we have angled faces, they're all unique unfortunately so use a vector of solids
  // and place individually
  G4VSolid* coilSolid = nullptr;
  std::vector<G4VSolid*>        coilsSolids;
  std::vector<G4LogicalVolume*> coilLVs;
  G4bool individualCoilsSolids = false;
  if (buildPole)
    {
      coilSolid = new G4Box(name + "_coil_solid",   // name
			    coilWidth*0.5,          // x half width
			    coilHeight*0.5,         // y half height
			    sLength*0.5 - lsl);     // z half length - same as yoke
      allSolids.push_back(coilSolid);
    }

  // Intersect and replace solids. Do it via replacmeent of the base class member G4VSolid*
  // as the intersection is only done if one of the angles is finite.

  // If we have angled faces, we need to also intersect the outer magnet container and
  // we therefore have an extra solid that should belong to that geometry component
  // and not the magnet outer - create a vector for holding any extra solids. These can
  // only be registered once the container object is created and that can only happen
  // once the solid and logical volumes are fully complete.
  std::vector<G4VSolid*> magnetContainerExtraSolids;
  
  if (BDS::IsFinite(angleIn) || BDS::IsFinite(angleOut))
    { 
      std::pair<G4ThreeVector,G4ThreeVector> faces = BDS::CalculateFaces(angleIn,angleOut);
      inputFaceNormal  = faces.first;
      outputFaceNormal = faces.second;
      
      // angled solid for magnet outer and coils
      G4VSolid* angledFaces = new G4CutTubs(name + "_angled_face_solid", // name
					    0,                           // inner radius
					    intersectionRadius,          // outer radius
					    length*0.5 - lsl,            // z half length
					    0,                           // start angle
					    CLHEP::twopi,                // sweep angle
					    inputFaceNormal,             // input face normal
					    outputFaceNormal);           // output face normal
      
      // angled solid for magnet outer container volume
      G4VSolid* angledFacesCont = new G4CutTubs(name + "_angled_face_cont_solid", // name
						0,                           // inner radius
						intersectionRadius,          // outer radius
						length*0.5 - lengthSafety,   // z half length
						0,                           // start angle
						CLHEP::twopi,                // sweep angle
						inputFaceNormal,             // input face normal
						outputFaceNormal);           // output face normal
      
      // angled solid for full magnet container
      G4VSolid* angledFacesMagCont = new G4CutTubs(name + "_angled_face_mag_cont_solid", // name
						   0,                    // inner radius
						   intersectionRadius,   // outer radius
						   containerLength*0.5,  // z half length
						   0,                    // start angle
						   CLHEP::twopi,         // sweep angle
						   inputFaceNormal,      // input face normal
						   outputFaceNormal);    // output face normal

      // register angled solids
      allSolids.push_back(angledFaces);
      allSolids.push_back(angledFacesCont);
      allSolids.push_back(angledFacesMagCont);
      
      // now do intersections overwriting existing pointers
      yokeSolid = new G4IntersectionSolid(name + "_yoke_solid", // name
					  yokeSolid,            // solid a - existing yoke
					  angledFaces);         // solid b
      
      containerSolid = new G4IntersectionSolid(name + "_outer_container_solid", // name
					       containerSolid,                  // solid a
					       angledFacesCont);                // solid b

      magnetContainerSolid = new G4IntersectionSolid(name + "_container_solid", // name
						     magnetContainerSolid,      // solid a
						     angledFacesMagCont);       // solid b
      
      individualCoilsSolids = true; // flag that we have individual coil solids
      if (buildPole)
        {
	  for (G4int i = 0; i < 4; i++)
            {
	      G4VSolid *coilS = new G4IntersectionSolid(name + "_pole_solid_" + std::to_string(i), // name
							angledFaces,                // solid a
							coilSolid,                  // solid b
							(G4RotationMatrix *) nullptr, // 0 rotation
							coilDisps[i]);              // translation
	      coilsSolids.push_back(coilS);
	      allSolids.push_back(coilS);
            }
        }
    }
  
  // logical volumes
  CreateLogicalVolumes(name, colour, material);
  // we only use one coil solid here so do that here
  G4LogicalVolume* coilLV = nullptr;
  if (buildPole)
    {
      G4Colour* coil = BDSColours::Instance()->GetColour("coil");
      G4VisAttributes* coilVis = new G4VisAttributes(*coil);
      coilVis->SetVisibility(true);
      allVisAttributes.push_back(coilVis);
      G4Material* coilMaterial = BDSMaterials::Instance()->GetMaterial("copper");
      if (individualCoilsSolids)
	{
	  for (G4int i = 0; i < 4; i++)
	    {
	      G4String theName =  name + "_coil_solid_" + std::to_string(i);
	      G4LogicalVolume* aCoilLV = new G4LogicalVolume(coilsSolids[i],
							     coilMaterial,
							     theName);
	      aCoilLV->SetVisAttributes(coilVis);
	      coilLVs.push_back(aCoilLV);
	      allLogicalVolumes.push_back(aCoilLV);
	    }
	}
      else
	{
	  coilLV = new G4LogicalVolume(coilSolid,
				       coilMaterial,
				       name + "_coil_lv");
	  coilLV->SetVisAttributes(coilVis);
	  allLogicalVolumes.push_back(coilLV);
	}
    }
  
  // placement
  // place yoke+pole (one solid) together
  G4bool buildPoleOriginal = buildPole;
  buildPole = false; // this will avoid PlaceComponents trying to place poles that don't exist
  PlaceComponents(name, 1); // places yoke
  buildPole = buildPoleOriginal;
  // place coils
  if (buildPole)
    {
      G4PVPlacement* aCoilPV = nullptr;
      // coils intersection was done such that the coordinates of the resultant solid
      // are the same as the plain coil G4Box. This allows us to programmatically place
      // the coils here irrespective if there's one or 4 unqiue coils.
      for (G4int i = 0; i < 4; i++)
	{
	  G4LogicalVolume* volToPlace = coilLV;
	  G4ThreeVector displacement = G4ThreeVector();
	  if (individualCoilsSolids)
	    {volToPlace = coilLVs[i];}
	  else
	    {displacement = coilDisps[i];} // with no intersection we have to displace it
	  G4String theName = name + "_coil_" + std::to_string(i) + "_pv";
	  G4RotationMatrix* rot = new G4RotationMatrix();
	  allRotationMatrices.push_back(rot);
	  if (buildVertically)
	    {
	      rot->rotateZ(CLHEP::halfpi);
	      displacement.transform(*rot);
	    }
	  aCoilPV = new G4PVPlacement(rot,             // no rotation
				      displacement,    // position
				      volToPlace,      // lv to be placed
				      theName,         // name
				      containerLV,     // mother lv to be placed in
				      false,           // no boolean operation
				      0,               // copy number
				      checkOverlaps);
	  allPhysicalVolumes.push_back(aCoilPV);
	}
    }

  // magnet outer instance
  CreateMagnetContainerComponent();
  magnetContainer->RegisterSolid(magnetContainerExtraSolids);

  // build the BDSMagnetOuter instance and return it
  BDSMagnetOuter* outer = new BDSMagnetOuter(containerSolid,
					     containerLV, ext,
					     magnetContainer);
  
  outer->RegisterSolid(allSolids);
  outer->RegisterLogicalVolume(allLogicalVolumes);
  outer->RegisterRotationMatrix(allRotationMatrices);
  outer->RegisterPhysicalVolume(allPhysicalVolumes);
  outer->RegisterVisAttributes(allVisAttributes);
  
  outer->RegisterSolid(yokeSolid);
  outer->RegisterLogicalVolume(yokeLV);

  outer->RegisterSensitiveVolume(yokeLV);

  // no need to proceed with end pieces if we didn't build poles - just return
  if (!buildPole)
    {return outer;}

  // continue with registration of objects and end piece construction  
  if (individualCoilsSolids)
    {outer->RegisterSensitiveVolume(coilLVs);}
  else
    {outer->RegisterSensitiveVolume(coilLV);}

  // skip rest of this construction if no end pieces required
  if (!buildEndPiece)
    {
      SetFaceNormals(outer); // would also update end pieces if they existed
      return outer;
    }
  
  // end pieces - note with bends both are likely to be different so build independently here

  // end piece coil solids - build in x-y plane and project in z as that's the
  // only way with an extruded solid then rotate so it's along z
  // curved sections are elliptical with major axis being the projected with of the coil
  // at an angle and the minor being the normal coil width.
  // x = x0 + acos(t)cos(phi) - bsin(t)sin(phi)
  // y = y0 + acos(t)sin(phi) + bsin(t)cos(phi)
  std::vector<G4TwoVector> inEPPoints;  // input face end piece points
  std::vector<G4TwoVector> outEPPoints; // output face end piece points

  G4double inXO = poleHalfWidth - lsl;

  // create an ellipse with no angle, then shear it to match the angle
  G4int nSegments = ceil((G4double)nSegmentsPerCircle / 4.0);
  G4double increment = CLHEP::halfpi/nSegments;
  for (G4double t = -CLHEP::pi; t <= -CLHEP::halfpi + 1e-9; t += increment)
    { // left side
      G4double x = -inXO + coilWidth*cos(t);
      G4double y = coilWidth*sin(t);
      inEPPoints.emplace_back(x,y);
    }
  
  for (G4double t = -CLHEP::halfpi; t <= 0 + 1e-9; t += increment)
    { // right side
      G4double x = inXO + coilWidth*cos(t);
      G4double y = coilWidth*sin(t);
      inEPPoints.emplace_back(x,y);
    }

  // shear it
  // (x')   (1        0) (x)
  // (y') = (tan(phi) 1) (y)
  // ->  x' = x; y' = tan(phi)*x + y
  // copy the points, flip in y and shear for output angle
  for (const auto point : inEPPoints) // copying this time
    {
      G4double outy = -1*(point.x()*tan(-angleOut) + point.y());
      outEPPoints.emplace_back(point.x(),outy);
    }
  for (auto& point : inEPPoints)  // modify in place for shearing original points
    {point.setY(point.x()*tan(-angleIn) + point.y());}

  G4TwoVector zOffsets(0,0); // the transverse offset of each plane from 0,0
  G4double zScale = 1;       // the scale at each end of the points = 1
  
  // these are projected by coilHeight in z as they'll be rotated later
  G4VSolid* endPieceSolidIn  = new G4ExtrudedSolid(name + "_end_coil_in_solid", // name
						   inEPPoints, // transverse 2d coordinates
						   coilHeight*0.5 - lsl, // z half length
						   zOffsets, zScale,
						   zOffsets, zScale);

  // these are projected by coilHeight in z as they'll be rotated later
  G4VSolid* endPieceSolidOut = new G4ExtrudedSolid(name + "_end_coil_out_solid", // name
						   outEPPoints, // transverse 2d coordinates
						   coilHeight*0.5 - lsl, // z half length
						   zOffsets, zScale,
						   zOffsets, zScale);
  
  // end piece container solids - simple extruded solid intersectd with cut tubs for
  // angled faces build about magnet zero so we get the coordinates right for general placement
  std::vector<G4TwoVector> contEPPoints;
  const G4double connector = 1*CLHEP::mm;
  G4double xmax = poleHalfWidth + coilWidth + connector;
  G4double ymax = poleHalfGap + coilHeight + cDY + connector;
  G4double yInn = poleHalfGap + cDY;
  contEPPoints.emplace_back(xmax + connector,  ymax);
  contEPPoints.emplace_back(-xmax,  ymax);
  contEPPoints.emplace_back(-xmax,  yInn);
  contEPPoints.emplace_back( xmax,  yInn);
  contEPPoints.emplace_back( xmax, -yInn);
  contEPPoints.emplace_back(-xmax, -yInn);
  contEPPoints.emplace_back(-xmax, -ymax);
  contEPPoints.emplace_back(xmax + connector, -ymax);

  G4double containerIntersectionRadius = std::hypot(ymax, xmax);
  
  if (buildVertically)
    {
      for (auto& point : contEPPoints)
	{point = BDS::Rotate(point, CLHEP::halfpi);}
    }
  
  // calculate length for each end piece container volume - could be different due to different angles
  G4double ePInLength  = coilWidth + lengthSafetyLarge; // adjustable for intersection
  G4double ePInLengthZ = ePInLength; // copy that will be final length of object
  if (BDS::IsFinite(angleIn))
    {
      G4double dzIn = tan(std::abs(angleIn)) * outerDiameter;
      ePInLength    = std::max(2*ePInLength, 2*dzIn); // 2x as long for intersection
    }

    BDSExtent ePExtOuter = BDSExtent(-xmax, xmax + connector,
                                     -ymax, ymax,
                                     -ePInLengthZ*0.5, ePInLengthZ*0.5);
    BDSExtent ePExtInner = BDSExtent(-xmax, xmax,
                                     -yInn, yInn,
                                     -ePInLengthZ*0.5, ePInLengthZ*0.5);

  G4VSolid* ePContSolidIn  = new G4ExtrudedSolid(name + "_end_coil_in_solid", // name
						 contEPPoints,   // transverse 2d coordinates
						 ePInLength*0.5, // z half length
						 zOffsets, zScale,
						 zOffsets, zScale);

  G4double ePOutLength  = coilWidth + lengthSafetyLarge; // adjustable for intersection
  G4double ePOutLengthZ = ePOutLength; // copy that will be final length of object
  if (BDS::IsFinite(angleOut))
    {
      G4double dzOut = tan(std::abs(angleOut)) * outerDiameter;
      ePOutLength    = std::max(2*ePOutLength, 2*dzOut); // 2x as long for intersection
    }

  G4VSolid* ePContSolidOut = new G4ExtrudedSolid(name + "_end_coil_out_solid", // name
						 contEPPoints,   // transverse 2d coordinates
						 ePOutLength*0.5, // z half length
						 zOffsets, zScale,
						 zOffsets, zScale);

  // If there is a finite input face angle, the extruded solid above will be longer and
  // we now make a cut tubs for the angled faces to intersect it with.
  if (BDS::IsFinite(angleIn))
    {
      
      G4ThreeVector inputfaceReversed = inputFaceNormal * -1;
      G4VSolid* ePContSolidInAng = new G4CutTubs(name + "_angled_face_solid", // name
						 0,                           // inner radius
						 containerIntersectionRadius, // outer radius
						 ePInLengthZ*0.5,             // z half length
						 0,                           // start angle
						 CLHEP::twopi,                // sweep angle
						 inputFaceNormal,             // input face normal
						 inputfaceReversed);          // output face normal

      // potential memory leak here with overwriting
      ePContSolidIn = new G4IntersectionSolid(name + "_end_coil_cont_solid", // name
					      ePContSolidIn,
					      ePContSolidInAng);
    }

  if (BDS::IsFinite(angleOut))
    {
      
      G4ThreeVector outputfaceReversed = outputFaceNormal * -1;
      G4VSolid* ePContSolidOutAng = new G4CutTubs(name + "_angled_face_solid", // name
						  0,                           // inner radius
						  containerIntersectionRadius, // outer radius
						  ePOutLengthZ*0.5,            // z half length
						  0,                           // start angle
						  CLHEP::twopi,                // sweep angle
						  outputfaceReversed,          // input face normal
						  outputFaceNormal);           // output face normal

      // potential memory leak here with overwriting
      ePContSolidOut = new G4IntersectionSolid(name + "_end_coil_cont_solid", // name
					       ePContSolidOut,
					       ePContSolidOutAng);
    }

  // end piece logical volumes
  // coil end piece lv
  G4Material* copper = BDSMaterials::Instance()->GetMaterial("copper");
  G4Material* emptyMaterial = BDSMaterials::Instance()->GetMaterial(BDSGlobalConstants::Instance()->EmptyMaterial());
  
  G4LogicalVolume* ePInLV = new G4LogicalVolume(endPieceSolidIn,           // solid
						copper,                    // material
						name + "_end_coil_in_lv"); // name
  
  G4LogicalVolume* ePContInLV = new G4LogicalVolume(ePContSolidIn,             // solid
						    emptyMaterial,             // material
						    name + "_end_cont_in_lv"); // name

  G4LogicalVolume* ePOutLV = new G4LogicalVolume(endPieceSolidOut,           // solid
						 copper,                     // material
						 name + "_end_coil_out_lv"); // name
  
  G4LogicalVolume* ePContOutLV = new G4LogicalVolume(ePContSolidOut,             // solid
						     emptyMaterial,              // material
						     name + "_end_cont_out_lv"); // name

  G4Colour* coilColour = BDSColours::Instance()->GetColour("coil");
  G4VisAttributes* coilVisIn  = new G4VisAttributes(*coilColour);
  coilVisIn->SetVisibility(true);
  G4VisAttributes* coilVisOut = new G4VisAttributes(*coilVisIn);

  ePInLV->SetVisAttributes(coilVisIn);
  ePOutLV->SetVisAttributes(coilVisOut);
  ePContInLV->SetVisAttributes(BDSGlobalConstants::Instance()->ContainerVisAttr());
  ePContOutLV->SetVisAttributes(BDSGlobalConstants::Instance()->ContainerVisAttr());

  auto defaultUserLimits = BDSGlobalConstants::Instance()->DefaultUserLimits();
  ePInLV->SetUserLimits(defaultUserLimits);
  ePOutLV->SetUserLimits(defaultUserLimits);
  ePContInLV->SetUserLimits(defaultUserLimits);
  ePContOutLV->SetUserLimits(defaultUserLimits);
  
  // placements
  G4RotationMatrix* endCoilInRM = new G4RotationMatrix();
  endCoilInRM->rotateX(-CLHEP::halfpi);
  if (buildVertically)
    {endCoilInRM->rotateY(CLHEP::halfpi);}
  G4ThreeVector endCoilTranslationInTop(0, coilDY-lsl, 0.5*coilWidth);
  G4ThreeVector endCoilTranslationInLow(0,-coilDY+lsl, 0.5*coilWidth);
  G4RotationMatrix* vertRot = new G4RotationMatrix();
  vertRot->rotateZ(CLHEP::halfpi);
  if (buildVertically)
    {
      endCoilTranslationInTop.transform(*vertRot);
      endCoilTranslationInLow.transform(*vertRot);
    }
  G4PVPlacement* ePInTopPv = new G4PVPlacement(endCoilInRM,        // rotation
					       endCoilTranslationInTop, // position
					       ePInLV,             // logical volume
					       name + "_end_piece_in_top_pv", // name
					       ePContInLV,         // mother lv to be placed in
					       false,              // no boolean operation
					       0,                  // copy number
					       checkOverlaps);     // check overlaps
  G4PVPlacement* ePInLowPv = new G4PVPlacement(endCoilInRM,        // rotation
					       endCoilTranslationInLow, // position
					       ePInLV,             // logical volume
					       name + "_end_piece_in_low_pv", // name
					       ePContInLV,         // mother lv to be placed in
					       false,              // no boolean operation
					       0,                  // copy number
					       checkOverlaps);     // check overlaps

  G4RotationMatrix* endCoilOutRM = new G4RotationMatrix(*endCoilInRM); // copy as independent objects
  
  G4ThreeVector endCoilTranslationOutTop(0, coilDY, -0.5*coilWidth);
  G4ThreeVector endCoilTranslationOutLow(0,-coilDY, -0.5*coilWidth);
  if (buildVertically)
    {
      endCoilTranslationOutTop.transform(*vertRot);
      endCoilTranslationOutLow.transform(*vertRot);
    }
  delete vertRot;
  G4PVPlacement* ePOutTopPv = new G4PVPlacement(endCoilOutRM,       // rotation
						endCoilTranslationOutTop, // position
						ePOutLV,            // logical volume
						name + "_end_piece_out_top_pv", // name
						ePContOutLV,        // mother lv to be placed in
						false,              // no boolean operation
						0,                  // copy number
						checkOverlaps);     // check overlaps
  G4PVPlacement* ePOutLowPv = new G4PVPlacement(endCoilOutRM,        // rotation
						endCoilTranslationOutLow, // position
						ePOutLV,             // logical volume
						name + "_end_piece_out_low_pv", // name
						ePContOutLV,        // mother lv to be placed in
						false,              // no boolean operation
						0,                  // copy number
						checkOverlaps);     // check overlaps
  
  // packaging - geometry component
  G4ThreeVector inputFaceNormalR = inputFaceNormal * -1; // just -1 as parallel but opposite
  BDSSimpleComponent* endPieceInSC = new BDSSimpleComponent(name + "_end_piece_in",
							    ePInLengthZ,
							    0/*angle*/,
							    ePContSolidIn,
							    ePContInLV,
							    inputFaceNormal,
							    inputFaceNormalR);
  
  endPieceInSC->RegisterPhysicalVolume(ePInTopPv);
  endPieceInSC->RegisterPhysicalVolume(ePInLowPv);
  endPieceInSC->RegisterRotationMatrix(endCoilInRM);
  endPieceInSC->RegisterVisAttributes(coilVisIn);
  endPieceInSC->RegisterLogicalVolume(ePInLV);
  endPieceInSC->RegisterSolid(endPieceSolidIn);
  endPieceInSC->RegisterSensitiveVolume(ePInLV);
  endPieceInSC->SetExtent(ePExtOuter);
  endPieceInSC->SetInnerExtent(ePExtInner);
  
  G4ThreeVector outputFaceNormalR = outputFaceNormal * -1; // just -1 as parallel but opposite
  BDSSimpleComponent* endPieceOutSC = new BDSSimpleComponent(name + "_end_piece_out",
							     ePOutLengthZ,
							     0/*angle*/,
							     ePContSolidOut,
							     ePContOutLV,
							     outputFaceNormalR,
							     outputFaceNormal);
  
  endPieceOutSC->RegisterPhysicalVolume(ePOutTopPv);
  endPieceOutSC->RegisterPhysicalVolume(ePOutLowPv);
  endPieceOutSC->RegisterRotationMatrix(endCoilOutRM);
  endPieceOutSC->RegisterVisAttributes(coilVisOut);
  endPieceOutSC->RegisterLogicalVolume(ePOutLV);
  endPieceOutSC->RegisterSolid(endPieceSolidOut);
  endPieceOutSC->RegisterSensitiveVolume(ePOutLV);
  endPieceOutSC->SetExtent(ePExtOuter);
  endPieceOutSC->SetInnerExtent(ePExtInner);

  // attach to the magnet outer
  outer->SetEndPieceBefore(endPieceInSC);
  outer->SetEndPieceAfter(endPieceOutSC);
  
  // update normals of outer and end pieces
  SetFaceNormals(outer);
  
  return outer;
}
