/* 
Beam Delivery Simulation (BDSIM) Copyright (C) BDSIM Collaboration, 2001 - 2026.

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
#include "BDSAcceleratorModel.hh"
#include "BDSCollimatorTipJaw.hh"
#include "BDSBeamPipeInfo.hh"
#include "BDSCollimatorJaw.hh"
#include "BDSColours.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSMaterials.hh"
#include "BDSSDType.hh"

#include "G4Box.hh"
#include "G4Para.hh"
#include "G4GenericTrap.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

#include <cmath>
#include <set>

BDSCollimatorTipJaw::BDSCollimatorTipJaw(const G4String&    nameIn,
                                         G4double    lengthIn,
                                         G4double    horizontalWidthIn,
                                         G4double    xHalfGapIn,
                                         G4double    yHalfHeightIn,
                                         G4double    xSizeLeftIn,
                                         G4double    xSizeRightIn,
                                         G4double    leftJawTiltIn,
                                         G4double    rightJawTiltIn,
                                         G4double    tipThicknessIn,
                                         G4bool      buildLeftJawIn,
                                         G4bool      buildRightJawIn,
                                         G4Material* collimatorMaterialIn,
                                         G4Material* collimatorTipMaterialIn,
                                         G4Material* vacuumMaterialIn,
                                         G4Colour*   colourIn,
                                         G4Colour*   tipColourIn):
BDSCollimatorJaw(nameIn, lengthIn, horizontalWidthIn, xHalfGapIn, yHalfHeightIn, xSizeLeftIn, xSizeRightIn,
                 leftJawTiltIn, rightJawTiltIn, buildLeftJawIn, buildRightJawIn, collimatorMaterialIn,
                 vacuumMaterialIn, colourIn, "jcoltip"),
  tipColour(tipColourIn),
  tipThickness(tipThicknessIn),
  collimatorTipMaterial(collimatorTipMaterialIn)
{;}

BDSCollimatorTipJaw::~BDSCollimatorTipJaw()
{;}

void BDSCollimatorTipJaw::Calculations()
{
  leftJawWidth -= (tipThickness + 2*lengthSafety);
  leftJawPos -= G4ThreeVector(leftJawWidth*0.5, 0, 0);
  rightJawWidth -= (tipThickness + 2*lengthSafety);
  rightJawPos += G4ThreeVector(rightJawWidth*0.5, 0, 0);

  leftJawTipPos = G4ThreeVector(leftJawHalfGap+0.5*tipThickness, 0, 0);
  rightJawTipPos = G4ThreeVector(rightJawHalfGap-0.5*tipThickness, 0, 0);
}

void BDSCollimatorTipJaw::CheckParameters()
{
  BDSCollimatorJaw::CheckParameters();

  if (!tipColour)
    {tipColour = BDSColours::Instance()->GetColour("collimatorTip");}

  // the remaining checks only apply to the jaw and vacuum geometry
  if (!buildAperture)
    {return;}

  // tip solids have a half width of tipThickness/2 - lengthSafety
  if (tipThickness * 0.5 - lengthSafety < 1e-3) // 1um minimum, could also be negative
    {throw BDSException(__METHOD_NAME__, "insufficient tipThickness for jcoltip \"" + name + "\"");}

  // bulk jaw solids have a half width of (jawWidth - tipThickness)/2 - lengthSafety
  if (buildLeftJaw && ((0.5 * horizontalWidth - leftJawHalfGap - tipThickness) * 0.5 - lengthSafety < 1e-3))
    {throw BDSException(__METHOD_NAME__, "left jaw of jcoltip \"" + name + "\" too thin given horizontalWidth, aperture and tipThickness");}
  if (buildRightJaw && ((0.5 * horizontalWidth - rightJawHalfGap - tipThickness) * 0.5 - lengthSafety < 1e-3))
    {throw BDSException(__METHOD_NAME__, "right jaw of jcoltip \"" + name + "\" too thin given horizontalWidth, aperture and tipThickness");}

  // shift of each jaw face at the ends of the element due to its tilt - tilt is ignored for a jaw that isn't built
  G4double tiltShiftLeft  = buildLeftJaw  ? std::tan(jawTiltLeft)  * chordLength * 0.5 : 0;
  G4double tiltShiftRight = buildRightJaw ? std::tan(jawTiltRight) * chordLength * 0.5 : 0;

  if (std::abs(tiltShiftLeft) > leftJawHalfGap)
    {throw BDSException(__METHOD_NAME__, "tilted left jaw not allowed to cross the mid-plane: \"" + name + "\"");}
  if (std::abs(tiltShiftRight) > rightJawHalfGap)
    {throw BDSException(__METHOD_NAME__, "tilted right jaw not allowed to cross the mid-plane: \"" + name + "\"");}

  // vacuum full width at each end of the element - see vacuum construction in Build()
  G4double vacuumWidthUpstream   = (leftJawHalfGap - tiltShiftLeft) + (rightJawHalfGap + tiltShiftRight);
  G4double vacuumWidthDownstream = (leftJawHalfGap + tiltShiftLeft) + (rightJawHalfGap - tiltShiftRight);
  if (std::min(vacuumWidthUpstream, vacuumWidthDownstream) * 0.5 - lengthSafety < 1e-3) // 1um minimum
    {throw BDSException(__METHOD_NAME__, "insufficient aperture between jaws for jcoltip \"" + name + "\"");}
}

void BDSCollimatorTipJaw::Build()
{
  BDSCollimatorJaw::Build();
  BuildTips();
}

void BDSCollimatorTipJaw::BuildTips()
{
  G4VisAttributes* tipVisAttr = new G4VisAttributes(*tipColour);
  RegisterVisAttributes(tipVisAttr);

  G4UserLimits* tipCollUserLimits = CollimatorUserLimits();
  RegisterUserLimits(tipCollUserLimits);

  G4VSolid* leftJawTipSolid = nullptr;
  if (buildLeftJaw && buildAperture)
    {
      G4double leftHalfLength = chordLength * 0.5 * std::cos(jawTiltLeft);
      leftJawTipSolid = new G4Para(name + "_leftjawtip_solid",
                                   tipThickness * 0.5 - lengthSafety,
                                   yHalfHeight - lengthSafety,
                                   leftHalfLength - lengthSafety,
                                   0, jawTiltLeft, 0);
    }
  else
    {
      leftJawTipSolid = new G4Box(name + "_leftjawtip_solid",
                                  tipThickness * 0.5 - lengthSafety,
                                  yHalfHeight - lengthSafety,
                                  chordLength * 0.5 - lengthSafety);
    }
  RegisterSolid(leftJawTipSolid);

  G4LogicalVolume* leftJawTipLV = new G4LogicalVolume(leftJawTipSolid, collimatorTipMaterial, name + "_leftjawtip_lv");
  leftJawTipLV->SetVisAttributes(tipVisAttr);
  leftJawTipLV->SetUserLimits(tipCollUserLimits);
  RegisterLogicalVolume(leftJawTipLV);
  BDSAcceleratorModel::Instance()->VolumeSet("collimators")->insert(leftJawTipLV);
  if (sensitiveOuter)
    {RegisterSensitiveVolume(leftJawTipLV, BDSSDType::collimatorcomplete);}

  // place the tip
  G4PVPlacement* leftJawTipPV = new G4PVPlacement(nullptr,
                                                  leftJawTipPos,
                                                  leftJawTipLV,
                                                  name + "_leftjawtip_pv",
                                                  containerLogicalVolume,
                                                  false,               // no boolean operation
                                                  0,                   // copy number
                                                  checkOverlaps);
  RegisterPhysicalVolume(leftJawTipPV);

  G4VSolid* rightJawTipSolid = nullptr;
  if (buildRightJaw && buildAperture)
    {
      G4double rightHalfLength = chordLength * 0.5 * std::cos(jawTiltRight);
      rightJawTipSolid = new G4Para(name + "_rightjawtip_solid",
                                    tipThickness * 0.5 - lengthSafety,
                                    yHalfHeight - lengthSafety,
                                    rightHalfLength - lengthSafety,
                                    0, jawTiltRight, 0);
    }
  else
    {
      rightJawTipSolid = new G4Box(name + "_rightjawtip_solid",
                                   tipThickness * 0.5 - lengthSafety,
                                   yHalfHeight - lengthSafety,
                                   chordLength * 0.5 - lengthSafety);
    }
  RegisterSolid(rightJawTipSolid);

  G4LogicalVolume* rightJawTipLV = new G4LogicalVolume(rightJawTipSolid, collimatorTipMaterial, name + "_rightjawtip_lv");
  rightJawTipLV->SetVisAttributes(tipVisAttr);
  rightJawTipLV->SetUserLimits(tipCollUserLimits);
  RegisterLogicalVolume(rightJawTipLV);
  BDSAcceleratorModel::Instance()->VolumeSet("collimators")->insert(rightJawTipLV);
  if (sensitiveOuter)
    {RegisterSensitiveVolume(rightJawTipLV, BDSSDType::collimatorcomplete);}

  // place the tip
  G4PVPlacement* rightJawTipPV = new G4PVPlacement(nullptr,
                                                   rightJawTipPos,
                                                   rightJawTipLV,
                                                   name + "_rightjawtip_pv",
                                                   containerLogicalVolume,
                                                   false,           // no boolean operation
                                                   0,               // copy number
                                                   checkOverlaps);
  RegisterPhysicalVolume(rightJawTipPV);
}