/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2020.

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
#include "BDSAperture.hh"
#include "BDSApertureFactory.hh"
#include "BDSBeamPipeFactoryGeneral.hh"
#include "BDSBeamPipe.hh"
#include "BDSBeamPipeInfo2.hh"
#include "BDSExtent.hh"

#include "globals.hh"                      // geant4 globals / types
#include "G4CutTubs.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"

#include <cmath>
#include <set>
#include <utility>               // for std::pair

BDSBeamPipeFactoryGeneral::BDSBeamPipeFactoryGeneral()
{;}

BDSBeamPipe* BDSBeamPipeFactoryGeneral::CreateBeamPipe(const G4String& name,
                                                G4double length,
                                                BDSBeamPipeInfo2* bpi)
{
  CleanUp();
  
  BDSAperture* apVacIn  = bpi->aperture;
  BDSAperture* apVacOut = bpi->apertureOut ? bpi->apertureOut : bpi->aperture;
  G4double containerThickness = bpi->vacuumOnly ? lengthSafetyLarge : bpi->beamPipeThickness + 2*lengthSafetyLarge;
  BDSAperture* apBpInnerIn  = apVacIn->Plus(lengthSafetyLarge);
  BDSAperture* apBpInnerOut = apVacOut->Plus(lengthSafetyLarge);
  BDSAperture* apBpOuterIn  = apBpInnerIn->Plus(bpi->beamPipeThickness);
  BDSAperture* apBpOuterOut = apBpInnerOut->Plus(bpi->beamPipeThickness);
  BDSAperture* apContIn     = apVacIn->Plus(containerThickness);
  BDSAperture* apContOut    = apVacOut->Plus(containerThickness);
  BDSAperture* apContSubIn  = apVacIn->Plus(containerThickness + lengthSafety);
  BDSAperture* apContSubOut = apVacOut->Plus(containerThickness + lengthSafety);
  
  BDSApertureFactory fac;
  vacuumSolid = fac.CreateSolid(name+"_vacuum", length - lengthSafety, apVacIn, apVacOut,
                                bpi->inputFaceNormal, bpi->outputFaceNormal);
  
  if (!bpi->vacuumOnly)
  {
    beamPipeSolid = fac.CreateSolidWithInner(name, length,
                                             apBpInnerIn, apBpOuterIn,
                                             apBpInnerOut, apBpOuterOut,
                                             bpi->inputFaceNormal, bpi->outputFaceNormal);
  }
  
  containerSolid = fac.CreateSolid(name + "_container_solid", length, apContIn, apContOut,
  bpi->inputFaceNormal, bpi->outputFaceNormal);
  
  containerSubtractionSolid = fac.CreateSolid(name+"_container_sub",
                                              length, apContSubIn, apContSubOut,
                                              nullptr, nullptr,
                                              0.2*length);
  
  BDSBeamPipeFactoryBase::CommonConstruction(name, bpi->vacuumMaterial, bpi->beamPipeMaterial, length);

  // record extents
  BDSExtent ext = std::max(apContIn->Extent(), apContOut->Extent());

  // true for containerIsGeneral - true for this factory
  return BDSBeamPipeFactoryBase::BuildBeamPipeAndRegisterVolumes(ext, /*TBC*/containerThickness, true);
}