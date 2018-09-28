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
#ifndef BDSCOLLIMATORJAW_H
#define BDSCOLLIMATORJAW_H

#include "BDSAcceleratorComponent.hh"

#include "globals.hh" // geant4 types / globals
#include "G4Material.hh"

class G4Colour;
class G4VSolid;

/**
 * @brief Collimator with only two jaw and no top bit.
 *
 * @author Laurie Nevay
 */

class BDSCollimatorJaw: public BDSAcceleratorComponent
{
public:
  BDSCollimatorJaw(G4String  nameIn,
		   G4double  lengthIn,
		   G4double  horizontalWidthIn,
		   G4double  jaw1OffsetIn,
		   G4double  jaw2OffsetIn,
		   G4double  xHalfGapIn,
		   G4double  yHalfHeightIn,
		   G4bool    buildLeftJawIn,
		   G4bool    buildRightJawIn,
		   G4Material* collimatorMaterialIn,
		   G4Material* vacuumMaterialIn,
		   G4Colour*   colourIn = nullptr);
  virtual ~BDSCollimatorJaw();

  /// Accessor for collimator material.
  virtual G4String Material() const {return collimatorMaterial->GetName();}

protected:
  virtual void Build();
  
  virtual void BuildContainerLogicalVolume();

  virtual void BuildJawCollimator();

  ///@{ Geometrical objects:
  G4VSolid* jawSolid;
  G4VSolid* vacuumSolid;
  ///@}

  G4double    horizontalWidth;    ///< Horizontal width.
  G4double    xsizeLeft;          ///< Offset of jaw 1
  G4double    xsizeRight;         ///< Offset of jaw 2
  G4double    xHalfGap;           ///< Half gap separation between jaws.
  G4double    jawHalfWidth;       ///< Half width of each jaw.
  G4double    yHalfHeight;        ///< Half height of each jaw.
  G4bool      buildLeftJaw;       ///< Build left jaw or not.
  G4bool      buildRightJaw;      ///< Build right jaw or not.
  G4bool      buildAperture;	  ///< Build aperture or not.
  G4Material* collimatorMaterial; ///< Material
  G4Material* vacuumMaterial;     ///< Vacuum material
  G4Colour*   colour;             ///< Colour of collimator

private:
  /// Private default constructor to force the use of the supplied one.
  BDSCollimatorJaw() = delete;

  /// @{ Assignment and copy constructor not implemented nor used
  BDSCollimatorJaw& operator=(const BDSCollimatorJaw&) = delete;
  BDSCollimatorJaw(BDSCollimatorJaw&) = delete;
  /// @}
};

#endif
