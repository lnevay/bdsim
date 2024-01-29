/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2024.

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
#ifndef BDSFIELDMAGMULTIPOLEOUTEROLD_H
#define BDSFIELDMAGMULTIPOLEOUTEROLD_H

#include "BDSFieldMag.hh"

#include "globals.hh" // geant4 types / globals
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

#include <vector>

class BDSMagnetStrength;

/**
 * @brief A simple paramaterisation of N-Pole outer yoke magnetic field.
 *
 * Field that roughly approximates the field in a yoke of an n-pole magnet.
 * The field is generated by the linear sum of infinitely long wire current
 * sources along +-z (local) placed at a radius of poleTipRadius according
 * to the order of the magnet. Note, these current sources are placed exaclty
 * in between the pole will be for each type of magnet. The field is normalised
 * to the field at the pole tip radius given innerFieldIn.
 * 
 * If a point close (< 6mm) to a pole is queried, the field is capped at the 
 * normalised value at the pole tip - assumed to be the highest value.
 *
 * This field of course does not accurately represent the field in the yoke
 * of a magnet as it neglects the permeability and the geometry, however, it
 * is a 1st order approximation.
 *
 * Does not own or delete the BDSFieldMag* innerFieldIn argument.
 *
 * @author Laurie Nevay
 */

class BDSFieldMagMultipoleOuterOld: public BDSFieldMag
{
public:
  BDSFieldMagMultipoleOuterOld(G4int              orderIn,
			    G4double           poleTipRadius,
			    const BDSFieldMag* innerFieldIn,
			    G4bool             kPositive,
                               G4double           arbitraryScaling = 1.0);

  virtual ~BDSFieldMagMultipoleOuterOld(){;}

  /// Access the field value.
  virtual G4ThreeVector GetField(const G4ThreeVector& position,
				 const double         t = 0) const;

private:
  const G4int       order;           ///< N-poles / 2.
  G4double          normalisation;   ///< Storage of the overal normalisation factor.
  G4bool            positiveField;   ///< Sign of magnetic field.
  G4double          poleTipRadius;   ///< Radius of transition between inner and outer fields.
  std::vector<G4TwoVector> currents; ///< Locations of inifite wire current sources.
};

#endif
