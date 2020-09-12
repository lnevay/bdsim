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
#ifndef BDSAPERTUREELLIPSE_H
#define BDSAPERTUREELLIPSE_H

#include "BDSAperture.hh"
#include "BDSTiltOffset.hh"

#include "G4Types.hh"

#include <algorithm>

class BDSExtent;
class BDSPolygon;

/**
 * @brief Elliptical aperture.
 *
 * @author Laurie Nevay
 */

class BDSApertureEllipse: public BDSAperture
{
public:
  BDSApertureEllipse() = delete;
  BDSApertureEllipse(G4double aIn,
                     G4double bIn,
                     G4int    nPointsIn = 0);
  virtual ~BDSApertureEllipse(){;}

  G4double a;
  G4double b;
  G4int    nPoints;

  G4bool    Equals(const BDSAperture* other) const override;
  void      CheckInfoOK()                    const override;

  /// Note this is true even if the ellipse is rotated.
  G4double  RadiusToEncompass()     const override {return std::max(a,b) + tiltOffset.Radius();}
  BDSExtent Extent()                const override;
  unsigned int MinimumNumberOfPoints() const override {return nPoints;}

  BDSApertureEllipse        operator+ (G4double number) const;
  const BDSApertureEllipse& operator+=(G4double number);
  BDSApertureEllipse        operator* (G4double number) const;
  const BDSApertureEllipse& operator*=(G4double number);
  
  BDSAperture* Plus(G4double number) const override;
  BDSAperture* Times(G4double number) const override;
  BDSAperture* Clone() const override;

  BDSPolygon Polygon(G4int nPointsIn = 0) const override;
};

#endif
