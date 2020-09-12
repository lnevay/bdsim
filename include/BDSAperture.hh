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
#ifndef BDSAPERTURE_H
#define BDSAPERTURE_H

#include "BDSApertureType.hh"
#include "BDSExtent.hh"
#include "BDSTiltOffset.hh"

#include "G4String.hh"
#include "G4TwoVector.hh"
#include "G4Types.hh"

#include <vector>

class BDSExtent;
class BDSPolygon;

/**
 * @brief A definition of an aperture shape.
 *
 * @author Laurie Nevay
 */

class BDSAperture
{
public:
  explicit BDSAperture(BDSApertureType apertureTypeIn);
  BDSAperture(BDSApertureType      apertureTypeIn,
	      const BDSTiltOffset& tiltOffsetIn);
  virtual ~BDSAperture(){;}

  inline void SetTiltOffset(const BDSTiltOffset& tiltOffsetIn) {tiltOffset = tiltOffsetIn;}
  G4bool FiniteTilt() const {return tiltOffset.HasFiniteTilt();}
  
  BDSApertureType apertureType;
  BDSTiltOffset   tiltOffset;

  friend bool operator == (const BDSAperture& lhs, const BDSAperture& rhs) {return (lhs.apertureType == rhs.apertureType) && lhs.tiltOffset == rhs.tiltOffset && lhs.Equals(&rhs);}
  friend bool operator != (const BDSAperture& lhs, const BDSAperture& rhs) {return !(lhs == rhs);}
  friend bool operator <  (const BDSAperture& lhs, const BDSAperture& rhs) {return lhs.LessThan(&rhs);}
  friend bool operator >  (const BDSAperture& lhs, const BDSAperture& rhs) {return rhs < lhs;}
  friend bool operator <= (const BDSAperture& lhs, const BDSAperture& rhs) {return !(lhs > rhs);}
  friend bool operator >= (const BDSAperture& lhs, const BDSAperture& rhs) {return !(lhs < rhs);}
  
  virtual BDSAperture* Plus(G4double number)  const = 0;
  virtual BDSAperture* Times(G4double number) const = 0;
  virtual BDSAperture* Clone() const = 0;
  
  virtual G4bool    Equals(const BDSAperture* other)   const = 0;
  virtual void      CheckInfoOK()                      const = 0;
  virtual G4double  RadiusToEncompass()                const = 0; // was IndicativeRadius RCWD
  virtual BDSExtent Extent()                           const = 0;
  
  /// Minimal number of 2D points to represent the shape, e.g. 3 for a triangle; 4 for a rectangle etc.
  virtual unsigned int MinimumNumberOfPoints() const = 0;

  /// Check flagged parameters (variable number) are non-zero and throw exception if not.
  void CheckRequiredParametersSet(G4double aper1,     G4bool aper1Set,
				  G4double aper2 = 0, G4bool aper2Set = false,
				  G4double aper3 = 0, G4bool aper3Set = false,
				  G4double aper4 = 0, G4bool aper4Set = false) const;

  /// Check value is >= 0 and throw exception with nameToWarn if not.
  void CheckParameterIsPositive(G4double        parameter,
				const G4String& nameToWarn) const;
  
  virtual G4bool LessThan(const BDSAperture* other) const;
  
  /// Generate a set of x,y points to represent the aperture. Number == 0 means minimal
  /// number of points which will be at least 3.
  BDSPolygon Polygon(unsigned int nPointsIn = 0) const;

protected:
  /// Protected interface - no protection over the value of nPoints used. This should therefore
  /// be called from the public Polygon() function that will check the value.
  virtual BDSPolygon PolygonNPoints(unsigned int nPointsIn) const = 0;
  
  /// Apply the members offsetX, offsetY and tilt to an extent. Utility
  /// function for Extent() implemented in derived classes.
  inline BDSExtent ExtentOffsetTilt(const BDSExtent& simpleExtent) const {return simpleExtent.TiltOffset(&tiltOffset);}

  /// Utility function to append points on an ellipse to a vector.
  static void AppendAngleEllipse(std::vector<G4TwoVector>& vec,
                                 G4double startAngle,
				 G4double finishAngle,
				 G4double radiusA,
				 G4double radiusB,
                                 unsigned int nPoints,
				 G4double xOffset,
				 G4double yOffset);
};

#endif
