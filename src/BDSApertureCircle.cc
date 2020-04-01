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
#include "BDSApertureCircle.hh"
#include "BDSApertureType.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSExtent.hh"
#include "BDSGlobalConstants.hh"
#include "BDSPolygon.hh"
#include "BDSTiltOffset.hh"
#include "BDSUtilities.hh"

#include "G4TwoVector.hh"
#include "G4Types.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <vector>

BDSApertureCircle::BDSApertureCircle(G4double radiusIn,
                                     G4int    nPointsIn):
  BDSAperture(BDSApertureType::circle),
  radius(radiusIn),
  nPoints(nPointsIn)
{
  if (nPoints == 0)
    {nPoints = BDSGlobalConstants::Instance()->NSegmentsPerCircle();}
}

G4bool BDSApertureCircle::Equals(const BDSAperture* other) const
{
  if (!other)
    {return false;}
  else if (other->apertureType != apertureType)
    {return false;}
  else if (other->tiltOffset != tiltOffset)
    {return false;}
  else
    {
      const BDSApertureCircle* oc = dynamic_cast<const BDSApertureCircle*>(other);
      return BDS::DoublesAreEqual(oc->radius, radius);
    }
}

void BDSApertureCircle::CheckInfoOK() const
{
  CheckRequiredParametersSet(radius, true);
  CheckParameterIsPositive(radius, "radius");
}

BDSExtent BDSApertureCircle::Extent() const
{
  BDSExtent simpleExtent(radius, radius, 0);
  return ExtentOffsetTilt(simpleExtent);
}

BDSApertureCircle BDSApertureCircle::operator+ (G4double number) const
{
  BDSApertureCircle result = BDSApertureCircle(*this);
  result += number;
  return result;
}

const BDSApertureCircle& BDSApertureCircle::operator+=(G4double number)
{
  radius += number;
  return *this;
}

BDSApertureCircle BDSApertureCircle::operator* (G4double number) const
{
  BDSApertureCircle result = BDSApertureCircle(*this);
  result *= number;
  return result;
}

const BDSApertureCircle& BDSApertureCircle::operator*=(G4double number)
{
  radius *= number;
  return *this;
}

BDSPolygon BDSApertureCircle::Polygon(G4int nPointsIn) const
{
  G4int np = nPointsIn == 0 ? nPoints : nPointsIn;
  if (np < 3 || np < MinimumNumberOfPoints())
    {throw BDSException(__METHOD_NAME__, "number of points for aperture specified < 3.");}

  std::vector<G4TwoVector> r;
  G4double dTheta = CLHEP::twopi / (G4double)np;
  for (G4int i = 0; i < np; i++)
    {r.push_back(G4TwoVector(radius * std::cos(i*dTheta), radius * std::sin(i*dTheta)));}

  return BDSPolygon(r).ApplyTiltOffset(tiltOffset);
}

