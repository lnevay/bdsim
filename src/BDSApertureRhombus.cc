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
#include "BDSAperture.hh"
#include "BDSApertureRhombus.hh"
#include "BDSApertureType.hh"
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSExtent.hh"
#include "BDSPolygon.hh"
#include "BDSTiltOffset.hh"
#include "BDSUtilities.hh"

#include "G4TwoVector.hh"
#include "G4Types.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

BDSApertureRhombus::BDSApertureRhombus(G4double xIn,
                                       G4double yIn,
                                       G4double cornerRadiusIn,
                                       unsigned int nPointsPerTwoPi):
  BDSApertureCurved(BDSApertureType::rhombus, nPointsPerTwoPi),
  x(xIn),
  y(yIn),
  cornerRadius(cornerRadiusIn),
  nPointsMinimum(0)
{
  nPointsMinimum = cornerRadius > 0 ? 24 : 4;
  if (!(cornerRadius > 0))
    {nPoints = 4;}
  else
    {
      if (nPoints < nPointsMinimum)
        {nPoints = nPointsMinimum;}
    }
}

G4bool BDSApertureRhombus::Equals(const BDSAperture* other) const
{
  if (!other)
    {return false;}
  else if (other->apertureType != apertureType)
    {return false;}
  else if (other->tiltOffset != tiltOffset)
    {return false;}
  else
    {
      const BDSApertureRhombus* oc = dynamic_cast<const BDSApertureRhombus*>(other);
      return BDS::DoublesAreEqual(oc->x, x) &&
        BDS::DoublesAreEqual(oc->y, y) &&
        BDS::DoublesAreEqual(oc->cornerRadius, cornerRadius);
    }
}

void BDSApertureRhombus::CheckInfoOK() const
{
  CheckRequiredParametersSet(x, true, y, true); // cornerRadius is not required
  CheckParameterIsPositive(x, "x");
  CheckParameterIsPositive(y, "y");
  if (BDS::IsFinite(cornerRadius))
    {
      if (cornerRadius < 0)
        {throw BDSException(__METHOD_NAME__, "parameter \"aper3\" is not >= 0");}
    }
}

G4double BDSApertureRhombus::RadiusToEncompass() const
{
  return std::max(x, y) + tiltOffset.Radius();
}

BDSExtent BDSApertureRhombus::Extent() const
{
  BDSExtent simpleExtent(x, y, 0);
  return ExtentOffsetTilt(simpleExtent);
}

unsigned int BDSApertureRhombus::MinimumNumberOfPoints() const
{
  return nPointsMinimum;
}

BDSApertureRhombus BDSApertureRhombus::operator+ (G4double number) const
{
  BDSApertureRhombus result = BDSApertureRhombus(*this);
  result += number;
  return result;
}

const BDSApertureRhombus& BDSApertureRhombus::operator+=(G4double number)
{
  x += number;
  y += number;
  if (cornerRadius > 0)
    {cornerRadius += number;}
  return *this;
}

BDSApertureRhombus BDSApertureRhombus::operator* (G4double number) const
{
  BDSApertureRhombus result = BDSApertureRhombus(*this);
  result *= number;
  return result;
}

const BDSApertureRhombus& BDSApertureRhombus::operator*=(G4double number)
{
  x *= number;
  y *= number;
  cornerRadius *= number;
  return *this;
}

BDSAperture* BDSApertureRhombus::Plus(G4double number) const
{
  BDSApertureRhombus result = (*this) + number;
  return new BDSApertureRhombus(result);
}

BDSAperture* BDSApertureRhombus::Times(G4double number) const
{
  BDSApertureRhombus result = (*this) * number;
  return new BDSApertureRhombus(result);
}

BDSAperture* BDSApertureRhombus::Clone() const
{
  return new BDSApertureRhombus(*this);
}

std::array<G4double,7> BDSApertureRhombus::ApertureNumbers() const
{
  return {x,y,cornerRadius,0,tiltOffset.OffsetX(),tiltOffset.OffsetY(),tiltOffset.Tilt()};
}

BDSPolygon BDSApertureRhombus::PolygonNPoints(unsigned int nPointsIn) const
{
  std::vector<G4TwoVector> vec;
  if (! (cornerRadius > 0))
    {
      vec.emplace_back( x,  0);
      vec.emplace_back( 0,  y);
      vec.emplace_back(-x,  0);
      vec.emplace_back( 0, -y);
      return BDSPolygon(vec);
    }
  else
    {
      nPointsIn = BDS::NextMultiple(nPointsIn, 8); // ensure multiple of 8
      // The rhombus may be asymmetric, in which case building a pi/2 range of curve
      // won't result in a smooth transition to the straight sections. We need to work
      // out what range of angle to cover to come to the right tangent.
      // consider top curved point -> alpha = angle between the y-axis and (0,y1) -> (x1,0) line
      G4double alpha = std::atan2(std::abs(x), std::abs(y));
      G4double halfAngle = CLHEP::halfpi - alpha;
      std::vector<G4TwoVector> topBit;
      G4double nPointsTopDouble = (2*halfAngle/CLHEP::twopi) * nPointsIn;
      G4int nPointsTop = std::max(4, (G4int)std::ceil(nPointsTopDouble)); // ensure at least 4 points
      G4double currentAngle = -halfAngle;
      G4double dAngle = 2*halfAngle / (G4double)nPointsTop;
      G4TwoVector rotationPointTop(0, y-cornerRadius);
      for (G4int i = 0; i <= nPointsTop; i++)
        {
          G4TwoVector r(0, cornerRadius);
          r.rotate(-currentAngle);
          topBit.push_back(rotationPointTop + r);
          currentAngle += dAngle;
        }

      std::vector<G4TwoVector> rightBit;
      G4double nPointsRightDouble = (2*alpha/CLHEP::twopi) * nPointsIn;
      G4int nPointsRight = std::max(4, (G4int)std::ceil(nPointsRightDouble)); // ensure at least 4 points
      currentAngle = -alpha;
      dAngle = 2*alpha / (G4double)nPointsRight;
      G4TwoVector rotationPointRight(x-cornerRadius, 0);
      for (G4int i = 0; i <= nPointsRight; i++)
        {
          G4TwoVector r(cornerRadius, 0);
          r.rotate(-currentAngle);
          topBit.push_back(rotationPointRight + r);
          currentAngle += dAngle;
        }

      std::vector<G4TwoVector> result;
      for (const auto& p : topBit)
        {result.push_back(p);}
      for (const auto& p : rightBit)
        {result.push_back(p);}
      for (const auto& p : topBit)
        {result.emplace_back(-p.x(), -p.y());}
      for (const auto& p : rightBit)
        {result.emplace_back(-p.x(), -p.y());}
      std::reverse(result.begin(), result.end());
      return BDSPolygon(result);
    }
}