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
#include "BDSDebug.hh"
#include "BDSException.hh"
#include "BDSExtent.hh"
#include "BDSGlobalConstants.hh"
#include "BDSPolygon.hh"
#include "BDSTiltOffset.hh"

#include "G4TwoVector.hh"
#include "G4Types.hh"

#include <algorithm>
#include <vector>

BDSPolygon::BDSPolygon(const std::vector<G4TwoVector>& pointsIn):
  points(pointsIn),
  extent(nullptr)
{
  if (points.size() < 3)
    {throw BDSException(__METHOD_NAME__, "polygon must have at least 3 x,y points to be valid.");}
}

BDSPolygon::~BDSPolygon()
{
  delete extent;
}

G4bool BDSPolygon::Inside(const G4TwoVector& point) const
{
  G4bool result = true;
  for (G4int i = 0; i < (G4int)size() - 1; i++)
    {
      G4TwoVector test = points[i] - point;
      G4TwoVector norm = (points[i+1] - points[i]).orthogonal();
      result &= test.dot(norm) > 0;
    }
  return result;
}

G4bool BDSPolygon::Inside(const BDSPolygon& other) const
{
  // test all segments don't intersect -> doesn't overlap
  // test one point is inside
  G4bool result = true;
  for (G4int i = 0; i < (G4int)size() - 1; i++)
    {
      for (G4int j = 0; j < (G4int)other.size() - 1; j++)
	{result &= !SegmentsIntersect(points[i], points[i+1], other.points[j], other.points[j+1]);}
    }
  if (!result) // don't bother with inside point check
    {return result;}

  // test one point is inside.
  result &= other.Inside(points[0]);
  return result;
}

{
  G4bool result = false;
  for (G4int i = 0; i < (G4int)size() - 1; i++)
    {
      for (G4int j = i + 1; j < (G4int)size() - 1; j++)
	{
	  result |= SegmentsIntersect(points[i], points[i+1], points[j], points[j+1]);
	  if (result)
	    {// for optional feedback
	      if (indexOfIntersectionA)
		{*indexOfIntersectionA = i;}
	      if (indexOfIntersectionB)
		{*indexOfIntersectionB = j;}
	    }
	}
    }
  return result;
}

BDSPolygon BDSPolygon::ApplyTiltOffset(const BDSTiltOffset& to) const
{
  BDSPolygon r(*this);
  
  G4double    t      = to.Tilt();
  G4TwoVector offset = to.GetOffset2D();
  for (auto& p : r)
    {
      p.rotate(t);
      p += offset;
    }
  return r;
}

BDSExtent BDSPolygon::Extent() const
{
  if (extent)
    {return BDSExtent(*extent);}

  G4double extXNeg = 0;
  G4double extXPos = 0;
  G4double extYNeg = 0;
  G4double extYPos = 0;
  for (const auto& p : points)
    {
      extXNeg = std::min(extXNeg, p.x());
      extXPos = std::max(extXPos, p.x());
      extYNeg = std::min(extYNeg, p.y());
      extYPos = std::max(extYPos, p.y());
    }

  G4double lsl = BDSGlobalConstants::Instance()->LengthSafetyLarge();
  extent = new BDSExtent(extXNeg, extXPos,
			 extYNeg, extYPos,
			 -lsl,    lsl);
  return BDSExtent(*extent);
}

G4bool BDSPolygon::SegmentsIntersect(const G4TwoVector& p1,
				     const G4TwoVector& p2,
				     const G4TwoVector& q1,
				     const G4TwoVector& q2) const
{return true;}
