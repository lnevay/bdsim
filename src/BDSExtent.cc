#include "BDSExtent.hh"
#include "BDSTiltOffset.hh"
#include "BDSUtilities.hh"

#include "globals.hh" // geant4 types / globals
#include "G4TwoVector.hh"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

BDSExtent::BDSExtent():
  extXNeg(0.0),
  extXPos(0.0),
  extYNeg(0.0),
  extYPos(0.0),
  extZNeg(0.0),
  extZPos(0.0)
{;}

BDSExtent::BDSExtent(const std::pair<G4double, G4double>& extXIn,
		     const std::pair<G4double, G4double>& extYIn,
		     const std::pair<G4double, G4double>& extZIn):
  extXNeg(extXIn.first),
  extXPos(extXIn.second),
  extYNeg(extYIn.first),
  extYPos(extYIn.second),
  extZNeg(extZIn.first),
  extZPos(extZIn.second)
{;}

BDSExtent::BDSExtent(G4double extXNegIn, G4double extXPosIn,
		     G4double extYNegIn, G4double extYPosIn,
		     G4double extZNegIn, G4double extZPosIn):
  extXNeg(extXNegIn),
  extXPos(extXPosIn),
  extYNeg(extYNegIn),
  extYPos(extYPosIn),
  extZNeg(extZNegIn),
  extZPos(extZPosIn)
{;}

BDSExtent::BDSExtent(G4double extXIn, G4double extYIn, G4double extZIn):
  extXNeg(-std::abs(extXIn)),
  extXPos( std::abs(extXIn)),
  extYNeg(-std::abs(extYIn)),
  extYPos( std::abs(extYIn)),
  extZNeg(-std::abs(extZIn)),
  extZPos( std::abs(extZIn))
{;}

BDSExtent::BDSExtent(G4ThreeVector extIn):
  BDSExtent(extIn.x(), extIn.y(), extIn.z())
{;}

BDSExtent::~BDSExtent()
{;}

BDSExtent BDSExtent::TiltOffset(const BDSTiltOffset* tiltOffset) const
{
  if (!tiltOffset)
    {return BDSExtent(*this);}
  BDSExtent tilted = Tilted(tiltOffset->GetTilt());
  BDSExtent offset = tilted.Offset(tiltOffset->GetOffset());
  return offset;
}

BDSExtent BDSExtent::Offset(G4ThreeVector offset) const
{
  return Offset(offset.x(), offset.y(), offset.z());
}

BDSExtent BDSExtent::Offset(G4double dx, G4double dy, G4double dz) const
{
  return BDSExtent(extXNeg + dx, extXPos + dx,
		   extYNeg + dy, extYPos + dy,
		   extZNeg + dz, extZPos + dz);
}

BDSExtent BDSExtent::Tilted(G4double angle) const
{
  if (!BDS::IsFinite(angle))
    {return BDSExtent(*this);}

  // rotate each vector (from origin to each corner) by angle
  // and check - safer than checking based on +ve / -ve angle
  G4TwoVector topRight = G4TwoVector(extXPos, extYPos);
  G4TwoVector botRight = G4TwoVector(extXPos, extYNeg);
  G4TwoVector botLeft  = G4TwoVector(extXNeg, extYNeg);
  G4TwoVector topLeft  = G4TwoVector(extXNeg, extYPos);

  topRight.rotate(angle);
  botRight.rotate(angle);
  botLeft.rotate(angle);
  topLeft.rotate(angle);
  
  G4double xMin = std::min(botLeft.x(),  topLeft.x());
  G4double xMax = std::max(topRight.x(), botRight.x());
  G4double yMin = std::min(botRight.y(), botLeft.y());
  G4double yMax = std::max(topRight.y(), topLeft.y());

  BDSExtent result = BDSExtent(xMin, xMax,
			       yMin, yMax,
			       extZNeg, extZPos);
  return result;
}

std::ostream& operator<< (std::ostream& out, BDSExtent const& ext)
{
  out << ext.extXNeg << " " << ext.extXPos << " ";
  out << ext.extYNeg << " " << ext.extYPos << " ";
  out << ext.extZNeg << " " << ext.extZPos;
  return out;
}

BDSExtent BDSExtent::Shift(G4double x, G4double y) const
{
  BDSExtent xShifted = ShiftX(x);
  return xShifted.ShiftY(y);
}

BDSExtent BDSExtent::ShiftX(G4double x) const
{
  BDSExtent result = BDSExtent(extXNeg + x, extXPos + x,
			       extYNeg, extYPos,
			       extZNeg, extZPos);
  return result;
}

BDSExtent BDSExtent::ShiftY(G4double y) const
{
  BDSExtent result = BDSExtent(extXNeg, extXPos,
			       extYNeg + y, extYPos + y,
			       extZNeg, extZPos);
  return result;
}

G4double BDSExtent::MaximumAbs() const
{
  std::vector<G4double> exts = {std::abs(extXNeg), extXPos,
				std::abs(extYNeg), extYPos,
				std::abs(extZNeg), extZPos};
  return *std::max_element(exts.begin(), exts.end());
}

G4double BDSExtent::MaximumAbsTransverse() const
{
  std::vector<G4double> exts = {std::abs(extXNeg), extXPos,
				std::abs(extYNeg), extYPos};
  return *std::max_element(exts.begin(), exts.end());
}

G4bool BDSExtent::Encompasses(G4ThreeVector point) const
{
  BDSExtent extentPoint = BDSExtent(point);
  return extentPoint < (*this);
}

G4bool BDSExtent::Encompasses(G4double x,
			      G4double y,
			      G4double z) const
{
  return Encompasses(G4ThreeVector(x,y,z));
}
