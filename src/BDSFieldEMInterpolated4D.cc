#include "BDSFieldEMInterpolated4D.hh"
#include "BDSInterpolator4D.hh"

#include "G4ThreeVector.hh"

#include <utility>

BDSFieldEMInterpolated4D::BDSFieldEMInterpolated4D(BDSInterpolator4D* eInterpolatorIn,
						   BDSInterpolator4D* bInterpolatorIn,
						   G4Transform3D      offset,
						   G4double           eScalingIn,
						   G4double           bScalingIn):
  BDSFieldEMInterpolated(offset, eScalingIn, bScalingIn),
  eInterpolator(eInterpolatorIn),
  bInterpolator(bInterpolatorIn)
{;}

BDSFieldEMInterpolated4D::~BDSFieldEMInterpolated4D()
{
  delete eInterpolator;
  delete bInterpolator;
}

std::pair<G4ThreeVector,G4ThreeVector> BDSFieldEMInterpolated4D::GetField(const G4ThreeVector& position,
									  const G4double       t) const
{
  G4ThreeVector e = eInterpolator->GetInterpolatedValue(position[0],position[1],position[2],t) * EScaling();
  G4ThreeVector b = bInterpolator->GetInterpolatedValue(position[0],position[1],position[2],t) * BScaling();
  return std::make_pair(b,e);
}
