#include "BDSFieldEMInterpolated1D.hh"
#include "BDSInterpolator1D.hh"

#include "G4ThreeVector.hh"

#include <utility>

BDSFieldEMInterpolated1D::BDSFieldEMInterpolated1D(BDSInterpolator1D* eInterpolatorIn,
						   BDSInterpolator1D* bInterpolatorIn,
						   G4Transform3D      offset,
						   G4double           eScalingIn,
						   G4double           bScalingIn):
  BDSFieldEMInterpolated(offset, eScalingIn, bScalingIn),
  eInterpolator(eInterpolatorIn),
  bInterpolator(bInterpolatorIn)
{;}

BDSFieldEMInterpolated1D::~BDSFieldEMInterpolated1D()
{
  delete eInterpolator;
  delete bInterpolator;
}

std::pair<G4ThreeVector,G4ThreeVector> BDSFieldEMInterpolated1D::GetField(const G4ThreeVector& position,
									  const G4double       /*t*/) const
{
  G4ThreeVector e = eInterpolator->GetInterpolatedValue(position[0]) * EScaling();
  G4ThreeVector b = bInterpolator->GetInterpolatedValue(position[0]) * BScaling();
  return std::make_pair(b,e);
}
