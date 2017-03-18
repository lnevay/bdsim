#ifndef BDSFIELDEMINTERPOLATED4D_H
#define BDSFIELDEMINTERPOLATED4D_H

#include "BDSFieldEMInterpolated.hh"

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include <utility>

class BDSInterpolator4D;

/**
 * @brief A 4D field from an interpolated array with any interpolation.
 *
 * Owns interpolators.
 *
 * This provides a simple interface for em fields to use a 4D
 * interpolator irrespective of which type of interpolator it is.
 *
 * @author Laurie Nevay
 */

class BDSFieldEMInterpolated4D: public BDSFieldEMInterpolated
{
public:
  BDSFieldEMInterpolated4D(BDSInterpolator4D* eInterpolatorIn,
			   BDSInterpolator4D* bInterpolatorIn,
  			   G4Transform3D      offset     = G4Transform3D::Identity,
			   G4double           eScalingIn = 1.0,
			   G4double           bScalingIn = 1.0);

  virtual ~BDSFieldEMInterpolated4D();

  /// Return the interpolated field value at a given point.
  virtual std::pair<G4ThreeVector,G4ThreeVector> GetField(const G4ThreeVector& position,
							  const G4double       t = 0) const;

  inline const BDSInterpolator4D* EInterpolator() const {return eInterpolator;}
  inline const BDSInterpolator4D* BInterpolator() const {return bInterpolator;}

private:
  /// Private default constructor to force use of provided one.
  BDSFieldEMInterpolated4D() = delete;

  BDSInterpolator4D* eInterpolator; ///< E Interpolator the field is based on.
  BDSInterpolator4D* bInterpolator; ///< B Interpolator the field is based on.
};

#endif
