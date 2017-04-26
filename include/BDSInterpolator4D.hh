#ifndef BDSINTERPOLATOR4D_H
#define BDSINTERPOLATOR4D_H

#include "BDSFieldValue.hh"

#include "G4Types.hh"
#include "G4ThreeVector.hh"

class BDSArray4DCoords;

/**
 * @brief Interface for all 4D interpolators.
 *
 * Does not own data - so multiple interpolators could be used on same data.
 *
 * @author Laurie Nevay
 */

class BDSInterpolator4D
{
public:
  explicit BDSInterpolator4D(BDSArray4DCoords* arrayIn);
  virtual ~BDSInterpolator4D(){;}

  /// Public interface to a 4D interpolator. Returns Geant4 type as that's what will be needed.
  G4ThreeVector GetInterpolatedValue(G4double x, G4double y, G4double z, G4double t) const;

  /// Alternative public interface to allow use of a 3-vector position.
  /// Returns Geant4 type as that's what will be needed.
  G4ThreeVector GetInterpolatedValue(G4ThreeVector position, G4double t) const;

  inline const BDSArray4DCoords* Array() const {return array;}
  
protected:
  /// Each derived class should implement this function. Note T suffix (was templated)
  /// to distinguish it from GetInterpolatedValue which returns Geant4 types and is
  /// the main interface to should have the clean name.
  virtual BDSFieldValue GetInterpolatedValueT(G4double x,
					      G4double y,
					      G4double z,
					      G4double t) const = 0;

  /// The field data.
  BDSArray4DCoords* array;
};

#endif
