#ifndef BDSINTERPOLATOR1D_H
#define BDSINTERPOLATOR1D_H

#include "BDSFieldValue.hh"

#include "G4Types.hh"
#include "G4ThreeVector.hh"

class BDSArray1DCoords;

/**
 * @brief Interface for all 1D interpolators.
 *
 * Does not own data - so multiple interpolators could be used on same data.
 *
 * @author Laurie Nevay
 */

class BDSInterpolator1D
{
public:
  explicit BDSInterpolator1D(BDSArray1DCoords* arrayIn);
  virtual ~BDSInterpolator1D(){;}

  /// Public interface to a 1D interpolator. Returns Geant4 type as that's what will be needed.
  G4ThreeVector GetInterpolatedValue(G4double x) const;

  /// Alternative public interface to allow use of a 3-vector position (ignores z).
  /// Returns Geant4 type as that's what will be needed.
  G4ThreeVector GetInterpolatedValue(G4ThreeVector position) const;

  inline const BDSArray1DCoords* Array() const {return array;}

protected:
  /// Each derived class should implement this function. Note T suffix (was templated)
  /// to distinguish it from GetInterpolatedValue which returns Geant4 types and is
  /// the main interface to should have the clean name.
  virtual BDSFieldValue GetInterpolatedValueT(G4double x) const = 0;
  
  /// The field data.
  BDSArray1DCoords* array;
};

#endif
