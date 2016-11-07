#ifndef BDSFIELDMAGOUTERMULTIPOLE_H
#define BDSFIELDMAGOUTERMULTIPOLE_H

#include "BDSFieldMag.hh"

#include "globals.hh" // geant4 types / globals
#include "G4ThreeVector.hh"

/**
 * @brief A simple parameterisation of N-Pole outer yoke magnetic field.
 *
 */

class BDSFieldMagOuterMultipole: public BDSFieldMag
{
public:
  BDSFieldMagOuterMultipole(const G4int    nPolesIn,
			    const G4double polField,
			    const G4double phiOffset);

  virtual ~BDSFieldMagOuterMultipole(){;}

  /// Access the field value.
  virtual G4ThreeVector GetField(const G4ThreeVector &position,
				 const double         t = 0) const;

private:
  const G4int    nPoles;
  const G4double fieldStrength;
  const G4double phiOffset;
  G4double itsSectorPhi;
};

#endif
