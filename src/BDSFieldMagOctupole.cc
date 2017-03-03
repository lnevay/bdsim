#include "BDSDebug.hh"
#include "BDSFieldMagOctupole.hh"
#include "BDSMagnetStrength.hh"

#include "G4ThreeVector.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

BDSFieldMagOctupole::BDSFieldMagOctupole(BDSMagnetStrength const* strength,
					 G4double          const   brho)
{
  // B''' = d^3By/dx^3 = Brho * (1/Brho d^3By/dx^3) = Brho * k3
  bTriplePrime = brho * (*strength)["k3"]  / (CLHEP::m3*CLHEP::m);
  bTPNormed    = bTriplePrime / 6.; 
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "B''' = " << bTriplePrime << G4endl;
#endif
}

G4ThreeVector BDSFieldMagOctupole::GetField(const G4ThreeVector &position,
					    const G4double       /*t*/) const
{
  // B_x = (3x^2y - y^3) * (B'''/3!)
  // B_y = (x^3-3xy^2) * (B'''/3!)
  // B_z = 0

  //shortcuts
  G4double x = position.x();
  G4double y = position.y();
  
  G4ThreeVector localField;
  localField[0] = (3 * std::pow(x,2) * y - std::pow(y,3)) * bTPNormed;
  localField[1] = (std::pow(x,3) - 3* x * std::pow(y,2)) * bTPNormed;
  localField[2] = 0;

  return localField;
}
