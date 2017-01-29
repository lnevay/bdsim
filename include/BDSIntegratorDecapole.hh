#ifndef BDSINTEGRATORDECAPOLE_H
#define BDSINTEGRATORDECAPOLE_H

#include "BDSIntegratorBase.hh"
#include "BDSMagnetStrength.hh"

#include "globals.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

/**
 * @brief Integrator for Decapolar field.
 * 
 * @author Laurie Nevay
 */

class BDSIntegratorDecapole: public BDSIntegratorBase
{
public:
  BDSIntegratorDecapole(BDSMagnetStrength const* strength,
			G4double                 brho,
			G4Mag_EqRhs*             eqOfMIn);
  
  virtual ~BDSIntegratorDecapole(){;}

  /// The stepper for integration. The stepsize is fixed, equal to h. The reason for this
  /// is so that intermediate steps can be calculated and therefore the error ascertained
  /// or distance from the chord.  Error calculation is not currently implemented.
  virtual void Stepper(const G4double y[],
		       const G4double dydx[],
		       const G4double h,
		       G4double yout[],
		       G4double yerr[]);

protected:
  /// Calculate the new particle coordinates. A first order Step along a solenoid inside the field.
  void AdvanceHelix(const G4double  yIn[],
		    G4ThreeVector Bfld,
		    G4double  h,
		    G4double  yDec[]);

private:
  /// Private default constructor to enforce use of supplied constructor
  BDSIntegratorDecapole();

  /// 4th derivative of magnetic field
  G4double bQuadruplePrime;

  /// Data stored in order to find the chord.
  G4ThreeVector yInitial, yMidPoint, yFinal;
};

#endif
