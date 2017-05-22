#include "BDSIntegratorKickerThin.hh"
#include "BDSMagnetStrength.hh"
#include "BDSStep.hh"
#include "BDSUtilities.hh"

#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

#include <cmath>

BDSIntegratorKickerThin::BDSIntegratorKickerThin(BDSMagnetStrength const* strength,
						 G4double                 brhoIn,
						 G4Mag_EqRhs*             eqOfMIn):
  BDSIntegratorMag(eqOfMIn, 6),
  eqOfM(eqOfMIn),
  hkick((*strength)["hkick"]),
  vkick((*strength)["vkick"]),
  brho(brhoIn)
{
  zeroStrength = (!BDS::IsFinite(hkick) && !BDS::IsFinite(vkick));
}

void BDSIntegratorKickerThin::Stepper(const G4double   yIn[],
				      const G4double[] /*dydx[]*/,
				      const G4double   h,
				      G4double         yOut[],
				      G4double         yErr[])
{
  G4ThreeVector pos    = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  G4ThreeVector mom    = G4ThreeVector(yIn[3], yIn[4], yIn[5]);

  const G4double fcof = eqOfM->FCof();
  if (fcof == 0 || zeroStrength)
    {// neutral particle - drift through
      AdvanceDriftMag(yIn, h, yOut);
      SetDistChord(0);
      return;
    }
  
  // global to local (thinElementLength is static base class member)
  // purposively control step length as we're dealing with a very short element
  BDSStep       localPosMom  = ConvertToLocal(pos, mom, h, false, thinElementLength);
  G4ThreeVector localPos     = localPosMom.PreStepPoint();
  G4ThreeVector localMom     = localPosMom.PostStepPoint();
  G4ThreeVector localMomUnit = localMom.unit();

  // only use for paraxial momenta, else advance particle as if in a drift
  if (std::abs(localMomUnit.z()) < 0.9)
    {
      AdvanceDriftMag(yIn, h, yOut);
      SetDistChord(0);
      return;
    }
  
  // Transverse position stays the same and the particle is advanced by h along local z.
  G4double x1 = localPos.x();
  G4double y1 = localPos.y();
  G4double z1 = localPos.z();

  G4double factor = 1; // for forwards / backwards action
  if (localMomUnit.z() < 0)
    {// backwards travelling
      z1 -= h;
      factor = -1;
    }
  else
    {z1 += h;} // forwards travelling

  G4double px = localMom.x();
  G4double py = localMom.y();
  G4double pz = localMom.z();

  // Scale the action for the current particle momentum w.r.t.
  // the design momentum. Even though a thin kicker, it represents
  // a thin dipole, which would affect different particles differently.
  G4double localMomMag = localMom.mag();
  G4double ratio = (fcof * brho ) / localMomMag;

  G4double dpx = ratio * hkick * localMomMag * factor;
  G4double dpy = ratio * vkick * localMomMag * factor;

  px += dpx;
  py += dpy;
  pz -= std::abs(dpx) + std::abs(dpy); // conservation of momentum
  
  G4ThreeVector localPosOut = G4ThreeVector(x1, y1, z1);
  G4ThreeVector localMomOut = G4ThreeVector(px, py, pz);
  
  ConvertToGlobal(localPosOut, localMomOut, yOut);

  for (G4int i = 0; i < nVariables; i++)
    {yErr[i] = 0;}
}
