/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "BDSDebug.hh"
#include "BDSFieldMagDipole.hh"
#include "BDSGlobalConstants.hh"
#include "BDSIntegratorDipoleFringe.hh"
#include "BDSMagnetStrength.hh"
#include "BDSStep.hh"
#include "BDSUtilities.hh"

#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

#include <cmath>

G4double BDSIntegratorDipoleFringe::thinElementLength = -1; // unphysical

BDSIntegratorDipoleFringe::BDSIntegratorDipoleFringe(BDSMagnetStrength const* strengthIn,
						     G4double                 brhoIn,
						     G4Mag_EqRhs*             eqOfMIn,
						     G4double                 minimumRadiusOfCurvatureIn,
						     const G4double&          tiltIn):
  BDSIntegratorDipoleRodrigues2(eqOfMIn, minimumRadiusOfCurvatureIn),
  polefaceCurvature((*strengthIn)["polefacecurv"]),
  rho(std::abs(brhoIn)/(*strengthIn)["field"]),
  fieldArcLength((*strengthIn)["length"]),
  fieldAngle((*strengthIn)["angle"]),
  tilt(tiltIn),
  finiteTilt(BDS::IsFinite(tiltIn))
{
  if (thinElementLength < 0)
    {thinElementLength = BDSGlobalConstants::Instance()->ThinElementLength();}

  if ((*strengthIn)["isentrance"])
    {
      polefaceAngle = (*strengthIn)["e1"];
      fringeCorr = BDS::FringeFieldCorrection(strengthIn, 1);
      secondFringeCorr = BDS::SecondFringeFieldCorrection(strengthIn, 1);
    }
  else  // must be exit face
    {
      polefaceAngle = (*strengthIn)["e2"];
      fringeCorr = BDS::FringeFieldCorrection(strengthIn, 0);
      secondFringeCorr = BDS::SecondFringeFieldCorrection(strengthIn, 0);
    }

  // thin sextupole strength for curved polefaces
  G4double thinSextStrength = (-polefaceCurvature / rho) * 1.0 / std::pow(std::cos(polefaceAngle),3);

  // prepare magnet strength object for thin sextupole needed for curved polefaces
  BDSMagnetStrength* sextStrength = new BDSMagnetStrength();
  std::vector<G4double> normalComponent;
  normalComponent.push_back(0); // zero strength quadrupolar component
  normalComponent.push_back(thinSextStrength); // sextupolar component
  std::vector<G4String> normKeys = sextStrength->NormalComponentKeys();
  auto kn = normalComponent.begin();
  auto nkey = normKeys.begin();
  for (; kn != normalComponent.end(); kn++, nkey++)
    {(*sextStrength)[*nkey] = (*kn);}
  // integrator for thin sextupole
  multipoleIntegrator = new BDSIntegratorMultipoleThin(sextStrength, brhoIn, eqOfMIn);
  delete sextStrength;
    
  zeroStrength = !BDS::IsFinite((*strengthIn)["field"]); // no fringe if no field
  BDSFieldMagDipole* dipoleField = new BDSFieldMagDipole(strengthIn);
  unitField = (dipoleField->FieldValue()).unit();
  delete dipoleField;
}

void BDSIntegratorDipoleFringe::Stepper(const G4double yIn[6],
					const G4double dydx[6],
					const G4double h,
					G4double       yOut[6],
					G4double       yErr[6])
{
  // unit normalisation
  const G4double fcof = eqOfM->FCof();

  BaseStepper(yIn, dydx, h, yOut, yErr, fcof, 1.0);
}

void BDSIntegratorDipoleFringe::BaseStepper(const G4double  yIn[6],
					    const G4double  dydx[6],
					    const G4double& h,
					    G4double        yOut[6],
					    G4double        yErr[6],
					    const G4double& fcof,
					    const G4double& momScaling)
{
  // Protect against neutral particles, and zero field: drift through.
  if (!BDS::IsFinite(fcof) || zeroStrength)
    {
      AdvanceDriftMag(yIn,h,yOut,yErr);
      FudgeDistChordToZero(); // see doxygen in header for explanation
      return;
    }

  // container for multipole kick output, used as dipole step input
  G4double yMultipoleOut[7];
  // copy input coords as initials as multipole kick method not called
  for (G4int i = 0; i < 3; i++)
    {
      yMultipoleOut[i]     = yIn[i];
      yMultipoleOut[i + 3] = yIn[i + 3];
    }
  // apply thin multipole kick if finite curvature. Does not step, stepping occurs in dipole integrator
  if (BDS::IsFinite(polefaceCurvature))
    {MultipoleStep(yIn, yMultipoleOut, h);}

  // container for dipole step output, used as fringe step input
  G4double yTemp[7];

  // do the dipole kick and step using base class
  BDSIntegratorDipoleRodrigues2::Stepper(yMultipoleOut, dydx, h, yTemp, yErr); // yErr is correct output variable

  // only apply the kick if we're taking a step longer than half the length of the item,
  // in which case, apply the full kick. This appears more robust than scaling the kick
  // by h / thinElementLength as the precise geometrical length depends on the geometry
  // ie if there's a beam pipe etc -> more length safetys.  The geometry layout should
  // prevent more than one step begin taken, but occasionally, a very small initial step
  // can be taken resulting in a double kick.
  G4double lengthFraction = h / thinElementLength;

  // don't do fringe kick if we're sampling the field for a long step
  // or if it's a half step inside the thin element apply the dipole
  // motion but not the one-off fringe kick
  if ((h > 1*CLHEP::cm) || (lengthFraction < 0.501))
    {
      // copy output from dipole kick output
      for (G4int i = 0; i < 3; i++)
        {
          yOut[i]     = yTemp[i];
          yOut[i + 3] = yTemp[i + 3];
        }
      return;
    }

  // position and momentum post dipole kick.
  G4ThreeVector pos = G4ThreeVector(yTemp[0], yTemp[1], yTemp[2]);
  G4ThreeVector mom = G4ThreeVector(yTemp[3], yTemp[4], yTemp[5]);

  BDSStep  localPosMom    = GlobalToCurvilinear(pos, mom, h, true);
  if (finiteTilt)
    {localPosMom = localPosMom.rotateZ(-tilt);}
  G4ThreeVector localPos  = localPosMom.PreStepPoint();
  G4ThreeVector localMom  = localPosMom.PostStepPoint();
  G4ThreeVector localMomU = localMom.unit();

  // check for forward going paraxial particles - only
  if (localMomU.z() < 0.9)
    {// copy output from dipole kick output
      for (G4int i = 0; i < 3; i++)
	{
	  yOut[i]     = yTemp[i];
	  yOut[i + 3] = yTemp[i + 3];
	}
      return; // note distchord comes from inherited BDSIntegratorDipoleRodrigues2
    }

  // calculate new position and momentum kick from fringe effect only
  G4ThreeVector localCLPosOut;
  G4ThreeVector localCLMomOut;

  // normalise to particle charge
  G4double charge = fcof / std::abs(fcof);
  G4double bendingRad = rho * charge;

  // normalise to momentum - passed in from derived class.
  bendingRad *= momScaling;

  // apply fringe field kick
  OneStep(localPos, localMomU, localCLPosOut, localCLMomOut, bendingRad);

  // convert to global coordinates for output
  if (finiteTilt)
    {
      //localCLPosOut = localCLPosOut.rotateZ(tilt); // we only need the momentum - ignore
      localCLMomOut = localCLMomOut.rotateZ(tilt);
    }
  
  G4ThreeVector globalMom = ConvertAxisToGlobal(localCLMomOut, true);
  globalMom *= mom.mag();
  // error along direction of travel really
  G4ThreeVector globalMomU = globalMom.unit();
  globalMomU *= 1e-8;

  // write out values and errors
  for (G4int i = 0; i < 3; i++)
    {
      yOut[i]     = pos[i];
      yOut[i + 3] = globalMom[i];
      yErr[i]     = globalMomU[i];
      yErr[i + 3] = 1e-40;
    }
}

void BDSIntegratorDipoleFringe::OneStep(const G4ThreeVector& posIn,
					const G4ThreeVector& momUIn,
					G4ThreeVector&       posOut,
					G4ThreeVector&       momOut,
					const G4double&      bendingRadius) const
{
  G4double x0 = posIn.x() / CLHEP::m;
  G4double y0 = posIn.y() / CLHEP::m;
  G4double s0 = posIn.z();
  G4double xp = momUIn.x();
  G4double yp = momUIn.y();
  G4double zp = momUIn.z();

  G4double x1 = x0;
  G4double y1 = y0;
  G4double s1 = s0;
  G4double xp1 = xp;
  G4double yp1 = yp;
  G4double zp1 = zp;

  G4double X11=0,X12=0,X21=0,X22 = 0;
  G4double Y11=0,Y12=0,Y21=0,Y22 = 0;

  // fringe field correction using normalised rho
  G4double fringeFieldCorrection = (fringeCorr/bendingRadius) * (1.0 - secondFringeCorr/bendingRadius);

  // calculate fringe field kick
  X11 = 1;
  X21 = std::tan(polefaceAngle) / (bendingRadius / CLHEP::m);
  X22 = 1;

  Y11 = 1;
  Y21 = -std::tan(polefaceAngle - fringeFieldCorrection) / (bendingRadius / CLHEP::m);
  Y22 = 1;

  x1  = X11*x0 + X12*xp;
  xp1 = X21*x0 + X22*xp;
  y1  = Y11*y0 + Y12*yp;
  yp1 = Y21*y0 + Y22*yp;

  // relies on normalised momenta otherwise this will be nan.
  zp1 = std::sqrt(1 - xp1*xp1 - yp1*yp1);
  if (std::isnan(zp1))
    {zp1 = zp;} // ensure not nan

  momOut = G4ThreeVector(xp1, yp1, zp1);
  posOut = G4ThreeVector(x1*CLHEP::m, y1*CLHEP::m, s1);
}

void BDSIntegratorDipoleFringe::MultipoleStep(const G4double  yIn[6],
                                              G4double        yMultipoleOut[7],
                                              const G4double& h)
{
  // convert to local curvilinear co-ordinates
  G4ThreeVector pos = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  G4ThreeVector mom = G4ThreeVector(yIn[3], yIn[4], yIn[5]);

  BDSStep localPosMom     = GlobalToCurvilinear(pos, mom, h, true);
  G4ThreeVector localPos  = localPosMom.PreStepPoint();
  G4ThreeVector localMom  = localPosMom.PostStepPoint();
  G4ThreeVector localMomU = localMom.unit();

  // initialise new position and momentum for curvature effect as input coords
  G4ThreeVector localCLPosOut = localPos;
  G4ThreeVector localCLMomOut = localMomU;

  // Zero step length as stepping has already been applied in base class stepper
  multipoleIntegrator->OneStep(localPos, localMomU, localMom.mag(), localCLPosOut, localCLMomOut, 0);

  G4ThreeVector globalMom = ConvertAxisToGlobal(localCLMomOut, true);
  globalMom *= mom.mag();

  // write out values
  for (G4int i = 0; i < 3; i++)
    {
      yMultipoleOut[i]     = pos[i];
      yMultipoleOut[i + 3] = globalMom[i];
    }
}

// fringe field correction terms are BDS namespace methods instead of class methods
// as the functions are called during BDSIntegratorDipoleFringe construction.
// The fringe field correction terms should be cached at construction rather than
// calculated every time the integrator is called.

G4double BDS::FringeFieldCorrection(BDSMagnetStrength const* strength,
                                    G4bool isEntrance)
{
  G4double fint;
  G4double pfAngle;
  if (isEntrance)
    {
      fint = (*strength)["fint"];
      pfAngle = (*strength)["e1"];
    }
  else // must be exit face
    {
      fint = (*strength)["fintx"];
      pfAngle = (*strength)["e2"];
    }
  G4double hgap = (*strength)["hgap"];
  G4double vertGap = 2 * hgap;
  G4double corrValue = fint * vertGap * (1.0 + std::pow(sin(pfAngle),2)) / cos(pfAngle);
  return corrValue;
}

G4double BDS::SecondFringeFieldCorrection(BDSMagnetStrength const* strength,
                                          G4bool isEntrance)
{
  G4double fint;
  G4double fintK2;
  G4double pfAngle;
  if (isEntrance)
    {
      fint    = (*strength)["fint"];
      fintK2  = (*strength)["fintk2"];
      pfAngle = (*strength)["e1"];
    }
  else // must be exit face
    {
      fint    = (*strength)["fintx"];
      fintK2  = (*strength)["fintxk2"];
      pfAngle = (*strength)["e2"];
    }
  G4double hgap    = (*strength)["hgap"];
  G4double vertGap = 2 * hgap;
  G4double corrValue = fint * fintK2 * vertGap * tan(pfAngle);
  return corrValue;
}

