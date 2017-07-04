#include "BDSIntegratorMultipoleThin.hh"
#include "BDSMagnetStrength.hh"
#include "BDSStep.hh"
#include "BDSUtilities.hh"

#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <list>
#include <vector>
#include <include/BDSGlobalConstants.hh>


BDSIntegratorMultipoleThin::BDSIntegratorMultipoleThin(BDSMagnetStrength const* strength,
						       G4double                 brhoIn,
						       G4Mag_EqRhs*             eqOfMIn):
  BDSIntegratorMag(eqOfMIn, 6),
  brho(brhoIn)
{
  b0l = (*strength)["field"] * brho;
  std::vector<G4String> normKeys = strength->NormalComponentKeys();
  std::vector<G4String> skewKeys = strength->SkewComponentKeys();
  std::vector<G4String>::iterator nkey = normKeys.begin();
  std::vector<G4String>::iterator skey = skewKeys.begin();
  for (G4double i = 0; i < normKeys.size(); i++, ++nkey, ++skey)
    {
      bnl.push_back((*strength)[*nkey] / std::pow(CLHEP::m,i+1));
      bsl.push_back((*strength)[*skey] / std::pow(CLHEP::m,i+1));
      nfact.push_back(Factorial(i));
    }
}

void BDSIntegratorMultipoleThin::Stepper(const G4double yIn[],
					 const G4double[] /*dydx[]*/,
					 const G4double h,
					 G4double       yOut[],
					 G4double       yErr[])
{
  G4ThreeVector pos    = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  G4ThreeVector mom    = G4ThreeVector(yIn[3], yIn[4], yIn[5]);
  G4double      momMag = mom.mag();
  
  // global to local
  BDSStep       localPosMom  = ConvertToLocal(pos, mom, h, false, thinElementLength);
  G4ThreeVector localPos     = localPosMom.PreStepPoint();
  G4ThreeVector localMom     = localPosMom.PostStepPoint();
  G4ThreeVector localMomUnit = localMom.unit();

  // only use for forward paraxial momenta, else advance particle as if in a drift
  if (localMomUnit.z() < 0.9)
    {
      AdvanceDriftMag(yIn, h, yOut);
      SetDistChord(0);
      return;
    }

  G4double x0  = localPos.x();
  G4double y0  = localPos.y();
  G4double z0  = localPos.z();
  G4double xp  = localMomUnit.x();
  G4double yp  = localMomUnit.y();
  G4double zp  = localMomUnit.z();

  // initialise output variables with input position as default
  G4double x1  = x0;
  G4double y1  = y0;
  G4double z1  = z0 + h; // new z position will be along z by step length h
  G4double xp1 = xp;
  G4double yp1 = yp;
  G4double zp1 = zp;

  // kicks come from pg 27 of mad-8 physical methods manual
  G4complex kick(0,0);
  G4complex position(x0, y0);
  G4complex result(0,0);
  // components of complex vector
  G4double knReal = 0;
  G4double knImag = 0;
  G4double momx;
  G4double momy;

  // normalise to momentum and charge
  G4double ratio = eqOfM->FCof() * brho / momMag;

  G4int n = 1;
  std::list<double>::iterator kn = bnl.begin();

  // sum higher order components into one kick
  for (; kn != bnl.end(); n++, ++kn)
    {
      momx = 0; //reset to zero
      momy = 0;
      knReal = (*kn) * ratio * std::pow(position,n).real() / nfact[n];
      knImag = (*kn) * ratio * std::pow(position,n).imag() / nfact[n];
      if (!std::isnan(knReal))
	{momx = knReal;}
      if (!std::isnan(knImag))
	{momy = knImag;}
      result = {momx,momy};
      kick += result;
    }

  // apply kick
  xp1 -= kick.real();
  yp1 += kick.imag();
  zp1 = std::sqrt(1 - std::pow(xp1,2) - std::pow(yp1,2));
  if (std::isnan(zp1))
    {zp1 = zp;}

  // reset n for skewed kicks.
  n=1;
  G4double ksReal = 0;
  G4double ksImag = 0;

  G4ThreeVector momOut = G4ThreeVector(xp1,yp1,zp1);
  G4complex skewkick(0,0);

  std::list<double>::iterator ks = bsl.begin();
  for (; ks != bsl.end(); n++, ++ks)
    {
      if (BDS::IsFinite(*ks))
        {
          //reset to zero
          momx = 0;
          momy = 0;
          ksReal = (*ks) * ratio * std::pow(position, n).real() / nfact[n];
          ksImag = (*ks) * ratio * std::pow(position, n).imag() / nfact[n];
          if (!std::isnan(ksReal))
            {momx = ksReal;}
          if (!std::isnan(ksImag))
            {momy = ksImag;}
          result = {momx,momy};
          skewkick += result;
        }
    }
  
  //apply kick
  xp1 -= skewkick.imag();
  yp1 += skewkick.real();
  zp1 = std::sqrt(1 - std::pow(xp1,2) - std::pow(yp1,2));
  if (std::isnan(zp1))
    {zp1 = zp;}

  // xp1 or yp1 may be > 1, so isnan check also needed for zp1.
  if (std::isnan(zp1) || (zp1 < 0.9))
    {
      AdvanceDriftMag(yIn, h, yOut);
      SetDistChord(0);
      return;
    }
  
  G4ThreeVector localPosOut     = G4ThreeVector(x1, y1, z1);
  G4ThreeVector localMomUnitOut = G4ThreeVector(xp1, yp1, zp1);
  ConvertToGlobal(localPosOut, localMomUnitOut, yOut, momMag);

  for (G4int i = 0; i < nVariables; i++)
    {yErr[i] = 0;}
}

G4int BDSIntegratorMultipoleThin::Factorial(G4int n)
{
  G4int result = 1;
  for (G4int i = 1; i <= n; i++)
    {result *= i;}
  return result;
}
