#include "BDSBunchCircle.hh"
#include "BDSDebug.hh"

#include "parser/options.h"

#include "Randomize.hh"
#include "CLHEP/Units/PhysicalConstants.h"

BDSBunchCircle::BDSBunchCircle():
  BDSBunch(), envelopeR(0.0), envelopeRp(0.0), envelopeT(0.0), envelopeE(0.0)
{
#ifdef BDSDEBUG 
  G4cout << __METHOD_NAME__ << G4endl;
#endif

  FlatGen  = new CLHEP::RandFlat(*CLHEP::HepRandom::getTheEngine());
}

BDSBunchCircle::~BDSBunchCircle() 
{
  delete FlatGen;
}

void BDSBunchCircle::SetOptions(const GMAD::Options& opt,
				G4Transform3D beamlineTransformIn)
{
  BDSBunch::SetOptions(opt, beamlineTransformIn);
  SetEnvelopeR(opt.envelopeR); 
  SetEnvelopeRp(opt.envelopeRp);
  SetEnvelopeT(opt.envelopeT);
  SetEnvelopeE(opt.envelopeE); 
}

void BDSBunchCircle::GetNextParticle(G4double& x0, G4double& y0, G4double& z0, 
				     G4double& xp, G4double& yp, G4double& zp,
				     G4double& t , G4double&  E, G4double& weight)
{
#ifdef BDSDEBUG 
  G4cout << __METHOD_NAME__ << G4endl;
#endif
  x0 = X0 * CLHEP::m;
  y0 = Y0 * CLHEP::m;
  z0 = Z0 * CLHEP::m;
  xp = Xp0 * CLHEP::rad;
  yp = Yp0 * CLHEP::rad;
  z0 = Z0 * CLHEP::m + (T0 - envelopeT * (1.-2.*FlatGen->shoot())) * CLHEP::c_light * CLHEP::s;

  G4double phiR  = FlatGen->shoot()*2*CLHEP::pi;
  G4double phiRp = FlatGen->shoot()*2*CLHEP::pi;
  G4double r     = FlatGen->shoot()*envelopeR;
  G4double rp    = FlatGen->shoot()*envelopeRp; 

  x0 += cos(phiR) * r * CLHEP::m;
  y0 += sin(phiR) * r * CLHEP::m;
  xp += cos(phiRp) * rp * CLHEP::rad;
  yp += sin(phiRp) * rp * CLHEP::rad;
  
  zp = CalculateZp(xp,yp,Zp0);

  ApplyTransform(x0,y0,z0,xp,yp,zp);
  
  t = 0.0;
  E = E0 * CLHEP::GeV * (1 + envelopeE * (1-2*FlatGen->shoot()));

  weight = 1.0; 
}
