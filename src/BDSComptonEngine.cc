/* BDSIM code.    Version 1.0
   Author: Grahame A. Blair, Royal Holloway, Univ. of London.
   Last modified 24.7.2002
   Copyright (c) 2002 by G.A.Blair.  ALL RIGHTS RESERVED. 
*/
//      ------------ BDSComptonEngine physics process --------
//                     by Grahame Blair, 19 October 2001

#include "BDSComptonEngine.hh"
#include "G4ios.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh" 
#include "CLHEP/Units/PhysicalConstants.h"

BDSComptonEngine::BDSComptonEngine(){}

 
BDSComptonEngine::BDSComptonEngine(G4LorentzVector InGam,
                                 G4LorentzVector InEl )
  : itsIncomingEl(InEl),itsIncomingGam(InGam)
{
 if(itsIncomingGam.e()<=0.)
   {G4Exception("BDSComptonEngine: Invalid Photon Energy", "-1", FatalException, "");}
} 
 
 
BDSComptonEngine::~BDSComptonEngine(){}


void BDSComptonEngine::PerformCompton()
{

  // Generate compton event; using method described in
  // H.Burkardt, SL/Note 93-73
  
  G4double phi=CLHEP::twopi * G4UniformRand() ;
  G4double sinphi=sin(phi);
  G4double cosphi=cos(phi);
  

  G4int ntry=0;
 
  G4LorentzVector GamInLab;
  G4double x;        
  G4ThreeVector Boost;
  G4double costh, costh2,sinth,sinth2;
  
  Boost = itsIncomingEl.boostVector();
   //   BoostToLab.boost(-Boost);
  G4LorentzRotation BoostToLab(-Boost);
  GamInLab=BoostToLab*(itsIncomingGam);
  G4ThreeVector LabGammaDir= GamInLab.vect().unit();
  
  G4double weight_CovT=0;  //ratio of Compton to Thompson cross sections:
  do {ntry++;
  // 1+cos^2 theta distribution
  if(G4UniformRand()>0.25){costh=2.*G4UniformRand()-1.;}  
  else
    {
      costh=G4UniformRand();
      G4double r1=G4UniformRand();
      G4double r2=G4UniformRand();
      if(r1>costh)costh=r1;
      if(r2>costh)costh=r2;
      if(G4UniformRand()<0.5)costh=-costh;
    }
  
  costh2=costh*costh;
  sinth2=1.-costh2;
  
  // x is ratio of scattered to unscattered photon energy:
  x = 1/(1+ GamInLab.e()*(1-costh)/CLHEP::electron_mass_c2);
  
  //calculate weight of Compton relative to Thompson cross sections:
  weight_CovT= x* x * (x+1/x-sinth2)/(1+costh2);
  } while(ntry<ntryMax && G4UniformRand()>weight_CovT);
  
  if(ntry==ntryMax)
    G4Exception("BDSComptonEngine:Max number of loops exceeded", "-1", FatalException, "");
  
  // G4LorentzVector ElInLab=BoostToLab*(itsIncomingEl);
  
  G4double Egam = x * GamInLab.e();
  
  // Generate final momenta:
  
  sinth=std::sqrt(sinth2);
  itsScatteredGam.setPx(Egam*sinth*cosphi);
  itsScatteredGam.setPy(Egam*sinth*sinphi);
  itsScatteredGam.setPz(Egam*costh);
  itsScatteredGam.setE(Egam);
  
  itsScatteredEl.setPx(-itsScatteredGam.px());
  itsScatteredEl.setPy(-itsScatteredGam.py());
  itsScatteredEl.setPz(GamInLab.e()-itsScatteredGam.pz());
  itsScatteredEl.setE( std::sqrt( std::pow(CLHEP::electron_mass_c2,2) +
			     std::pow(itsScatteredEl.px(),2)+
			     std::pow(itsScatteredEl.py(),2)+
			     std::pow(itsScatteredEl.pz(),2) ) );
  
  // Rotates the reference frame from Uz to newUz (unit vector).
  itsScatteredGam.rotateUz(LabGammaDir);
  itsScatteredEl.rotateUz(LabGammaDir);
  
  // now boost back to the original frame:
  G4LorentzRotation BoostFromLab(Boost);
  itsScatteredGam *= BoostFromLab;
  itsScatteredEl *= BoostFromLab;
}
