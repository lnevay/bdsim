#include "BDSIntegratorEulerOld.hh"
#include "BDSUtilities.hh"

#include "globals.hh"
#include "G4AffineTransform.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <stdexcept>

BDSIntegratorEulerOld::BDSIntegratorEulerOld(G4Mag_EqRhs* eqOfMIn):
  BDSIntegratorMag(eqOfMIn, 6)
{;}

void BDSIntegratorEulerOld::Stepper(const G4double yIn[],
				    const G4double dydx[],
				    const G4double h,
				    G4double       yOut[],
				    G4double       yErr[])
{
  G4double yTemp[7];

  G4ThreeVector pos     = G4ThreeVector(yIn[0], yIn[1], yIn[2]);
  G4ThreeVector mom     = G4ThreeVector(yIn[3], yIn[4], yIn[5]);
  G4ThreeVector momUnit = mom.unit();

  auxNavigator->LocateGlobalPointAndSetup(pos);
  G4AffineTransform GlobalAffine = auxNavigator->GetGlobalToLocalTransform();
  G4ThreeVector     localMomUnit = GlobalAffine.TransformAxis(momUnit);
  
  if (localMomUnit.z() < 0.9 || mom.mag() < 40.0)
    {
      backupStepper->Stepper(yIn, dydx, h, yOut, yErr);
      SetDistChord(backupStepper->DistChord());
      return;
    }

  try
    {
      // do two half steps
      AdvanceHelix(yIn,   0.5*h, yTemp);
      AdvanceHelix(yTemp, 0.5*h, yOut); 
      
      // do a full step
      AdvanceHelix(yIn, h, yTemp);

      // Error is difference between two half steps and full step
      for(G4int i = 0; i < nVariables; i++)
	{yErr[i] = yOut[i] - yTemp[i];}
    }
  catch (const std::exception&)
    {// thrown in case we should use a runge kutta
      backupStepper->Stepper(yIn, dydx, h, yOut, yErr);
      SetDistChord(backupStepper->DistChord());
    }
}

void BDSIntegratorEulerOld::AdvanceChord(const G4double       h,
					 G4ThreeVector&       localPos,
					 G4ThreeVector&       localMom,
					 const G4ThreeVector& localA)
{
  // determine effective curvature
  G4double localAMag = localA.mag();
  if(BDS::IsFinite(localAMag))
    {
      // chord distance (simple quadratic approx)
      G4double h2 = h*h;
      G4double dc = h2*localAMag/8;
      SetDistChord(dc);
      
      G4double dx = localMom.x()*h + localA.x()*h2/2;
      G4double dy = localMom.y()*h + localA.y()*h2/2;

      // this can go negative for very long step queries that
      // presumably cause a very large deflection. This results in nan
      // and bad tracking from Geant4 / a crash. Throw exception that'll
      // be caught higher up and the backup stepper will be used.
      G4double dz = std::sqrt(h2*(1.-h2*localAMag*localAMag/12)-dx*dx-dy*dy);
      if (std::isnan(dz))
	{throw std::out_of_range("non-paraxial in old euler method");}
      
      // check for precision problems
      G4double ScaleFac=(dx*dx+dy*dy+dz*dz)/h2;
      if(ScaleFac>1.0000001)
	{
	  ScaleFac=std::sqrt(ScaleFac);
	  dx/=ScaleFac;
	  dy/=ScaleFac;
	  dz/=ScaleFac;
	}
      
      localPos.setX(localPos.x()+dx);
      localPos.setY(localPos.y()+dy);
      localPos.setZ(localPos.z()+dz);
      
      localMom = localMom + h*localA;
    }
  else
    {
      localPos += h*localMom;
      SetDistChord(0);
    }
}
