#include "TRKDipole.hh"

#include <cmath>
#include <cstdlib>
#include "vector3.hh"
#include "vector6.hh"

TRKDipole::TRKDipole(double strengthIn, double bFieldIn, TRKTrackingElement::TRKType typeIn, int trackingStepsIn, std::string nameIn, double lengthIn, TRKAperture *apertureIn, TRKPlacement *placementIn):
  TRKTrackingElement(typeIn, trackingStepsIn,nameIn,lengthIn,apertureIn,placementIn), strength(strengthIn), bField(bFieldIn), drift(NULL)
{
  drift = new TRKDrift(typeIn, trackingStepsIn,nameIn+"Drift",lengthIn,apertureIn,placementIn);
}

TRKDipole::~TRKDipole() {
  delete drift;
}

void TRKDipole::ThinTrack(const double /*vIn*/[], double /*vOut*/[], double /*h*/) {  

}

void TRKDipole::HybridTrack(const double vIn[], double vOut[], double h) {  
  if (std::abs(strength)<=1e-12) {
    return drift->Track(vIn,vOut,h);
  }

  double x0 = vIn[0];
  double y0 = vIn[1];
  double z0 = vIn[2];
  double xp = vIn[3];
  double yp = vIn[4];
  double zp = vIn[5];

  // adapted from BDSDipoleStepper.cc
  // 1) bfield
  // 2) kicker strength , overwrites first calculation?

  vector3 yhat(0.,1.,0.);

  vector3 vhat(vIn3[3],vIn[4],vIn[5]);
  
  vector3 vnorm = vhat.cross(yhat);
  // TODO: double initMag = vhat.mag(); // not correct, need to be included in bfield in Factory
  double R = DBL_MAX;

  //  if(bField!=0)
    //    R=-(InitMag/CLHEP::GeV)/(0.299792458 * itsBField/CLHEP::tesla) * CLHEP::m;
  //  else
  //    R=DBL_MAX;

 // include the sign of the charge of the particles

  // if(  fPtrMagEqOfMot->FCof()<0) R*=-1.;
  // else if ( fPtrMagEqOfMot->FCof()==0) R=DBL_MAX;
  
  double Theta   = h/R;

  double CosT_ov_2, SinT_ov_2, CosT, SinT;
  CosT_ov_2=cos(Theta/2);
  SinT_ov_2=sin(Theta/2);
  
  CosT=(CosT_ov_2*CosT_ov_2)- (SinT_ov_2*SinT_ov_2);
  SinT=2*CosT_ov_2*SinT_ov_2;
  
  //BDSLocalRadiusOfCurvature=R;
  
  //  itsDist=fabs(R)*(1.-CosT_ov_2);
  
  vector3 dPos=R*(SinT*vhat + (1-CosT)*vnorm);

  // strength
  
  if(fabs(strength)<std::numeric_limits<double>::epsilon()) {  // no gradient
  
    vOut[0] = x0+dPos.X();
    vOut[1] = y0+dPos.Y();
    vOut[2] = z0+dPos.Z();
    
    vOut[3] = cosT*vhat.X() + sinT*vnorm.X();
    vOut[4] = cosT*vhat.Y() + sinT*vnorm.Y();
    vOut[5] = cosT*vhat.Z() + sinT*vnorm.Z();
  }

  G4double x1,x1p,y1,y1p,z1p;
  //G4double z1;
  
  G4double NomEnergy = BDSGlobalConstants::Instance()->GetBeamTotalEnergy();
  G4double NomR = -(NomEnergy/CLHEP::GeV)/(0.299792458 * itsBField/CLHEP::tesla) * CLHEP::m;

  G4double NominalPath = sqrt(NomR*NomR - LocalR.z()*LocalR.z()) - fabs(NomR)*cos(itsAngle/2);
  
  G4double EndNomPath = sqrt(NomR*NomR - itsFinalPoint.z()*itsFinalPoint.z()) - fabs(NomR)*cos(itsAngle/2);

  if(R<0)
    {
      NominalPath*=-1;
      EndNomPath*=-1;
    }

  G4double x0=LocalR.x() - NominalPath;
  G4double y0=LocalR.y();
  G4double z0=LocalR.z();

  G4double theta_in = asin(z0/NomR);
  
  LocalRp.rotateY(-theta_in);

  G4double xp=LocalRp.x();
  G4double yp=LocalRp.y();
  G4double zp=LocalRp.z();

  // Save for Synchrotron Radiation calculations:
  BDSLocalRadiusOfCurvature=R;
  
  G4double rootK=sqrt(fabs(kappa*zp));
  G4double rootKh=rootK*h*zp;
  G4double X11,X12,X21,X22;
  G4double Y11,Y12,Y21,Y22;

  if (kappa>0)
    {
      X11= cos(rootKh);
      X12= sin(rootKh)/rootK;
      X21=-fabs(kappa)*X12;
      X22= X11;
      
      Y11= cosh(rootKh);
      Y12= sinh(rootKh)/rootK;
      Y21= fabs(kappa)*Y12;
      Y22= Y11;
    }
  else // if (kappa<0)
    {
      X11= cosh(rootKh);
      X12= sinh(rootKh)/rootK;
      X21= fabs(kappa)*X12;
      X22= X11;
      
      Y11= cos(rootKh);
      Y12= sin(rootKh)/rootK;
      Y21= -fabs(kappa)*Y12;
      Y22= Y11;
    }
  // else // should not happen as already returned in that case
  //   {
  //     X11 = 1;
  //     X12 = 0;
  //     X21 = 0;
  //     X22 = 1;
  //     Y11 = 1;
  //     Y12 = 0;
  //     Y21 = 0;
  //     Y22 = 1;
  //   }

  x1      = X11*x0 + X12*xp;    
  x1p     = X21*x0 + X22*xp;

  y1      = Y11*y0 + Y12*yp;    
  y1p     = Y21*y0 + Y22*yp;

  z1p=sqrt(1 - x1p*x1p -y1p*y1p);
  /* 
  x1 -=(kappa/ (24*R) ) * h2*h2;
  x1p-=(kappa/ (6*R) ) * h*h2;
  */
  G4double dx=x1-x0;
  G4double dy=y1-y0;
  // Linear chord length
  
  LocalR.setX(dx +itsInitialR.x() + EndNomPath - NominalPath);
  LocalR.setY(dy + itsInitialR.y());
  LocalR.setZ(itsFinalPoint.z());
  

  LocalRp.setX(x1p);
  LocalRp.setY(y1p);
  LocalRp.setZ(z1p);
  LocalRp.rotateY(theta_in);
  
  GlobalPosition=LocalAffine.TransformPoint(LocalR); 
  
  LocalRp.rotateY(-h/R);
  GlobalTangent=LocalAffine.TransformAxis(LocalRp);


  GlobalTangent*=InitMag;
  
  // gab: replace += with =
  yOut[0] = GlobalPosition.x(); 
  yOut[1] = GlobalPosition.y(); 
  yOut[2] = GlobalPosition.z(); 
  
  yOut[3] = GlobalTangent.x();
  yOut[4] = GlobalTangent.y();
  yOut[5] = GlobalTangent.z();




}

void TRKDipole::ThickTrack(const double /*vIn*/[], double /*vOut*/[], double /*h*/) {

}

std::ostream& operator<< (std::ostream &out, const TRKDipole &element) {
  return out << element.strength << " ";
}
