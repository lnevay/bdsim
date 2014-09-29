#include "TRKThinSymplectic.hh"


#include "TRKBunch.hh"
#include "TRKDipole.hh"
#include "TRKParticle.hh"

TRKThinSymplectic::TRKThinSymplectic(int trackingStepsIn) : 
  TRKThin(trackingStepsIn) {
}

TRKThinSymplectic::~TRKThinSymplectic() { 
}

void TRKThinSymplectic::Track(TRKDrift* el, TRKBunch* bunch) { 

  const double h = el->GetLength()/trackingSteps;

  TRKBunchIter iter = bunch->begin();
  TRKBunchIter end = bunch->end();
  
  for (;iter!=end;++iter) {
    TRKParticle* part = *iter;
    
    for (int i=0; i<trackingSteps; i++) {

      // double x0 = part->X();
      // double y0 = part->Y();
      // double z0 = part->Z();
      double xp = part->Xp();
      double yp = part->Yp();
      // double zp = part->Zp();

      // TODO, E0, P0 beta0 get from BDSParticle and GConstants
      double beta0=0,E0=0,P0=0, momentum0=0,E=0;
      /// from Sixtrack Physics Manual 3.1 Drift
      double delta = (part->mom().mag() - momentum0) / momentum0;

      double pt = (E-E0)/P0/TRK::SOL;

      double pz = std::sqrt((1+delta)*(1+delta)-xp*xp-yp*yp);
      double betaz=pz/(1/beta0+pt);

      vector3 dv(xp/pz*h,yp/pz*h,(1-beta0/betaz)*h);
      part->pluspos(dv);
    }
  }	
}

void TRKThinSymplectic::Track(TRKDipole* el, TRKBunch* bunch) { 
  /// from Sixtrack Physics Manual 3.2.2 Thin Dipole exact Hamiltonian

  // in order to represent a dipole of length L the map is combined with two surrounding drift spaces.

  const double strength = el->GetStrength();
  if (std::abs(strength)<=1e-12) {
    return Track((TRKDrift*)el,bunch);
  }

  const double h = el->GetLength()/trackingSteps;

  for (int i=0; i<trackingSteps; i++) {
      
    TRKBunchIter iter = bunch->begin();
    TRKBunchIter end = bunch->end();
  
    // TODO how to use half? add option?
    // half drift
    Track((TRKDrift*)el,bunch);
    
    for (;iter!=end;++iter) {
      TRKParticle* part = *iter;
	
      double x0 = part->X();
      double y0 = part->Y();
      // double z0 = part->Z();
      double xp = part->Xp();
      // double yp = part->Yp();
      double zp = part->Zp();
      
      // 3 steps , T11(L/2) T1 (L),T11(L/2)

      double v[6];

      // todo 
      double hx;

      //      TRKParticle* v = new TRKParticle(part);
      v[0] = x0 + hx*h/2* x0 * xp / zp;
      v[1] = y0 + hx*h/2* x0 * xp / zp;

    }

    // TODO how to use half? add option?
    // 2nd half drift
    Track((TRKDrift*)el,bunch);
  }
}
