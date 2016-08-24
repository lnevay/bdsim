#ifndef BDSBUNCHCOMPOSITE_H
#define BDSBUNCHCOMPOSITE_H 

#include "globals.hh"

#include "BDSBunchInterface.hh"

namespace GMAD {
  class Options;
}

/**
 * @brief A distribution that allows mixing of three different 
 * distributions in each primary dimension.
 * 
 * @author Stewart Boogert
 */

class BDSBunchComposite: public BDSBunchInterface
{ 
protected:
  BDSBunchInterface* xBunch; 
  BDSBunchInterface* yBunch;
  BDSBunchInterface* zBunch;

public:
  BDSBunchComposite(); 
  virtual ~BDSBunchComposite();
  void SetOptions(const GMAD::Options& opt); 
  void GetNextParticle(G4double& x0, G4double& y0, G4double& z0, 
		       G4double& xp, G4double& yp, G4double& zp,
		       G4double& t , G4double&  E, G4double& weight);
  
};

#endif
