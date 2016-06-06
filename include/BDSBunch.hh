#ifndef BDSBUNCH_H
#define BDSBUNCH_H

#include <string>
#include "BDSBunchInterface.hh"

#include "globals.hh"
#include "parser/options.h"

/**
 * @brief A delegator class for the BDSIM particle distribution generator.
 * 
 * The user should instantiate this class only and set the options from the parser.
 * This class generates new particle coordinates for each event based on the 
 * GMAD::Options structure passed from the parser. 
 *
 * @author Stewart Boogert
 */

// can become a singleton? - JS

class BDSBunch
{
protected:
  /// Particle distribution
  std::string        distribType;
  /// BDSBunch holds the distribution information and provides the next particle
  BDSBunchInterface *bdsBunch; 
  
public:
  BDSBunch(); 
  ~BDSBunch();
  /// Copy bunch parameters from parser options.
  void SetOptions(const GMAD::Options& opt);

  /// Interface to supply coordinates for a single particle.
  void GetNextParticle(G4double& x0, G4double& y0, G4double& z0, 
		       G4double& xp, G4double& yp, G4double& zp,
		       G4double& t , G4double&  E, G4double& weight); 

  /// Return the distribution name.
  std::string        GetDistributionType()  {return distribType;}

  /// Interface.
  BDSBunchInterface* GetBDSBunchInterface() {return bdsBunch;}
};

#endif
