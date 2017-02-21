#ifndef BDSCAVITYRF_H
#define BDSCAVITYRF_H

#include "BDSCavity.hh"

#include "globals.hh" // geant4 globals / types

class BDSFieldInfo;

/** 
 * @brief RF Cavity class
 *  
 * @author Stuart Walker
 */

class BDSCavityRF: public BDSCavity
{
public:
  BDSCavityRF(G4String            name,
	      G4double            length,
	      const BDSFieldInfo* vacuumFieldIn);

  virtual ~BDSCavityRF(){;}
    
  /// Creates field objects - doesn't nothing by default and derived classes can override.
  virtual void BuildField();

private:
  /// Private default constructor to force the use of the supplied one.
  BDSCavityRF();
  
  /// @{ Assignment and copy constructor not implemented nor used
  BDSCavityRF& operator=(const BDSCavityRF&) = delete;
  BDSCavityRF(BDSCavityRF&) = delete;
  /// @}
};
#endif
