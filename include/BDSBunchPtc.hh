#ifndef BDSBUNCHPTC_H
#define BDSBUNCHPTC_H

#include "BDSBunch.hh"

#include "globals.hh"

#include <vector>

/**
 * @brief A bunch distribution that reads a PTC inrays file.
 * 
 * @author Stewart Boogert
 */

class BDSBunchPtc: public BDSBunch
{  
public: 
  BDSBunchPtc();
  virtual ~BDSBunchPtc(); 
  virtual void SetOptions(const GMAD::Options& opt,
			  G4Transform3D beamlineTransformIn = G4Transform3D::Identity);
  virtual void GetNextParticle(G4double& x0, G4double& y0, G4double& z0, 
                               G4double& xp, G4double& yp, G4double& zp,
                               G4double& t , G4double&  E, G4double& weight);
  
private:
  /// Load the PTC file into memory
  void LoadPtcFile();

  /// Assign the distribution file by finding the full path of it.
  void SetDistribFile(G4String distribFileNameIn);
  
  G4int    nRays;     ///< Number of rays in file (1 counting).
  G4String fileName;  ///< File name.
  G4int    iRay;      ///< Iterator counter for current ray.
  std::vector<double*> ptcData; ///< Data.

  G4bool loopedOver;  ///< Wether we've reset to loop over the file again.
};

#endif
