#ifndef BDSFIELDMAGINTERPOLATED_H
#define BDSFIELDMAGINTERPOLATED_H

#include "BDSFieldMag.hh"

#include "globals.hh"
#include "G4Transform3D.hh"

/**
 * @brief Class to provide scaling and a base class pointer for interpolator fields.
 *
 * @author Laurie Nevay
 */

class BDSFieldMagInterpolated: public BDSFieldMag
{
public:
  explicit BDSFieldMagInterpolated(G4Transform3D offset,
				   G4double      scalingIn = 1.0);

  virtual ~BDSFieldMagInterpolated(){;}

  inline G4double Scaling() const {return scaling;}
  inline void     SetScaling(G4double scalingIn) {scaling = scalingIn;}
  
private:
  G4double scaling; ///< Field value scaling.
};

#endif
