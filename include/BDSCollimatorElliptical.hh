#ifndef BDSCOLLIMATORELLIPTICAL_H
#define BDSCOLLIMATORELLIPTICAL_H

#include "BDSCollimatorBase.hh"

#include "globals.hh" // geant4 globals / types

/**
 * @brief A class for an elliptical collimator.
 */

class BDSCollimatorElliptical: public BDSCollimatorBase
{
public:
  BDSCollimatorElliptical(G4String name,
			  G4double length,
			  G4double outerDiameter,
			  G4double xApertureIn        = 0,
			  G4double yApertureIn        = 0,
			  G4double xOutApertureIn     = 0,
			  G4double yOutApertureIn     = 0,
			  G4String collimatorMaterial = "copper",
			  G4String vacuumMaterial     = "vacuum",
			  G4String colour             = "collimator");
  virtual ~BDSCollimatorElliptical(){;};

  virtual void BuildInnerCollimator();

private:
  /// Private default constructor to force the use of the supplied one.
  BDSCollimatorElliptical();

  /// @{ Assignment and copy constructor not implemented nor used
  BDSCollimatorElliptical& operator=(const BDSCollimatorElliptical&) = delete;
  BDSCollimatorElliptical(BDSCollimatorElliptical&) = delete;
  ///@}
};

#endif
