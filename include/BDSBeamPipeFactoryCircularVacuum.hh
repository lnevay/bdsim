#ifndef BDSBEAMPIPEFACTORYCIRCULARVACUUM_H
#define BDSBEAMPIPEFACTORYCIRCULARVACUUM_H

#include "BDSBeamPipeFactoryBase.hh"
#include "BDSBeamPipe.hh"

/**
 * @brief Factory for vacuum-only circular volumes.
 *
 * No metal beam pipe is made surrounding the vacuum, there is only vacuum.
 * Singleton pattern.
 * 
 * @author Laurie Nevay
 */

class BDSBeamPipeFactoryCircularVacuum: public BDSBeamPipeFactoryBase
{
public:
  static BDSBeamPipeFactoryCircularVacuum* Instance(); ///< Singleton accessor.
  
  virtual ~BDSBeamPipeFactoryCircularVacuum();

  virtual BDSBeamPipe* CreateBeamPipe(G4String    nameIn,
				      G4double    lengthIn,
				      G4double    aper1               = 0,
				      G4double    aper2               = 0,
				      G4double    aper3               = 0,
				      G4double    aper4               = 0,
				      G4Material* vacuumMaterialIn    = nullptr,
				      G4double    beamPipeThicknessIn = 0,
				      G4Material* beamPipeMaterialIn  = nullptr);
  
  virtual BDSBeamPipe* CreateBeamPipe(G4String      nameIn,
				      G4double      lengthIn,
				      G4ThreeVector inputFaceNormalIn,
				      G4ThreeVector outputFaceNormalIn,
				      G4double      aper1               = 0,
				      G4double      aper2               = 0,
				      G4double      aper3               = 0,
				      G4double      aper4               = 0,
				      G4Material*   vacuumMaterialIn    = nullptr,
				      G4double      beamPipeThicknessIn = 0,
				      G4Material*   beamPipeMaterialIn  = nullptr);

private:
  /// Private default constructor - singleton pattern.
  BDSBeamPipeFactoryCircularVacuum();

  /// Singleton instance.
  static BDSBeamPipeFactoryCircularVacuum* instance;

  //abstract common build features to one function
  //use member variables unique to this factory to pass them around

  BDSBeamPipe* CommonFinalConstruction(G4String    nameIn,
				       G4Material* vacuumMaterialIn,
				       G4double    lengthIn,
				       G4double    containerRadiusIn);
};
  
#endif
