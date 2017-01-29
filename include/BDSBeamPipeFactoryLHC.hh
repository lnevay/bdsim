#ifndef BDSBEAMPIPEFACTORYLHC_H
#define BDSBEAMPIPEFACTORYLHC_H

#include "BDSBeamPipeFactoryBase.hh"
#include "BDSBeamPipe.hh"

/**
 * @brief Factory for simple lhc aperture model beam pipes.
 * 
 * Singleton pattern.
 * 
 * @author Laurie Nevay
 */

class BDSBeamPipeFactoryLHC: public BDSBeamPipeFactoryBase
{
public:
  static BDSBeamPipeFactoryLHC* Instance(); ///< Singleton accessor.
  
  virtual ~BDSBeamPipeFactoryLHC();

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
  BDSBeamPipeFactoryLHC(); ///< Private default constructor - singleton pattern.
  static BDSBeamPipeFactoryLHC* instance; ///< Singleton instance.

  //abstract common build features to one function
  //use member variables unique to this factory to pass them around

  /// only the solids are unique, once we have those, the logical volumes and placement in the
  /// container are the same.  group all this functionality together
  BDSBeamPipe* CommonFinalConstruction(G4String    nameIn,
				       G4Material* vacuumMaterialIn,
				       G4Material* beamPipeMaterialIn,
				       G4double    lengthIn,
				       G4double    widthIn,
				       G4double    heightIn);
  void CreateGeneralAngledSolids(G4String      nameIn,
				 G4double      lengthIn,
				 G4double      aper1In,
				 G4double      aper2In,
				 G4double      aper3In,
				 G4double      beamPipeThicknessIn,
				 G4ThreeVector inputfaceIn,
				 G4ThreeVector outputfaceIn);

  /// function to create the container subtraction solid and avoid passing aper1,2,3 around further
  void CreateContainerSubtractionSolid(G4String& nameIn,
				       G4double& lengthIn,
				       G4double& beamPipeThicknessIn,
				       G4double& aper1In,
				       G4double& aper2In,
				       G4double& aper3In);

};
  
#endif
