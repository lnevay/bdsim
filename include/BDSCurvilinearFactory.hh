#ifndef BDSCURVILINEARFACTORY_H
#define BDSCURVILINEARFACTORY_H

#include "globals.hh"
#include "G4ThreeVector.hh"

class BDSSimpleComponent;
class BDSTiltOffset;

class G4VSolid;

/**
 * @brief Factory to create curvilinear geometry for parallel world.
 * 
 * Creates cylinders for straight or angled sections that can be used
 * in a parallel world (only!) for curvilinear coordinates. These have
 * logical volumes with nullptr for material, so can therefore only be
 * used in a parallel world.
 * 
 * The length supplied should be the full length desired. The length of the
 * object returned will be this length - 2x length safety.
 * 
 * @author Laurie Nevay
 */

class BDSCurvilinearFactory
{
public:
  BDSCurvilinearFactory();
  ~BDSCurvilinearFactory();

  /// Build a straight curvilinear volume.
  BDSSimpleComponent* CreateCurvilinearVolume(const G4String name,
					      const G4double chordLength,
					      const G4double radius);

  /// Build a straight section of curvilinear geometry, but with angled
  /// faces - typically used for a bend section.  Note this is straight though.
  /// The normal vectors are w.r.t. the chord of the volume to be built.
  BDSSimpleComponent* CreateCurvilinearVolume(const G4String       name,
					      const G4double       arcLength,
					      const G4double       chordLength,
					      const G4double       radius,
					      const G4double       angle,
					      const G4ThreeVector  inputFaceNormal,
					      const G4ThreeVector  outputFaceNormal,
					      const BDSTiltOffset* tiltOffset = nullptr);

  /// Convenience method to use the angled method above but with
  /// the angle evently spit between each face.
  BDSSimpleComponent* CreateCurvilinearVolume(const G4String       name,
					      const G4double       arcLength,
					      const G4double       chordLength,
					      const G4double       radius,
					      const G4double       angle,
					      const BDSTiltOffset* tiltOffset = nullptr);

private:
  /// Common construction.  G4VSolid* can't be const as G4LogicalVolume won't take a
  /// const solid.
  BDSSimpleComponent* CommonConstruction(const G4String      name,
					 const G4double      arcLength,
					 const G4double      chordLength,
					 const G4double      radius,
					 G4VSolid*           solid,
					 const G4double      angle);
  
  /// Cache of length safety from BDSGlobalConstants.
  const G4double lengthSafety;
};

#endif
