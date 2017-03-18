#ifndef BDSFIELDFACTORY_H
#define BDSFIELDFACTORY_H

#include "BDSFieldType.hh"
#include "BDSMagnetType.hh"

#include "globals.hh" // geant4 globals / types

#include <map>
#include <vector>

namespace GMAD
{
  class Field;
}

class BDSField;
class BDSFieldInfo;
class BDSFieldObjects;
class BDSMagnetStrength;

class G4EquationOfMotion;
class G4MagIntegratorStepper;
class G4Mag_EqRhs;

/**
 * @brief Factory that produces fields and their associated objects.
 * 
 * Field objects are created according to a BDSFieldType and the associated
 * and required Geant4 objects to properly implement a field. These are packaged
 * together in one object. This factory does not own any of its products. Construction
 * follows in this order:
 *
 * 1 field
 * 2 equation of motion (based on field object)
 * 3 integrator
 * 4 chord finder
 * 5 field manager
 * 6 package it up
 *
 * This also makes use of BDSParser singleton class to create a series of BDSFieldInfo 
 * field specfications as defined by the parser.
 * 
 * @author Laurie Nevay
 */

class BDSFieldFactory
{
public:
  /// Public accessor method for singleton pattern.
  static BDSFieldFactory* Instance();
  
  ~BDSFieldFactory();

  /// Main interface to field factory.
  BDSFieldObjects* CreateField(const BDSFieldInfo&      info,
			       const BDSMagnetStrength* scalingStrength = nullptr,
			       const G4String           scalingKey      = "none");

  /// Return a BDSFieldInfo instance from the parser definitions. Will
  /// exit if no matching field definition found.
  BDSFieldInfo* GetDefinition(G4String name) const;

private:
  /// Create a purely magnetic field.
  BDSFieldObjects* CreateFieldMag(const BDSFieldInfo&      info,
				  const BDSMagnetStrength* scalingStrength = nullptr,
				  const G4String           scalingKey      = "none");

  /// Create a general EM field.
  BDSFieldObjects* CreateFieldEM(const BDSFieldInfo& info);

  /// Create an electric field.
  BDSFieldObjects* CreateFieldE(const BDSFieldInfo& info);

  /// Create an irregular (special) field.
  BDSFieldObjects* CreateFieldIrregular(const BDSFieldInfo& info);

  /// Create a purely magnetic integrator. As it's purely magnetic, this
  /// requires a G4Mag_EqRhs* equation of motion instance.
  G4MagIntegratorStepper* CreateIntegratorMag(const BDSFieldInfo&      info,
					      G4Mag_EqRhs*             eqOfM,
					      const BDSMagnetStrength* strength);

  /// Create an integrator for a general EM field. As it's a general field,
  /// this takes a G4EquationOfMotion* equation of motion instance.
  G4MagIntegratorStepper* CreateIntegratorEM(const BDSFieldInfo& info,
					     G4EquationOfMotion* eqOfM);

  /// Create an integrator for a general E field. 
  /// Same ones as EM but keep this method for clarity as Geant4 unclear - only based
  /// on their examples. examples/extended/field/field02/src/F02ElectricFieldSetup.cc
  G4MagIntegratorStepper* CreateIntegratorE(const BDSFieldInfo& info,
					    G4EquationOfMotion* eqOfM);

  /// Create only a local field object
  BDSField* CreateFieldMagLocal(const BDSFieldType       type,
				const BDSMagnetStrength* strength,
				const G4double           brho);

  /// Create a pure magnetic field as described by an equation, such as a quadupole or
  /// dipole field.  All associated objects are created and packaged together.
  BDSFieldObjects* CreateFieldEquation(const BDSFieldType       type,
				       const BDSMagnetStrength* strength,
				       const G4double           brho);

  /// Create a pure magnetic field as described by an equation, such as a quadupole or
  /// dipole field.  All associated objects are created and packaged together.
  BDSFieldObjects* CreateFieldMagEquation(const BDSMagnetType      type,
					  const BDSMagnetStrength* strength,
					  const G4double           brho);

  /// Create a special teleporter 'field' that shifts particles at the end of rings to
  /// match up correctly.
  BDSFieldObjects* CreateTeleporter(const BDSFieldInfo& info);
  
  /// Private default constructor as singleton class.
  BDSFieldFactory();

  /// Instance - singleton pattern.
  static BDSFieldFactory* instance;
  
  /// Splits the G4String member variable formatAndName on the ":" character.
  /// Whatever is before is taken as the fromat string and whatever is after is
  /// taken as the filepath.
  void ParseFormatAndFilename();

  /// Reset all pointers to nullptr that are temporarily used during construction
  /// to avoid mistaken contamination between uses of the factory
  void CleanUp();

  /// Prepare all required definitions that can be used dynamically.
  void PrepareFieldDefinitions(const std::vector<GMAD::Field>& definitions,
			       const G4double defaultBRho);

  /// BDSFieldInfo definitions prepare from parser vector of definitions.
  std::map<G4String, BDSFieldInfo*> parserDefinitions;
};
#endif
