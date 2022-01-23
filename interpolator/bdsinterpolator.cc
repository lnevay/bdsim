/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2022.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "BDSException.hh"
#include "BDSExecOptions.hh"
#include "BDSFieldFactory.hh"
#include "BDSFieldInfo.hh"
#include "BDSFieldObjects.hh"
#include "BDSFieldQueryInfo.hh"
#include "BDSFieldQueryRaw.hh"
#include "BDSDetectorConstruction.hh"
#include "BDSGlobalConstants.hh"
#include "BDSParser.hh"
#include "BDSPhysicsUtilities.hh"
#include "BDSWarning.hh"

#include "globals.hh"      // geant4 types / globals
#include "G4Field.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"

#include "parser/beam.h"

#include <stdexcept>
#include <string>
#include <vector>


int main(int argc, char** argv)
{
  /// Print header & program information
  G4cout<<"bdsinterpolator : version @BDSIM_VERSION@"<<G4endl;
  G4cout<<"                  (C) 2001-@CURRENT_YEAR@ Royal Holloway University London"<<G4endl;
  G4cout<<"                  http://www.pp.rhul.ac.uk/bdsim"<<G4endl;
  G4cout<<G4endl;

   /// Initialize executable command line options reader object
  const BDSExecOptions* execOptions = new BDSExecOptions(argc,argv);

  /// Parse lattice file
  auto fileName = execOptions->InputFileName();
  G4cout << __FUNCTION__ << "> Using input file : "<< fileName << G4endl;
  BDSParser::Instance(fileName);

  /// Update options generated by parser with those from executable options.
  BDSParser::Instance()->AmalgamateOptions(execOptions->Options());
  delete execOptions; /// No longer needed.

  /// Force construction of global constants
  BDSGlobalConstants::Instance();

  BDSParticleDefinition* designParticle = nullptr;
  BDSParticleDefinition* beamParticle = nullptr;
  G4bool beamDifferentFromDesignParticle = false;
  auto beamDefinition = BDSParser::Instance()->GetBeam();
  if (!beamDefinition.particle.empty())
    {
      BDS::ConstructDesignAndBeamParticle(beamDefinition,
					  BDSGlobalConstants::Instance()->FFact(),
					  designParticle,
					  beamParticle,
					  beamDifferentFromDesignParticle);
      BDSFieldFactory::SetDesignParticle(designParticle);
    }

  BDSFieldQueryRaw querier;

  try
  {
    std::vector<BDSFieldQueryInfo*> queries = BDSDetectorConstruction::PrepareFieldQueries(nullptr);
    for (const auto& query: queries)
      {
	G4cout << "Query: \"" << query->name << "\"" << G4endl;
        if (query->fieldObject.empty())
	  {
	    G4String msg = "\"fieldObject\" variable is empty in query definition \"" + query->name;
	    msg += "\" - it must have a value\nContinuing to next query...";
	    BDS::Warning(msg);
	    continue;
	  }
	BDSFieldInfo* recipe = BDSFieldFactory::Instance()->GetDefinition(query->fieldObject);
	recipe->SetProvideGlobalTransform(false);
	BDSFieldObjects* completeField = BDSFieldFactory::Instance()->CreateField(*recipe);
	G4Field* field = nullptr;
	if (completeField)
	  {field = completeField->GetField();}
	if (!field)
	  {G4cout << "No field constructed - skipping" << G4endl; continue;}
	querier.QueryFieldRaw(field, query);
      }
  }
  catch (BDSException& e)
    {G4cout << e.what();}
  catch (std::exception& e)
    {G4cout << e.what();}
  
  delete BDSGlobalConstants::Instance();
  delete BDSParser::Instance();

  return 0;
}
