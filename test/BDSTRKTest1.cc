/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2021.

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
#include "parser/beam.h"
#include "parser/options.h"

#include "tracker/TRKBunch.hh"
#include "tracker/TRKLine.hh"
#include "tracker/TRKFactory.hh"
#include "tracker/TRKStrategy.hh"
#include "tracker/TRKTracker.hh"

#include "BDSExecOptions.hh"
#include "BDSGlobalConstants.hh"
#include "BDSOutput.hh"
#include "BDSOutputFactory.hh"
#include "BDSParser.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"
#include "BDSRandom.hh" // for random number generator from CLHEP

int main (int argc, char** argv)
{
  /// Initialize executable command line options reader object
  const BDSExecOptions* execOptions = new BDSExecOptions(argc,argv);
  execOptions->Print();

  /// Parse lattice file
  BDSParser* parser = BDSParser::Instance(execOptions->InputFileName());
  /// Update options generated by parser with those from executable options.
  parser->AmalgamateOptions(execOptions->Options());
  parser->AmalgamateBeam(execOptions->Beam(), execOptions->Options().recreate);
  /// Check options for consistency
  parser->CheckOptions();

  delete execOptions;

  // Force construction of global constants after parser has been initialised
  /// (requires materials first). This uses the options from BDSParser.
  const BDSGlobalConstants* globalConstants = BDSGlobalConstants::Instance();
  
  /// Local shortcut to options
  const GMAD::Options& options = BDSParser::Instance()->GetOptions();
  const GMAD::Beam&    beam    = BDSParser::Instance()->GetBeam();
  
  /// Initialise output
  BDSOutput* output = BDSOutputFactory::CreateOutput(globalConstants->OutputFormat(),
						     globalConstants->OutputFileName());

  /// Initialise random numbers
  BDSRandom::CreateRandomNumberGenerator();
  BDSRandom::SetSeed(); // set the seed from options or from exec options

  /// Build bunch
  BDSParticleDefinition* designParticle = nullptr;
  BDSParticleDefinition* beamParticle   = nullptr;
  bool beamDifferentFromDesignParticle = false;
  long int nGenerate = globalConstants->NGenerate();
  BDS::ConstructDesignAndBeamParticle(BDSParser::Instance()->GetBeam(),
				      globalConstants->FFact(),
				      designParticle,
				      beamParticle,
				      beamDifferentFromDesignParticle);

  TRKBunch* bunch = new TRKBunch(beam, beamParticle, nGenerate);
  /// Write primaries to output file
  output->WriteTrackerBunch("primaries",bunch,true);
  
  /// Build beamline
  TRKFactory* factory   = new TRKFactory(options, designParticle, output);
  TRKLine* beamline     = factory->CreateLine(BDSParser::Instance()->GetBeamline());
  TRKStrategy* strategy = factory->CreateStrategy();

  /// Build tracker
  TRKTracker tracker(beamline,strategy,options,output);

  //run tracking - all output through bdsim / samplers
  tracker.Track(bunch);
  
  // free memory (good code test)
  delete designParticle;
  delete beamParticle;
  delete factory;
  delete bunch;
  delete output;

  //done
  return 0;
}
