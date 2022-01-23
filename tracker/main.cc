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
#include "parser/beam.h"
#include "parser/options.h"

#include "tracker/EventOutput.hh"
#include "tracker/TRKBunch.hh"
#include "tracker/TRKLine.hh"
#include "tracker/TRKFactory.hh"
#include "tracker/TRKStrategy.hh"
#include "tracker/FixedStepTracker.hh"
#include "tracker/VariableStepTracker.hh"
#include "tracker/TrackerUtils.hh"

#include "BDSExecOptions.hh"
#include "BDSGlobalConstants.hh"
#include "BDSOutput.hh"
#include "BDSParser.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"
#include "BDSRandom.hh" // for random number generator from CLHEP


#include "Options.hh"
#include "OpticsAccumulator.hh"
#include <stdexcept>

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

  auto eventOutput =
      std::make_shared<trk::EventOutput>(globalConstants->OutputFileName());

  std::shared_ptr<trk::OpticsAccumulator> opticsAccumulator;
  if (options.trackerOptics)
  {
      opticsAccumulator = std::make_shared<trk::OpticsAccumulator>();
  }

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

//  if (beamDifferentFromDesignParticle) {
//    throw std::runtime_error(
//        "Tracker beam particle must match design particle.");
//  }

  TRKBunch* bunch = new TRKBunch(beam, designParticle, beamParticle, nGenerate);

  /// Build beamline
  TRKFactory *factory =
      new TRKFactory(options, designParticle, eventOutput, opticsAccumulator);
  TRKLine* beamline     = factory->CreateLine(BDSParser::Instance()->GetBeamline());
  TRKStrategy* strategy = factory->CreateStrategy();

  /// Get the mass of the particle
  std::cout << "Mass: " << bunch->mass << std::endl;

  eventOutput->NewFile();
  eventOutput->WritePrimaries(bunch);

  /// Build tracker
  auto tracker = trk::CreateTracker(beamline, strategy, options, eventOutput);

  //run tracking - all output through bdsim / samplers
  tracker->Track(bunch);
  eventOutput->WriteEvents();
  eventOutput->CloseFile();

  if(opticsAccumulator)
  {
      opticsAccumulator->Write(globalConstants->OutputFileName());
  }

  // // free memory (good code test)
  delete designParticle;
  delete beamParticle;
  delete factory;
  delete bunch;

  //done
  return 0;
}
