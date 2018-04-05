/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2018.

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
#include "optionsBase.h"

// include git commit version.
#include "version.h"

#include <iostream>

using namespace GMAD;

OptionsBase::OptionsBase()
{
  // Default Values for Options

  // executable options
  inputFileName         = "optics.mad";
  visMacroFileName      = "";
  visDebug              = false;
  outputFileName        = "output";
  outputFormat          = "rootevent";
  survey                = false;
  surveyFileName        = "survey.dat";
  batch                 = false;
  verbose               = false;
  verboseEvent          = false;
  verboseStep           = false;
  verboseEventNumber    = -1;
  verboseRunLevel       = 0;
  verboseEventLevel     = 0;
  verboseTrackingLevel  = 0;
  verboseSteppingLevel  = 0;
  circular              = false;
  seed                  = -1;
  nGenerate             = 1;
  recreate              = false;
  recreateFileName      = "";
  startFromEvent        = 0;
  writeSeedState        = false;
  useASCIISeedState     = false;
  seedStateFileName     = "";
  generatePrimariesOnly = false;
  exportGeometry        = false;
  exportType            = "gdml";
  exportFileName        = "geometry";
  bdsimPath             = "";
  
  // very important options
  physicsList           = ""; //default - only transportation
  physicsVerbose        = 0;
  
  eventOffset           = 0;
  recreateSeedState     = true;

  elossHistoBinWidth      = 1.0; // m
  
  ffact                 = 1.0;

  beamlineX         = 0;
  beamlineY         = 0;
  beamlineZ         = 0;
  beamlinePhi       = 0;
  beamlineTheta     = 0;
  beamlinePsi       = 0;
  beamlineAxisX     = 0;
  beamlineAxisY     = 0;
  beamlineAxisZ     = 0;
  beamlineAngle     = 0;
  beamlineAxisAngle = false;

  eventNumberOffset       = 0;

  // general geometrical prameters
  checkOverlaps           = false;
  xsize=0.0, ysize=0.0;

  // magnet geometry
  magnetGeometryType   = "polessquare";
  outerMaterialName    = "iron";
  outerDiameter        = 0.6;
  thinElementLength    = 1e-7;
  hStyle               = false;
  vhRatio              = -1; // vhRatio < 0 as signal to use geometry factory default
  coilWidthFraction    = -1;
  coilHeightFraction   = -1;

  preprocessGDML = true;
  
  // geometry debugging
  // always split sbends into smaller chunks by default
  dontSplitSBends      = false;
  includeFringeFields  = true;

  yokeFields           = true;
  sensitiveBeamlineComponents = true;

  // beam pipe / aperture
  beampipeThickness    = 0.0025;
  apertureType         = "circular";
  aper1                = 0.025; // also beampipeRadius
  aper2                = 0.025;
  aper3                = 0.025;
  aper4                = 0.025;
  beampipeMaterial     = "StainlessSteel";
  ignoreLocalAperture  = false;
  
  vacMaterial          = "Vacuum";
  emptyMaterial        = "G4_Galactic";
  worldMaterial        = "air";

  worldVolumeMargin = 5; //m

  vacuumPressure       = 1e-12;
  sensitiveBeamPipe    = true;

  // tunnel options
  buildTunnel         = false;
  buildTunnelStraight = false;
  tunnelType          = "circular";
  tunnelThickness     = 0.1;
  tunnelSoilThickness = 1.0;
  tunnelMaterial      = "concrete";
  soilMaterial        = "soil";
  buildTunnelFloor    = true;
  tunnelFloorOffset   = 1.0; // m
  tunnelAper1         = 2.0; // m
  tunnelAper2         = 2.0; // m
  tunnelSensitive     = true;
  tunnelVisible       = true;
  tunnelOffsetX       = 0;
  tunnelOffsetY       = 0;

  removeTemporaryFiles = true;
  
  // samplers
  samplerDiameter     = 5; // m

  // beam loss monitors geometry
  blmRad                   = 0.05;
  blmLength                = 0.18;
  sensitiveBLMs            = true;

  // physics processes
  turnOnOpticalAbsorption  = true;
  turnOnMieScattering      = true;
  turnOnRayleighScattering = true;
  turnOnOpticalSurface     = true;
  scintYieldFactor         = 1.0;
  maximumPhotonsPerStep    = -1;  ///< -1 -> no action taken (could want 0)
  maximumBetaChangePerStep = 10;
  maximumTracksPerEvent    = 0;   ///< 0 -> no action taken
  minimumKineticEnergy     = 0;
  minimumRange             = 0;
  defaultRangeCut          = 1e-3;
  prodCutPhotons           = 1e-3;
  prodCutElectrons         = 1e-3;
  prodCutPositrons         = 1e-3;
  prodCutProtons           = 1e-3;
  neutronTimeLimit         = 1e-6;
  neutronKineticEnergyLimit = 0;
  useLENDGammaNuclear      = false;
  
  // biasing options
  defaultBiasVacuum        = "";
  defaultBiasMaterial      = "";

  // tracking options
  integratorSet            = "bdsimmatrix";
  lengthSafety             = 1e-9;   // be very careful adjusting this as it affects all the geometry
  maximumTrackingTime      = -1;      // s, nonsensical - used for testing
  maximumStepLength        = 20;      // m, quite big
  maximumTrackLength       = 1e90;    // m, no limit but smaller than DBL_MAX for safe *CLHEP::m
  chordStepMinimum         = 1e-9;    // m // minimum step in a field for an integrator
  chordStepMinimumYoke     = 1e-6;
  deltaIntersection        = 1e-8;    // m - should be greater than lengthSafety!

  // default value in Geant4, old value 0 - error must be greater than this
  minimumEpsilonStep       = 5e-25;
  maximumEpsilonStep       = 1e-7;    // default value in Geant4, old value 1e-7
  deltaOneStep             = 1e-6;    // maximum allowed spatial error in position (1um)
  stopSecondaries          = false;
  killNeutrinos            = true;
  minimumRadiusOfCurvature = 0.05; // 5cm - typical aperture

  // output / analysis options
  numberOfEventsPerNtuple  = 0;
  trajCutGTZ               = 1e99;  // minimum z position, so large default value
  trajCutLTR               = 0.0;   // maximum radius in mm, so small default value
  trajConnect              = false; // connect disconnected trajectory trees
  trajNoTransportation     = false;
  storeElossLinks          = false;
  storeElossLocal          = false;
  storeElossGlobal         = false;
  storeElossTime           = false;
  storeTrajectory          = false;
  storeTrajectoryDepth     = 0;
  storeTrajectoryParticle   = "";
  storeTrajectoryParticleID = "";
  storeTrajectoryEnergyThreshold = -1.0;
  storeSamplerCharge       = false;
  storeSamplerMass         = false;
  storeSamplerRigidity     = false;
  storeSamplerIon          = false;
  writePrimaries           = true;

  // circular options
  nturns                   = 1;

  printFractionEvents = 0.1;
  printFractionTurns  = 0.2;

  // visualisation
  nSegmentsPerCircle       = 50;

  // scoring map
  nbinsx = 1;
  nbinsy = 1;
  nbinsz = 1;
  xmin   = 0;
  xmax   = 0;
  ymin   = 0;
  ymax   = 0;
  zmin   = 0;
  zmax   = 0;
  useScoringMap = false;

  // tracker options
  trackingType             = "hybrid";
  trackingSteps            = 1;
  dontUseAperture          = false;
}


void OptionsBase::print() const
{
  std::cout<<"Options                 " << std::endl;
  std::cout<<"n particles           : " << nGenerate                << std::endl;
  std::cout<<"BV sign               : " << ffact                    << std::endl;
  std::cout<<"Optical absorption on : " << turnOnOpticalAbsorption  << std::endl;
  std::cout<<"Mie scattering on     : " << turnOnMieScattering      << std::endl;
  std::cout<<"Rayleigh scatering on : " << turnOnRayleighScattering << std::endl;
  std::cout<<"Optical surface on    : " << turnOnOpticalSurface     << std::endl;
}

