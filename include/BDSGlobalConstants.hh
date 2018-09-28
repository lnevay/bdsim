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
#ifndef BDSGLOBALCONSTANTS_H
#define BDSGLOBALCONSTANTS_H 

#include "BDSIntegratorSetType.hh"
#include "BDSMagnetGeometryType.hh"
#include "BDSOutputType.hh"
#include "BDSParticleDefinition.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Version.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "parser/beam.h"
#include "parser/options.h"

#include <map>

class G4UserLimits;
class G4VisAttributes;

class BDSBeamPipeInfo;
class BDSTunnelInfo;

/**
 * @brief A class that holds global options and constants.
 *
 * This wraps (but does not inherit) an instance of GMAD::Options and also converts
 * from the std C++ types used in the options to the Geant4 types
 * solely used in BDSIM. Additionally, some more refined information 
 * is provided - for example, simple strings are converted to enum types
 * for outputs and geometry, as well as G4Materials.
 *
 * This is almost entirely constant apart from a very few variables 
 * that by necessity can be updated. Generally, this is not used
 * as a way to dynamically pass around information at a global level
 * but purely as constants.
 * 
 * Singleton pattern. The (private) constructor requires a GMAD::Options
 * instance, but to maintain the singleton pattern this is accessed through
 * BDSParser singleton instance.
 *
 * Note, some options that are used by BDSBunchInterface derived classes 
 * don't have accessors here as they're only used directly there.
 */

class BDSGlobalConstants
{
protected:
  /// Protected constructor based on a set of gmad options.
  explicit BDSGlobalConstants(const GMAD::Options& opt,
			      GMAD::Beam&          beamIn);

private:
  /// Singleton instance
  static BDSGlobalConstants* instance;

  /// Options instance that this is largely based on and extends
  const GMAD::Options& options;

  /// Copy of beam definition. Can't be const as we may update the contents.
  GMAD::Beam& beam;

  ///@{ Unused default constructors
  BDSGlobalConstants() = delete;
  BDSGlobalConstants(const BDSGlobalConstants&) = delete;
  BDSGlobalConstants& operator=(const BDSGlobalConstants&) = delete;
  ///@}
public:
   /// Access method 
  static BDSGlobalConstants* Instance();
  ~BDSGlobalConstants();

  // Options that access GMAD::options instance
  // Executable options
  inline G4String VisMacroFileName()       const {return G4String(options.visMacroFileName);}
  inline G4bool   VisDebug()               const {return G4bool  (options.visDebug);}
  inline G4String OutputFileName()         const {return G4String(options.outputFileName);}
  inline G4bool   OutputFileNameSet()      const {return G4bool  (options.HasBeenSet("outputFileName"));}
  inline BDSOutputType OutputFormat()      const {return outputType;}
  inline G4bool   Survey()                 const {return G4bool  (options.survey);}
  inline G4String SurveyFileName()         const {return G4String(options.surveyFileName);}
  inline G4bool   Batch()                  const {return G4bool  (options.batch);}
  inline G4bool   Verbose()                const {return G4bool  (options.verbose);}
  inline G4bool   VerboseEvent()           const {return G4bool  (options.verboseEvent);}
  inline G4bool   VerboseStep()            const {return G4bool  (options.verboseStep);}
  inline G4int    VerboseEventNumber()     const {return G4int   (options.verboseEventNumber);}
  inline G4int    VerboseRunLevel()        const {return G4int   (options.verboseRunLevel);}
  inline G4int    VerboseEventLevel()      const {return G4int   (options.verboseEventLevel);}
  inline G4int    VerboseTrackingLevel()   const {return G4int   (options.verboseTrackingLevel);}
  inline G4int    VerboseSteppingLevel()   const {return G4int   (options.verboseSteppingLevel);}
  inline G4bool   Circular()               const {return G4bool  (options.circular);}
  inline G4int    Seed()                   const {return G4int   (options.seed);}
  inline G4bool   SeedSet()                const {return G4bool  (options.HasBeenSet("seed"));}
  inline G4bool   Recreate()               const {return G4bool  (options.recreate);}
  inline G4String RecreateFileName()       const {return G4String(options.recreateFileName);}
  inline G4int    StartFromEvent()         const {return G4int   (options.startFromEvent);}
  inline G4bool   WriteSeedState()         const {return G4bool  (options.writeSeedState);}
  inline G4bool   UseASCIISeedState()      const {return G4bool  (options.useASCIISeedState);}
  inline G4String SeedStateFileName()      const {return G4String(options.seedStateFileName);}
  inline G4String BDSIMPath()              const {return G4String(options.bdsimPath);}
  inline G4int    NGenerate()              const {return numberToGenerate;}
  inline G4bool   GeneratePrimariesOnly()  const {return G4bool  (options.generatePrimariesOnly);}
  inline G4bool   ExportGeometry()         const {return G4bool  (options.exportGeometry);}
  inline G4String ExportType()             const {return G4String(options.exportType);}
  inline G4String ExportFileName()         const {return G4String(options.exportFileName);}  

  // regular options from here on
  G4int PrintModuloEvents()                  const;
  G4int PrintModuloTurns()                   const;
  inline G4bool   PhysicsVerbose()           const {return G4bool  (options.physicsVerbose);}
  inline G4double PrintFractionEvents()      const {return G4double(options.printFractionEvents);}
  inline G4double PrintFractionTurns()       const {return G4double(options.printFractionTurns);}
  inline G4double LengthSafety()             const {return G4double(options.lengthSafety*CLHEP::m);}
  inline G4double HorizontalWidth()          const {return G4double(options.horizontalWidth)*CLHEP::m;}
  inline G4bool   IgnoreLocalAperture()      const {return G4bool  (options.ignoreLocalAperture);}
  inline G4String OuterMaterialName()        const {return G4String(options.outerMaterialName);}
  inline G4bool   DontSplitSBends()          const {return G4bool  (options.dontSplitSBends);}
  inline G4bool   BuildTunnel()              const {return G4bool  (options.buildTunnel);}
  inline G4bool   BuildTunnelStraight()      const {return G4bool  (options.buildTunnelStraight);}
  inline G4double TunnelOffsetX()            const {return G4double(options.tunnelOffsetX)*CLHEP::m;}
  inline G4double TunnelOffsetY()            const {return G4double(options.tunnelOffsetY)*CLHEP::m;}
  inline G4double ElossHistoBinWidth()       const {return G4double(options.elossHistoBinWidth)*CLHEP::m;}
  inline G4double BeamTotalEnergy()          const {return G4double(beam.beamEnergy)*CLHEP::GeV;}
  inline G4double BlmRad()                   const {return G4double(options.blmRad)*CLHEP::m;}
  inline G4double BlmLength()                const {return G4double(options.blmLength)*CLHEP::m;}
  inline G4double DefaultRangeCut()          const {return G4double(options.defaultRangeCut)*CLHEP::m;}
  inline G4double ProdCutPhotons()           const {return G4double(options.prodCutPhotons)*CLHEP::m;}
  inline G4double ProdCutElectrons()         const {return G4double(options.prodCutElectrons)*CLHEP::m;}
  inline G4double ProdCutPositrons()         const {return G4double(options.prodCutPositrons)*CLHEP::m;}
  inline G4double ProdCutProtons()           const {return G4double(options.prodCutProtons)*CLHEP::m;}
  inline G4bool   ProdCutPhotonsSet()        const {return G4bool  (options.HasBeenSet("prodCutPhotons"));}
  inline G4bool   ProdCutElectronsSet()      const {return G4bool  (options.HasBeenSet("prodCutElectrons"));}
  inline G4bool   ProdCutPositronsSet()      const {return G4bool  (options.HasBeenSet("prodCutPositrons"));}
  inline G4bool   ProdCutProtonsSet()        const {return G4bool  (options.HasBeenSet("prodCutProtons"));}
  inline G4double NeutronTimeLimit()         const {return G4double(options.neutronTimeLimit)*CLHEP::s;}
  inline G4double NeutronKineticEnergyLimit()const {return G4double(options.neutronKineticEnergyLimit)*CLHEP::GeV;}
  inline G4bool   UseLENDGammaNuclear()      const {return G4bool  (options.useLENDGammaNuclear);}
  inline G4bool   UseElectroNuclear()        const {return G4bool  (options.useElectroNuclear);}
  inline G4bool   UseMuonNuclear()           const {return G4bool  (options.useMuonNuclear);}
  inline G4bool   UseGammaToMuMu()           const {return G4bool  (options.useGammaToMuMu);}
  inline G4bool   UsePositronToMuMu()        const {return G4bool  (options.usePositronToMuMu);}
  inline G4bool   UsePositronToHadrons()     const {return G4bool  (options.usePositronToHadrons);}

  inline G4double DeltaIntersection()        const {return G4double(options.deltaIntersection)*CLHEP::m;}
  inline G4double ChordStepMinimum()         const {return G4double(options.chordStepMinimum)*CLHEP::m;}
  inline G4double ChordStepMinimumYoke()     const {return G4double(options.chordStepMinimumYoke)*CLHEP::m;}
  inline G4double DeltaOneStep()             const {return G4double(options.deltaOneStep)*CLHEP::m;}
  inline G4double MinimumEpsilonStep()       const {return G4double(options.minimumEpsilonStep);}
  inline G4double MaximumEpsilonStep()       const {return G4double(options.maximumEpsilonStep);}
  inline G4double MaxTime()                  const {return G4double(options.maximumTrackingTime)*CLHEP::s;}
  inline G4double MaxStepLength()            const {return G4double(options.maximumStepLength)*CLHEP::m;}
  inline G4double MaxTrackLength()           const {return G4double(options.maximumTrackLength)*CLHEP::m;}
  inline G4int    TurnsToTake()              const {return G4int   (options.nturns);}
  inline G4double FFact()                    const {return G4double(options.ffact);}
  inline G4double BeamlineS()                const {return G4double(options.beamlineS*CLHEP::m);}
  inline G4double ParticleTotalEnergy()      const {return G4double(beam.E0)*CLHEP::GeV;}
  inline G4bool   SensitiveBeamPipe()        const {return G4bool  (options.sensitiveBeamPipe);}
  inline G4bool   SensitiveOuter()           const {return G4bool  (options.sensitiveOuter);}
  inline G4bool   SensitiveVacuum()          const {return G4bool  (options.sensitiveVacuum);}
  inline G4bool   SensitiveBLMs()            const {return G4bool  (options.sensitiveBLMs);}
#if G4VERSION_NUMBER != 1030
  inline G4bool   CheckOverlaps()            const {return G4bool  (options.checkOverlaps);}
#else
  // Overlap checking slow in 10.3.0 therefore disabled,
  // see https://bitbucket.org/jairhul/bdsim/issues/151/overlap-checking-in-103-gives-warnings-and
  inline G4bool   CheckOverlaps()            const {return false;}
#endif
  inline G4int    EventNumberOffset()        const {return G4int   (options.eventNumberOffset);}
  inline G4bool   WritePrimaries()           const {return G4bool  (options.writePrimaries);}
  inline G4bool   StoreELossTurn()           const {return G4bool  (options.storeElossTurn || options.circular);}
  inline G4bool   StoreELossLinks()          const {return G4bool  (options.storeElossLinks);}
  inline G4bool   StoreELossLocal()          const {return G4bool  (options.storeElossLocal);}
  inline G4bool   StoreELossGlobal()         const {return G4bool  (options.storeElossGlobal);}
  inline G4bool   StoreElossTime()           const {return G4bool  (options.storeElossTime);}
  inline G4bool   StoreElossStepLength()     const {return G4bool  (options.storeElossStepLength);}
  inline G4bool   StoreElossPreStepKineticEnergy() const {return G4bool  (options.storeElossPreStepKineticEnergy);}
  inline G4bool   StoreTrajectory()          const {return G4bool  (options.storeTrajectory);}
  inline G4int    StoreTrajectoryDepth()     const {return G4int   (options.storeTrajectoryDepth);}
  inline G4String StoreTrajectoryParticle()  const {return G4String(options.storeTrajectoryParticle);}
  inline G4String StoreTrajectoryParticleID()const {return G4String(options.storeTrajectoryParticleID);}
  inline G4double StoreTrajectoryEnergyThreshold() const {return G4double (options.storeTrajectoryEnergyThreshold*CLHEP::GeV);}
  std::vector<int>                           StoreTrajectorySamplerIDs();
  std::vector<std::pair<double,double>>      StoreTrajectoryELossSRange();
  inline G4bool   StoreSamplerCharge()       const {return G4bool  (options.storeSamplerCharge);}
  inline G4bool   StoreSamplerKineticEnergy() const {return G4bool  (options.storeSamplerKineticEnergy);}
  inline G4bool   StoreSamplerMass()         const {return G4bool  (options.storeSamplerMass);}
  inline G4bool   StoreSamplerRigidity()     const {return G4bool  (options.storeSamplerRigidity);}
  inline G4bool   StoreSamplerIon()          const {return G4bool  (options.storeSamplerIon);}
  inline G4bool   StoreModel()               const {return G4bool  (options.storeModel);}
  inline G4bool   TrajConnect()              const {return G4bool  (options.trajConnect);}
  inline G4double TrajCutGTZ()               const {return G4double(options.trajCutGTZ*CLHEP::m);}
  inline G4double TrajCutLTR()               const {return G4double(options.trajCutLTR*CLHEP::m);}
  inline G4bool   TrajNoTransportation()     const {return G4bool  (options.trajNoTransportation);}
  inline G4bool   StopSecondaries()          const {return G4bool  (options.stopSecondaries);}
  inline G4bool   KillNeutrinos()            const {return G4bool  (options.killNeutrinos);}
  inline G4double MinimumRadiusOfCurvature() const {return G4double(options.minimumRadiusOfCurvature*CLHEP::m);}
  inline G4double ScintYieldFactor()         const {return G4double(options.scintYieldFactor);}
  inline G4int    MaximumPhotonsPerStep()    const {return G4int   (options.maximumPhotonsPerStep);}
  inline G4int    MaximumBetaChangePerStep() const {return G4int   (options.maximumBetaChangePerStep);}
  inline G4long   MaximumTracksPerEvent()    const {return G4long  (options.maximumTracksPerEvent);}
  inline G4double MinimumKineticEnergy()     const {return G4double(options.minimumKineticEnergy*CLHEP::GeV);}
  inline G4double MinimumRange()             const {return G4double(options.minimumRange*CLHEP::m);}
  inline G4String VacuumMaterial()           const {return G4String(options.vacMaterial);}
  inline G4String EmptyMaterial()            const {return G4String(options.emptyMaterial);}
  inline G4String WorldMaterial()            const {return G4String(options.worldMaterial);}
  inline G4double WorldVolumeMargin()        const {return G4double(options.worldVolumeMargin*CLHEP::m);}
  inline G4bool   YokeFields()               const {return G4bool  (options.yokeFields);}
  inline G4bool   TurnOnOpticalAbsorption()  const {return G4bool  (options.turnOnOpticalAbsorption);}
  inline G4bool   TurnOnRayleighScattering() const {return G4bool  (options.turnOnRayleighScattering);}
  inline G4bool   TurnOnMieScattering()      const {return G4bool  (options.turnOnMieScattering);}
  inline G4bool   TurnOnOpticalSurface()     const {return G4bool  (options.turnOnOpticalSurface);}
  inline G4int    NumberOfEventsPerNtuple()  const {return G4int   (options.numberOfEventsPerNtuple);}
  inline G4bool   IncludeFringeFields()      const {return G4bool  (options.includeFringeFields);}
  inline G4int    NSegmentsPerCircle()       const {return G4int   (options.nSegmentsPerCircle);}
  inline G4double ThinElementLength()        const {return G4double(options.thinElementLength*CLHEP::m);}
  inline G4bool   HStyle()                   const {return G4bool  (options.hStyle);}
  inline G4double VHRatio()                  const {return G4double(options.vhRatio);}
  inline G4double CoilWidthFraction()        const {return G4double(options.coilWidthFraction);}
  inline G4double CoilHeightFraction()       const {return G4double(options.coilHeightFraction);}
  inline G4bool   PreprocessGDML()           const {return G4bool  (options.preprocessGDML);}
  inline G4int    NBinsX()                   const {return G4int   (options.nbinsx);}
  inline G4int    NBinsY()                   const {return G4int   (options.nbinsy);}
  inline G4int    NBinsZ()                   const {return G4int   (options.nbinsz);}
  inline G4double XMin()                     const {return G4double(options.xmin) * CLHEP::m;}
  inline G4double YMin()                     const {return G4double(options.ymin) * CLHEP::m;}
  inline G4double ZMin()                     const {return G4double(options.zmin) * CLHEP::m;}
  inline G4double XMax()                     const {return G4double(options.xmax) * CLHEP::m;}
  inline G4double YMax()                     const {return G4double(options.ymax) * CLHEP::m;}
  inline G4double ZMax()                     const {return G4double(options.zmax) * CLHEP::m;}
  inline G4bool   UseScoringMap()            const {return G4bool  (options.useScoringMap);}
  inline G4bool   MatchDistrFileLength()     const {return G4bool  (beam.matchDistrFileLength);}
  inline G4bool   RemoveTemporaryFiles()     const {return G4bool  (options.removeTemporaryFiles);}
  inline G4bool   SampleElementsWithPoleface() const {return G4bool  (options.sampleElementsWithPoleface);}
  inline G4double NominalMatrixRelativeMomCut() const {return G4double (options.nominalMatrixRelativeMomCut);}
  inline G4bool   TeleporterFullTransform()  const {return G4bool  (options.teleporterFullTransform);}
  
  // options that require members in this class (for value checking or because they're from another class)
  inline G4int                 TurnsTaken()              const {return turnsTaken;}
  inline G4String              ParticleName()            const {return beam.particleName;}
  inline G4double              SamplerDiameter()         const {return samplerDiameter;}
  inline BDSParticleDefinition* BeamParticleDefinition() const {return beamParticleDefinition;}
  inline BDSBeamPipeInfo*      DefaultBeamPipeModel()    const {return defaultBeamPipeModel;}
  inline BDSMagnetGeometryType MagnetGeometryType()      const {return magnetGeometryType;}
  inline BDSTunnelInfo*        TunnelInfo()              const {return tunnelInfo;}
  inline G4VisAttributes*      GetInvisibleVisAttr()     const {return invisibleVisAttr;}
  inline G4VisAttributes*      VisibleDebugVisAttr()     const {return visibleDebugVisAttr;}
  inline G4VisAttributes*      ContainerVisAttr()        const {return options.visDebug ? visibleDebugVisAttr : invisibleVisAttr;}
  inline G4UserLimits*         DefaultUserLimits()       const {return defaultUserLimits;}
  inline BDSIntegratorSetType  IntegratorSet()           const {return integratorSet;}
  inline G4Transform3D         BeamlineTransform()       const {return beamlineTransform;}

  /// @{ Setter
  inline void SetSamplerDiameter(const G4double& samplerDiameterIn) {samplerDiameter = samplerDiameterIn;}
  inline void SetBeamParticleDefinition(BDSParticleDefinition* particleDefinitionIn);
  inline void IncrementTurnNumber()  {turnsTaken += 1;}
  inline void ResetTurnNumber()      {turnsTaken = 1;}
  inline void SetNumberToGenerate(G4int number) {numberToGenerate = number;}
  /// @}

  // laserwire stuff that probably shouldn't be in global constants
  inline G4double      GetLaserwireWavelength()     const {return itsLaserwireWavelength;}
  inline G4ThreeVector GetLaserwireDir()            const {return itsLaserwireDir;}
  inline G4bool        GetLaserwireTrackPhotons()   const {return itsLaserwireTrackPhotons;}
  inline G4bool        GetLaserwireTrackElectrons() const {return itsLaserwireTrackElectrons;}
  inline G4double      GetLaserwireWavelength(G4String aName) const {return lwWavelength.at(aName);}
  inline G4ThreeVector GetLaserwireDir(G4String aName)        const {return lwDirection.at(aName);}
  void SetLaserwireWavelength(G4String aName, G4double aWavelength);
  void SetLaserwireDir(G4String aName, G4ThreeVector aDirection);

private:
  BDSParticleDefinition* beamParticleDefinition; ///< Initial bunch parameters

  /// Particle name
  G4String particleName;
  
  /// Number of particles to generate can be set from outside (by e.g. BDSBunchPtc)
  G4int numberToGenerate;

  /// Cache of sampler diameter in this class so it can be updated.
  G4double samplerDiameter;

  ///@{ Magnet geometry
  BDSMagnetGeometryType magnetGeometryType;
  ///@}

  /// Default beam pipe model information
  BDSBeamPipeInfo* defaultBeamPipeModel;
  
  /// Tunnel model
  BDSTunnelInfo* tunnelInfo;
  
  // test map container for laserwire parameters - Steve
  std::map<const G4String, G4double> lwWavelength;
  std::map<const G4String, G4ThreeVector> lwDirection;
  G4double      itsLaserwireWavelength;
  G4ThreeVector itsLaserwireDir;
  G4bool        itsLaserwireTrackPhotons;
  G4bool        itsLaserwireTrackElectrons;

  /// Prepare the G4Transform3D instance from the options for the initial beam line transform.
  void InitialiseBeamlineTransform();
  
  void InitVisAttributes();
  G4VisAttributes* invisibleVisAttr;
  G4VisAttributes* visibleDebugVisAttr;

  void InitDefaultUserLimits();
  G4UserLimits* defaultUserLimits;
  
  /// Turn Control
  G4int turnsTaken;

  BDSOutputType        outputType;         ///< Output type enum for output format to be used.
  BDSIntegratorSetType integratorSet;      ///< Integrator type enum for integrator set to be used.
  G4Transform3D        beamlineTransform;  ///< Transform for start of beam line.
};

inline void BDSGlobalConstants::SetBeamParticleDefinition(BDSParticleDefinition* particleDefinitionIn)
{beamParticleDefinition = particleDefinitionIn;}

inline void BDSGlobalConstants::SetLaserwireWavelength(G4String aName, G4double aWavelength)
{lwWavelength[aName]=aWavelength;}

inline void BDSGlobalConstants::SetLaserwireDir(G4String aName, G4ThreeVector aDirection)
{lwDirection[aName]=aDirection;}


#endif
