#ifndef OPTIONSBASE_H
#define OPTIONSBASE_H

#include <string>

namespace GMAD
{
  /**
   * @brief Basic options class independent of Geant4.
   * 
   * This is the basic options structure that can be saved to output
   * and should only contain basic built in C++ types. No Geant4 types
   * are to be used. See also the Options class which inherits this one.
   */
  
  class OptionsBase
  {
  public:
    OptionsBase();

    std::string gitVersion; ///< The version of BDSIM as given by the git repository.

    std::string inputFileName;    ///< input filename
    std::string visMacroFileName; ///< visualisation filename
    bool        visDebug;         ///< flag for visualisation debug
  
    ///@{ Parameter for output format
    std::string outputFileName;
    std::string outputFormat;
    ///@}
  
    ///@{ Parameter for survey
    bool        survey;
    std::string surveyFileName;
    ///@}
  
    ///@{ Parameter for gflash shower parameterisation
    bool   gflash;
    double gflashemax;
    double gflashemin;
    ///@}
  
    bool batch; ///< Flag for batch / interactive mode
    
    ///@{ Geant4 verbose levels
    bool verbose;
    bool verboseEvent;
    bool verboseStep;
    int  verboseEventNumber;
    
    int  verboseRunLevel;
    int  verboseEventLevel;
    int  verboseTrackingLevel;
    int  verboseSteppingLevel;
    ///@}
  
    bool circular;                 ///< Flag for circular machine
    int  seed;                     ///< The seed value for the random number generator
    int  nGenerate;                ///< The number of primary events to simulate
    bool recreate;                 ///< Whether to recreate from a file or not.
    std::string recreateFileName;  ///< The file path to recreate a run from.
    int  startFromEvent;           ///< Event to start from when recreating.
    bool writeSeedState;           ///< Write the seed state each event to a text file.
    bool useASCIISeedState;        ///< Whether to use the seed state from an ASCII file.
    std::string seedStateFileName; ///< Seed state file path.
    
    /// Whether to only generate primary coordinates and quit, or not.
    bool generatePrimariesOnly; 
    
    ///@{ Parameter for controlling geometry export
    bool        exportGeometry;
    std::string exportType;
    std::string exportFileName;
    ///@}

    /// String that points to path where files are searched;
    /// based on environment variable or else input filename
    std::string bdsimPath;
    
    /// list of physics processes
    std::string physicsList;

    ///@{ beam parameters
    std::string particleName;
    std::string distribType;
    std::string xDistribType;
    std::string yDistribType;
    std::string zDistribType;
    std::string distribFile;
    std::string distribFileFormat;
    ///@}
    
    int nlinesIgnore; ///< ignore first lines in the input bunch file
    int eventOffset;  ///< Event number to start from when recreating from a root file.
    bool recreateSeedState; ///< Load seed state when recreating events.
    
    double elossHistoBinWidth;
    double elossHistoTransBinWidth;
    
    /// magnetic field flip (+1 default, -1: flip sign)
    double ffact;
    double beamEnergy;

    ///@{ initial beam centroid
    double X0, Y0, Z0, S0;
    double Xp0, Yp0, Zp0;
    double T0;
    double E0;
    ///@}
    
    /// bunch length
    double sigmaT;

    ///@{ initial twiss parameters
    double betx, bety, alfx, alfy, emitx, emity, dispx, dispy, dispxp, dispyp;
    ///@}
    
    ///@{ for the gaussian beam distribution
    double sigmaX, sigmaXp, sigmaY, sigmaYp;
    ///@}
    
    ///@{ for the circle/square beam distribution
    double envelopeX, envelopeXp, envelopeY, envelopeYp, envelopeT, envelopeE;
    double envelopeR, envelopeRp;
    ///@}
    
    ///@{ for the gaussian sigma matrix distribution
    double sigma11, sigma12, sigma13, sigma14, sigma15, sigma16;
    double          sigma22, sigma23, sigma24, sigma25, sigma26;
    double                   sigma33, sigma34, sigma35, sigma36;
    double                            sigma44, sigma45, sigma46;
    double                                     sigma55, sigma56;
    double                                              sigma66;
    ///@}
    
    ///@{ for the elliptic shell distribution
    double shellX, shellXp, shellY, shellYp;
    double shellXWidth, shellXpWidth, shellYWidth, shellYpWidth;
    ///@}
    
    ///@{ for the ring beam distribution
    double Rmin, Rmax;
    ///@}
    
    ///@{ for the halo distribution
    double      haloEmitX;
    double      haloEmitY;
    double      haloEnvelopeEmitX;
    double      haloEnvelopeEmitY;
    double      haloEnvelopeCollMinX;
    double      haloEnvelopeCollMaxX;
    double      haloEnvelopeCollMinXp;
    double      haloEnvelopeCollMaxXp;
    double      haloEnvelopeCollMinY;
    double      haloEnvelopeCollMaxY;
    double      haloEnvelopeCollMinYp;
    double      haloEnvelopeCollMaxYp;
    double      haloPSWeightParameter;
    std::string haloPSWeightFunction;
    ///@}
    
    /// for the gaussian, elliptic shell, ring distributions
    double sigmaE;

    int    eventNumberOffset;
    
    /// bdsim options
    int       checkOverlaps;
    /// for element specification
    double xsize, ysize;

    /// default magnet geometry parameters
    std::string magnetGeometryType;
    std::string outerMaterialName;
    double      outerDiameter;

    /// geometry debug, don't split bends into multiple segments
    bool      dontSplitSBends;
    
    bool      includeIronMagFields;
    bool      sensitiveBeamlineComponents;

    ///@{ default beampipe parameters
    double      beampipeThickness;
    std::string apertureType;
    double      aper1;
    double      aper2;
    double      aper3;
    double      aper4;
    std::string beampipeMaterial;
    ///@}
    
    /// vacuum material
    std::string vacMaterial;
    /// world volume
    std::string emptyMaterial;

    double    vacuumPressure;
    bool      sensitiveBeamPipe;
    
    ///@{ tunnel geometry parameters
    bool        buildTunnel;
    bool        buildTunnelStraight;
    std::string tunnelType;
    double      tunnelThickness;
    double      tunnelSoilThickness;
    std::string tunnelMaterial;
    std::string soilMaterial;
    bool        buildTunnelFloor;
    double      tunnelFloorOffset;
    double      tunnelAper1;
    double      tunnelAper2;
    bool        tunnelSensitive;
    bool        tunnelVisible;
    double      tunnelOffsetX;
    double      tunnelOffsetY;
    ///@}
    
    /// Sampler
    double   samplerDiameter;

    ///BLM geometry
    double   blmRad;
    double   blmLength;
    bool     sensitiveBLMs;

    bool     useEMLPB;
    bool     useHadLPB;
    bool     doPlanckScattering;    
    double   planckScatterFe;
    bool     turnOnCerenkov;
    bool     turnOnOpticalAbsorption;
    bool     turnOnMieScattering;
    bool     turnOnRayleighScattering;
    bool     turnOnOpticalSurface;
    bool     turnOnBirksSaturation;

    ///Cross section biasing parameters
    double   scintYieldFactor;
    double   LPBFraction;
    double   thresholdCutCharged;
    double   thresholdCutPhotons;
    double   defaultRangeCut;
    double   prodCutPhotons;
    double   prodCutPhotonsP;
    double   prodCutPhotonsA;
    double   prodCutElectrons;
    double   prodCutElectronsP;
    double   prodCutElectronsA;
    double   prodCutPositrons;
    double   prodCutPositronsP;
    double   prodCutPositronsA;
    double   prodCutProtons;
    double   prodCutProtonsP;
    double   prodCutProtonsA;
    
    /// Biasing options
    std::string defaultBiasVacuum;
    std::string defaultBiasMaterial;

    /// Tracking related parameters
    double   lengthSafety;
    double   maximumTrackingTime; ///< maximum tracking time per volume [s]
    double   deltaChord;
    double   chordStepMinimum;
    double   deltaIntersection;
    double   minimumEpsilonStep;
    double   maximumEpsilonStep;
    double   deltaOneStep;
    bool     stopTracks;    
    bool     stopSecondaries;
    bool     killNeutrinos;

    int         numberOfEventsPerNtuple;
    double      trajCutGTZ;
    double      trajCutLTR;
    bool        storeTrajectory;
    int         storeTrajectoryDepth;
    std::string storeTrajectoryParticle;
    double      storeTrajectoryEnergyThreshold;

    /// Ring parameters
    int      nturns;

    double   printModuloFraction;

    /// Visualisation
    int nSegmentsPerCircle; ///< Number of facets per 2pi in visualisation

    std::string trackingType;
    int         trackingSteps;
    bool        dontUseAperture;

    /// print some properties
    void print() const;
  };
}


#endif //__OPTIONSBASE_H
