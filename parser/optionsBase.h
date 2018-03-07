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
    
    int eventOffset;  ///< Event number to start from when recreating from a root file.
    bool recreateSeedState; ///< Load seed state when recreating events.
    
    double elossHistoBinWidth;
    
    /// magnetic field flip (+1 default, -1: flip sign)
    double ffact;

    ///@{ Intial beam line transform w.r.t. the world coordinate frame.
    double beamlineX;
    double beamlineY;
    double beamlineZ;
    double beamlinePhi;
    double beamlineTheta;
    double beamlinePsi;
    double beamlineAxisX;
    double beamlineAxisY;
    double beamlineAxisZ;
    double beamlineAngle;
    bool   beamlineAxisAngle;
    ///@}   

    int    eventNumberOffset;
    
    /// bdsim options
    bool       checkOverlaps;
    /// for element specification
    double xsize, ysize;

    /// default magnet geometry parameters
    std::string magnetGeometryType;
    std::string outerMaterialName;
    double      outerDiameter;
    double      thinElementLength;
    bool        hStyle; ///< H Style dipoles (if not, C Style).
    double      vhRatio;
    double      coilWidthFraction;
    double      coilHeightFraction;

    /// geometry control
    bool preprocessGDML;

    /// geometry debug, don't split bends into multiple segments
    bool      dontSplitSBends;
    
    bool      includeIronMagFields;
    bool      sensitiveBeamlineComponents;

    bool        includeFringeFields;

    ///@{ default beampipe parameters
    double      beampipeThickness;
    std::string apertureType;
    double      aper1;
    double      aper2;
    double      aper3;
    double      aper4;
    std::string beampipeMaterial;
    bool        ignoreLocalAperture;
    ///@}
    
    std::string vacMaterial;   ///< vacuum material
    std::string emptyMaterial; ///< material in container volumes
    std::string worldMaterial;

    double    worldVolumeMargin; ///< Padding margin for world volume size.

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

    bool removeTemporaryFiles;
    
    /// Sampler
    double   samplerDiameter;

    ///BLM geometry
    double   blmRad;
    double   blmLength;
    bool     sensitiveBLMs;

    ///@{ Physics processes
    bool     turnOnOpticalAbsorption;
    bool     turnOnMieScattering;
    bool     turnOnRayleighScattering;
    bool     turnOnOpticalSurface;
    ///@}

    ///@{Cross section biasing parameters
    double   scintYieldFactor;
    int      maximumPhotonsPerStep;
    int      maximumBetaChangePerStep;
    long     maximumTracksPerEvent;
    double   minimumKineticEnergy;
    double   minimumRange;
    double   defaultRangeCut;
    double   prodCutPhotons;
    double   prodCutElectrons;
    double   prodCutPositrons;
    double   prodCutProtons;
    ///@}

    /// Biasing options
    std::string defaultBiasVacuum;
    std::string defaultBiasMaterial;

    /// Tracking related parameters
    std::string integratorSet;
    double   lengthSafety;
    double   maximumTrackingTime; ///< Maximum tracking time per track [s].
    double   maximumStepLength;   ///< Maximum permitted step length in any volume.
    double   maximumTrackLength;  ///< Maximum permitted track length [m].
    double   chordStepMinimum;
    double   deltaIntersection;
    double   minimumEpsilonStep;
    double   maximumEpsilonStep;
    double   deltaOneStep;
    bool     stopSecondaries;
    bool     killNeutrinos;
    double   minimumRadiusOfCurvature; ///< Minimum allowed radius of curvature. 

    int         numberOfEventsPerNtuple;

    bool        storeElossLinks;
    bool        storeElossLocal;
    bool        storeElossGlobal;

    bool        storeTrajectory;
    int         storeTrajectoryDepth;
    std::string storeTrajectoryParticle;
    std::string storeTrajectoryParticleID;
    double      storeTrajectoryEnergyThreshold;

    double      trajCutGTZ;
    double      trajCutLTR;
    bool        trajConnect;
    bool        trajNoTransportation;

    bool        writePrimaries;

    /// Ring parameters
    int      nturns;

    double   printModuloFraction;

    /// Visualisation
    int nSegmentsPerCircle; ///< Number of facets per 2pi in visualisation

    /// Scoring Map
    int    nbinsx;
    int    nbinsy;
    int    nbinsz;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    bool   useScoringMap;

    std::string trackingType;
    int         trackingSteps;
    bool        dontUseAperture;

    /// print some properties
    void print() const;
  };
}


#endif //__OPTIONSBASE_H
