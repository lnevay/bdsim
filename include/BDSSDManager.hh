#ifndef BDSSDMANAGER_H
#define BDSSDMANAGER_H

class BDSEnergyCounterSD;
class BDSSamplerSD;
class BDSTerminatorSD;

/**
 * @brief A singleton class that holds all required sensitive
 * detector class instances.  
 * 
 * Each sensitive detector class
 * need only be instantiated once and attached to the relevant
 * volume. This instantiates all necessary SDs and holds them.
 * 
 * @author Laurie Nevay
 */

class BDSSDManager
{
public:
  static BDSSDManager* Instance(); /// Singleton accessor.

  ~BDSSDManager();

  /// SD for samplers (plane type)
  inline BDSSamplerSD*       GetSamplerPlaneSD() const
  {return samplerPlane;}

  /// SD for samplers (cylinder type)
  inline BDSSamplerSD*       GetSamplerCylinderSD() const
  {return samplerCylinder;}

  /// SD for measuring turns around circular machine and terminating
  /// particles appropriately.
  inline BDSTerminatorSD*    GetTerminatorSD() const
  {return terminator;}
  
  /// SD for general energy counter.
  inline BDSEnergyCounterSD* GetEnergyCounterSD() const
  {return eCounter;}

  /// SD for tunnel energy counter
  inline BDSEnergyCounterSD* GetEnergyCounterTunnelSD() const
  {return tunnelECounter;}

private:
  /// Private default constructor for singleton.
  BDSSDManager();

  BDSSDManager(const BDSSDManager&) = delete;
  BDSSDManager& operator=(const BDSSDManager&) = delete;
 
  static BDSSDManager* _instance;

  /// @{ SD instance.
  BDSSamplerSD*       samplerPlane;
  BDSSamplerSD*       samplerCylinder;
  BDSTerminatorSD*    terminator;
  BDSEnergyCounterSD* eCounter;
  BDSEnergyCounterSD* tunnelECounter;
  /// @}
};

#endif
