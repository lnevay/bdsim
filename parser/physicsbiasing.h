#ifndef __PHYSICSBIASING_H
#define __PHYSICSBIASING_H

#include <string>
#include <list>
#include <vector>

namespace GMAD {
  struct Array;
}

namespace GMAD {
  enum class PhysicsBiasingType {ALL=1,PRIMARY=2,SECONDARY=3};
  
  /**
   * @brief Physics biasing class for parser
   * 
   * @author Jochem Snuverink <Jochem.Snuverink@rhul.ac.uk> 
   */
  class PhysicsBiasing {
    
  public:
    std::string name; ///< name
    /// particle name
    std::string particle;
    /// geant4 process: single string, but can have multiple processes separated with a space
    std::string process;
    /// factors corresponding to process
    std::list<double> factor;
    /// flag which particles are biased
    std::list<PhysicsBiasingType> flag;
    /// logical volumes where biasing is attached to
    std::string logicalVolumes;
    
    /// constructor
    PhysicsBiasing();
    /// reset
    void clear();
    /// print some properties
    void print()const;

    /// set methods by property name, numeric values
    void set_value(std::string property, double value);
    /// set methods by property name, list of numeric values
    void set_value(std::string property, Array* value);
    /// set methods by property name, string values
    void set_value(std::string property, std::string value);
  };

}

#endif
