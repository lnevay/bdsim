#ifndef OPTIONS_H
#define OPTIONS_H

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "published.h"
#include "optionsBase.h"

namespace GMAD
{
  /**
   * @brief Options class
   * 
   * Options passed with option and beam command. This inherits
   * the OptionsBase class which contains all the members and 
   * provides templated filling functions. This separation allows
   * the OptionsBase class to be more easily written out to ROOT files
   * or other formats for strong reproducibility in a BDSIM run.
   *
   * @author I. Agapov, J. Snuverink
   */
  class Options: public Published<OptionsBase>, public OptionsBase
  {
  public:
    Options();
    explicit Options(const GMAD::OptionsBase& options);
    
    /// set methods by property name
    template<typename T>
    void set_value(std::string name, T value);
    
    /// get method (only for doubles)
    double get_value(std::string name) const;

    /// Take another instance of options and copy the values that have
    /// been set (through setKeys, which although private each instance
    /// has access to as C++ treats encapsulation at the class level).
    /// If override is true, the input option will override the existing
    /// one in this instance.
    void Amalgamate(const Options& optionsIn, bool override);

    /// Whether a parameter has been set using the set_value method or not.
    bool HasBeenSet(std::string name) const;
    
  private:
    /// publish members so these can be looked up from parser
    void PublishMembers();

    /// A list of all the keys that have been set in this instance.
    std::vector<std::string> setKeys;
  };

  template<typename T>
  void Options::set_value(std::string name, T value)
  {
#ifdef BDSDEBUG
    std::cout << "options> Setting value " << std::setw(25) << std::left << name << value << std::endl;
#endif
    // member method can throw runtime_error, catch and exit gracefully
    try
      {
	set(this, name, value);
	setKeys.push_back(name);
      }
    catch (std::runtime_error)
    {
      std::cerr << "Error: options> unknown option \"" << name << "\" with value " << value << std::endl;
      exit(1);
    }
  }
}

#endif
