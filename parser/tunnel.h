#ifndef TUNNEL_H
#define TUNNEL_H

#include <iomanip>
#include <iostream>
#include <string>

#include "published.h"

namespace GMAD
{
  /**
   * @brief Tunnel class for parser
   * 
   * @author Jochem Snuverink
   */
  class Tunnel: public Published<Tunnel> {

  public:
    std::string name; ///< name
    /// geometry type
    std::string type;

    /// radius, aperture parameters
    double   aper1;
    double   aper2;
    /// offset x and y
    double   offsetX, offsetY;
    /// tunnel geometry parameters
    double   thickness;
    double   soilThickness;
    double   floorOffset;

    /// visibility
    int      visible;
  
    /// material
    std::string material;
    /// soil
    std::string soilMaterial;

    /// start and end element by name
    std::string startElement, endElement;
  
    /// constructor
    Tunnel();
    /// reset
    void clear();
    /// print some properties
    void print()const;
    /// set methods by property name and value
    template <typename T>
      void set_value(std::string property, T value);

  private:
    /// Publish members
    void PublishMembers();
  };
  
  template <typename T>
    void Tunnel::set_value(std::string property, T value)
    {
#ifdef BDSDEBUG
      std::cout << "tunnel> Setting value " << std::setw(25) << std::left << property << value << std::endl;
#endif
      // member method can throw runtime_error, catch and exit gracefully
      try {
        set(this,property,value);
      }
      catch(std::runtime_error) {
        std::cerr << "Error: tunnel> unknown option \"" << property << "\" with value " << value  << std::endl;
        exit(1);
      }
    }
}

#endif
