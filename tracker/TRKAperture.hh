#ifndef TRKAperture_h
#define TRKAperture_h 1

#include "TRK.hh"

#include "TRKParticle.hh"
/**
 * @brief Aperture base class
 */

class TRKAperture {
public : 
  TRKAperture(){};
  virtual ~TRKAperture(){}; 

  /// checks if aperture is hit
  virtual bool OutsideAperture(const TRKParticle& particle)=0;
  //for now literally checks if outside the aperture
  //in future will check including the rigidity and strength
  //if we have to give this check function the strength, this could
  //be complicated because each element is different in type..
  //we could therefore base the heuristic on the relativistic gamma
  //and the length of the element - the gamma will give us a rough measure
  //of rigidity and combined with length we can set an approximate fraction
  //of the aperture to multiply by for the test 
  //TBC!!

  //Should also provide a << cout stream implementation for each and declare here

};


#endif 
