#ifndef TRKQuadrupole_h
#define TRKQuadrupole_h

#include "TRKElement.h"

class TRKQuadrupole : public TRKElement, TRKTrackingInterface {
public : 
  TRKQuadrupole();  
  TRKQuadrupole(double length, double b)
  ~TRKQuadrupole();
  
private: 
  
};

#endif
