#include "TRKQuadrupole.hh"

TRKQuadrupole::TRKQuadrupole() {}

TRKQuadrupole::TRKQuadrupole(std::string name, double length, double size_x, double size_y, double str, TRKAperture *aperture, TRKPlacement *placement):
  TRKTrackingElement(TRKTrackingElement::thin, 10,name,length,size_x,size_y,aperture,placement),
  strength(str)
{
  
}
TRKQuadrupole::~TRKQuadrupole() {}

void TRKQuadrupole::Track(const double vIn[], double vOut[], double h){}
void TRKQuadrupole::Track(const double vIn[], double vOut[]){}

