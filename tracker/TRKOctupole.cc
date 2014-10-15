#include "TRKOctupole.hh"
#include "TRKStrategy.hh"

/*
#include <cmath>
#include <complex>
#include <cstdlib>
#include "vector3.hh"
#include "vector6.hh"
*/

TRKOctupole::TRKOctupole(double strengthIn, std::string nameIn, double lengthIn, TRKAperture *apertureIn, TRKPlacement *placementIn):
  TRKElement(nameIn,lengthIn,apertureIn,placementIn), strength(strengthIn)
{
}

TRKOctupole::~TRKOctupole() {
}

void TRKOctupole::Track(TRKBunch* bunch, TRKStrategy* strategy) {
  strategy->Track(this,bunch);
}

void TRKOctupole::Print(std::ostream &out)const {
  TRKElement::Print(out);
  out << "; Strength: " << strength << "T/m^3";
}
