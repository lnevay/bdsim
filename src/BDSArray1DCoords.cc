#include "BDSArray1DCoords.hh"

#include <ostream>

#include "globals.hh"

BDSArray1DCoords::BDSArray1DCoords(G4int nXIn,
				   G4double xMinIn, G4double xMaxIn):
  BDSArray2DCoords(nXIn,1,
		   xMinIn,xMaxIn,
		   0,   1)
{;}

std::ostream& operator<< (std::ostream& out, BDSArray1DCoords const &a)
{
  return a.Print(out);
}
