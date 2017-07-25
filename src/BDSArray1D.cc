#include "BDSArray1D.hh"

#include "G4Types.hh"

#include <ostream>

BDSArray1D::BDSArray1D(G4int nXIn):
  BDSArray2D(nXIn, 1)
{;}

std::ostream& operator<< (std::ostream& out, BDSArray1D const &a)
{
  return a.Print(out);
}
