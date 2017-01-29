#ifndef BDSFIELDFORMAT_H
#define BDSFIELDFORMAT_H

#include "BDSTypeSafeEnum.hh"
#include "globals.hh"         // geant4 globals / types

/**
 * @brief Type definition for field file formats.
 * 
 * @author Laurie Nevay
 */

struct fieldformats_def
{
  enum type
    {
      none,
      bdsim1d, bdsim2d, bdsim3d, bdsim4d,
      poisson2d, poisson2dquad, poisson2ddipole
    };
};

typedef BDSTypeSafeEnum<fieldformats_def,int> BDSFieldFormat;

namespace BDS
{
  /// Function that gives corresponding enum value for string (case-insensitive)
  BDSFieldFormat DetermineFieldFormat(G4String fieldformat);
}

#endif


