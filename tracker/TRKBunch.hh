/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2021.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef TRKBunch_h
#define TRKBunch_h

#include <ostream>
#include <vector>
#include <utility>
#include <string>

class BDSParticleDefinition;

namespace GMAD
{
  class Beam;
}

#include "TRKParticle.hh"


/**
 * @brief Vector of particles
 */

class TRKBunch
{
public:
  TRKBunch();
  TRKBunch(const GMAD::Beam& beam,
	   BDSParticleDefinition* designParticle,
	   BDSParticleDefinition* beamParticle,
	   long int nGenerate = 1);
  TRKBunch(const std::vector<TRKParticle>& particleVectorIn);
  TRKBunch(std::vector<TRKParticle>&& particleVectorIn);
  ~TRKBunch();

  typedef std::vector<TRKParticle>::iterator iterator;
  typedef std::vector<TRKParticle>::const_iterator const_iterator;
  iterator       begin()        {return bunch.begin();}
  iterator       end()          {return bunch.end();}
  const_iterator begin()  const {return bunch.begin();}
  const_iterator end()    const {return bunch.end();}
  int            size()   const {return bunch.size();}
  bool           empty()  const {return bunch.empty();}


  /// erase method, returns iterator to element
  iterator Erase(iterator iter) {return bunch.erase(iter);}
  iterator Erase(iterator start, iterator finish) {return bunch.erase(start,finish);}

  /// output stream
  friend std::ostream& operator<< (std::ostream &out, const TRKBunch &beam);

  double mass;
  double charge;
  double totalEnergy;
  double kineticEnergy;

  
private:
  long int population;
  std::vector <TRKParticle> bunch;

  /// populate particles using options and random number generator
  void Populate(const GMAD::Beam& beam,
		BDSParticleDefinition* designParticle,
		BDSParticleDefinition* beamParticle);
  
  std::pair<double,int> GetParticleMassAndCharge(std::string particlename);

};

#endif
