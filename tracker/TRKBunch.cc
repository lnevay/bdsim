/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2020.

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
#include <cstdlib>
#include <string>
#include <iostream>
#include <utility>

#include "TRK.hh"
#include "TRKBunch.hh"
#include "TRKParticle.hh"
#include "TRKParticleDefinition.hh"

#include "BDSBunch.hh"
#include "BDSBunchFactory.hh"
#include "BDSDebug.hh"
#include "BDSParticleCoordsFullGlobal.hh"
#include "BDSParticleDefinition.hh"

#include "parser/beam.h"

#include "CLHEP/Units/SystemOfUnits.h"

TRKBunch::TRKBunch(const GMAD::Beam& beam,
		   BDSParticleDefinition* particle,
		   long int nGenerate):
  population(nGenerate)
{
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << "Initialisation" << std::endl;
#endif
  mass   = particle->Mass();
  charge = particle->Charge();
  totalEnergy   = particle->TotalEnergy();
  kineticEnergy = particle->KineticEnergy();

  //populate particles using options & random number generator
  Populate(beam, particle);
}

TRKBunch::TRKBunch(const std::vector<TRKParticle>& particleVectorIn)
{
  bunch = particleVectorIn;
}

TRKBunch::~TRKBunch()
{;}

void TRKBunch::Populate(const GMAD::Beam& beam,
			BDSParticleDefinition* particle)
{
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << "Generating particles" << std::endl;
#endif
  //particle definition / generation
  //just now this allows every particle to have a different
  //mass and charge
  //if we assume all the same, then it could be done on a global
  //basis, which would save around 20% memory on each particle...

  // Initialise bunch
  BDSBunch* bdsbunch = BDSBunchFactory::CreateBunch(particle, beam);

  //must have positive number of particles
  if (population < 0)
    {population = std::abs(population);}
  //must have at least 1 particle
  if (population == 0)
    {population = 1;}
  //initialise the vector of particles
  bunch.reserve(population);

  BDSParticleDefinition const *referenceParticle =
      bdsbunch->ParticleDefinition();

  auto p0 = referenceParticle->Momentum();
  auto gamma0 = referenceParticle->Gamma();
  auto beta0 = referenceParticle->Beta();

  for (int i = 0; i < population; i++)
    {
      // bdsbunch generates values in CLHEP mm standard.
      BDSParticleCoordsFullGlobal c = bdsbunch->GetNextParticle();
      double energy = c.local.totalEnergy;
      double p = std::sqrt(energy*energy - mass*mass);

      // Momenta from BDSBunch are px/|p|, but we want px/p0
      // (i.e. normalised w.r.t reference momentum, not particle momentum).
      auto px = c.local.xp * p / p0;
      auto py = c.local.yp * p / p0;

      auto z = 0; // ???  Not sure what to do with this one yet.
      auto pz = energy / p0 - 1/beta0;

      bunch.emplace_back(c.local.x, px, c.local.y, py,
			 z, pz,  // longitudinal coordinates
			 beta0, gamma0 // reference relativistic
			 );

      //weight not required - maybe should be kept though to pass on to bdsim
    }
}

std::ostream& operator<< (std::ostream &out, const TRKBunch &beam)
{
  std::vector<TRKParticle>::const_iterator iter = beam.bunch.begin();
  std::vector<TRKParticle>::const_iterator end  = beam.bunch.end();

  for (;iter!=end;++iter)
    {out << *iter << std::endl;}
  return out;
}
