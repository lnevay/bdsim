/*
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway,
University of London 2001 - 2022.

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
#ifndef BDSBUNCHRADIOACTIVEDECAYSOURCE_H
#define BDSBUNCHRADIOACTIVEDECAYSOURCE_H

#include "BDSBunch.hh"

#include "TH3D.h"

#include "G4Transform3D.hh"

namespace CLHEP
{
    class RandFlat;
}

/**
 * @brief A ...
 *
 * @author Eliott Ramoisiaux
 */

class BDSBunchRadioactiveDecaySource: public BDSBunch
{
public:
    BDSBunchRadioactiveDecaySource();
    virtual ~BDSBunchRadioactiveDecaySource();
    virtual void SetOptions(const BDSParticleDefinition* beamParticle,
                            const GMAD::Beam& beam,
                            const BDSBunchType& distrType,
                            G4Transform3D beamlineTransformIn = G4Transform3D::Identity,
                            const G4double beamlineS = 0);
    virtual BDSParticleCoordsFull GetNextParticleLocal();

protected:
    std::string radioactiveDecaySourceFile;
    std::vector<G4String> radioactiveDecaySourceName;
    CLHEP::RandFlat* flatGen;
    std::vector<TH3D*> histVector;
    TH1D* activity;
    std::vector<CLHEP::HepRotation> rotation;
    std::vector<G4ThreeVector> translation, bins;
    double x,y,z;
};

#endif
