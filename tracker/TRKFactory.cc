/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2019.

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
#include <iostream>

#include "BDSDebug.hh"
#include "BDSParticleDefinition.hh"
#include "BDSSamplerRegistry.hh"

#include "TRKAperture.hh"
#include "TRKApertureCircular.hh"
#include "TRKApertureRectangular.hh"
#include "TRKApertureEllipsoidal.hh"
#include "TRKFactory.hh"
#include "TRKDipole.hh"
#include "TRKElement.hh"
#include "TRKDecapole.hh"
#include "TRKLine.hh"
#include "TRKOctupole.hh"
#include "TRKQuadrupole.hh"
#include "TRKRBend.hh"
#include "TRKSampler.hh"
#include "TRKSextupole.hh"
#include "TRKSolenoid.hh"
#include "TRKSBend.hh"

//tracking strategies / routines
#include "TRKStrategy.hh"
#include "TRK.hh"
#include "TRKParticleDefinition.hh"
#include "TRKThin.hh"
#include "TRKThinSymplectic.hh"
#include "TRKThick.hh"
#include "TRKHybrid.hh"

#include "parser/beam.h"
#include "parser/element.h"
#include "parser/elementtype.h"
#include "parser/options.h"

#include "CLHEP/Units/SystemOfUnits.h"

TRKFactory::TRKFactory(const GMAD::Options&   options,
		       BDSParticleDefinition* particle,
		       BDSOutput*             outputIn):
  output(outputIn)
{
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << "Initialisation" << std::endl;
#endif
  // define charge and momentum from options
  charge   = particle->Charge();
  momentum = particle->Momentum();
  energy   = particle->TotalEnergy();
  brho     = particle->BRho();
  
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << "Rigidity (Brho) : " << brho/(CLHEP::tesla*CLHEP::m) << " T*m" << std::endl;
  std::cout << "ffact    : " << options.ffact << std::endl;
  std::cout << "momentum : " << momentum << " GeV" << std::endl;
  std::cout << "charge   : " << charge << std::endl;
#endif
  
  /// start placement, could be updated after every new element
  placement = nullptr;//new TRKPlacement();

  //pull out info from options
  strategy        = SetStrategyEnum(options.trackingType);
  aperturetype    = SetApertureEnum(options.apertureType);
  beampiperadius  = options.aper1;
  trackingsteps   = options.trackingSteps;
  defaultaperture = new TRKApertureCircular(beampiperadius);
  dontuseaperture = options.dontUseAperture;
  if (options.nturns > 1)
    {circular = true;}
  else
    {circular = false;}
}

std::ostream& operator<< (std::ostream& out, const TRKFactory& factory)
{
  out << "TRKFactory Instance - details:"                      << std::endl
      << "Tracking Strategy (enum): " << factory.strategy         << std::endl
      << "Aperture Type:            " << factory.aperturetype     << std::endl
      << "Beam Pipe Radius (m):     " << factory.beampiperadius   << std::endl
      << "Tracking Steps:           " << factory.trackingsteps    << std::endl
      << "Default Aperture Type:    " << *factory.defaultaperture << std::endl
      << "Don't Use Aperture:       " << factory.dontuseaperture;
  return out;
}

TRK::Strategy TRKFactory::SetStrategyEnum(std::string sIn)
{
  if (sIn == "thin") {return TRK::THIN;}
  else if (sIn == "thinsymplectic") {return TRK::THINSYMPLECTIC;}
  else if (sIn == "thick") {return TRK::THICK;}
  else if (sIn == "hybrid") {return TRK::HYBRID;}
  else
    {
      std::cout << "Unknown strategy " << sIn << std::endl;
      exit(1);
    }
}

TRKStrategy* TRKFactory::CreateStrategy()
{
  TRKStrategy* result = nullptr;
  switch(strategy)
    {
    case TRK::THIN:
      {result = new TRKThin(trackingsteps); break;}
    case TRK::THINSYMPLECTIC:
      {result = new TRKThinSymplectic(trackingsteps); break;}
    case TRK::THICK:
      {result = new TRKThick(trackingsteps); break;}
    case TRK::HYBRID:
      {result = new TRKHybrid(trackingsteps); break;}
    default:
      {break;}
    }
  if (result)
    {result->SetMomentumAndEnergy(momentum, energy);}
  return result;
}

TRK::Aperture TRKFactory::SetApertureEnum(std::string aIn)
{
  if      (aIn == "circular")    {return TRK::CIRCULAR;}
  else if (aIn == "rectangular") {return TRK::RECTANGULAR;}
  else if (aIn == "ellipsoidal") {return TRK::ELLIPSOIDAL;}
  else
    {
      std::cout << "Unknown aperture type " << aIn << std::endl;
      exit(1);
    }
}

TRKAperture* TRKFactory::CreateAperture(GMAD::Element& element)
{
  //this is annoyingly complex because of the poor implementation of aperture
  //in bdsim.  this is made to allow both general form and individual elements
  //to have their own aperture definitions. individual aperture type is not 
  //possible. it will just default to the general kind.
  //default case = aperturetype
  
  if (dontuseaperture)
    {return NULL;} //no aperture at all - will never be check with this setting so no seg fault
  else
    {return new TRKApertureCircular(element.aper1);}
  /*
    // THIS SHOULD BE FIXED GIVEN THE NEW APERTURE MODELS AVAILABLE
  else if ((element.aperX != 0) && (element.aperY !=0)) {
    //must have been specified - now check whether one of the asymmetric aperture types is specified
    if (aperturetype == TRK::RECTANGULAR) {
      return new TRKApertureRectangular(element.aperX,element.aperY);}
    else {
      //in effect the default if two apertures are specified for an element
      return new TRKApertureEllipsoidal(element.aperX,element.aperY);}
  }
  else if (aperturetype == TRK::RECTANGULAR) {
    //no specified x and y aper, but told its rectangular -> square
    return new TRKApertureRectangular(beampiperadius,beampiperadius);}
  else {
    //must be circular then
    return defaultaperture;}
  */
}

TRKLine* TRKFactory::CreateLine(const GMAD::FastList<GMAD::Element>& beamline_list)
{
  TRKLine* line = new TRKLine("beamline",circular);
  // Quick fix to get Primary branch filled.  index=-1 being key here.
  line->AddElement(new TRKSampler("primaries", -1, output));
  for (auto it : beamline_list)
    {
      TRKElement* element = CreateElement(it);
      if (element)
	{
	  line->AddElement(element);
	  // if the element is flagged as a sampler, we create
	  // a sampler element and put it in the beam line
	  if (it.samplerType != "none")
	    {
	      TRKElement* sampler = CreateSampler(it);
	      line->AddElement(sampler);
	    }
	}
    }
  return line;
}


TRKElement* TRKFactory::CreateElement(GMAD::Element& element)
{
  TRKElement* trkelement = NULL;
  switch (element.type)
    {
    case GMAD::ElementType::_LINE:
      trkelement = CreateLine(element);
      break;
    case GMAD::ElementType::_DRIFT:
      trkelement = CreateDrift(element);
      break;
    case GMAD::ElementType::_SBEND:
      trkelement = CreateSBend(element);
      break;
    case GMAD::ElementType::_RBEND:
      trkelement = CreateRBend(element);;
      break;
    case GMAD::ElementType::_QUAD:
      trkelement = CreateQuadrupole(element);
      break;
    case GMAD::ElementType::_SEXTUPOLE:
      trkelement = CreateSextupole(element);
      break;
    case GMAD::ElementType::_OCTUPOLE:
      trkelement = CreateOctupole(element);
      break;
    // case GMAD::ElementType::_DECAPOLE:
    //   trkelement = CreateDecapole(element);
    //    break;
    case GMAD::ElementType::_SOLENOID:
      trkelement = CreateSolenoid(element);
      break;
    case GMAD::ElementType::_MULT:
      //TEMPORARY
      trkelement = CreateDrift(element);
      break;
    case GMAD::ElementType::_ELEMENT:
      //TEMPORARY
      trkelement = CreateDrift(element);
      break;
    case GMAD::ElementType::_RCOL:
      trkelement = CreateDrift(element);
      break;
    case GMAD::ElementType::_TRANSFORM3D:
      //TEMPORARY
      trkelement = NULL;
      break;
    case GMAD::ElementType::_VKICKER:
    case GMAD::ElementType::_HKICKER:
      trkelement = CreateDipole(element);
      break;
    default:
      trkelement = NULL;
      break;
  }
  // TBC - thin hkicker and vkicker
  // TBC - implement sampler here based on element.samplerType - str - defualt 'none'

  if (trkelement)
    {AddCommonProperties(trkelement,element);}

  return trkelement;
}

void TRKFactory::AddCommonProperties(TRKElement* trkelement, GMAD::Element& element)
{
  // offset and tilt
  //  if (element.phi!=0.0 || element.theta!=0.0 || element.psi!=0.0) {
  // only tilt for now
  if (element.tilt!=0.0)
    {trkelement->SetTilt(element.tilt,0,0);} // todo check rotation correct!
  if (element.xdir!=0.0 || element.ydir!=0.0)
    {trkelement->SetOffset(element.xdir, element.ydir);}
}

TRKElement* TRKFactory::CreateLine(GMAD::Element& /*element*/)
{
  return nullptr;
}

TRKElement* TRKFactory::CreateDrift(GMAD::Element& element)
{
  TRKAperture* aperture = CreateAperture(element);
  return new TRKDrift(element.name,
		      element.l,
		      aperture,
		      placement);
}

TRKElement* TRKFactory::CreateDipole(GMAD::Element& /*element*/)
{
  // bfield , see componentfactory and bdskicker.cc
  // strength (bprime)
  //TRKAperture* aperture = CreateAperture(element);
  return NULL;
}

TRKElement* TRKFactory::CreateQuadrupole(GMAD::Element& element)
{
  TRKAperture* aperture = CreateAperture(element);
  return new TRKQuadrupole(element.k1,
			   element.name,
			   element.l,
			   aperture,
			   placement);
}

TRKElement* TRKFactory::CreateSextupole(GMAD::Element& element)
{
  double bPrime = - brho * (element.k2 / CLHEP::m3); // to be checked
  TRKAperture* aperture = CreateAperture(element);
  return new TRKSextupole(bPrime,
			  element.name,
			  element.l,
			  aperture,
			  placement);
}

TRKElement* TRKFactory::CreateOctupole(GMAD::Element& element)
{
  double bPrime = - brho * (element.k3 / CLHEP::m2 / CLHEP::m2); // to be checked
  TRKAperture* aperture = CreateAperture(element);
  return new TRKOctupole(bPrime,
			 element.name,
			 element.l,
			 aperture,
			 placement);
}

TRKElement* TRKFactory::CreateDecapole(GMAD::Element& /*element*/)
{
  //TRKAperture* aperture = CreateAperture(element);
  return NULL;
}

TRKElement* TRKFactory::CreateSolenoid(GMAD::Element& element)
{
  //
  // magnetic field
  //
  // B = B/Brho * Brho = ks * Brho
  // brho is in Geant4 units, but ks is not -> multiply ks by m^-1
  double bField;
  if(element.B != 0){
    bField = element.B * CLHEP::tesla;
    //    element.ks  = (bField/brho) / CLHEP::m;
  }
  else{
    bField = (element.ks/CLHEP::m) * brho;
    //    element.B = bField/CLHEP::tesla;
  }

  TRKAperture* aperture = CreateAperture(element);
  return new TRKSolenoid(bField,
			 element.name,
			 element.l,
			 aperture,
			 placement);
}

TRKElement* TRKFactory::CreateSBend(GMAD::Element& element)
{
  double angle;
  // B Field used to define bend: theta = qBl / p 
  if(element.B != 0)
  {    
    angle = charge * (element.B * CLHEP::tesla) * element.l / momentum;
  }
  else
  {
    angle = element.angle;
  }
  
  TRKAperture* aperture = CreateAperture(element);
  return new TRKSBend(angle,
      element.name,
      element.l,
      aperture,
      placement);
}

TRKElement* TRKFactory::CreateRBend(GMAD::Element& element)
{
  TRKAperture* aperture = CreateAperture(element);
  return new TRKRBend(element.angle,
      element.name,
      element.l,
      aperture,
      placement);
}

TRKElement* TRKFactory::CreateSampler(GMAD::Element& element)
{
  std::string name = element.name;
  int samplerIndex = BDSSamplerRegistry::Instance()->RegisterSampler(name, nullptr);
  TRKElement* result = new TRKSampler(name, samplerIndex, output);
  return result;
}
