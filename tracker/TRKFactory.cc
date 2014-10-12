#include <cstdlib>
#include <iostream>

#include "TRKFactory.hh"

//debug
#include "BDSDebug.hh"

//individual beam elements
#include "TRKLine.hh"
#include "TRKDipole.hh"
#include "TRKDecapole.hh"
#include "TRKQuadrupole.hh"
#include "TRKOctupole.hh"
#include "TRKSextupole.hh"
#include "TRKSampler.hh"

//tracking strategies / routines
#include "TRKStrategy.hh"
#include "TRK.hh"
#include "TRKThin.hh"
#include "TRKThinSymplectic.hh"
#include "TRKThick.hh"
#include "TRKHybrid.hh"

//parser info
#include "parser/element.h"
#include "parser/elementlist.h"
#include "parser/enums.h"
#include "parser/options.h"

//aperture
#include "TRKAperture.hh"
#include "TRKApertureCircular.hh"
#include "TRKApertureRectangular.hh"
#include "TRKApertureEllipsoidal.hh"

#include "CLHEP/Units/SystemOfUnits.h"

TRKFactory::TRKFactory(Options& options) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << "Initialisation" << std::endl;
#endif
  // define charge and momentum from options
  charge = 1;
  momentum = 0;
  // magnetic rigidity
  brho = options.ffact * momentum / 0.299792458 / charge;
  
  /// start placement, could be updated after every new element
  placement = new TRKPlacement();
  /// circular flag? //not available from Options, as command line option!
  // could use GC.nturns instead!
  //  bool circular = options.

  //pull out info from options
  strategy        = setStrategyEnum(options.trackingType);
  aperturetype    = setApertureEnum(options.apertureType);
  beampiperadius  = options.beampipeRadius;
  trackingsteps   = options.trackingSteps;
  defaultaperture = new TRKApertureCircular(beampiperadius);
  dontuseaperture = options.dontUseAperture;
  circular        = options.circular;
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

TRK::Strategy TRKFactory::setStrategyEnum(std::string sIn)
{
  if (sIn == "thin") {return TRK::THIN;}
  else if (sIn == "thinsymplectic") {return TRK::THINSYMPLECTIC;}
  else if (sIn == "thick") {return TRK::THICK;}
  else if (sIn == "hybrid") {return TRK::HYBRID;}
  else {
    std::cout << "Unknown strategy " << sIn << std::endl;
    exit(1);
  }
}

TRKStrategy* TRKFactory::createStrategy() {
  switch(strategy) {

  case TRK::THIN:
    return new TRKThin(trackingsteps);
  case TRK::THINSYMPLECTIC:
    return new TRKThinSymplectic(trackingsteps);
  case TRK::THICK:
    return new TRKThick(trackingsteps);
  case TRK::HYBRID:
    return new TRKHybrid(trackingsteps);
  default:
    return NULL;
  }
}

TRK::Aperture TRKFactory::setApertureEnum(std::string aIn)
{
  if (aIn == "circular")         {return TRK::CIRCULAR;}
  else if (aIn == "rectangular") {return TRK::RECTANGULAR;}
  else if (aIn == "ellipsoidal") {return TRK::ELLIPSOIDAL;}
  else {
    std::cout << "Unknown aperture type " << aIn << std::endl;
    exit(1);
  }
}

TRKAperture* TRKFactory::createAperture(Element& element) {
  //this is annoyingly complex because of the poor implementation of aperture
  //in bdsim.  this is made to allow both general form and individual elements
  //to have their own aperture definitions. individual aperture type is not 
  //possible. it will just default to the general kind.
  //default case = aperturetype
  
  if (dontuseaperture) {
    return NULL; //no aperture at all - will never be check with this setting so no seg fault
  }
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
}

TRKLine* TRKFactory::createLine(ElementList& beamline_list) {
  TRKLine* line = new TRKLine("beamline",circular);
  std::list<struct Element>::iterator it = beamline_list.begin();

  for(;it!=beamline_list.end();it++){
    TRKElement* element = createElement(*it);
    if (element) {
      line->AddElement(element);
#ifdef TRKDEBUG
      std::cout << __METHOD_NAME__ << "element created: " << *element << std::endl;
#endif
      // update placement
    }
  }
  return line;
}


TRKElement* TRKFactory::createElement(Element& element) {
  switch (element.type) {

  case _LINE:
    return createLine(element);
  case _DRIFT:
    return createDrift(element);
  case _PCLDRIFT:
    return createDrift(element);
  case _SBEND:
    //TEMPORARY
    return createDrift(element);
  case _RBEND:
    //TEMPORARY
    return createDrift(element);;
  case _QUAD:
    return createQuadrupole(element);
  case _SEXTUPOLE:
    return createSextupole(element);
  case _OCTUPOLE:
    return createOctupole(element);
    //  case _DECAPOLE:
    //TEMPORARY
    //return createDrift(element);
  case _SOLENOID:
    //TEMPORARY
    return createDrift(element);
  case _MULT:
    //TEMPORARY
    return createDrift(element);
  case _ELEMENT:
    //TEMPORARY
    return createDrift(element);
  case _SAMPLER:
    return createSampler(element);
  case _TRANSFORM3D:
    //TEMPORARY
    return NULL;
  case _VKICK:
    //TEMPORARY
    return createDrift(element);
  case _HKICK:
    return createDipole(element);
  default:
    return NULL;
  }
  /* types not specified - see parser/enums.h for types
  _SCREEN
  _CSAMPLER
  */
}

TRKElement* TRKFactory::createLine(Element& /*element*/) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  // method not needed?
  return NULL;
}

TRKElement* TRKFactory::createDrift(Element& element) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  TRKAperture* aperture = createAperture(element);
  return new TRKDrift(element.name,
		      element.l,
		      aperture,
		      placement);
}

TRKElement* TRKFactory::createDipole(Element& /*element*/) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  // bfield , see componentfactory and bdskicker.cc
  // strength (bprime)
  //TRKAperture* aperture = createAperture(element);
  return NULL;
}

TRKElement* TRKFactory::createQuadrupole(Element& element) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  double bPrime = - brho * (element.k1 / CLHEP::m2);
  TRKAperture* aperture = createAperture(element);
  return new TRKQuadrupole(bPrime,
			   element.name,
			   element.l,
			   aperture,
			   placement);
}

TRKElement* TRKFactory::createSextupole(Element& element) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  double bPrime = - brho * (element.k2 / CLHEP::m3); // to be checked
  TRKAperture* aperture = createAperture(element);
  return new TRKSextupole(bPrime,
			  element.name,
			  element.l,
			  aperture,
			  placement);
}

TRKElement* TRKFactory::createOctupole(Element& element) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  double bPrime = - brho * (element.k3 / CLHEP::m2 / CLHEP::m2); // to be checked
  TRKAperture* aperture = createAperture(element);
  return new TRKOctupole(bPrime,
			 element.name,
			 element.l,
			 aperture,
			 placement);
}

TRKElement* TRKFactory::createDecapole(Element& /*element*/) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  //TRKAperture* aperture = createAperture(element);
  return NULL;
}

TRKElement* TRKFactory::createSampler(Element& element) {
#ifdef TRKDEBUG
  std::cout << __METHOD_NAME__ << std::endl;
#endif
  return new TRKSampler(element.name);
}
