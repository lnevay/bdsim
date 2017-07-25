#include "element.h"
#include "elementtype.h"
#include "parameters.h"
#include "parser.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace GMAD;

namespace {
  // helper method
  void print(std::list<Element> l, int ident=0)
  {
    if(ident == 0) std::cout << "using line " << Parser::Instance()->current_line << std::endl;
  
    for(std::list<Element>::iterator it=l.begin();it!=l.end();++it)
      {
	(*it).print(ident);
      }
  }
}

Element::Element():
  lst(nullptr)
{
  flush();

  PublishMembers();
}

void Element::PublishMembers()
{
  publish("l",    &Element::l);
  publish("ks",   &Element::ks);
  publish("k1",   &Element::k1);
  publish("k2",   &Element::k2);
  publish("k3",   &Element::k3);
  publish("k4",   &Element::k4);
  publish("angle",&Element::angle);
  publish("B",    &Element::B);
  publish("e1",   &Element::e1);
  publish("e2",   &Element::e2);
  publish("fint", &Element::fint);
  publish("fintx",&Element::fintx);
  publish("hgap", &Element::hgap);
  publish("kick", &Element::kick);
  publish("hkick",&Element::hkick);
  publish("vkick",&Element::vkick);
  publish("knl",  &Element::knl);
  publish("ksl",  &Element::ksl);
  publish("gradient",  &Element::gradient);
  publish("E",         &Element::E);
  publish("frequency", &Element::frequency);
  publish("phase",     &Element::phase);
  publish("tOffset",   &Element::tOffset);
  
  publish("beampipeThickness",&Element::beampipeThickness);
  publish("aper",&Element::aper1);
  alternativeNames["aper"] = "aper1";
  publish("aperture", &Element::aper1);
  alternativeNames["aperture"] = "aper1";
  publish("aper1",    &Element::aper1);
  publish("aperture1",&Element::aper1);
  alternativeNames["aperture1"] = "aper1";
  publish("aper2",    &Element::aper2);
  publish("aperture2",&Element::aper2);
  alternativeNames["aperture2"] = "aper2";
  publish("aper3",    &Element::aper3);
  publish("aperture3",&Element::aper3);
  alternativeNames["aperture3"] = "aper3";
  publish("aper4",    &Element::aper4);
  publish("aperture4",&Element::aper4);
  alternativeNames["aperture4"] = "aper4";
  publish("outerDiameter",&Element::outerDiameter);
  publish("xsize",&Element::xsize);
  publish("ysize",&Element::ysize);
  publish("xsizeOut",&Element::xsizeOut);
  publish("ysizeOut",&Element::ysizeOut);
  publish("tilt",&Element::tilt);

  publish("offsetX",&Element::offsetX);
  publish("offsetY",&Element::offsetY);
  publish("x",&Element::xdir);
  alternativeNames["x"] = "xdir";
  publish("y",&Element::ydir);
  alternativeNames["y"] = "ydir";
  publish("z",&Element::zdir);
  alternativeNames["z"] = "zdir";
  publish("xdir",&Element::xdir);
  publish("ydir",&Element::ydir);
  publish("zdir",&Element::zdir);
  publish("phi",&Element::phi);
  publish("theta",&Element::theta);
  publish("psi",&Element::psi);
  publish("region",&Element::region);
  publish("fieldOuter",  &Element::fieldOuter);
  publish("fieldVacuum", &Element::fieldVacuum);
  publish("fieldAll",    &Element::fieldAll);
  publish("waveLength",&Element::waveLength);

  // screen
  publish("tscint",&Element::tscint);
  publish("twindow",&Element::twindow);
  publish("tmount",&Element::tmount);
  publish("windowScreenGap",&Element::windowScreenGap);
  publish("screenXSize",&Element::screenXSize);
  publish("screenYSize",&Element::screenYSize);
  publish("layerThicknesses",&Element::layerThicknesses);
  publish("layerMaterials",&Element::layerMaterials);
  publish("layerIsSampler",&Element::layerIsSampler);

  // for AWAKE spectrometer
  publish("screenPSize",&Element::screenPSize);
  publish("screenEndZ",&Element::screenEndZ);
  publish("poleStartZ",&Element::poleStartZ);
  publish("screenWidth",&Element::screenWidth);
  publish("awakeMagnetOffsetX", &Element::awakeMagnetOffsetX);
  
  publish("numberWedges",&Element::numberWedges);
  publish("wedgeLength",&Element::wedgeLength);
  publish("degraderHeight",&Element::degraderHeight);
  publish("materialThickness",&Element::materialThickness);
  publish("degraderOffset",&Element::degraderOffset);

  publish("geometry",&Element::geometryFile);
  publish("bmap",    &Element::fieldAll);
  alternativeNames["bmap"] = "fieldAll";
  publish("outerMaterial",&Element::outerMaterial);
  publish("material",&Element::material);
  publish("yokeOnInside", &Element::yokeOnInside);
  publish("apertureType",&Element::apertureType);
  publish("magnetGeometryType",&Element::magnetGeometryType);
  publish("beampipeMaterial",&Element::beampipeMaterial);
  publish("vacuumMaterial",&Element::vacuumMaterial);
  publish("scintmaterial",&Element::scintmaterial);
  publish("windowmaterial",&Element::windowmaterial);
  publish("mountmaterial",&Element::mountmaterial);
  publish("spec",&Element::spec);
  publish("cavityModel",&Element::cavityModel);
  publish("bias",&Element::bias);
  publish("biasMaterial",&Element::biasMaterial);
  publish("biasVacuum",&Element::biasVacuum);
  publish("samplerName",&Element::samplerName);
  publish("samplerType",&Element::samplerType);
  publish("r",&Element::samplerRadius); // historic
  publish("samplerRadius",&Element::samplerRadius);
  alternativeNames["r"] ="samplerRadius";
  
  publish("blmLocZ",&Element::blmLocZ);
  publish("blmLocTheta",&Element::blmLocTheta);

  publish("colour", &Element::colour);
}

std::string Element::getPublishedName(std::string nameIn)const
{
  auto it = alternativeNames.find(nameIn);
  if (it != alternativeNames.end()) {
    return it->second;
  }
  // if not found return name
  return nameIn;
}

bool Element::isSpecial()const {
  bool isSpecial = false;

  if (type == ElementType::_TRANSFORM3D ||
      type == ElementType::_MARKER ||
      type == ElementType::_LINE ||
      type == ElementType::_REV_LINE )
    {isSpecial = true;}

  return isSpecial;
}

void Element::print(int ident)const{
  for(int i=0;i<ident;i++)
    {std::cout << "--";}

  std::cout << name << " : " << type << std::endl;
  if (l>0.0)
    {std::cout << "l     = " << l << "m" << std::endl;}
  if (samplerType != "none")
    {std::cout << "samplerType = " << samplerType << std::endl;}

  switch(type) {
  case ElementType::_DRIFT:
    break;
  case ElementType::_SBEND:
  case ElementType::_RBEND:
    std::cout << "B     = " << B     << std::endl
	      << "angle = " << angle << std::endl
	      << "k1    = " << k1    << std::endl;
    break;
  case ElementType::_QUAD:
    std::cout << "k1    = " << k1    << std::endl;
    break;
  case ElementType::_SEXTUPOLE:
    std::cout << "k2    = " << k2    << std::endl;
    break;
  case ElementType::_OCTUPOLE:
    std::cout << "k3    = " << k3    << std::endl;
    break;
  case ElementType::_DECAPOLE:
    std::cout << "k4    = " << k4 << std::endl;
    break;
  case ElementType::_SOLENOID:
    std::cout << "ks    = " << ks << std::endl;
    break;

  case ElementType::_MULT:
  case ElementType::_THINMULT:
    std::cout << " , knl={";
    for(auto it=knl.begin();it!=knl.end();++it)
      {std::cout << (*it);}
    std::cout << "},  ksl={";
    for(auto it=ksl.begin();it!=ksl.end();++it)
      {std::cout << (*it);}
    std::cout << "}" << std::endl;
    break;
    
  case ElementType::_ECOL:
  case ElementType::_RCOL:
    std::cout << "x half aperture = " << xsize <<" m" << std::endl
	      << "y half aperture = " << ysize <<" m" << std::endl
	      << "material = \""      << material << "\"" << std::endl;
    break;

  case ElementType::_ELEMENT:
    std::cout << "outerDiameter = "  << outerDiameter << "m" << std::endl
	      << "precision region " << region << std::endl
	      << "Geometry file : "  << geometryFile << std::endl
	      << "Field object  : "  << fieldAll << std::endl;
    break;
    
  case ElementType::_AWAKESCREEN:
    std::cout << "twindow         = " << twindow*1e6         << " um" << std::endl
	      << "tscint          = " << tscint*1e6          << " um" << std::endl
	      << "windowScreenGap = " << windowScreenGap*1e6 << " um" << std::endl
	      << "windowmaterial  = " << windowmaterial      << std::endl
	      << "scintmaterial   = " << scintmaterial       << std::endl;
      break;

  case ElementType::_AWAKESPECTROMETER:
    std::cout << "twindow         = " << twindow*1e6         << " um" << std::endl
	      << "tscint          = " << tscint*1e6          << " um" << std::endl
	      << "screenPSize     = " << screenPSize*1e6     << " um" << std::endl
	      << "windowScreenGap = " << windowScreenGap*1e6 << " um" << std::endl
	      << "windowmaterial  = " << windowmaterial      << std::endl
	      << "tmount          = " << tmount*1e6          << " um" << std::endl
	      << "mountmaterial   = " << mountmaterial       << std::endl
	      << "scintmaterial   = " << scintmaterial       << std::endl;
    break;

  case ElementType::_LASER:
    std::cout << "lambda = " << waveLength << "m" << std::endl
	      << "xSigma = " << xsize << "m" << std::endl
	      << "ySigma = " << ysize << "m" << std::endl
	      << "xdir = "   << xdir << std::endl
	      << "ydir = "   << ydir << std::endl
	      << "zdir = "   << zdir << std::endl;
    break;

  case ElementType::_SCREEN:
    std::cout << "angle=" << angle <<"rad" << std::endl
	      << "precision region " << region << std::endl;
    break;
    
  case ElementType::_TRANSFORM3D:
    std::cout << "xdir= "  << xdir    << "m" << std::endl
	      << "ydir= "  << ydir    << "m" << std::endl
	      << "zdir= "  << zdir    << "m" << std::endl
	      << "phi= "   << phi   << "rad" << std::endl
	      << "theta= " << theta << "rad" << std::endl
	      << "psi= "   << psi   << "rad" << std::endl;
    break;
  default:
    break;
  }
  
  if(lst != nullptr)
    {
      ::print(*lst,++ident);
    }
}

void Element::flush()
{
  type = ElementType::_NONE;
  name = "";
  l = 0;
  ks = 0;
  k1 = 0;
  k2 = 0;
  k3 = 0;
  k4 = 0;
  angle = 0;
  B = 0;
  e1 = 0;
  e2 = 0;
  fint = 0;
  fintx = 0;
  hgap  = 0;
  kick  = 0;
  hkick = 0;
  vkick = 0;
  knl.clear();
  ksl.clear();
  gradient  = 0;
  E         = 0;
  frequency = 0;
  phase     = 0;
  tOffset   = 0;
  
  // degrader
  numberWedges = 1;
  wedgeLength = 0;
  degraderHeight = 0;
  materialThickness = 0;
  degraderOffset = 0;

  // new aperture model
  beampipeThickness = 0;
  aper1 = 0;
  aper2 = 0;
  aper3 = 0;
  aper4 = 0;
  apertureType = "";
  beampipeMaterial = "";
  vacuumMaterial = "";

  // magnet geometry
  magnetGeometryType  = "";
  outerMaterial = "";
  outerDiameter = 0;
  yokeOnInside  = true;
  
  tilt = 0;
  xsize = 0;
  ysize = 0;
  xsizeOut = 0;
  ysizeOut = 0;
  offsetX = 0;
  offsetY = 0;
  tscint = 0.0003;
  twindow = 0;
  tmount = 0;
  windowScreenGap = 0;
  screenXSize = 0;
  screenYSize = 0;

  screenPSize        = 0;
  screenEndZ         = 0;
  poleStartZ         = 0;
  screenWidth        = 0;
  awakeMagnetOffsetX = 0.13;

  layerThicknesses.clear();
  layerMaterials.clear();
  layerIsSampler.clear();
  
  xdir = 0;
  ydir = 0;
  zdir = 0;
  waveLength = 0;
  gradient = 0;
  phi = 0;
  theta = 0;
  psi = 0;
  
  blmLocZ.clear();
  blmLocTheta.clear();

  bias = ""; biasMaterial=""; biasVacuum="";
  biasMaterialList.clear();
  biasVacuumList.clear();

  samplerName = "";
  samplerType = "none"; // allowed "none", "plane", "cylinder"
  samplerRadius = 0;
  
  region      = "";
  fieldOuter  = "";
  fieldVacuum = "";
  fieldAll    = "";

  geometryFile ="";
  material="";  
  windowmaterial = "vacuum";
  mountmaterial="";
  scintmaterial = "";
  spec = "";
  cavityModel = "";

  colour = "";

  angleSet = false;
}

double Element::property_lookup(std::string property_name)const{
  double value;
  try {
    value = get<double>(this,property_name);
  }
  catch (std::runtime_error) {
    std::cerr << "element.cc> Error: unknown property \"" << property_name << "\" (only works on numerical properties)" << std::endl; 
    exit(1);
  }
  return value;
}

void Element::set(const Parameters& params,std::string nameIn, ElementType typeIn)
{
  // common parameters for all elements
  type = typeIn;
  name = nameIn;
  
  set(params);

  if (params.setMap.at("angle"))
    {angleSet = true;}
}

void Element::set(const Parameters& params)
{
  for (auto& i : params.setMap)
    {
      if(i.second == true)
	{
	  std::string property = i.first;

	  // method can in theory throw runtime_error (shouldn't happen), catch and exit gracefully
	  try {
	    Published<Element>::set(this,(Element*)&params,property);
	  }
	  catch(std::runtime_error) {
	    std::cerr << "Error: parser> unknown property \"" << property << "\" for element " << name  << std::endl;
	    exit(1);
	  }

	  // split bias into tokens and add to both material and vacuum
	  if (property == "bias")
	    {
	      std::stringstream ss(bias);
	      std::string tok;
	      while(ss >> tok)
		{
		  biasMaterialList.push_back(tok);
		  biasVacuumList.push_back(tok);
		}
	    }
	  else if (property == "biasMaterial")
	    {
	      std::stringstream ss(biasMaterial);
	      std::string tok;
	      while(ss >> tok) {biasMaterialList.push_back(tok);}
	    }
	  else if (property == "biasVacuum")
	    {
	      std::stringstream ss(biasVacuum);
	      std::string tok;
	      while(ss >> tok) {biasVacuumList.push_back(tok);}
	    }
	}
    }
}

void Element::setSamplerInfo(std::string samplerTypeIn, std::string samplerNameIn, double samplerRadiusIn)
{
  if (samplerType != "none") {
    std::cout << "WARNING: overwriting already defined sampler info for element: " << name << std::endl;
  }

  samplerType   = samplerTypeIn;
  samplerName   = samplerNameIn;
  samplerRadius = samplerRadiusIn;
}
