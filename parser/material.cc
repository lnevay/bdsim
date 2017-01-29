#include "material.h"

#include "array.h"

using namespace GMAD;

Material::Material()
{
  clear();
  PublishMembers();
}

void Material::clear()
{
  A = 0;
  Z = 0;
  density = 0;      //g*cm-3
  temper = 300;     //kelvin
  pressure = 0;     //atm
  state = "solid";  //allowed values: "solid", "liquid", "gas"

  components.clear();
  componentsFractions.clear();
  componentsWeights.clear();
}

void Material::PublishMembers()
{
  publish("name",   &Material::name);
  publish("A",      &Material::A);
  publish("Z",      &Material::Z);
  publish("density",&Material::density);
  publish("T"      ,&Material::temper);
  publish("P"      ,&Material::pressure);
  publish("state",  &Material::state);
  publish("components"         ,&Material::components);
  publish("componentsWeights"  ,&Material::componentsWeights);
  publish("componentsFractions",&Material::componentsFractions);
}

void Material::print()const
{
  std::cout << "material: " << name << " ";

  if (components.empty()) {
    std::cout << "A= " << A << "g/mole "
	      << "Z= " << Z << " ";
  } else {
    std::cout << "ncomponents= " << components.size() << " ";
  }
  std::cout << "density= " << density  << "g/cm3 "
	    << "state= "   << state    << " "
	    << "T= "       << temper   << "K "
	    << "P= "       << pressure << "atm "
	    << std::endl;
}

// template specialisation for Array pointers, to be merged into templated function
void Material::set_value(std::string property, Array* value)
{
#ifdef BDSDEBUG
  std::cout << "parser> Setting value " << std::setw(25) << std::left << property << std::endl;
#endif
  if(property=="components")
    {
      value->set_vector(components);
    } 
  else if(property=="componentsWeights")
    {
      value->set_vector(componentsWeights);
    }
  else if(property=="componentsFractions")
    {
      value->set_vector(componentsFractions);
    }
  else
    {
      std::cerr << "Error: parser> unknown material option \"" << property << "\", or doesn't expect vector type" << std::endl;
      exit(1);
    }
}
