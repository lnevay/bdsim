#include "BDSFieldType.hh"
#include "BDSDebug.hh"
#include "globals.hh"

#include <map>
#include <string>

// dictionary for BDSFieldType for reflexivity
template<>
std::map<BDSFieldType, std::string>* BDSFieldType::dictionary =
  new std::map<BDSFieldType, std::string> ({
      {BDSFieldType::none,             "none"},
      {BDSFieldType::bzero,            "bzero"},
      {BDSFieldType::ezero,            "ezero"},
      {BDSFieldType::ebzero,           "ebzero"},
      {BDSFieldType::teleporter,       "teleporter"},
      {BDSFieldType::bmap1d,           "bmap1d"},
      {BDSFieldType::bmap2d,           "bmap2d"},
      {BDSFieldType::bmap3d,           "bmap3d"},
      {BDSFieldType::bmap4d,           "bmap4d"},
      {BDSFieldType::emap1d,           "emap1d"},
      {BDSFieldType::emap2d,           "emap2d"},
      {BDSFieldType::emap3d,           "emap3d"},
      {BDSFieldType::emap4d,           "emap4d"},
      {BDSFieldType::ebmap1d,          "ebmap1d"},
      {BDSFieldType::ebmap2d,          "ebmap2d"},
      {BDSFieldType::ebmap3d,          "ebmap3d"},
      {BDSFieldType::ebmap4d,          "ebmap4d"},
      {BDSFieldType::mokka,            "mokka"},
      {BDSFieldType::solenoid,         "solenoid"},
      {BDSFieldType::dipole,           "dipole"},
      {BDSFieldType::quadrupole,       "quadrupole"},
      {BDSFieldType::dipolequadrupole, "dipolequadrupole"},
      {BDSFieldType::sextupole,        "sextupole"},
      {BDSFieldType::octupole,         "octupole"},
      {BDSFieldType::decapole,         "decapole"},
      {BDSFieldType::multipole,        "multipole"},
      {BDSFieldType::muonspoiler,      "muonspoiler"},
      {BDSFieldType::skewquadrupole,   "skewquadrupole"},
      {BDSFieldType::skewsextupole,    "skewsextupole"},
      {BDSFieldType::skewoctupole,     "skewoctupole"},
      {BDSFieldType::skewdecapole,     "skewdecapole"},
      {BDSFieldType::rfcavity,         "rfcavity"},
      {BDSFieldType::rf,               "rf"},
      {BDSFieldType::dipole3d,         "dipole3d"}
});	

BDSFieldType BDS::DetermineFieldType(G4String bType)
{
  std::map<G4String, BDSFieldType> types;
  types["none"]             = BDSFieldType::none;
  types["bzero"]            = BDSFieldType::bzero;
  types["ezero"]            = BDSFieldType::ezero;
  types["ebzero"]           = BDSFieldType::ebzero;
  types["teleporter"]       = BDSFieldType::teleporter;
  types["bmap1d"]           = BDSFieldType::bmap1d;
  types["bmap2d"]           = BDSFieldType::bmap2d;
  types["bmap3d"]           = BDSFieldType::bmap3d;
  types["bmap4d"]           = BDSFieldType::bmap4d;
  types["emap1d"]           = BDSFieldType::emap1d;
  types["emap2d"]           = BDSFieldType::emap2d;
  types["emap3d"]           = BDSFieldType::emap3d;
  types["emap4d"]           = BDSFieldType::emap4d;
  types["ebmap1d"]          = BDSFieldType::ebmap1d;
  types["ebmap2d"]          = BDSFieldType::ebmap2d;
  types["ebmap3d"]          = BDSFieldType::ebmap3d;
  types["ebmap4d"]          = BDSFieldType::ebmap4d;
  types["mokka"]            = BDSFieldType::mokka;
  types["solenoid"]         = BDSFieldType::solenoid;
  types["dipole"]           = BDSFieldType::dipole;
  types["quadrupole"]       = BDSFieldType::quadrupole;
  types["dipolequadrupole"] = BDSFieldType::dipolequadrupole;
  types["sextupole"]        = BDSFieldType::sextupole;
  types["octupole"]         = BDSFieldType::octupole;
  types["decapole"]         = BDSFieldType::decapole;
  types["multipole"]        = BDSFieldType::multipole;
  types["muonspoiler"]      = BDSFieldType::muonspoiler;
  types["skewquadrupole"]   = BDSFieldType::skewquadrupole;
  types["skewsextupole"]    = BDSFieldType::skewsextupole;
  types["skewoctupole"]     = BDSFieldType::skewoctupole;
  types["skewdecapole"]     = BDSFieldType::skewdecapole;
  types["rfcavity"]         = BDSFieldType::rfcavity;
  types["rf"]               = BDSFieldType::rf;
  types["dipole3d"]         = BDSFieldType::dipole3d;
  
  bType.toLower();

  auto result = types.find(bType);
  if (result == types.end())
    {
      // it's not a valid key
      G4cerr << __METHOD_NAME__ << bType << " is not a valid field type" << G4endl;

      G4cout << "Available field types are:" << G4endl;
      for (auto it : types)
	{G4cout << "\"" << it.first << "\"" << G4endl;}
      exit(1);
    }

#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "determined field type to be " << result->second << G4endl;
#endif
  return result->second;
}
