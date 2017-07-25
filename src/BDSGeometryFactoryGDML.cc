#ifdef USE_GDML
#include "BDSAcceleratorModel.hh"
#include "BDSDebug.hh"
#include "BDSGeometryExternal.hh"
#include "BDSGeometryFactoryGDML.hh"
#include "BDSGeometryInspector.hh"
#include "BDSGlobalConstants.hh"

#include "globals.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

class G4VSolid;

BDSGeometryFactoryGDML* BDSGeometryFactoryGDML::instance = nullptr;

BDSGeometryFactoryGDML::BDSGeometryFactoryGDML()
{;}

BDSGeometryFactoryGDML::~BDSGeometryFactoryGDML()
{
  instance = nullptr;
}

BDSGeometryFactoryGDML* BDSGeometryFactoryGDML::Instance()
{
  if (!instance)
    {instance = new BDSGeometryFactoryGDML();}
  return instance;
}

BDSGeometryExternal* BDSGeometryFactoryGDML::Build(G4String componentName,
						   G4String fileName,
						   std::map<G4String, G4Colour*>* mapping,
						   G4double /*suggestedLength*/,
						   G4double /*suggestedOuterDiameter*/)
{
  CleanUp();

  G4String tempFileName = ReplaceStringInFile(componentName, fileName,
					      "PREPEND", componentName);
  BDSAcceleratorModel::Instance()->RegisterTemporaryFile(tempFileName);
  
  G4GDMLParser* parser = new G4GDMLParser();
  parser->Read(tempFileName, /*validate=*/true);

  G4VPhysicalVolume* containerPV = parser->GetWorldVolume();
  G4LogicalVolume*   containerLV = containerPV->GetLogicalVolume();
  G4VSolid*       containerSolid = containerLV->GetSolid();
  G4ThreeVector gdmlWorldOrigin(0,0,0);
  if (containerPV->GetName() == "world_volume_lv_PV")
    {
      gdmlWorldOrigin = parser->GetPosition("PygdmlOrigin"); // TODO check if Pygdml geometry
      gdmlWorldOrigin[2] = 0.0;
    }

  // record all pvs and lvs used in this loaded geometry
  std::vector<G4VPhysicalVolume*> pvsGDML;
  std::vector<G4LogicalVolume*>   lvsGDML;
  GetAllLogicalAndPhysical(containerPV, pvsGDML, lvsGDML);

  auto visesGDML = ApplyColourMapping(lvsGDML, mapping);

  /// Now overwrite container lv vis attributes
  containerLV->SetVisAttributes(BDSGlobalConstants::Instance()->GetContainerVisAttr());

  std::pair<BDSExtent, BDSExtent> outerInner = BDS::DetermineExtents(containerSolid);
  
  BDSGeometryExternal* result = new BDSGeometryExternal(containerSolid, containerLV,
                                                        outerInner.first,  /*outer*/
                                                        outerInner.second, /*inner*/
                                                        gdmlWorldOrigin);
  result->RegisterLogicalVolume(lvsGDML);
  result->RegisterPhysicalVolume(pvsGDML);
  result->RegisterVisAttributes(visesGDML);

  delete parser;
  return result;
}

void BDSGeometryFactoryGDML::GetAllLogicalAndPhysical(const G4VPhysicalVolume*         volume,
						      std::vector<G4VPhysicalVolume*>& pvsIn,
						      std::vector<G4LogicalVolume*>&   lvsIn)
{
  const auto& lv = volume->GetLogicalVolume();
  lvs.push_back(lv);
  for (G4int i = 0; i < lv->GetNoDaughters(); i++)
    {
      const auto& pv = lv->GetDaughter(i);
      pvsIn.push_back(pv);
      lvsIn.push_back(lv);
      GetAllLogicalAndPhysical(pv, pvsIn, lvsIn); // recurse into daughter
    }
}


G4String BDSGeometryFactoryGDML::ReplaceStringInFile(G4String componentName,
						     G4String fileName,
						     G4String key,
						     G4String replacement)
{
  // open input file in read mode
  std::ifstream ifs(fileName);

  // verify file open.
  if (!ifs.is_open())
    {
      G4cerr << __METHOD_NAME__ << "Cannot open file \"" << fileName << "\"" << G4endl;
      exit(1);
    }
  
  // temporary file name
  std::string tempFileName = "temp_" + componentName + ".gdml";
  //std::string tempFileName = std::string(std::tmpnam(nullptr)) + "_.gdml";
  std::ofstream fout(tempFileName);

#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "Original file:  " << fileName     << G4endl;
  G4cout << __METHOD_NAME__ << "Temporary file: " << tempFileName << G4endl;
#endif

  int lenOfKey = key.size();
  
  // loop over and replace
  std::string buffer;
  while (std::getline(ifs, buffer))
    {// if we find key, replace it
      int f = buffer.find(key);    
      if (f != -1)
	{
	  std::string outputString = std::string(buffer);
	  outputString.replace(f, lenOfKey, replacement);
	  fout << outputString << "\n"; // getline strips \n
	}
      else // copy line to temp file as is
	{fout << buffer << "\n";}
    }

  // clean up
  ifs.close();
  fout.close();
  
  return G4String(tempFileName);
}


#endif
