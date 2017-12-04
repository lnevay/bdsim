/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2017.

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
#ifndef BDSGDMLPREPROCESSOR_H
#define BDSGDMLPREPROCESSOR_H

#include "G4String.hh"

#include <vector>
#include <map>

#include "xercesc/dom/DOM.hpp"
#include "xercesc/dom/DOMNodeList.hpp"

namespace BDS
{
  G4String PreprocessGDML(const G4String& file,
			  const G4String& prefix);
}

class BDSGDMLPreprocessor
{
public:
  BDSGDMLPreprocessor();
  ~BDSGDMLPreprocessor();

  G4String PreprocessFile(const G4String& file,
			  const G4String& prefix);
private:

  void ReadNodes(xercesc::DOMNodeList* nodes);
  void ReadNode(xercesc::DOMNode* node);
  void ReadAttributes(xercesc::DOMNamedNodeMap* attributeMap);

  std::vector<std::string> names;
  std::map<std::string, int> count;  

};

#endif
