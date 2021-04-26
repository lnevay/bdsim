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

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "BDSDicomFileMgr.hh"

#include "BDSDicomFileCT.hh"

#include "dcmtk/dcmdata/dcdeftag.h"
#include "G4tgrFileIn.hh"
#include "G4UIcommand.hh"

BDSDicomFileMgr* BDSDicomFileMgr::theInstance = nullptr;
int BDSDicomFileMgr::verbose = 1;

BDSDicomFileMgr* BDSDicomFileMgr::GetInstance()
{
    if( !theInstance ) {
        theInstance = new BDSDicomFileMgr;
    }
    return theInstance;
}

BDSDicomFileMgr::BDSDicomFileMgr()
{
    fCompression = 1.;
    theCTFileAll = nullptr;
    theStructureNCheck = 4;
    theStructureNMaxROI = 100;
}

void BDSDicomFileMgr::Convert(G4String filePath, G4String fileName)
{
    G4tgrFileIn fin = G4tgrFileIn::GetInstance(filePath + fileName);
    std::vector<G4String> wl;
    // Read each file in file list
    theFileOutName = "test.g4dcm";
    int ii;
    for( ii = 0;; ii++) {
        if( ! fin.GetWordsInLine(wl) ) break;
        if( wl[0] == ":COMPRESSION" ) {
            CheckNColumns(wl,2);
            SetCompression(wl[1]);
        } else if( wl[0] == ":FILE" ) {
            CheckNColumns(wl,2);
            G4cout << "@@@@@@@ Reading FILE: " << wl[1]  << G4endl;
            AddFile(filePath + wl[1]);
        } else if( wl[0] == ":FILE_OUT" ) {
            CheckNColumns(wl,2);
            theFileOutName = wl[1];
        } else if( wl[0] == ":MATE_DENS" ) {
            CheckNColumns(wl,3);
            AddMaterialDensity(wl);
        } else if( wl[0] == ":MATE" ) {
            CheckNColumns(wl,3);
            AddMaterial(wl);
        } else if( wl[0] == ":CT2D" ) {
            CheckNColumns(wl,3);
            AddCT2Density(wl);
        } else {
            G4Exception("DICOM2G4",
                        "Wrong argument",
                        FatalErrorInArgument,
                        G4String("UNKNOWN TAG IN FILE  "+wl[0]).c_str());
        }

    }
    //@@@@@@ Process files
    ProcessFiles();
}

void BDSDicomFileMgr::CheckNColumns(std::vector<G4String> wl, size_t vsizeTh )
{
    if( wl.size() != vsizeTh ) {
        G4cerr << " Reading line " << G4endl;
        for( size_t ii = 0; ii < wl.size(); ii++){
            G4cerr << wl[ii] << " ";
        }
        G4cerr << G4endl;
        G4Exception("DICOM2G4",
                    "D2G0010",
                    FatalErrorInArgument,
                    ("Wrong number of columns in line " + std::to_string(wl.size()) + " <> "
                     + std::to_string(vsizeTh)).c_str());
    }
}

void BDSDicomFileMgr::SetCompression( G4String fComp )
{
    fCompression = G4UIcommand::ConvertToDouble(fComp);
}

void BDSDicomFileMgr::AddFile( G4String fileName )
{
    DcmFileFormat dfile;
    if( ! (dfile.loadFile(fileName.c_str())).good() ) {
        G4Exception("DicomHandler::ReadFile",
                    "dfile.loadFile",
                    FatalErrorInArgument,
                    ("Error reading file " + fileName).c_str());
    }
    DcmDataset* dset = dfile.getDataset();

    OFString dModality;
    if( !dset->findAndGetOFString(DCM_Modality,dModality).good() ) {
        G4Exception("DicomHandler::ReadData ",
                    "",
                    FatalException,
                    " Have not read Modality");
    }

    if( dModality == "CT"  || dModality == "OT") {
        BDSDicomFileCT* df = new BDSDicomFileCT(dset);
        df->ReadData();
        df->SetFileName( fileName );
        // reorder by location
        theCTFiles[df->GetMaxZ()] = df;
        G4cout << "Number of voxels: " << df->GetNoVoxels() << G4endl;
    } else if( dModality == "RTSTRUCT" ) {
        BDSDicomFileStructure* df = new BDSDicomFileStructure(dset);
        df->ReadData();
        df->SetFileName( fileName );
        //    theFiles.insert(msd::value_type(dModality,df));
        theStructFiles.push_back(df);
    } else if( dModality == "RTPLAN" ) {
        BDSDicomFilePlan *df = new BDSDicomFilePlan(dset);
        df->ReadData();
        df->SetFileName(fileName);
        //    theFiles.insert(msd::value_type(dModality,df));
        thePlanFiles.push_back(df);
    } else {
        G4Exception("DicomFileMgr::AddFIle()",
                    "DFM001",
                    FatalErrorInArgument,
                    (G4String("File is not of type CT or RTSTRUCT or RTPLAN, but: ")
                     + dModality).c_str());
    }

}

void BDSDicomFileMgr::AddMaterial( std::vector<G4String> wl )
{
    if( theMaterials.size() > 0 && bMaterialsDensity ) {
        G4Exception("DicomFileMgr::AddMaterial",
                    "DFM005",
                    FatalException,
                    "Trying to add a Material with :MATE and another with :MATE_DENS, check your input file");
    }
    bMaterialsDensity = false;
    theMaterials[G4UIcommand::ConvertToDouble(wl[2])] = wl[1];
}

void BDSDicomFileMgr::AddMaterialDensity( std::vector<G4String> wl )
{
    if( theMaterialsDensity.size() > 0 && !bMaterialsDensity ) {
        G4Exception("DicomFileMgr::AddMaterial",
                    "DFM005",
                    FatalException,
                    "Trying to add a Material with :MATE and another with :MATE_DENS, check your input file");
    }
    bMaterialsDensity = true;
    theMaterialsDensity[G4UIcommand::ConvertToDouble(wl[2])] = wl[1];
}

void BDSDicomFileMgr::AddCT2Density( std::vector<G4String> wl)
{
    theCT2Density[G4UIcommand::ConvertToInt(wl[1])] = G4UIcommand::ConvertToDouble(wl[2]);
    G4cout << this << " AddCT2density " << theCT2Density.size() << G4endl;//GDEB

}

G4double BDSDicomFileMgr::Hounsfield2density(Uint32 Hval)
{
    if( theCT2Density.size() == 0 ) {
        G4Exception("Hounsfield2density",
                    "DCM004",
                    FatalException,
                    "No :CT2D line in input file");
    }
    std::map<G4int,G4double>::const_iterator ite = theCT2Density.begin();
    G4int minHval = (*ite).first;
    if( G4int(Hval) < minHval ) {
        G4Exception("DicomHandler::Hounsfield2density",
                    "",
                    FatalException,
                    ("Hval value too small, change input file "+std::to_string(Hval) + " < "
                     + std::to_string(minHval)).c_str());
    }

    ite = theCT2Density.end(); ite--;
    G4int maxHval = (*ite).first;
    if( G4int(Hval) > maxHval ) {
        G4Exception("DicomHandler::Hval2density",
                    "",
                    FatalException,
                    ("Hval value too big, change CT2Density.dat file "+std::to_string(Hval) + " > "
                     + std::to_string(maxHval)).c_str());
    }

    G4float density = -1.;
    G4double deltaCT = 0;
    G4double deltaDensity = 0;

    ite = theCT2Density.upper_bound(Hval);
    std::map<G4int,G4double>::const_iterator itePrev = ite; itePrev--;

    deltaCT = (*ite).first - (*itePrev).first;
    deltaDensity = (*ite).second - (*itePrev).second;

    // interpolating linearly
    density = (*ite).second - (((*ite).first-Hval)*deltaDensity/deltaCT );

    if(density < 0.) {
        G4Exception("DicomFileMgr::Hounsfiled2Density",
                    "DFM002",
                    FatalException,
                    G4String("@@@ Error negative density = " + std::to_string(density) + " from HV = "
                             + std::to_string(Hval)).c_str());
    }

    //  G4cout << " Hval2density " << Hval << " -> " << density << G4endl;//GDEB
    return density;

}

size_t BDSDicomFileMgr::GetMaterialIndex( G4double Hval)
{
    std::map<G4double,G4String>::iterator ite = theMaterials.upper_bound(Hval);
    if( ite == theMaterials.end() ) {
        ite--;
        G4Exception("DicomFileMgr::GetMaterialIndex",
                    "DFM004",
                    FatalException,
                    ("Hounsfiled value too big, change input file "+std::to_string(Hval) + " > "
                     + std::to_string((*ite).first)).c_str());
    }

    size_t dist = std::distance( theMaterials.begin(), ite );

    return  dist;

}

size_t BDSDicomFileMgr::GetMaterialIndexByDensity( G4double density )
{
    std::map<G4double,G4String>::iterator ite = theMaterialsDensity.upper_bound(density);
    if( ite == theMaterialsDensity.end() ) {
        ite--;
        G4Exception("DicomFileMgr::GetMaterialIndexByDensity",
                    "DFM003",
                    FatalException,
                    ("Density too big, change input file "+std::to_string(density) + " > "
                     + std::to_string((*ite).first)).c_str());
    }

    size_t dist = std::distance( theMaterialsDensity.begin(), ite );

    return  dist;

}

void BDSDicomFileMgr::ProcessFiles()
{
    if( theCTFiles.size() == 0 ) {
        G4Exception("CheckCTSlices",
                    "DCM004",
                    JustWarning,
                    "No :FILE of type CT in input file");
    } else {
        CheckCTSlices();
        BuildCTMaterials();
        MergeCTFiles();
    }
    DumpToTextFile();
}

void BDSDicomFileMgr::CheckCTSlices()
{
    size_t nSlices = theCTFiles.size();
    G4cout << " DicomFileMgr::Checking CT slices: " << nSlices << G4endl;

    G4bool uniformSliceThickness = true;

    if(nSlices > 1) {
        if(nSlices == 2) {
            mdct::const_iterator ite = theCTFiles.begin();
            BDSDicomFileCT* one = (*ite).second;
            ite++;
            BDSDicomFileCT* two = (*ite).second;

            G4double real_distance = (two->GetLocation()-one->GetLocation())/2.;

            if(one->GetMaxZ() != two->GetMinZ()) {
                one->SetMaxZ(one->GetLocation()+real_distance);
                two->SetMinZ(two->GetLocation()-real_distance);
                //one->SetMinZ(one->GetLocation()-real_distance);
                //two->SetMaxZ(two->GetLocation()+real_distance);
                if(uniformSliceThickness) {
                    one->SetMinZ(one->GetLocation()-real_distance);
                    two->SetMaxZ(two->GetLocation()+real_distance);
                }
            }
        } else {
            mdct::iterator ite0 = theCTFiles.begin();
            mdct::iterator ite1 = ite0; ite1++;
            mdct::iterator ite2 = ite1; ite2++;
            for(; ite2 != theCTFiles.end(); ++ite0, ++ite1, ++ite2) {
                BDSDicomFileCT* prev = (BDSDicomFileCT*)((*ite0).second);
                BDSDicomFileCT* slice = (BDSDicomFileCT*)((*ite1).second);
                BDSDicomFileCT* next = (BDSDicomFileCT*)((*ite2).second);
                G4double real_up_distance = (next->GetLocation() -
                                             slice->GetLocation())/2.;
                G4double real_down_distance = (slice->GetLocation() -
                                               prev->GetLocation())/2.;
                G4double real_distance = real_up_distance + real_down_distance;
                G4double stated_distance = slice->GetMaxZ()-slice->GetMinZ();

                if(std::fabs(real_distance - stated_distance) > 1.E-9) {
                    unsigned int sliceNum = std::distance(theCTFiles.begin(),ite1);
                    G4cerr << "\tDicomFileMgr::CheckCTSlices - Slice Distance Error in slice [" << sliceNum
                           << "]: Distance between this slice and slices up and down = "
                           << real_distance
                           << " <> Slice width = " << stated_distance
                           << " Slice locations " <<prev->GetLocation() << " : " << slice->GetLocation()
                           << " : " << next->GetLocation()
                           << " DIFFERENCE= " << real_distance - stated_distance
                           << G4endl;
                    G4cerr << "!! WARNING: Geant4 will reset slice width so that it extends between "
                           << "lower and upper slice " << G4endl;

                    slice->SetMinZ(slice->GetLocation()-real_down_distance);
                    slice->SetMaxZ(slice->GetLocation()+real_up_distance);

                    if(ite0 == theCTFiles.begin()) {
                        prev->SetMaxZ(slice->GetMinZ());
                        // Using below would make all slice same thickness
                        //prev->SetMinZ(prev->GetLocation()-real_min_distance);
                        if(uniformSliceThickness) {
                            prev->SetMinZ(prev->GetLocation()-real_down_distance);
                        }
                    }
                    if(static_cast<unsigned int>(std::distance(theCTFiles.begin(),ite2)+1)==
                       nSlices) {
                        next->SetMinZ(slice->GetMaxZ());
                        // Using below would make all slice same thickness
                        //next->SetMaxZ(next->GetLocation()+real_max_distance);
                        if(uniformSliceThickness) {
                            next->SetMaxZ(next->GetLocation()+real_up_distance); }
                    }
                }
            }
        }
    }
}

void BDSDicomFileMgr::BuildCTMaterials()
{
    G4cout << " DicomFileMgr::Building Materials " << theCTFiles.size() << G4endl;//GDEB
    mdct::const_iterator ite = theCTFiles.begin();
    for( ; ite != theCTFiles.end(); ite++ ) {
        (*ite).second->BuildMaterials();
    }
}

void BDSDicomFileMgr::MergeCTFiles()
{
    G4cout << " DicomFileMgr::Merging CT Files " << theCTFiles.size() << G4endl;//GDEB
    mdct::const_iterator ite = theCTFiles.begin();
    theCTFileAll = new BDSDicomFileCT( *((*ite).second) );
    ite++;
    for( ; ite != theCTFiles.end(); ite++ ) {
        (*theCTFileAll) +=  *((*ite).second);
    }
}

void BDSDicomFileMgr::DumpToTextFile()
{
    G4cout << " DicomFileMgr::Dumping To Text File " << G4endl; //GDEB
    if( theCTFiles.size() != 0 ) {
        std::ofstream fout(theFileOutName);

        if( !bMaterialsDensity ) {
            fout << theMaterials.size() << std::endl;
            std::map<G4double,G4String>::const_iterator ite;
            G4int ii = 0;
            for(ite = theMaterials.begin(); ite != theMaterials.end(); ite++, ii++){
                fout << ii << " \"" << (*ite).second << "\"" << std::endl;
            }
        } else {
            fout << theMaterialsDensity.size() << std::endl;
            std::map<G4double,G4String>::const_iterator ite;
            G4int ii = 0;
            for(ite = theMaterialsDensity.begin(); ite != theMaterialsDensity.end(); ite++, ii++){
                fout << ii << " \"" << (*ite).second << "\"" << std::endl;
            }
        }

        theCTFileAll->DumpHeaderToTextFile(fout);
        for( mdct::const_iterator itect = theCTFiles.begin(); itect != theCTFiles.end(); itect++ ) {
            (*itect).second->DumpMateIDsToTextFile(fout);
        }
        for( mdct::const_iterator itect = theCTFiles.begin(); itect != theCTFiles.end(); itect++ ) {
            (*itect).second->DumpDensitiesToTextFile(fout);
        }
        for( mdct::const_iterator itect = theCTFiles.begin(); itect != theCTFiles.end(); itect++ ) {
            (*itect).second->BuildStructureIDs();
            (*itect).second->DumpStructureIDsToTextFile(fout);
        }

        std::vector<BDSDicomFileStructure*> dfs = GetStructFiles();
        for( size_t i1 = 0; i1 < dfs.size(); i1++ ){
            std::vector<BDSDicomROI*> rois = dfs[i1]->GetROIs();
            for( size_t i2 = 0; i2 < rois.size(); i2++ ){
                fout << rois[i2]->GetNumber()+1 << " \"" << rois[i2]->GetName() << "\"" <<G4endl;
            }
        }
    }

    for( size_t i1 = 0; i1 < thePlanFiles.size(); i1++ ){
        thePlanFiles[i1]->DumpToFile();
    }
}
