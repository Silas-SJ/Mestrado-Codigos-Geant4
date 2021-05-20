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
// 
/// \file mpdDetectorConstruction.cc
/// \brief Implementation of the mpdDetectorConstruction class

#include "mpdDetectorConstruction.hh"
#include "mpdCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* mpdDetectorConstruction::fMagFieldMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdDetectorConstruction::mpdDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true),
   fNofLayers(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdDetectorConstruction::~mpdDetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* mpdDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_CONCRETE");
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* mpdDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 10;
 // G4double absoThickness = 10.*mm;
 // G4double gapThickness =  5.*mm;
 // G4double calorSizeXY  = 10.*cm;

  G4double absoThickness = 500.*mm;
  G4double gapThickness =  3500.*mm;
  calorSizeX  = 230.*m;
  calorSizeY  = 25.*m;
 
  G4double scintThickness = 20.*mm;
  G4double scintSizeXY  = 40.*cm;
    
  G4double detabsThickness = 5.*mm;
  G4double detabsSizeXY  = 40.*cm;
  alturapredio = (absoThickness + gapThickness) * (fNofLayers + 1); // gambiarra (usar no primary)
 
  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = fNofLayers * layerThickness;
//  auto worldSizeX = 1.2 * calorSizeX;
//  auto worldSizeY = 1.2 * calorSizeY;
  auto worldSizeX = 500*m;
  auto worldSizeY = 200*m;
  auto worldSizeZ  = 1.2 * calorThickness; 
  
  // Get materials
//  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
//  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto absorberMaterial = G4Material::GetMaterial("G4_CONCRETE");
//  auto gapMaterial = G4Material::GetMaterial("liquidArgon");
  auto gapMaterial = G4Material::GetMaterial("G4_AIR");
  auto scintMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto detabsMaterial = G4Material::GetMaterial("G4_Pb");
  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("mpdDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeX/2, calorSizeY/2, calorThickness/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
  //                                 
  // Layer
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeX/2, calorSizeY/2, layerThickness/2); //its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name
  auto layerLV2
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer2");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., calorThickness/2+layerThickness/2), // its position
                  layerLV2,       // its logical volume
                  "Layer2",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,        // number of replica
                 layerThickness);  // witdth of replica
  
  //                               
  // Building concrete
  //
  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeX/2, calorSizeY/2, absoThickness/2); // its size
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2), // its position
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV,       // its logical volume
                  "Abso",           // its name
                  layerLV2,          // its mother  volume
                  false,            // no boolean operation
                  1,                // copy number
                  fCheckOverlaps);  // checking overlaps
  /*
   //New               
  auto absorberS2 
    = new G4Box("Abso2",            // its name
                 calorSizeX/2, calorSizeY/2, absoThickness/2); // its size
                         
  auto absorberLV2
    = new G4LogicalVolume(
                 absorberS2,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV2");        // its name
                                   
   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2), // its position
                 absorberLV2,       // its logical volume                         
                 "Abso2",           // its name
                 layerLV2,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV2,       // its logical volume
                  "Abso2",           // its name
                  layerLV2,          // its mother  volume
                  false,            // no boolean operation
                  1,                // copy number
                  fCheckOverlaps);  // checking overlaps
  */
  //                               
  // Gap
  //
  auto gapS 
    = new G4Box("Gap",             // its name
                 calorSizeX/2, calorSizeY/2, gapThickness/2); // its size
                         
  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name
  auto gapLV2
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV2");         // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV2,            // its logical volume
                   "Gap2",            // its name
                   layerLV2,          // its mother  volume
                   false,            // no boolean operation
                   1,                // copy number
                   fCheckOverlaps);  // checking overlaps
    //
    // Scintillator
    //
    auto scintS
      = new G4Box("Scint",             // its name
                   scintSizeXY/2, scintSizeXY/2, scintThickness/2); // its size
                           
    auto scintLV
      = new G4LogicalVolume(
                   scintS,             // its solid
                   scintMaterial,      // its material
                   "ScintLV");         // its name
    //G4double detposX = 5*m-calorSizeX/2;
   // G4double detposY = 2*m-calorSizeY/2;
     detposX = 5*m-calorSizeX/2;
     detposY = 2*m-calorSizeY/2;
     detposZ = scintThickness/2+gapThickness/2-70*cm;
    new G4PVPlacement(                // top scintillator
                   0,                // no rotation
                   //G4ThreeVector(detposX, detposY, scintThickness/2+gapThickness/2-70*cm), // its position
                   G4ThreeVector(detposX, detposY, detposZ),
                   scintLV,            // its logical volume
                   "Scint",            // its name
                   gapLV2,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps
    
    new G4PVPlacement(               // bottom scintillator
                   0,                // no rotation
                   G4ThreeVector(detposX, detposY, scintThickness/2+gapThickness/2-detabsThickness/2-74*cm), // its position
                   scintLV,            // its logical volume
                   "Scint",            // its name
                   gapLV2,          // its mother  volume
                   false,            // no boolean operation
                   1,                // copy number
                   fCheckOverlaps);  // checking overlaps

    
       //
       // Detector Absorber
       //
       auto detabsS
         = new G4Box("Detabs",             // its name
                      detabsSizeXY/2, detabsSizeXY/2, detabsThickness/2); // its size
                              
       auto detabsLV
         = new G4LogicalVolume(
                      detabsS,             // its solid
                      detabsMaterial,      // its material
                      "DetabsLV");         // its name
                                        
       new G4PVPlacement(
                      0,                // no rotation
                      G4ThreeVector(detposX, detposY, scintThickness/2+gapThickness/2-detabsThickness/2-72*cm), // its position
                      detabsLV,            // its logical volume
                      "Detabs",            // its name
                      gapLV2,          // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      fCheckOverlaps);  // checking overlaps
  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdDetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //

  auto absoSD
    = new mpdCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto gapSD 
    = new mpdCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);

  auto scintSD
    = new mpdCalorimeterSD("ScintSD", "ScintHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintSD);
  SetSensitiveDetector("ScintLV",scintSD);

  auto detabsSD
    = new mpdCalorimeterSD("DetabsSD", "DetabsHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(detabsSD);
  SetSensitiveDetector("DetabsLV",detabsSD);
    
  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
