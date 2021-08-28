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
//  nistManager->FindOrBuildMaterial("G4_Pb");
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
  fNofLayers = 1;
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
  alturapredio = (absoThickness + gapThickness) * (9 + 1);
//  alturapredio = (absoThickness + gapThickness) * (fNofLayers + 1); // gambiarra (usar no primary)
 
  auto layerThickness = absoThickness + gapThickness; // 4m
  auto calorThickness = 9 * layerThickness;   // for 10 layers = 40m
//  auto calorThickness = fNofLayers * layerThickness;
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
  /*                       
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name
 */
  auto layerLV3
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer3");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., calorThickness/2+layerThickness/2), // its position
                  layerLV3,       // its logical volume
                  "Layer3",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
  auto layerLV4
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer4");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)), // its position
                  layerLV4,       // its logical volume
                  "Layer4",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  auto layerLV5
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer5");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*2), // its position
                  layerLV5,       // its logical volume
                  "Layer5",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
  
 auto layerLV6
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer6");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*3), // its position
                  layerLV6,       // its logical volume
                  "Layer6",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
 
  auto layerLV7
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer7");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*4), // its position
                  layerLV7,       // its logical volume
                  "Layer7",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
                 
  auto layerLV8
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer8");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*5), // its position
                  layerLV8,       // its logical volume
                  "Layer8",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
 
  auto layerLV9
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer9");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*6), // its position
                  layerLV9,       // its logical volume
                  "Layer9",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
                 
  auto layerLV10
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer10");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*7), // its position
                  layerLV10,       // its logical volume
                  "Layer10",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
                         
 auto layerLV11
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer11");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*8), // its position
                  layerLV11,       // its logical volume
                  "Layer11",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
                  
  auto layerLV12
    = new G4LogicalVolume(
                   layerS,           // its solid
                   defaultMaterial,  // its material
                   "Layer12");         // its name

    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., (calorThickness/2+layerThickness/2) - (layerThickness)*9), // its position
                  layerLV12,       // its logical volume
                  "Layer12",           // its name
                  worldLV,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps                
                  

/*
  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,        // number of replica
                 layerThickness);  // witdth of replica
*/
  //                               
  // Building concrete
  //
  auto absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeX/2, calorSizeY/2, absoThickness/2); // its size
   /*                      
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
   */          
 auto absorberLV3
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV3");         
  new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV3,       // its logical volume
                  "Abso3",           // its name
                  layerLV3,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
 auto absorberLV4
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV4"); 
  new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV4,       // its logical volume
                  "Abso4",           // its name
                  layerLV4,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
                  
 auto absorberLV5
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV5"); 
  new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV5,       // its logical volume
                  "Abso5",           // its name
                  layerLV5,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
                  
 auto absorberLV6
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV6");   
  new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV6,       // its logical volume
                  "Abso6",           // its name
                  layerLV6,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
                  
  auto absorberLV7
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV7");                 
  new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV7,       // its logical volume
                  "Abso7",           // its name
                  layerLV7,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
                  
 auto absorberLV8
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV8");        
   new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV8,       // its logical volume
                  "Abso8",           // its name
                  layerLV8,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps
                  
  auto absorberLV9
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV9"); 
   new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV9,       // its logical volume
                  "Abso9",           // its name
                  layerLV9,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
auto absorberLV10
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV10");                
   new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV10,       // its logical volume
                  "Abso10",           // its name
                  layerLV10,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
                  
 auto absorberLV11
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV11");                
    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV11,       // its logical volume
                  "Abso11",           // its name
                  layerLV11,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
                  
 auto absorberLV12
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV12");                 
    new G4PVPlacement(
                  0,                // no rotation
                  G4ThreeVector(0., 0., -gapThickness/2), // its position
                  absorberLV12,       // its logical volume
                  "Abso12",           // its name
                  layerLV12,          // its mother  volume
                  false,            // no boolean operation
                  0,                // copy number
                  fCheckOverlaps);  // checking overlaps 
 
         
  //                               
  // Gap
  //
  auto gapS 
    = new G4Box("Gap",             // its name
                 calorSizeX/2, calorSizeY/2, gapThickness/2); // its size
  /*                       
  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name
  */
  auto gapLV3
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV3");         // its name
                               
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), // its position
                 gapLV3,            // its logical volume                         
                 "Gap3",            // its name
                 layerLV3,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
auto gapLV4
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV4");         // its name                                               
  new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV4,            // its logical volume
                   "Gap4",            // its name
                   layerLV4,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps
 auto gapLV5
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV5");         // its name
                                                 
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV5,            // its logical volume
                   "Gap5",            // its name
                   layerLV5,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps  
                               
 auto gapLV6
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV6");         // its name
                                                 
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV6,            // its logical volume
                   "Gap6",            // its name
                   layerLV6,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps  
                   
 auto gapLV7
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV7");         // its name
                                                 
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV7,            // its logical volume
                   "Gap7",            // its name
                   layerLV7,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps  
 
 auto gapLV8
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV8");         // its name
                                               
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV8,            // its logical volume
                   "Gap8",            // its name
                   layerLV8,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps  

auto gapLV9
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV9");         // its name
                               
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV9,            // its logical volume
                   "Gap9",            // its name
                   layerLV9,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps 

auto gapLV10
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV10");         // its name
                                                  
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV10,            // its logical volume
                   "Gap10",            // its name
                   layerLV10,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps 
 
auto gapLV11
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV11");         // its name
                                                  
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV11,            // its logical volume
                   "Gap11",            // its name
                   layerLV11,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps   
 
 
auto gapLV12
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV12");         // its name
                                              
new G4PVPlacement(
                   0,                // no rotation
                   G4ThreeVector(0., 0., absoThickness/2), // its position
                   gapLV12,            // its logical volume
                   "Gap12",            // its name
                   layerLV12,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
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
                 //  gapLV3,          // its mother  volume
                     gapLV3,          // its mother  volume
                   false,            // no boolean operation
                   0,                // copy number
                   fCheckOverlaps);  // checking overlaps
    
    new G4PVPlacement(               // bottom scintillator
                   0,                // no rotation
                   G4ThreeVector(detposX, detposY, scintThickness/2+gapThickness/2-detabsThickness/2-74*cm), // its position
                   scintLV,            // its logical volume
                   "Scint",            // its name
                   gapLV3,          // its mother  volume
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
                      gapLV3,          // its mother  volume
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
 ///   -----------------------------ABSOSD---------------------------------------------
  auto abso3SD
    = new mpdCalorimeterSD("Absorber3SD", "Absorber3HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso3SD);
  SetSensitiveDetector("AbsoLV3",abso3SD);
  auto abso4SD
    = new mpdCalorimeterSD("Absorber4SD", "Absorber4HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso4SD);
  SetSensitiveDetector("AbsoLV4",abso4SD);
  auto abso5SD
    = new mpdCalorimeterSD("Absorber5SD", "Absorber5HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso5SD);
  SetSensitiveDetector("AbsoLV5",abso5SD);
  auto abso6SD
    = new mpdCalorimeterSD("Absorber6SD", "Absorber6HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso6SD);
  SetSensitiveDetector("AbsoLV6",abso6SD);
  auto abso7SD
    = new mpdCalorimeterSD("Absorber7SD", "Absorber7HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso7SD);
  SetSensitiveDetector("AbsoLV7",abso7SD);
  auto abso8SD
    = new mpdCalorimeterSD("Absorber8SD", "Absorber8HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso8SD);
  SetSensitiveDetector("AbsoLV8",abso8SD);
  auto abso9SD
    = new mpdCalorimeterSD("Absorber9SD", "Absorber9HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso9SD);
  SetSensitiveDetector("AbsoLV9",abso9SD);
  auto abso10SD
    = new mpdCalorimeterSD("Absorber10SD", "Absorber10HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso10SD);
  SetSensitiveDetector("AbsoLV10",abso10SD);
  auto abso11SD
    = new mpdCalorimeterSD("Absorber11SD", "Absorber11HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso11SD);
  SetSensitiveDetector("AbsoLV11",abso11SD);
  auto abso12SD
    = new mpdCalorimeterSD("Absorber12SD", "Absorber12HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(abso12SD);
  SetSensitiveDetector("AbsoLV12",abso12SD);


 ///   -----------------------------GAPSD---------------------------------------------
  /*
  auto gapSD 
    = new mpdCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);
  */
  
  auto gap3SD 
    = new mpdCalorimeterSD("Gap3SD", "Gap3HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap3SD);
  SetSensitiveDetector("GapLV3",gap3SD);  // gap of detector layer

  auto gap4SD 
    = new mpdCalorimeterSD("Gap4SD", "Gap4HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap4SD);
  SetSensitiveDetector("GapLV4",gap4SD);  
  
 auto gap5SD 
    = new mpdCalorimeterSD("Gap5SD", "Gap5HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap5SD);
  SetSensitiveDetector("GapLV5",gap5SD);  

 auto gap6SD 
    = new mpdCalorimeterSD("Gap6SD", "Gap6HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap6SD);
  SetSensitiveDetector("GapLV6",gap6SD);  

 auto gap7SD 
    = new mpdCalorimeterSD("Gap7SD", "Gap7HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap7SD);
  SetSensitiveDetector("GapLV7",gap7SD);  

 auto gap8SD 
    = new mpdCalorimeterSD("Gap8SD", "Gap8HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap8SD);
  SetSensitiveDetector("GapLV8",gap8SD);  


 auto gap9SD 
    = new mpdCalorimeterSD("Gap9SD", "Gap9HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap9SD);
  SetSensitiveDetector("GapLV9",gap9SD);  

 auto gap10SD 
    = new mpdCalorimeterSD("Gap10SD", "Gap10HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap10SD);
  SetSensitiveDetector("GapLV10",gap10SD);  

 auto gap11SD 
    = new mpdCalorimeterSD("Gap11SD", "Gap11HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap11SD);
  SetSensitiveDetector("GapLV11",gap11SD);  
  
 auto gap12SD 
    = new mpdCalorimeterSD("Gap12SD", "Gap12HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gap12SD);
  SetSensitiveDetector("GapLV12",gap12SD);  

 ///   -----------------------------LAYERSD---------------------------------------------
 
  auto layer3SD 
    = new mpdCalorimeterSD("Layer3SD", "layer3HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer3SD);
  SetSensitiveDetector("Layer3",layer3SD);

   auto layer4SD 
    = new mpdCalorimeterSD("Layer4SD", "layer4HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer4SD);
  SetSensitiveDetector("Layer4",layer4SD);
  
   auto layer5SD 
    = new mpdCalorimeterSD("Layer5SD", "layer5HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer5SD);
  SetSensitiveDetector("Layer5",layer5SD);
  
   auto layer6SD 
    = new mpdCalorimeterSD("Layer6SD", "layer6HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer6SD);
  SetSensitiveDetector("Layer6",layer6SD);
  
   auto layer7SD 
    = new mpdCalorimeterSD("Layer7SD", "layer7HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer7SD);
  SetSensitiveDetector("Layer7",layer7SD);
  
   auto layer8SD 
    = new mpdCalorimeterSD("Layer8SD", "layer8HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer8SD);
  SetSensitiveDetector("Layer8",layer8SD);
  
   auto layer9SD 
    = new mpdCalorimeterSD("Layer9SD", "layer9HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer9SD);
  SetSensitiveDetector("Layer9",layer9SD);
  
   auto layer10SD 
    = new mpdCalorimeterSD("Layer10SD", "layer10HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer10SD);
  SetSensitiveDetector("Layer10",layer10SD);
  
   auto layer11SD 
    = new mpdCalorimeterSD("Layer11SD", "layer11HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer11SD);
  SetSensitiveDetector("Layer11",layer11SD);
  
   auto layer12SD 
    = new mpdCalorimeterSD("Layer12SD", "layer12HitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(layer12SD);
  SetSensitiveDetector("Layer12",layer12SD);
  
   ///   -----------------------------SCINTSD---------------------------------------------
  
  auto scintSD
    = new mpdCalorimeterSD("ScintSD", "ScintHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintSD);
  SetSensitiveDetector("ScintLV",scintSD);

 ///   -----------------------------DETABSSD---------------------------------------------

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
