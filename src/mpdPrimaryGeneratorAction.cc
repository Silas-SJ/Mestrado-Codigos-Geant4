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
/*
 
 spec-pion.dat
 TF1 * f1 = new TF1("f1","[0]*exp([1]*log([2]*x))",0,4000);
 ****************************************
 Minimizer is Minuit / Migrad
 Chi2                      =  2.80559e-18
 NDf                       =           17
 Edm                       =  2.12624e-26
 NCalls                    =           57
 p0                        =  8.11901e-06   +/-   4.90648e-07
 p1                        =     -1.34047   +/-   0.00919871
 p2                        =     0.378673   +/-   0.016582
 
 */

// 
/// \file mpdPrimaryGeneratorAction.cc
/// \brief Implementation of the mpdPrimaryGeneratorAction class

#include "mpdPrimaryGeneratorAction.hh"
#include "mpdDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath> 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorAction::mpdPrimaryGeneratorAction(mpdDetectorConstruction* da)
//: G4VUserPrimaryGeneratorAction(), DetConst(da),
 : fParticleGun(nullptr),DetConst(da)
{
  G4int nofParticles = 100;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  fParticleGun->SetParticleDefinition(particleDefinition);
  //fParticleGun->SetParticlePosition(G4ThreeVector(100., 30., -38.));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(2.*GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorAction::~mpdPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("mpdPrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  
//  G4double taxa = 450.;
 // G4double worldX = 200.*m; //276
 // G4double worldY = 20.*m;  // 30
 // G4double worldZ = 38.*m;  //43.2
  //G4double r = std::sqrt(pow(worldX,2) + pow(worldY,2) + pow(worldZ,2)); 
 // G4double r = worldZHalfLength + 25.*m;
//   for(int i = 0; i<taxa;i++){ 
//    G4double costheta =  2.*G4UniformRand()-1.0;
    G4double costheta =  G4UniformRand();
    G4double sintheta = std::sqrt(1. - pow(costheta,2));
//    G4double phi = -CLHEP::pi + CLHEP::twopi*G4UniformRand();
    G4double phi = CLHEP::twopi*G4UniformRand();
    G4double px = sintheta*cos(phi); // 
    G4double py = sintheta*sin(phi);
    G4double pz = costheta;

 std::cout << "px = " << px << "  py = " << py << "  pz = " << pz << std::endl;   
 
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    //fParticleGun->SetParticlePosition(G4ThreeVector(x, y, -z));
//   }
   
   //fParticleGun->SetParticlePosition(G4ThreeVector(x, y, -z)); //(old)
  
 /* 
  G4double cos_theta = zenith_PDF->GetRandom();
    G4double phi = -CLHEP::pi + CLHEP::twopi*G4UniformRand(); //phi uniform in [-pi, pi];
    G4double sin_theta = std::sqrt(1 - cos_theta*cos_theta);
    G4double x = sin_theta*cos(phi);
    G4double y = sin_theta*sin(phi);
    G4double z = cos_theta;
  */
  
  
  //Set gun position  (old)

//   G4double gunX = -15*m-DetConst->GetcalorSizeX()/2;
//   G4double gunY = -2*m-DetConst->GetcalorSizeY()/2;


   G4double gunX = -DetConst->GetcalorSizeX()/2-50*(2*G4UniformRand()-1)*m;
   G4double gunY = -DetConst->GetcalorSizeY()/2-50*(2*G4UniformRand()-1)*m;

   std::cout << "gunX = " << gunX << "  gunY = " << gunY << std::endl; //(old)
  
  fParticleGun
   ->SetParticlePosition(G4ThreeVector(gunX, gunY, 1*m - worldZHalfLength));
  //  ->SetParticlePosition(G4ThreeVector(gunX, gunY , 40.5*m/2));
  //  ->SetParticlePosition(G4ThreeVector(2*G4UniformRand()*worldBox->GetXHalfLength()-1, 2*G4UniformRand()*worldBox->GetYHalfLength()-1, 0.));
 //   ->SetParticlePosition(G4ThreeVector(0., 0., 0.));
 
 
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

