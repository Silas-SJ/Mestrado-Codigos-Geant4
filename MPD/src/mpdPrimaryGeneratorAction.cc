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

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdPrimaryGeneratorAction::mpdPrimaryGeneratorAction(mpdDetectorConstruction* da)
//: G4VUserPrimaryGeneratorAction(), DetConst(da),
 : fParticleGun(nullptr),DetConst(da)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleEnergy(0.04*GeV);
  //fParticleGun->SetParticlePosition(G4ThreeVector(100., 30., -38.));

/* Create muon (pion) angular spectrum */
  for (int n = 0; n < 150; n++)
  { 
    thetabin[n] = cos((M_PI/2.*n)/200.)*cos((M_PI/2.*n)/200.);
  }

/* Load cosmic pion spectrum */
 std::ifstream pispecFile;
 pispecFile.open("spec-pion.dat");
 int it = 0;
 while(pispecFile >> energybin[it] >> fluxbin[it]){ 
 it++;
 }
 pispecFile.close();
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

/* Use pion spectrum distribution */
  CLHEP::RandGeneral piEnergyDist(fluxbin,20);
//  fParticleGun->SetParticleEnergy(energybin[0]+piEnergyDist.shoot()*2023.6555*MeV);
//  fParticleGun->SetParticleEnergy(0.03*GeV);
  
//    G4double costheta =  G4UniformRand();
    CLHEP::RandGeneral muonDist(thetabin,150);
    Theta = muonDist.shoot()*M_PI/2;
    G4double costheta = cos(Theta); 
    G4double sintheta = std::sqrt(1. - pow(costheta,2));
//    Theta = acos(costheta);
    Phi = CLHEP::twopi*G4UniformRand();
    pX = sintheta*cos(Phi); // 
    pY = sintheta*sin(Phi);
    pZ = costheta;
//   std::cout << "Theta = " << Theta << "  Phi = " << Phi <<  std::endl;
//   std::cout << "px = " << pX << "  py = " << pY << "  pz = " << pZ << std::endl;   
 
//    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pX, pY, pZ));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

   X = -DetConst->GetcalorSizeX()/2-50*(2*G4UniformRand()-1)*m;
   Y = -DetConst->GetcalorSizeY()/2-50*(2*G4UniformRand()-1)*m;
//   X = 50*(2*G4UniformRand()-1)*m;
//   Y = 50*(2*G4UniformRand()-1)*m;
   Z = 1*m - worldZHalfLength; 
//   std::cout << "gunX = " << X << "  gunY = " << Y << std::endl; //(old)
   
   gunenergy = fParticleGun->GetParticleEnergy();
  
//  fParticleGun ->SetParticlePosition(G4ThreeVector(X, Y, Z));

 
  fParticleGun->SetParticlePosition(G4ThreeVector(-110.*m, -10.5*m, 18.5*m/2));
//  fParticleGun->SetParticlePosition(G4ThreeVector(-110.*m, -10.5*m, Z));

  fParticleGun->GeneratePrimaryVertex(anEvent);
 // fParticleGun->GeneratePrimaryVertex(anEvent);
 // fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

