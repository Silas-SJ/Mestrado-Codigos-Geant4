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
/// \file mpdPrimaryGeneratorAction.hh
/// \brief Definition of the mpdPrimaryGeneratorAction class

#ifndef mpdPrimaryGeneratorAction_h
#define mpdPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"


class G4ParticleGun;
class G4Event;
class mpdDetectorConstruction;

/// The primary generator action class with particle gun.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class mpdPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  mpdPrimaryGeneratorAction(mpdDetectorConstruction* detconst);
  virtual ~mpdPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  
  // set methods
  void SetRandomFlag(G4bool value);

   G4double GetEnergyPrimary() { return gunenergy; }
   G4double GetMomentumX() { return pX; }
   G4double GetMomentumY() { return pY; }
   G4double GetMomentumZ() { return pZ; }
   G4double GetPositionX() { return X; }
   G4double GetPositionY() { return Y; }
   G4double GetPositionZ() { return Z; }
   G4double GetTheta() { return Theta; }
   G4double GetPhi()   { return Phi; }
private:
  G4ParticleGun*  fParticleGun; // G4 particle gunß
    mpdDetectorConstruction* DetConst;
    G4double thetabin[200];
    G4double energybin[20];
    G4double fluxbin[20];
    G4double gunenergy;
    G4double pX;
    G4double pY;
    G4double pZ;
    G4double X;
    G4double Y;
    G4double Z;
    G4double Theta;
    G4double Phi;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
