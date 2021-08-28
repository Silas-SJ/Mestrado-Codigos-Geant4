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
/// \file mpdEventAction.hh
/// \brief Definition of the mpdEventAction class

#ifndef mpdEventAction_h
#define mpdEventAction_h 1

#include "G4UserEventAction.hh"

#include "mpdCalorHit.hh"

#include "globals.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in the hits collections.
class mpdPrimaryGeneratorAction;
class mpdEventAction : public G4UserEventAction
{
public:
  mpdEventAction(mpdPrimaryGeneratorAction* ga);
  virtual ~mpdEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
    
private:
  // methods
  mpdCalorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
  void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                            G4double gapEdep, G4double gapTrackLength) const;
 
  mpdPrimaryGeneratorAction* genaction; 
  // data members                   
  G4int  fAbs3HCID;
  G4int  fAbs4HCID;
  G4int  fAbs5HCID;
  G4int  fAbs6HCID;
  G4int  fAbs7HCID;
  G4int  fAbs8HCID;
  G4int  fAbs9HCID;
  G4int  fAbs10HCID;
  G4int  fAbs11HCID;
  G4int  fAbs12HCID;

  G4int  fGap3HCID;
  G4int  fGap4HCID; 
  G4int  fGap5HCID; 
  G4int  fGap6HCID; 
  G4int  fGap7HCID; 
  G4int  fGap8HCID; 
  G4int  fGap9HCID; 
  G4int  fGap10HCID; 
  G4int  fGap11HCID; 
  G4int  fGap12HCID; 
  
  G4int flayer3HCID;
  G4int flayer4HCID;
  G4int flayer5HCID;
  G4int flayer6HCID;
  G4int flayer7HCID;
  G4int flayer8HCID;
  G4int flayer9HCID;
  G4int flayer10HCID;
  G4int flayer11HCID;
  G4int flayer12HCID;
  
  G4int  fScintHCID;
  G4int fDetabsHCID;

};
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
