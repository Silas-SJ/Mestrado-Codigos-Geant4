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
/// \file mpdEventAction.cc
/// \brief Implementation of the mpdEventAction class

#include "mpdEventAction.hh"
#include "mpdCalorimeterSD.hh"
#include "mpdCalorHit.hh"
#include "mpdAnalysis.hh"
#include "mpdPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdEventAction::mpdEventAction(mpdPrimaryGeneratorAction* ga)
 : G4UserEventAction(),genaction(ga),
  /* fAbsHCID(-1),
   fGapHCID(-1),
  // fGap2HCID(-1), // gap of detector layer
   flayer3HCID(-1),
   fScintHCID(-1),
   fDetabsHCID(-1)
   */
  fAbs3HCID(-1), fAbs4HCID(-1), fAbs5HCID(-1), fAbs6HCID(-1), fAbs7HCID(-1), fAbs8HCID(-1),
  fAbs9HCID(-1), fAbs10HCID(-1), fAbs11HCID(-1), fAbs12HCID(-1),
  fGap3HCID(-1), fGap4HCID(-1), fGap5HCID(-1), fGap6HCID(-1), fGap7HCID(-1), fGap8HCID(-1),
  fGap9HCID(-1), fGap10HCID(-1), fGap11HCID(-1), fGap12HCID(-1),
  flayer3HCID(-1), flayer4HCID(-1), flayer5HCID(-1), flayer6HCID(-1), flayer7HCID(-1),
  flayer8HCID(-1), flayer9HCID(-1), flayer10HCID(-1), flayer11HCID(-1), flayer12HCID(-1),
  fScintHCID(-1),
  fDetabsHCID(-1)
  
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdEventAction::~mpdEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorHitsCollection*
mpdEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<mpdCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("mpdEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: " 
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: " 
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: " 
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdEventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections IDs (only once)
  if ( fAbs3HCID == -1 ) {
    fAbs3HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber3HitsCollection");
    fAbs4HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber4HitsCollection");  
    fAbs5HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber5HitsCollection");   
    fAbs6HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber6HitsCollection");   
    fAbs7HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber7HitsCollection");  
    fAbs8HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber8HitsCollection");
   fAbs9HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber9HitsCollection"); 
   fAbs10HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber10HitsCollection");   
   fAbs11HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber11HitsCollection");   
   fAbs12HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Absorber12HitsCollection");     
          
    fGap3HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap3HitsCollection");
   fGap4HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap4HitsCollection");
   fGap5HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap5HitsCollection");        
   fGap6HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap6HitsCollection");
   fGap7HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap7HitsCollection");        
   fGap8HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap8HitsCollection");    
  fGap9HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap9HitsCollection");     
   fGap10HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap10HitsCollection");    
   fGap11HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap11HitsCollection");   
   fGap12HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("Gap12HitsCollection"); 
                  
    flayer3HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer3HitsCollection");
    flayer4HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer4HitsCollection");   
    flayer5HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer5HitsCollection");
    flayer6HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer6HitsCollection");   
   flayer7HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer7HitsCollection");   
   flayer8HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer8HitsCollection");   
   flayer9HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer9HitsCollection");    
   flayer10HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer10HitsCollection");    
   flayer11HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer11HitsCollection");    
   flayer12HCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("layer12HitsCollection");   
      
    fScintHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("ScintHitsCollection");
    fDetabsHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("DetabsHitsCollection");
  }

  // Get hits collections

  auto abso3HC = GetHitsCollection(fAbs3HCID, event);
  auto abso4HC = GetHitsCollection(fAbs4HCID, event);
  auto abso5HC = GetHitsCollection(fAbs5HCID, event);
  auto abso6HC = GetHitsCollection(fAbs6HCID, event);
  auto abso7HC = GetHitsCollection(fAbs7HCID, event);
  auto abso8HC = GetHitsCollection(fAbs8HCID, event);
  auto abso9HC = GetHitsCollection(fAbs9HCID, event);
  auto abso10HC = GetHitsCollection(fAbs10HCID, event);
  auto abso11HC = GetHitsCollection(fAbs11HCID, event);
  auto abso12HC = GetHitsCollection(fAbs12HCID, event);
  
  auto gap3HC = GetHitsCollection(fGap3HCID, event);
  auto gap4HC = GetHitsCollection(fGap4HCID, event);
  auto gap5HC = GetHitsCollection(fGap5HCID, event);
  auto gap6HC = GetHitsCollection(fGap6HCID, event);
  auto gap7HC = GetHitsCollection(fGap7HCID, event);
  auto gap8HC = GetHitsCollection(fGap8HCID, event);
  auto gap9HC = GetHitsCollection(fGap9HCID, event);
  auto gap10HC = GetHitsCollection(fGap10HCID, event);
  auto gap11HC = GetHitsCollection(fGap11HCID, event);
  auto gap12HC = GetHitsCollection(fGap12HCID, event);

  auto layer3HC = GetHitsCollection(flayer3HCID, event);
  auto layer4HC = GetHitsCollection(flayer4HCID, event);
  auto layer5HC = GetHitsCollection(flayer5HCID, event);
  auto layer6HC = GetHitsCollection(flayer6HCID, event);
  auto layer7HC = GetHitsCollection(flayer7HCID, event);
  auto layer8HC = GetHitsCollection(flayer8HCID, event);
  auto layer9HC = GetHitsCollection(flayer9HCID, event);
  auto layer10HC = GetHitsCollection(flayer10HCID, event);
  auto layer11HC = GetHitsCollection(flayer11HCID, event);
  auto layer12HC = GetHitsCollection(flayer12HCID, event);
  
  auto scintHC = GetHitsCollection(fScintHCID, event);
  auto detabsHC = GetHitsCollection(fDetabsHCID, event);
  
  
  
  // Get hit with total values
  auto abso3Hit = (*abso3HC)[abso3HC->entries()-1];
  auto gap3Hit = (*gap3HC)[gap3HC->entries()-1];
 // auto scintHit = (*scintHC)[abso3HC->entries()-1];
//  auto detabsHit = (*detabsHC)[gap3HC->entries()-1];
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

  PrintEventStatistics(
      abso3Hit->GetEdep(), abso3Hit->GetTrackLength(),
      gap3Hit->GetEdep(), gap3Hit->GetTrackLength());
      
  }  
  
  // Fill histograms, ntuple
  //

  // get analysis manager
  
  auto analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
  
  analysisManager->FillH1(0, abso3Hit->GetEdep());
  analysisManager->FillH1(1, gap3Hit->GetEdep());
  analysisManager->FillH1(2, abso3Hit->GetTrackLength());
  analysisManager->FillH1(3, gap3Hit->GetTrackLength());
  
  // fill ntuple
  for (uint i=0;i<(abso3HC->entries()-1);i++){
   auto abso3Hitlayer = (*abso3HC)[i];
   auto gap3Hitlayer = (*gap3HC)[i];
   auto layer3Hitlayer = (*layer3HC)[i];
   auto scintHitlayer = (*scintHC)[i];
   auto detabsHitlayer = (*detabsHC)[i];
   
   analysisManager->FillNtupleDColumn(0, abso3Hitlayer->GetEdep());
   analysisManager->FillNtupleDColumn(1, gap3Hitlayer->GetEdep());
   analysisManager->FillNtupleDColumn(2, abso3Hitlayer->GetTrackLength());
   analysisManager->FillNtupleDColumn(3, gap3Hitlayer->GetTrackLength());
   analysisManager->FillNtupleIColumn(4, scintHitlayer->GetPionDecay());
   analysisManager->FillNtupleIColumn(5, scintHitlayer->GetEdep());
   analysisManager->FillNtupleIColumn(6, detabsHitlayer->GetPionDecay());
   analysisManager->FillNtupleIColumn(7, detabsHitlayer->GetEdep());
   analysisManager->FillNtupleIColumn(8, scintHitlayer->GetMuonDecay());
   analysisManager->FillNtupleIColumn(9, detabsHitlayer->GetMuonDecay());
 //  analysisManager->FillNtupleDColumn(10, abso3Hitlayer->GetLayerID());
   analysisManager->FillNtupleDColumn(10, layer3Hitlayer->GetPionPassed());
   analysisManager->FillNtupleDColumn(11, eventID);
   analysisManager->FillNtupleDColumn(12, abso3Hitlayer->GetPionDecay());
   analysisManager->FillNtupleDColumn(13, gap3Hitlayer->GetPionDecay());
   analysisManager->FillNtupleDColumn(14, genaction->GetTheta());
   analysisManager->FillNtupleDColumn(15, genaction->GetPhi());
   analysisManager->FillNtupleDColumn(16, genaction->GetEnergyPrimary());
   analysisManager->FillNtupleDColumn(17, genaction->GetMomentumX());
   analysisManager->FillNtupleDColumn(18, genaction->GetMomentumY());
   analysisManager->FillNtupleDColumn(19, genaction->GetMomentumZ());
   analysisManager->FillNtupleDColumn(20, genaction->GetPositionX());
   analysisManager->FillNtupleDColumn(21, genaction->GetPositionY());
   analysisManager->FillNtupleDColumn(22, genaction->GetPositionZ());
   analysisManager->FillNtupleDColumn(23, scintHitlayer->GetPionPassed());
   analysisManager->FillNtupleDColumn(24, detabsHitlayer->GetPionPassed());
   analysisManager->FillNtupleDColumn(25, gap3Hitlayer->GetPionPassed());
   analysisManager->FillNtupleDColumn(26, scintHitlayer->GetPionCapture());
   analysisManager->FillNtupleDColumn(27, detabsHitlayer->GetPionCapture());
   analysisManager->FillNtupleDColumn(28, gap3Hitlayer->GetPionDecay());
   analysisManager->AddNtupleRow();
  }
}  



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
