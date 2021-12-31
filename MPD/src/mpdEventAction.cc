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
#include "mpdPrimaryGeneratorAction.hh"
#include "mpdEventAction.hh"
#include "mpdCalorimeterSD.hh"
#include "mpdCalorHit.hh"
#include "mpdDigi.hh"
#include "mpdDigitizer.hh"
#include "mpdAnalysis.hh"
#include "mpdTreeManager.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4DigiManager.hh"
#include "G4EventManager.hh"
#include "G4VHitsCollection.hh"
#include "Randomize.hh"
#include <iomanip>
#include <set>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdEventAction::mpdEventAction(mpdPrimaryGeneratorAction* ga,mpdTreeManager * tree)
 : G4UserEventAction(),genaction(ga),
   fAbsHCID(-1),
   fGapHCID(-1),
   fScintHCID(-1),
   fDetabsHCID(-1),
   fPmtHCID(-1),
   fPmtHSCID(-1),
   fTreeManager(tree)
{
    mpdDigitizer * mpdDM = new mpdDigitizer( "mpdDigitizer" );
    G4DigiManager::GetDMpointer()->AddNewModule(mpdDM);
}

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
//  G4double pmtTime[10000];
//  G4double pmtCharge[10000];
//  G4double pmtTotalCharge1[12];
//  G4double pmtTotalCharge2[12];
//  G4double pmtAvgTime1[12];
//  G4double pmtAvgTime2[12];
//  G4double ScintCharge[12];
//  G4int PMTNumber; //numero do cintilador
//  G4int ntime1;
//  G4int ntime2;
    
  // Get hits collections IDs (only once)
  if ( fAbsHCID == -1 ) {
    fAbsHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
    fScintHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("ScintHitsCollection");
    fDetabsHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("DetabsHitsCollection");
    fPmtHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("pmtHitCollection");
    fPmtHSCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("pmtSignalHitCollection");
  }
  /*
    G4cout << "fAbsHCID: " << fAbsHCID << G4endl;
    G4cout << "fGapHCID: " << fGapHCID << G4endl;
    G4cout << "fScintHCID: " << fScintHCID << G4endl;
    G4cout << "fDetabsHCID: " << fDetabsHCID << G4endl;
    G4cout << "fPmtHCID: " << fPmtHCID << G4endl;
    G4cout << "fPmtHSCID: " << fPmtHSCID << G4endl;
*/
    G4DigiManager * fDM = G4DigiManager::GetDMpointer();
    
  // Get hits collections
  auto absHC = GetHitsCollection(fAbsHCID, event);
  auto gapHC = GetHitsCollection(fGapHCID, event);
  auto scintHC = GetHitsCollection(fScintHCID, event);
  auto detabsHC = GetHitsCollection(fDetabsHCID, event);
//  auto pmtHC = GetHitsCollection(fPmtHCID, event);
//  auto pmtSHC = GetHitsCollection(fPmtHSCID, event);
  // Get hit with total values
//  auto absoHit = (*absoHC)[absoHC->entries()-1];
//  auto gapHit = (*gapHC)[gapHC->entries()-1];
  //auto scintHit = (*scintHC)[absoHC->entries()-1];
  //auto detabsHit = (*detabsHC)[gapHC->entries()-1];
  // Print per event (modulo n)
  //
  std::map<G4int,G4int> ScintPionDecay;
  std::map<G4int,G4int> ScintMuonDecay;
  std::map<G4int,std::set<G4int>> ScintPionPassed;
  std::map<G4int,std::set<G4int>> ScintMuonPassed;
  std::set<G4int> ScintPionID;
  std::set<G4int> ScintMuonID;
  std::map<G4int,G4int> DetabsPionDecay;
  std::map<G4int,G4int> DetabsMuonDecay;
  std::map<G4int,std::set<G4int>> DetabsPionPassed;
  std::map<G4int,std::set<G4int>> DetabsMuonPassed;
  std::set<G4int> DetabsPionID;
  std::set<G4int> DetabsMuonID;
  auto eventID = event->GetEventID();

  for (uint i=0;i<(scintHC->entries());i++){
   auto scintHitlayer = (*scintHC)[i];     
/*	   G4cout << "Pion TrackID: " << scintHitlayer->GetPionTrackID() << "  DetNo: " << scintHitlayer->GetLayerID() << "  Edep: " << scintHitlayer->GetEdep() << "  TrackLength: " << scintHitlayer->GetTrackLength() << G4endl;
	   G4cout << "Muon TrackID: " << scintHitlayer->GetMuonTrackID() << "  DetNo: " << scintHitlayer->GetLayerID() << "  Edep: " << scintHitlayer->GetEdep() << "  TrackLength: " << scintHitlayer->GetTrackLength() << G4endl;
	   G4cout << "Pion Decay: " << scintHitlayer->GetPionDecay() << "  DetNo: " << scintHitlayer->GetLayerID() << "  Edep: " << scintHitlayer->GetEdep() << "  TrackLength: " << scintHitlayer->GetTrackLength() << G4endl;
	   G4cout << "Muon Decay: " << scintHitlayer->GetMuonDecay() << "  DetNo: " << scintHitlayer->GetLayerID() << "  Edep: " << scintHitlayer->GetEdep() << "  TrackLength: " << scintHitlayer->GetTrackLength() << G4endl;
*/
   if(scintHitlayer->GetPionTrackID()!=-111) {
    ScintPionID.insert(scintHitlayer->GetPionTrackID());
    ScintPionPassed[scintHitlayer->GetLayerID()]=ScintPionID;        
   }
   if(scintHitlayer->GetMuonTrackID()!=-111) {
    ScintMuonID.insert(scintHitlayer->GetMuonTrackID());
    ScintMuonPassed[scintHitlayer->GetLayerID()]=ScintMuonID;        
   }
   ScintPionDecay[scintHitlayer->GetLayerID()]+=scintHitlayer->GetPionDecay();
   ScintMuonDecay[scintHitlayer->GetLayerID()]+=scintHitlayer->GetMuonDecay();
  }

  for (uint i=0;i<(detabsHC->entries());i++){
   auto detabsHitlayer = (*detabsHC)[i];     
   if(detabsHitlayer->GetPionTrackID()!=-111) {
    DetabsPionID.insert(detabsHitlayer->GetPionTrackID());
    DetabsPionPassed[detabsHitlayer->GetLayerID()]=ScintPionID;        
   }
   if(detabsHitlayer->GetMuonTrackID()!=-111) {
    DetabsMuonID.insert(detabsHitlayer->GetMuonTrackID());
    DetabsMuonPassed[detabsHitlayer->GetLayerID()]=DetabsMuonID;        
   }
   DetabsPionDecay[detabsHitlayer->GetLayerID()]+=detabsHitlayer->GetPionDecay();
   DetabsMuonDecay[detabsHitlayer->GetLayerID()]+=detabsHitlayer->GetMuonDecay();
  }
  std::map<G4int,G4int> ScintPionPassou; //PionPassed
  for(auto i=ScintPionPassed.begin();i!=ScintPionPassed.end();++i){
  ScintPionPassou[i->first]=i->second.size();
  // G4cout << "Number of pions: " << i->second.size() << " on detector: " << i->first << G4endl;    
  }
  std::map<G4int,G4int> ScintMuonPassou; //MuonPassed
  for(auto i=ScintMuonPassed.begin();i!=ScintMuonPassed.end();++i){
   ScintMuonPassou[i->first]=i->second.size();
  // G4cout << "Number of muons: " << i->second.size() << " on detector: " << i->first << G4endl;    
  }
  /*
  for(auto i=ScintPionDecay.begin();i!=ScintPionDecay.end();++i){
   G4cout << "Number of pion decay: " << i->second << " on detector: " << i->first << G4endl;    
  }
  for(auto i=ScintMuonDecay.begin();i!=ScintMuonDecay.end();++i){
   G4cout << "Number of muon decay: " << i->second << " on detector: " << i->first << G4endl;    
  }
  */
  std::map<G4int,G4int> DetabsPionPassou; //PionPassed
  for(auto i=DetabsPionPassed.begin();i!=DetabsPionPassed.end();++i){
  DetabsPionPassou[i->first]=i->second.size();
 //  G4cout << "Number of pions: " << i->second.size() << " on detabsorber: " << i->first << G4endl;    
  }
  /*
  for(auto i=DetabsMuonPassed.begin();i!=DetabsMuonPassed.end();++i){
   G4cout << "Number of muons: " << i->second.size() << " on detabsorber: " << i->first << G4endl;    
  }
  for(auto i=DetabsPionDecay.begin();i!=DetabsPionDecay.end();++i){
   G4cout << "Number of pion decay: " << i->second << " on detabsorber: " << i->first << G4endl;    
  }
  for(auto i=DetabsMuonDecay.begin();i!=DetabsMuonDecay.end();++i){
   G4cout << "Number of muon decay: " << i->second << " on detabsorber: " << i->first << G4endl;    
  }
*/
  std::map<G4int,G4int> GapPionDecay;
  std::map<G4int,G4int> GapMuonDecay;
  std::map<G4int,std::set<G4int>> GapPionPassed;
  std::map<G4int,std::set<G4int>> GapMuonPassed;
  std::set<G4int> GapPionID;
  std::set<G4int> GapMuonID;
  std::map<G4int,G4int> AbsPionDecay;
  std::map<G4int,G4int> AbsMuonDecay;
  std::map<G4int,std::set<G4int>> AbsPionPassed;
  std::map<G4int,std::set<G4int>> AbsMuonPassed;
  std::set<G4int> AbsPionID;
  std::set<G4int> AbsMuonID;

  for (uint i=0;i<(gapHC->entries());i++){
   auto gapHitlayer = (*gapHC)[i];
/*         G4cout << "Pion TrackID: " << gapHitlayer->GetPionTrackID() << "  DetNo: " << gapHitlayer->GetLayerID() << "  Edep: " << gapHitlayer->GetEdep() << "  TrackLength: " << gapHitlayer->GetTrackLength() << G4endl;
           G4cout << "Muon TrackID: " << gapHitlayer->GetMuonTrackID() << "  DetNo: " << gapHitlayer->GetLayerID() << "  Edep: " << gapHitlayer->GetEdep() << "  TrackLength: " << gapHitlayer->GetTrackLength() << G4endl;
           G4cout << "Pion Decay: " << gapHitlayer->GetPionDecay() << "  DetNo: " << gapHitlayer->GetLayerID() << "  Edep: " << gapHitlayer->GetEdep() << "  TrackLength: " << gapHitlayer->GetTrackLength() << G4endl;
           G4cout << "Muon Decay: " << gapHitlayer->GetMuonDecay() << "  DetNo: " << gapHitlayer->GetLayerID() << "  Edep: " << gapHitlayer->GetEdep() << "  TrackLength: " << gapHitlayer->GetTrackLength() << G4endl;
*/
   if(gapHitlayer->GetPionTrackID()!=-111) {
    GapPionID.insert(gapHitlayer->GetPionTrackID());
    GapPionPassed[gapHitlayer->GetLayerID()]=GapPionID;
   }
   if(gapHitlayer->GetMuonTrackID()!=-111) {
    GapMuonID.insert(gapHitlayer->GetMuonTrackID());
    GapMuonPassed[gapHitlayer->GetLayerID()]=GapMuonID;
   }
   GapPionDecay[gapHitlayer->GetLayerID()]+=gapHitlayer->GetPionDecay();
   GapMuonDecay[gapHitlayer->GetLayerID()]+=gapHitlayer->GetMuonDecay();
  }

  for (uint i=0;i<(absHC->entries());i++){
   auto absHitlayer = (*absHC)[i];
   if(absHitlayer->GetPionTrackID()!=-111) {
    AbsPionID.insert(absHitlayer->GetPionTrackID());
    AbsPionPassed[absHitlayer->GetLayerID()]=GapPionID;
   }
   if(absHitlayer->GetMuonTrackID()!=-111) {
    AbsMuonID.insert(absHitlayer->GetMuonTrackID());
    AbsMuonPassed[absHitlayer->GetLayerID()]=AbsMuonID;
   }
   AbsPionDecay[absHitlayer->GetLayerID()]+=absHitlayer->GetPionDecay();
   AbsMuonDecay[absHitlayer->GetLayerID()]+=absHitlayer->GetMuonDecay();
  }
  std::map<G4int,G4int> GapPionPassou; //PionPassed
  for(auto i=GapPionPassed.begin();i!=GapPionPassed.end();++i){
   GapPionPassou[i->first]=i->second.size();
  // G4cout << "Number of pions: " << i->second.size() << " on gap: " << i->first << G4endl;
  }
  /*
  for(auto i=GapMuonPassed.begin();i!=GapMuonPassed.end();++i){
   G4cout << "Number of muons: " << i->second.size() << " on gap: " << i->first << G4endl;
  }
  for(auto i=GapPionDecay.begin();i!=GapPionDecay.end();++i){
   G4cout << "Number of pion decay: " << i->second << " on gap: " << i->first << G4endl;
  }
  for(auto i=GapMuonDecay.begin();i!=GapMuonDecay.end();++i){
   G4cout << "Number of muon decay: " << i->second << " on gap: " << i->first << G4endl;
  }
  for(auto i=AbsPionPassed.begin();i!=AbsPionPassed.end();++i){
   G4cout << "Number of pions: " << i->second.size() << " on absorber: " << i->first << G4endl;   
  }
  for(auto i=AbsMuonPassed.begin();i!=AbsMuonPassed.end();++i){
   G4cout << "Number of muons: " << i->second.size() << " on absorber: " << i->first << G4endl;  
  }
  for(auto i=AbsPionDecay.begin();i!=AbsPionDecay.end();++i){
   G4cout << "Number of pion decay: " << i->second << " on absorber: " << i->first << G4endl;
  }
  for(auto i=AbsMuonDecay.begin();i!=AbsMuonDecay.end();++i){
   G4cout << "Number of muon decay: " << i->second << " on absorber: " << i->first << G4endl;
  }
*/
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

//    PrintEventStatistics(
//      absoHit->GetEdep(), absoHit->GetTrackLength(),
//      gapHit->GetEdep(), gapHit->GetTrackLength());
  }  

 // Digitization
 mpdDigitizer * mpdDM =
 (mpdDigitizer*)fDM->FindDigitizerModule( "mpdDigitizer" );
    if (!mpdDM){
          G4ExceptionDescription msg;
          msg << "Cannot find mpdDigitizer";
          G4Exception("mpdEventAction::EndOfEventAction()",
            "", FatalException, msg);
    }
    
 mpdDM->Digitize();

 G4int DCID = fDM->GetDigiCollectionID("mpdDigitsCollection");
  //  G4cout << "DCID: " << DCID << G4endl;
 mpdDigitsCollection * DC =
 (mpdDigitsCollection*)fDM->GetDigiCollection(DCID);
   // G4cout << "DC->entries() " << DC->entries() << G4endl;
/*
 if(DC) {
  G4int n_digi =  DC->entries();
  for (G4int i=0;i<n_digi;i++) {
      G4cout << "PMT: " << (*DC)[i]->GetPMT() << " Carga: " << (*DC)[i]->GetPhoton_eletron() << " Energia:  " << (*DC)[i]->GetOpticalPhotonEnergy() << " Photons Abs:  " <<
          (*DC)[i]->GetOpticalPhotonAbs() << G4endl;
  }
 }
*/
//std::vector <G4double> tv0;
//std::vector <G4double> tv1;
//std::map<int,std::vector<G4double>> tv;
std::vector <G4int> otpv;
std::vector <G4double> ottv;
 if(DC) {
  G4int n_digi =  DC->entries();
  for (G4int i=0;i<n_digi;i++) {
//   tv0 = (*DC)[i]->GetOverThresholdTimeVector0();
//   tv1 = (*DC)[i]->GetOverThresholdTimeVector1();
//   G4cout << "TimeVector0 Size: " << tv0.size() << G4endl;
//   G4cout << "TimeVector1 Size: " << tv1.size() << G4endl;
//   G4cout << "TimeVector Size: " << tv.size() << G4endl;
   otpv = (*DC)[i]->GetOverThresholdPMTVec();
   ottv = (*DC)[i]->GetOverThresholdTimeVec();
  }
 }
  // Fill histograms, ntuple
  //
fTreeManager->FillNtuple(
genaction->GetTheta(), genaction->GetPhi(), genaction->GetEnergyPrimary(), genaction->GetMomentumX(), genaction->GetMomentumY(), genaction->GetMomentumZ(), genaction->GetPositionX(), genaction->GetPositionY(), genaction->GetPositionZ(), ScintPionDecay, ScintPionPassou, ScintMuonDecay, ScintMuonPassou, DetabsPionDecay, DetabsPionPassou, GapPionDecay, GapPionPassou,  otpv, ottv);
//genaction->GetTheta(), genaction->GetPhi(), genaction->GetEnergyPrimary(), genaction->GetMomentumX(), genaction->GetMomentumY(), genaction->GetMomentumZ(), genaction->GetPositionX(), genaction->GetPositionY(), genaction->GetPositionZ(), tv0, tv1);
/*
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
 
  // fill histograms
  analysisManager->FillH1(0, absoHit->GetEdep());
  analysisManager->FillH1(1, gapHit->GetEdep());
  analysisManager->FillH1(2, absoHit->GetTrackLength());
  analysisManager->FillH1(3, gapHit->GetTrackLength());
  
  // fill ntuple
  for (uint i=0;i<(absoHC->entries()-1);i++){
   auto absoHitlayer = (*absoHC)[i];
   auto gapHitlayer = (*gapHC)[i];
   auto scintHitlayer = (*scintHC)[i];
   auto detabsHitlayer = (*detabsHC)[i];
      G4cout << "Edep: " << absoHitlayer->GetEdep() << G4endl;
   analysisManager->FillNtupleDColumn(0, absoHitlayer->GetEdep());
   analysisManager->FillNtupleDColumn(1, gapHitlayer->GetEdep());
   analysisManager->FillNtupleDColumn(2, absoHitlayer->GetTrackLength());
   analysisManager->FillNtupleDColumn(3, gapHitlayer->GetTrackLength());
   analysisManager->FillNtupleIColumn(4, scintHitlayer->GetPionDecay());
   analysisManager->FillNtupleIColumn(5, detabsHitlayer->GetPionDecay());
   analysisManager->FillNtupleIColumn(6, scintHitlayer->GetMuonDecay());
   analysisManager->FillNtupleIColumn(7, detabsHitlayer->GetMuonDecay());
   analysisManager->FillNtupleDColumn(8, absoHitlayer->GetLayerID());
   analysisManager->FillNtupleDColumn(9, eventID);

   analysisManager->AddNtupleRow();

  }
*/
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
