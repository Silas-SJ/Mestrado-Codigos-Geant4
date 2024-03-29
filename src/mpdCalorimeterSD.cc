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
/// \file mpdCalorimeterSD.cc
/// \brief Implementation of the mpdCalorimeterSD class

#include "mpdCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "mpdPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"

#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 mpdCalorimeterSD::mpdCalorimeterSD(
                                    const G4String& name, 
                                    const G4String& hitsCollectionName,
                                    G4int nofCells
                                    )
 : G4VSensitiveDetector(name), 
   fHitsCollection(nullptr), 
   fNofCells(nofCells)


{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdCalorimeterSD::~mpdCalorimeterSD()
{ 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new mpdCalorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofCells for cells + one more for total sums 
  for (G4int i=0; i<fNofCells+1; i++ ) {
    fHitsCollection->insert(new mpdCalorHit());
  }
  
   G4RunManager *rm = G4RunManager::GetRunManager();
   eventID= rm->GetCurrentEvent()->GetEventID();
   oldID = -9999;
 //std::cout << "Valor TrackID " << oldID << std::endl;
}

 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool mpdCalorimeterSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

 // if ( edep==0. && stepLength == 0. ) return false; // It was inactive

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Get calorimeter cell id 
  auto layerNumber = touchable->GetReplicaNumber(1);

  // Get hit accounting data for this cell
  auto hit = (*fHitsCollection)[layerNumber];
    
    if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber; 
    G4Exception("mpdCalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
    }

  // Get hit for total accounting
  auto hitTotal 
    = (*fHitsCollection)[fHitsCollection->entries()-1];
  auto piondecay = false;
  auto muondecay = false;
  auto pionpassed= false; 
  auto pioncapture= false; 
  
  ParentID = step->GetTrack()->GetParentID();
              // ::::::: Number of pions that captured at rest :::::::::: //
              
     if (abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211)
     {
       if (step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="hBertiniCaptureAtRest")
       {
          pioncapture=true;
      hit->AddPionCapture();
       }
     }
 
              // ::::::: Number of different pions that passed :::::::::: //

     if (abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211)
     {
       if (step->GetTrack()->GetTrackID() != oldID && ParentID != oldID)
       {
         pionpassed=true;
         hit->AddPionPassed();
    /*     
       if (step->GetPreStepPoint()->GetProcessDefinedStep() && step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() && step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName())
      {
     if (step->GetPreStepPoint()->GetPhysicalVolume() && step->GetPreStepPoint()-> GetPhysicalVolume()->GetName() != "Scint" && step->GetPostStepPoint()->GetPhysicalVolume() && step->GetPostStepPoint()->GetPhysicalVolume()->GetName() != "Scint")
    {
    std::cout << "Pion dont passed:" << step->GetPostStepPoint()-> GetPhysicalVolume()->GetName() << std::endl;
    }
   }  
    */
      }
       
      // std::cout << "Pion passed:" << "1" << std::endl;
      // std::cout << "TrackID_if:" << step->GetTrack()->GetTrackID() << std::endl;
      // std::cout << "oldID_if:" << oldID << std::endl;
      // std::cout << "Volume:" << step->GetPostStepPoint()-> GetPhysicalVolume() -> GetName() << std::endl;
      /*
     if (step->GetPostStepPoint()->GetPhysicalVolume() && step->GetPostStepPoint()-> GetPhysicalVolume() -> GetName())
   {
     std::cout << "Volume:" << step->GetPostStepPoint()-> GetPhysicalVolume() -> GetName() << std::endl;
    }  
    
    if (step->GetPreStepPoint()->GetProcessDefinedStep() && step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName())
   {
     std::cout << "PreStep Process_if: " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
    }
    if (step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName())
    {
      std::cout << "PostStep Process_if: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
   
    }    
     */  
      oldID = step->GetTrack()->GetTrackID();
   }
  
     
     
   /* 
            // ::::::: Identify an event:::::::::: //
        
   if (eventID == 202 || eventID == 285 || eventID == 402 || eventID == 547 || eventID == 977)
  {
   // std::cout << "Valor EventID: " << eventID << std::endl;
   std::cout << "Particle name: " << step->GetTrack()->GetDefinition()->GetPDGEncoding() << std::endl;
    std::cout << "Valor TrackID: " << step->GetTrack()->GetTrackID() << std::endl;
    if (step->GetPreStepPoint()->GetProcessDefinedStep() && step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName())
    {
    // std::cout << "Particle name: " << step->GetTrack()->GetDefinition()->GetPDGEncoding() << std::endl;
      std::cout << "PreStep Process: " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
      
     std::cout << "PreStep Position: " << step->GetPreStepPoint()->GetPosition() << std::endl;
    }
   
   // if (step->GetPreStepPoint()->GetPhysicalVolume() && step->GetPreStepPoint()-> GetPhysicalVolume()->GetName())
   // {
  //  std::cout << "volume_Pre:" << step->GetPostStepPoint()-> GetPhysicalVolume()->GetName() << std::endl;
   // }
    
    if (step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName())
    {
     std::cout << "PostStep Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
     
     std::cout << "PostStep Position: " << step->GetPostStepPoint()->GetPosition() << std::endl;
     
    }
    if (step->GetPostStepPoint()->GetPhysicalVolume() && step->GetPostStepPoint()->GetPhysicalVolume()->GetName())
    {
    std::cout << "volume_Pos:" << step->GetPostStepPoint()-> GetPhysicalVolume()->GetName() << std::endl;
    }
  }       
   */

    hit->Add(edep, stepLength, layerNumber);
    hitTotal->Add(edep, stepLength, layerNumber); 
    
    
 // :::::::Number of pions and muons that decay :::::::::: //
 
      if (abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211)
  {
      if (step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay")
     {
      piondecay=true;
      hit->AddPionDecay();
 // std::cout << "pion decay at layer " << layerNumber << std::endl;
//  std::cout << "pion decay:" << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
     }
  }
      else if (abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 13)
  {
      if (step->GetPostStepPoint()->GetProcessDefinedStep() && step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay")
 {
      muondecay=true;
     hit->AddMuonDecay();
 //(teste)     std::cout << "muon decay at layer " << layerNumber << std::endl;
   //std::cout << "muon decay:" << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
    }
  }
    
    
    
/*
  // Add values
 // if ( edep==0. && stepLength == 0. ){
    hit->Add(edep, stepLength, layerNumber);
    hitTotal->Add(edep, stepLength, layerNumber);
    if(piondecay) hit->AddPionDecay();
    if(muondecay) hit->AddMuonDecay();
    if(pionpassed) hit->AddPionPassed();
    if(pioncapture) hit->AddPionCapture();
   //}
   */
 
  /* 
   // :::::::Process of pion interation:::::::::: //
   if (step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211){
    std::cout << "PreStep Process Pion: " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
    }
    
   if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() && abs(step->GetTrack()->GetDefinition()->GetPDGEncoding()) == 211){

 std::cout << "PostStep Process Pion: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;  
    }
*/

// std::cout << "Particle: " << step->GetTrack()->GetDefinition()->GetParticleName()<< std::endl;
// std::cout << "PreStep Process: " << step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
 //std::cout << "PostStep Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
 
 /*
   if (step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Decay"){
  std::cout << "decay at layer " << layerNumber << "  particle: " << step->GetTrack()->GetDefinition()->GetParticleName() << std::endl;
    }
 */
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
