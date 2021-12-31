//
//  mpdPMTSD.cc
//
//
//  Created by Helio Nogima on 7/11/21.
//

#include "mpdPMTSD.hh"
#include "mpdPMTHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh" 
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"


mpdPMTSD::mpdPMTSD(G4String name)
  :G4VSensitiveDetector(name),pmtHitCollection(0),pmtSignalHitCollection(0), pmtPositionsX(0),pmtPositionsY(0),pmtPositionsZ(0)
{
  collectionName.insert("pmtHitCollection");
  collectionName.insert("pmtSignalHitCollection");
//  G4cout << "Passa no construtor do mpdPMTSD" << G4endl;
}

//=============================================================================
mpdPMTSD::~mpdPMTSD()
{}

//=============================================================================
void mpdPMTSD::Initialize(G4HCofThisEvent* HCE){
  pmtHitCollection = new mpdPMTHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
//
  pmtSignalHitCollection = new mpdPMTHitsCollection
                            (SensitiveDetectorName, collectionName[1]);
//

  static G4int HCID = -1;
  if(HCID<0){
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection( HCID, pmtHitCollection );

//
  static G4int HSCID = -1 ;
  if(HSCID<0){
     HSCID = GetCollectionID(1);
  }
  HCE->AddHitsCollection( HSCID, pmtSignalHitCollection );
//
}

//=============================================================================
G4bool mpdPMTSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist){
  return ProcessHits_constStep(aStep,ROhist);
}

//=============================================================================
G4bool mpdPMTSD::ProcessHits_constStep(const G4Step* aStep,
				       G4TouchableHistory* ){

  if(aStep->GetTrack()->GetDefinition() 
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

//G4VPhysicalVolume* physVol=aStep->GetPostStepPoint()->GetTouchable()->GetVolume();
G4int pmtNumber1=aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);

//G4cout << "  pmtNumber1: " << pmtNumber1 << G4endl;   
// mpdPMTHit* hit3 = new mpdPMTHit();
 mpdPMTHit* hit4 = new mpdPMTHit();

//       hit3->IncPhotonCount(); 
//       hit3->SetPMTNumber(pmtNumber1);

G4String particleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();

G4String volume = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
if (volume=="photocathode")
{ 

   G4int pmtNumber2=aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(2);
   G4double edep = aStep->GetTotalEnergyDeposit();
   G4double pmtCharge = 0;
   G4double pmtTime = -111;
//    hit3->SetEdep(edep);
//   G4cout << "edep = " << edep<< G4endl;

   if (edep>1.*CLHEP::eV)
     {
//       hit3->IncPhotonAbs();
//       hit3->SetPMTNumber(pmtNumber2);

       pmtCharge = CLHEP::RandGauss::shoot(1.62e-13,3.24e-14); // Fator de 10^6 de ganho da PMT
       pmtTime=aStep->GetPostStepPoint()->GetGlobalTime();
       hit4->SetEdep(edep);
       hit4->SetPMTNumber(pmtNumber2);
       hit4->SetPMTTime(pmtTime);
       hit4->SetPMTCharge(pmtCharge);
       pmtSignalHitCollection->insert(hit4);
//G4cout << "pmtNumber1: " << pmtNumber1 << "  pmtNumber2: " << pmtNumber2 << G4endl;   
//       G4cout << "pmtTime = " << pmtTime/ns << G4endl ;

//       G4cout << "pmtCharge = " << pmtCharge << G4endl ;
    
     }

}

//  pmtHitCollection->insert(hit3);
  return true;
}

//=============================================================================
void mpdPMTSD::EndOfEvent(G4HCofThisEvent* )
{}

//=============================================================================
void mpdPMTSD::clear(){
}

//=============================================================================
void mpdPMTSD::DrawAll(){
} 

//=============================================================================
void mpdPMTSD::PrintAll(){
} 


