//
//  mpdPMTSD.hh
//
//
//  Created by Helio Nogima on 7/11/21.
//

#ifndef mpdPMTSD_h
#define mpdPMTSD_h 1

#include "G4DataVector.hh"
#include "G4VSensitiveDetector.hh"
#include "mpdPMTHit.hh"

class G4Step;
class G4HCofThisEvent;

class mpdPMTSD : public G4VSensitiveDetector
{

public:
  mpdPMTSD(G4String name);
  ~mpdPMTSD();
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  

  G4bool ProcessHits_constStep(const G4Step* aStep,
			       G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  

  inline void InitPMTs(G4int nPMTs){
    if(pmtPositionsX)delete pmtPositionsX;
    if(pmtPositionsY)delete pmtPositionsY;
    if(pmtPositionsZ)delete pmtPositionsZ;
    pmtPositionsX=new G4DataVector(nPMTs);
    pmtPositionsY=new G4DataVector(nPMTs);
    pmtPositionsZ=new G4DataVector(nPMTs);
  }


  inline void SetPMTPos(G4int n,G4double x,G4double y,G4double z){
    if(pmtPositionsX)pmtPositionsX->insertAt(n,x);
    if(pmtPositionsY)pmtPositionsY->insertAt(n,y);
    if(pmtPositionsZ)pmtPositionsZ->insertAt(n,z);
  }
  
private:
  mpdPMTHitsCollection* pmtHitCollection;
  mpdPMTHitsCollection* pmtSignalHitCollection;

  
  G4DataVector* pmtPositionsX;
  G4DataVector* pmtPositionsY;
  G4DataVector* pmtPositionsZ;
};

#endif
