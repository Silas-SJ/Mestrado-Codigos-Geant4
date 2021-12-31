//
//  mpdPMTHit.cc
//
//
//  Created by Helio Nogima on 7/11/21.
//

#include "mpdPMTHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"  
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

G4Allocator<mpdPMTHit> mpdPMTHitAllocator;

//=============================================================================
mpdPMTHit::mpdPMTHit()
:G4VHit(),pmtNumber(-1),photons(0),photons_abs(0),physVol(0),drawit(false)
{}

//=============================================================================
mpdPMTHit::~mpdPMTHit()
{}

//=============================================================================
mpdPMTHit::mpdPMTHit(const mpdPMTHit &right)
  : G4VHit()
{
  pmtNumber=right.pmtNumber;
  photons=right.photons;
  photons_abs=right.photons_abs;
  physVol=right.physVol;
  edep=right.edep;
  pmtTime=right.pmtTime;
  pmtCharge=right.pmtCharge;
  drawit=right.drawit;
}

//============================================================================= 
const mpdPMTHit& mpdPMTHit::operator=(const mpdPMTHit &right){
  pmtNumber= right.pmtNumber;
  photons=right.photons;
  photons_abs=right.photons_abs;
  physVol=right.physVol;
  edep=right.edep;
  pmtTime=right.pmtTime;
  pmtCharge=right.pmtCharge;
  drawit=right.drawit;
  return *this;
}

//=============================================================================
G4int mpdPMTHit::operator==(const mpdPMTHit &right) const{
  return (pmtNumber==right.pmtNumber);
}

//=============================================================================
void mpdPMTHit::Draw(){
  if(drawit&&physVol){ 
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager){
      G4VisAttributes attribs(G4Colour(1.,0.,0.));
      attribs.SetForceSolid(true);
      G4RotationMatrix rot;
      if(physVol->GetRotation()) rot=*(physVol->GetRotation());
      G4Transform3D trans(rot,physVol->GetTranslation());
      pVVisManager->Draw(*physVol,attribs,trans);
    }
  }
}

//=============================================================================
void mpdPMTHit::Print(){
}


