//
//  mpdPMTHit.hh
//
//
//  Created by Helio Nogima on 7/11/21.
//

#ifndef mpdPMTHit_h
#define mpdPMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh" 
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Threading.hh"

class G4VTouchable;

class mpdPMTHit : public G4VHit
{
public:
  
  mpdPMTHit();
  ~mpdPMTHit();
  mpdPMTHit(const mpdPMTHit &right);

  const mpdPMTHit& operator=(const mpdPMTHit &right);
  G4int operator==(const mpdPMTHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();

//    inline void SetDrawit(G4bool b){drawit=b;}
//    inline G4bool GetDrawit(){return drawit;}
    
    void IncPhotonCount();
    void IncPhotonAbs();
    void SetEdep(G4double de);
    void AddEdep(G4double de);
    void SetPMTNumber(G4int n);
    void SetPMTPhysVol(G4VPhysicalVolume* physV);
    void SetPMTPos(G4double x,G4double y,G4double z);
    void SetPMTTime (G4double sctime);
    void SetPMTCharge (G4double sccharge);
    void AddPMTCharge (G4double sccharge);
    
    G4int GetPhotonCount() const;
    G4int GetPhotonAbs() const;
    G4double GetEdep() const;
    G4int GetPMTNumber() const;
    G4VPhysicalVolume* GetPMTPhysVol() const;
    G4ThreeVector GetPMTPos() const;
    G4double  GetPMTTime () const;
    G4double GetPMTCharge () const;
    
private:
  G4int pmtNumber;
  G4int photons;
  G4int photons_abs;
  G4ThreeVector pos;
  G4VPhysicalVolume* physVol;
  G4bool drawit;
  G4double edep;
  G4double pmtTime;
  G4double pmtCharge;
};

typedef G4THitsCollection<mpdPMTHit> mpdPMTHitsCollection;
extern G4Allocator<mpdPMTHit> mpdPMTHitAllocator;

inline void* mpdPMTHit::operator new(size_t){
  void *aHit;
  aHit = (void *) mpdPMTHitAllocator.MallocSingle();
  return aHit;
}

inline void mpdPMTHit::operator delete(void *aHit){
  mpdPMTHitAllocator.FreeSingle((mpdPMTHit*) aHit);
}

inline void mpdPMTHit::IncPhotonCount(){photons++;}
inline G4int mpdPMTHit::GetPhotonCount() const {return photons;}

inline void mpdPMTHit::IncPhotonAbs(){photons_abs++;}
inline G4int mpdPMTHit::GetPhotonAbs() const {return photons_abs;}

inline void mpdPMTHit::SetEdep(G4double de){ edep = de; }
inline void mpdPMTHit::AddEdep(G4double de){ edep += de; }
inline G4double mpdPMTHit::GetEdep() const { return edep; }

inline void mpdPMTHit::SetPMTNumber(G4int n) { pmtNumber = n; }
inline G4int mpdPMTHit::GetPMTNumber() const { return pmtNumber; }

inline void mpdPMTHit::SetPMTPhysVol(G4VPhysicalVolume* physV){this->physVol=physV;}
inline G4VPhysicalVolume* mpdPMTHit::GetPMTPhysVol() const {return physVol;}

inline void mpdPMTHit::SetPMTPos(G4double x,G4double y,G4double z){
  pos=G4ThreeVector(x,y,z);
}

inline G4ThreeVector mpdPMTHit::GetPMTPos() const {return pos;}

inline void mpdPMTHit::SetPMTTime (G4double sctime) { pmtTime = sctime; }
inline G4double  mpdPMTHit::GetPMTTime() const { return pmtTime ;}

inline void mpdPMTHit::SetPMTCharge (G4double sccharge) { pmtCharge = sccharge; }
inline void mpdPMTHit::AddPMTCharge (G4double sccharge) { pmtCharge+= sccharge; }
inline G4double mpdPMTHit::GetPMTCharge() const { return pmtCharge; }

#endif


