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
/// \file mpdCalorHit.hh
/// \brief Definition of the mpdCalorHit class

#ifndef mpdCalorHit_h
#define mpdCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class mpdCalorHit : public G4VHit
{
  public:
    mpdCalorHit();
    mpdCalorHit(const mpdCalorHit&);
    virtual ~mpdCalorHit();

    // operators
    const mpdCalorHit& operator=(const mpdCalorHit&);
    G4bool operator==(const mpdCalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double de, G4double dl, G4int layer);
    void AddPionDecay();
    void AddMuonDecay();
    void AddPionPassed();
    void AddPionCapture();

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;
    G4int GetPionDecay() const;
    G4int GetMuonDecay() const;
    G4int GetPionPassed() const; 
    G4int GetPionCapture() const;
    G4int GetLayerID() const;
    
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength; ///< Track length in the  sensitive volume
    G4int fPionDecay;
    G4int fPionPassed; 
    G4int fPionCapture; 
    G4int fMuonDecay;
    G4int fLayerID;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using mpdCalorHitsCollection = G4THitsCollection<mpdCalorHit>;

extern G4ThreadLocal G4Allocator<mpdCalorHit>* mpdCalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* mpdCalorHit::operator new(size_t)
{
  if (!mpdCalorHitAllocator) {
    mpdCalorHitAllocator = new G4Allocator<mpdCalorHit>;
  }
  void *hit;
  hit = (void *) mpdCalorHitAllocator->MallocSingle();
  return hit;
}

inline void mpdCalorHit::operator delete(void *hit)
{
  if (!mpdCalorHitAllocator) {
    mpdCalorHitAllocator = new G4Allocator<mpdCalorHit>;
  }
  mpdCalorHitAllocator->FreeSingle((mpdCalorHit*) hit);
}

inline void mpdCalorHit::Add(G4double de, G4double dl, G4int layer) {
  fEdep += de;
  fTrackLength += dl;
  fLayerID = layer;
}
inline void mpdCalorHit::AddPionDecay() {
  fPionDecay++;
}
inline void mpdCalorHit::AddMuonDecay() {
  fMuonDecay++;
}  
inline void mpdCalorHit::AddPionPassed() {
  fPionPassed++; 
}
inline void mpdCalorHit::AddPionCapture() {
  fPionCapture++; 
}
inline G4double mpdCalorHit::GetEdep() const {
  return fEdep;
}
inline G4double mpdCalorHit::GetTrackLength() const {
  return fTrackLength;
}
inline G4int mpdCalorHit::GetPionDecay() const {
  return fPionDecay;
}
inline G4int mpdCalorHit::GetPionPassed() const {
  return fPionPassed; 
}
inline G4int mpdCalorHit::GetPionCapture() const {
  return fPionCapture; 
}
inline G4int mpdCalorHit::GetMuonDecay() const {
  return fMuonDecay;
}
inline G4int mpdCalorHit::GetLayerID() const {
 return fLayerID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
