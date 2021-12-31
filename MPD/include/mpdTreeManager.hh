//
//  mpdTreeManager.hh
//  
//
//  Created by Helio Nogima on 7/11/21.
//

#ifndef mpdTreeManager_hh
#define mpdTreeManager_hh

#include <stdio.h>
#include "globals.hh"
#include <vector>

class TFile;
class TTree;

class mpdTreeManager
{
  public:
    mpdTreeManager();
    ~mpdTreeManager();

    void Book();
    void Save();
//    void FillNtuple(G4double energyAbs, G4double energyGap, G4double trackLAbs, G4double trackLGap, G4int Npe, G4double PMTCharge1[12], G4double PMTCharge2[12], G4double PMTMediumTime1[12], G4double PMTMediumTime2[12], G4double ScintCharge[12]);
//    void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::vector<G4double> ottv0, std::vector<G4double> ottv1);
//   void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,std::vector<G4double>> ottv);  
   void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,G4int> ScintPionDecay, std::map<G4int,G4int> ScintPionPassou, std::map<G4int,G4int> ScintMuonDecay, std::map<G4int,G4int> ScintMuonPassou, std::map<G4int,G4int> DetabsPionDecay, std::map<G4int,G4int> DetabsPionPassou, std::map<G4int,G4int> GapPionDecay, std::map<G4int,G4int> GapPionPassou, std::vector<G4int> ovp, std::vector<G4double> ott);  
//   void FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::vector<std::pair<int,double>> otv);  
  private:
    TFile*   fRootFile;
    TTree*   fTree1;

    G4double fTheta;
    G4double fPhi;
    G4double fEnergy;
    G4double fPx, fPy, fPz;
    G4double fX, fY, fZ;
    
    G4double fEabs;
    G4double fEgap;
    G4double fLabs;
    G4double fLgap;
    
    G4int fPionDecay;
    G4int fMuonDecay;
    G4int fPionPassed;
    G4int fMuonPassed;
              // for Scintillator ///
    std::map<G4int,G4int> fScintPionDecay;
    std::map<G4int,G4int> fScintPionPassou;
    std::map<G4int,G4int> fScintMuonDecay;
    std::map<G4int,G4int> fScintMuonPassou;
            // for Detabs ///
    std::map<G4int,G4int> fDetabsPionDecay;
    std::map<G4int,G4int> fDetabsPionPassou;  
             // for Gap ///
    std::map<G4int,G4int> fGapPionDecay;
    std::map<G4int,G4int> fGapPionPassou;
    
    std::vector <G4int> fOTpmt;
    std::vector <G4double> fOTtime;
//    std::vector <std::pair<G4int,G4double>> fOTV;
//    std::vector <G4double> fOTTimeVector0;
//    std::vector <G4double> fOTTimeVector1;
//    std::map<G4int,std::vector<G4double>> fOTTimeVector;
  
    G4int    fNpe;
//    G4double fPMTCharge1[12], fPMTCharge2[12];
//    G4double fPMTMediumTime1[12], fPMTMediumTime2[12];
//    G4double fScintCharge[12];
};

#endif /* mpdTreeManager_hh */
