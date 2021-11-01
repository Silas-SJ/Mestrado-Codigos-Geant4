//
//  TreeManager.cc
//  
//
//  Created by Helio Nogima on 7/11/21.
//

#include "mpdTreeManager.hh"
#include "G4UnitsTable.hh"
#include <TTree.h>
#include <TFile.h>

mpdTreeManager::mpdTreeManager()
:fRootFile(0),
 fTree1(0),
 fEabs(0.), fEgap(0.) ,fLabs(0.), fLgap(0.)
{
  fTree1 = 0;
}

mpdTreeManager::~mpdTreeManager()
{
  if (fRootFile) delete fRootFile;
}

void mpdTreeManager::Book()
{
    G4String fileName = "mpdTree.root";
    fRootFile = new TFile(fileName,"RECREATE");
    if (! fRootFile) {
      G4cout << " mpdTreeManager::Book :"
             << " problem creating the ROOT TFile "
             << G4endl;
      return;
    }
    fTree1 = new TTree("Ntuple1", "Edep");
    
    fTree1->Branch("Theta", &fTheta, "Theta/D");
    fTree1->Branch("Phi", &fPhi, "Phi/D");
    fTree1->Branch("Energy",&fEnergy,"Energy/D");
    fTree1->Branch("Px",&fPx,"Px/D");
    fTree1->Branch("Py",&fPy,"Py/D");
    fTree1->Branch("Pz",&fPz,"Pz/D");
    fTree1->Branch("X",&fX,"X/D");
    fTree1->Branch("Y",&fY,"Y/D");
    fTree1->Branch("Z",&fZ,"Z/D");

    fTree1->Branch("Eabs", &fEabs, "Eabs/D");
    fTree1->Branch("Egap", &fEgap, "Egap/D");
    fTree1->Branch("Labs", &fLabs, "Labs/D");
    fTree1->Branch("Lgap", &fLgap, "Lgap/D");
    
    fTree1->Branch("PionDecay", &fPionDecay, "PionDecay/I");
    fTree1->Branch("MuonDecay", &fMuonDecay, "MuonDecay/I");
    fTree1->Branch("PionPassed", &fPionPassed, "PionPassed/I");
    fTree1->Branch("MuonPassed", &fMuonPassed, "MuonPassed/I");
    
    fTree1->Branch("Npe",&fNpe,"Npe/I");
//    fTree1->Branch("PMTCharge1",&fPMTCharge1,"PMTCharge1[12]/D");
//    fTree1->Branch("PMTCharge2",&fPMTCharge2,"PMTCharge2[12]/D");
//    fTree1->Branch("PMTMediumTime1",&fPMTMediumTime1,"PMTMediumTime1[12]/D");
//    fTree1->Branch("PMTMediumTime2",&fPMTMediumTime2,"PMTMediumTime2[12]/D");
//    fTree1->Branch("ScintCharge",&fScintCharge,"ScintCharge[12]/D");
//    fTree1->Branch("OTTimeVector0","std::vector<G4double>",&fOTTimeVector0);
//    fTree1->Branch("OTTimeVector1","std::vector<G4double>",&fOTTimeVector1);
//    fTree1->Branch("OTTimeVector","std::map<G4int,std::vector<G4double>>",&fOTTimeVector);
    fTree1->Branch("OTpmt","std::vector<int>",&fOTpmt);
    fTree1->Branch("OTtime","std::vector<double>",&fOTtime);
//    fTree1->Branch("OTV","std::vector<pair<G4int,G4double>>",&fOTV);
}

void mpdTreeManager::Save()
{
  if (! fRootFile) return;
  
  fRootFile->Write();       // Writing to the file
  fRootFile->Close();       // and closing the tree (and the file)
  
  G4cout << "\n----> Tree saved\n" << G4endl;
}

//void mpdTreeManager::FillNtuple(G4double energyAbs, G4double energyGap, G4double trackLAbs, G4double trackLGap, G4int Npe, G4double PMTCharge1[12], G4double PMTCharge2[12], G4double PMTMediumTime1[12], G4double PMTMediumTime2[12], G4double ScintCharge[12])
//void mpdTreeManager::FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::vector<G4double> ottv0, std::vector<G4double> ottv1)
//void mpdTreeManager::FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,std::vector<G4double>> ottv)
void mpdTreeManager::FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::vector<G4int> otp, std::vector<G4double> ott)
//void mpdTreeManager::FillNtuple(G4double theta, G4double phi, G4double energy, G4double px, G4double py, G4double pz, G4double x, G4double y, G4double z, std::map<G4int,std::vector<G4double>> ottv)
{
fTheta = theta;
fPhi = phi;
fEnergy = energy;
fPx = px;
fPy = py;
fPz = pz;
fX = x;
fY = y;
fZ = z;
//fOTTimeVector0 = ottv0;
//fOTTimeVector1 = ottv1;
//fOTTimeVector = ottv;
fOTpmt=otp;
fOTtime=ott;
/*
    fEabs = energyAbs;
    fEgap = energyGap;
    fLabs = trackLAbs;
    fLgap = trackLGap;
    fNpe =Npe;
    for (int i = 0; i<12 ; i++){
        fPMTCharge1[i]=PMTCharge1[i];
        fPMTCharge2[i]=PMTCharge2[i];
        fPMTMediumTime1[i]=PMTMediumTime1[i];
        fPMTMediumTime2[i]=PMTMediumTime2[i];
        fScintCharge[i]=ScintCharge[i];
    }
*/
//    fPMTCharge1 = PMTCharge1;
//    fPMTCharge2 = PMTCharge2;
//    fPMTMediumTime1 = PMTMediumTime1;
//    fPMTMediumTime2 = PMTMediumTime2;
//    fScintCharge = ScintCharge;
 /*
    fTree1->SetBranchAddress("Theta", &fTheta);
    fTree1->SetBranchAddress("Theta", &fPhi);
    fTree1->SetBranchAddress("Energy",&fEnergy);
    fTree1->SetBranchAddress("Px",&fPx);
    fTree1->SetBranchAddress("Py",&fPy);
    fTree1->SetBranchAddress("Pz",&fPz);
    fTree1->SetBranchAddress("X",&fX);
    fTree1->SetBranchAddress("Y",&fY);
    fTree1->SetBranchAddress("Z",&fZ);
    
    fTree1->SetBranchAddress("Eabs", &fEabs);
    fTree1->SetBranchAddress("Egap", &fEgap);
    fTree1->SetBranchAddress("Labs", &fLabs);
    fTree1->SetBranchAddress("Lgap", &fLgap);
    fTree1->SetBranchAddress("PionPassed", &fPionPassed);
    fTree1->SetBranchAddress("MuonPassed", &fMuonPassed);
    
    fTree1->SetBranchAddress("PionDecay", &fPionDecay);
    fTree1->SetBranchAddress("MuonDecay", &fMuonDecay);
    fTree1->SetBranchAddress("PionPassed", &fPionPassed);
    fTree1->SetBranchAddress("MuonPassed", &fMuonPassed);
    
    fTree1->SetBranchAddress("Npe",&fNpe);
    fTree1->SetBranchAddress("PMTCharge1",&PMTCharge1);
    fTree1->SetBranchAddress("PMTCharge2",&PMTCharge2);
    fTree1->SetBranchAddress("PMTMediumTime1",&PMTMediumTime1);
    fTree1->SetBranchAddress("PMTMediumTime2",&PMTMediumTime2);
    fTree1->SetBranchAddress("ScintCharge",&ScintCharge);
*/
  if (fTree1) fTree1->Fill();
}
