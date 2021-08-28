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
/// \file mpdRunAction.cc
/// \brief Implementation of the mpdRunAction class

#include "mpdRunAction.hh"
#include "mpdAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdRunAction::mpdRunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in mpdAnalysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms
  analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
  analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
  analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
  analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("mpd", "Edep and TrackL");
  analysisManager->CreateNtupleDColumn("Eabs"); // 0
  analysisManager->CreateNtupleDColumn("Egap"); // 1
  analysisManager->CreateNtupleDColumn("Labs"); // 2
  analysisManager->CreateNtupleDColumn("Lgap"); // 3
  analysisManager->CreateNtupleIColumn("scintPiondecay");  // 4
  analysisManager->CreateNtupleIColumn("scintPionEnergy"); // 5
  analysisManager->CreateNtupleIColumn("detabsPiondecay"); // 6
  analysisManager->CreateNtupleIColumn("detabsPionEnergy"); // 7
  analysisManager->CreateNtupleIColumn("scintMuondecay");  // 8
  analysisManager->CreateNtupleIColumn("detabsMuondecay");  // 9
  analysisManager->CreateNtupleDColumn("Layer");  // 10
  analysisManager->CreateNtupleDColumn("EvID");  // 11
  analysisManager->CreateNtupleDColumn("absopiondecay");  // 12
  analysisManager->CreateNtupleDColumn("gappiondecay");  // 13
  analysisManager->CreateNtupleDColumn("theta"); // 14
  analysisManager->CreateNtupleDColumn("phi");  // 15
  analysisManager->CreateNtupleDColumn("energy");  // 16
  analysisManager->CreateNtupleDColumn("px");  // 17
  analysisManager->CreateNtupleDColumn("py");  // 18
  analysisManager->CreateNtupleDColumn("pz");  // 19
  analysisManager->CreateNtupleDColumn("x");   // 20
  analysisManager->CreateNtupleDColumn("y");   // 21
  analysisManager->CreateNtupleDColumn("z");   // 22
  analysisManager->CreateNtupleDColumn("scintPionPassed");  // 23
  analysisManager->CreateNtupleDColumn("detabsPionPassed");  // 24
  analysisManager->CreateNtupleDColumn("gapPionPassed");  // 25
  analysisManager->CreateNtupleDColumn("scintPionCapture"); //26
  analysisManager->CreateNtupleDColumn("detabsPionCapture"); // 27
  analysisManager->CreateNtupleDColumn("gap2Piondecay"); // 28
  analysisManager->FinishNtuple();
  
  // Open an output file
  //
  G4String fileName = "mpd";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdRunAction::~mpdRunAction()
{
  // For run end keep information
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdRunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
/*
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "mpd";
  analysisManager->OpenFile(fileName);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {
    G4cout << G4endl << " ----> print histograms statistic ";
    if(isMaster) {
      G4cout << "for the entire run " << G4endl << G4endl; 
    }
    else {
      G4cout << "for the local thread " << G4endl << G4endl; 
    }
    
    G4cout << " EAbs : mean = " 
       << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
    
    G4cout << " EGap : mean = " 
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
       << " rms = " 
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
    G4cout << " LAbs : mean = " 
      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;

    G4cout << " LGap : mean = " 
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
      << " rms = " 
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
  }

  // save histograms & ntuple
  // For run and don't keep information
 /* 
  analysisManager->Write();
  analysisManager->CloseFile();
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
