//
//  mpdDigitizer.hh
//
//
//  Created by Helio Nogima on 7/11/21.
//

#ifndef mpdDigitizer_h
#define mpdDigitizer_h 1

#include "G4VDigitizerModule.hh"
#include "mpdDigi.hh"
#include "globals.hh"
//#include "g4std/vector"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class TH1D;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class mpdDigitizer : public G4VDigitizerModule
{
public:
  
  mpdDigitizer(G4String name);
  ~mpdDigitizer();
  
  void Digitize();
  void SetThreshold(G4double val) { Threshold = val;}
  
private:
  
  mpdDigitsCollection*  DigitsCollection;
//  mpdDigitsCollection2* DigitsColleciont2;
//  mpdDigi2* mpdDigi2;
  G4double Threshold; // for TKR digi
  G4double OpPhotonsTotalEnergy; // for CAL analysis
  G4double probList[200];
};

#endif








