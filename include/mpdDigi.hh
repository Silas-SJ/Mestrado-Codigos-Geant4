

#ifndef mpdDigi_h
#define mpdDigi_h 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class mpdDigi : public G4VDigi
{

public:
  
  mpdDigi();
  ~mpdDigi();
  mpdDigi(const mpdDigi&);
  const mpdDigi& operator=(const mpdDigi&);
  int operator==(const mpdDigi&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();

private:
  
  G4double Photon_eletron;
  G4double Energy;
  G4double Stime;
  G4double Scharge;
  G4int pmt;
  G4int Nphotons;
//  std::vector <G4double> Tvec0;
//  std::vector <G4double> Tvec1;  
//  std::map <G4int,std::vector<G4double>> Tvec;
  std::vector <G4int> PMTvec;
  std::vector <G4double> Timevec;
public:
  
  inline void SetPhoton_eletron(G4double Photon_e)   {Photon_eletron = Photon_e;}
  inline void SetOpticalPhotonEnergy(G4double Ene)  {Energy = Ene;}
  inline void SetOpticalPhotonAbs(G4int nOpPh)  {Nphotons = nOpPh;}
  inline void SetPMT(G4double PMT)   {pmt = PMT;}
  inline void SetSignalTime(G4double stime)   {Stime = stime;}
  inline void SetSignalCharge(G4double scharge)   {Stime = scharge;}
//  inline void SetOverThresholdTimeVector0(std::vector <G4double> otv) {Tvec0=otv;}
//  inline void SetOverThresholdTimeVector1(std::vector <G4double> otv) {Tvec1=otv;}
//  inline void SetOverThresholdTimeVector(std::map<int,std::vector<G4double>> otv) {Tvec=otv;}
  inline void SetOverThresholdPMTVec(std::vector <G4int> otp) {PMTvec=otp;}
  inline void SetOverThresholdTimeVec(std::vector <G4double> otv) {Timevec=otv;}

  inline G4double GetPhoton_eletron() {return Photon_eletron;}
  inline G4double GetOpticalPhotonEnergy()  {return Energy;}
  inline G4int GetOpticalPhotonAbs()  {return Nphotons;}
  inline G4int GetPMT() {return pmt;};
  inline G4double GetSignalTime() {return Stime;}
  inline G4double GetSignalCharge() {return Scharge;}
//  inline std::vector <G4double> GetOverThresholdTimeVector0() {return Tvec0;}
//  inline std::vector <G4double> GetOverThresholdTimeVector1() {return Tvec1;}
//  inline std::map<G4int,std::vector<G4double>> GetOverThresholdTimeVector() {return Tvec;}
  inline std::vector <G4int> GetOverThresholdPMTVec() {return PMTvec;}
  inline std::vector <G4double> GetOverThresholdTimeVec() {return Timevec;}
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

typedef G4TDigiCollection<mpdDigi> mpdDigitsCollection;

extern G4Allocator<mpdDigi> mpdDigiAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void* mpdDigi::operator new(size_t)
{
  void* aDigi;
  aDigi = (void*) mpdDigiAllocator.MallocSingle();
  return aDigi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void mpdDigi::operator delete(void* aDigi)
{
  mpdDigiAllocator.FreeSingle((mpdDigi*) aDigi);
}

#endif
