//============================
//Modificado em 08/02/2011
//============================

#include "mpdDigi.hh"

G4Allocator<mpdDigi> mpdDigiAllocator;

//=============================================================================

mpdDigi::mpdDigi()
:Nphotons(0)
{
  Photon_eletron = 0;
  Energy =0.;
 // G4cout << "Passa no construtor do mpdDigi" << G4endl;
}

//=============================================================================

mpdDigi::~mpdDigi()
{}

//=============================================================================

mpdDigi::mpdDigi(const mpdDigi& right)
  :G4VDigi()
{
  Photon_eletron = right.Photon_eletron;
  Energy= right.Energy;
  Nphotons= right.Nphotons;
  Stime = right.Stime;
  Scharge = right.Scharge;
//  Tvec0 = right.Tvec0;
//  Tvec1 = right.Tvec1;
  PMTvec = right.PMTvec;
  Timevec = right.Timevec;
//  Tvec = right.Tvec;
}

//=============================================================================

const mpdDigi& mpdDigi::operator=(const mpdDigi& right)
{
  Photon_eletron = right.Photon_eletron;
  Energy= right.Energy;
  Nphotons= right.Nphotons;
  Stime = right.Stime;
  Scharge = right.Scharge;
//  Tvec0 = right.Tvec0;
//  Tvec1 = right.Tvec1;
  PMTvec = right.PMTvec;
  Timevec = right.Timevec;
//  Tvec = right.Tvec;
  return *this;
}

//=============================================================================

int mpdDigi::operator==(const mpdDigi& right) const
{ 
// return ((Photon_eletron==right.Photon_eletron)&&(Energy==right.Energy)&&(Stime==right.Stime)&&(Scharge==right.Scharge)&&(Nphotons==right.Nphotons)&&(Tvec0==right.Tvec0)&&(Tvec1==right.Tvec1)&&(PMTvec==right.PMTvec)&&(Timevec==right.Timevec)&&(Tvec==right.Tvec));
 return ((Photon_eletron==right.Photon_eletron)&&(Energy==right.Energy)&&(Stime==right.Stime)&&(Scharge==right.Scharge)&&(Nphotons==right.Nphotons)&&(PMTvec==right.PMTvec)&&(Timevec==right.Timevec));
}

//=============================================================================

void mpdDigi::Draw()
{;}

//=============================================================================

void mpdDigi::Print()
{;}


