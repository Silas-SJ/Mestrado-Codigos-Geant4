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
/// \file mpdActionInitialization.cc
/// \brief Implementation of the mpdActionInitialization class

#include "mpdActionInitialization.hh"
#include "mpdPrimaryGeneratorAction.hh"
#include "mpdDetectorConstruction.hh"
#include "mpdRunAction.hh"
#include "mpdEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdActionInitialization::mpdActionInitialization(mpdDetectorConstruction* da)
//mpdActionInitialization::mpdActionInitialization(mpdPrimaryGeneratorAction* ga)
: G4VUserActionInitialization(), DetConst(da)
//: PrimGen(ga)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mpdActionInitialization::~mpdActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdActionInitialization::BuildForMaster() const
{
  SetUserAction(new mpdRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mpdActionInitialization::Build() const
{
//  SetUserAction(new mpdPrimaryGeneratorAction(DetConst));
 // SetUserAction(PrimGen);
  auto PrimGen = new mpdPrimaryGeneratorAction(DetConst);
  SetUserAction(PrimGen);
  SetUserAction(new mpdRunAction);
  SetUserAction(new mpdEventAction(PrimGen));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
