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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4BetaMinusDecay.hh                                               //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   25 October 2014                                                   //
//  Description: performs beta- decay of radioactive nuclei, and returns      //
//               daughter particles in rest frame of parent nucleus           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4BMDNeutronDecay_h
#define G4BMDNeutronDecay_h 1

#include "G4NuclearDecay.hh"
#include "G4BetaDecayType.hh"
#include "Randomize.hh"
#include "CoulombFunctions.hh"


class G4BMDNeutronDecay : public G4NuclearDecay
{
public:
    G4BMDNeutronDecay(const G4ParticleDefinition* theParentNucleus,
                     const G4double& theBR, const G4double& endpointE,
                     const G4double& ex, const G4Ions::G4FloatLevelBase& flb,
                     const G4BetaDecayType& type,
                     const G4double aResonance,
                     const G4double aWidth,
                     const G4int aL,
                     const G4double aEndpointExcitation);

    virtual ~G4BMDNeutronDecay();

    virtual G4DecayProducts* DecayIt(G4double);

    virtual G4double GetWeight(G4double transferLevel);

    virtual G4double GetConfiguration();

    virtual void DumpNuclearInfo();

private:
    void SetUpBetaSpectrumSampler(const G4int& parentZ, const G4int& parentA,
                                  const G4BetaDecayType& type);

    G4double endpointEnergy;
    const G4double resonance;
    const G4double width;
    const G4int l;
    const G4double endpointExcitation;
    G4RandGeneral* spectrumSampler;
    G4double origTransferMass;
    G4int transferA;
    G4int transferZ;
    G4double nucleusMass;
    G4double parentMass;
    G4double totalQ;
    G4double eMass;
    G4double neutronMass;
    G4double minLevel;
    G4double maxLevel;
    G4double maxWeight;
    G4int daughterZ;
    G4int daughterA;
    const G4BetaDecayType betaType;
    CoulombFunctions* cfNeutron;
};
#endif

