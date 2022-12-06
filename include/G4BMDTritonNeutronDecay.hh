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

#ifndef G4BMDTritonNeutronDecay_h
#define G4BMDTritonNeutronDecay_h 1

#include "G4NuclearDecay.hh"
#include "G4BetaDecayType.hh"
#include "Randomize.hh"
#include "CoulombFunctions.hh"


class G4BMDTritonNeutronDecay : public G4NuclearDecay
{
public:
    G4BMDTritonNeutronDecay(const G4ParticleDefinition* theParentNucleus,
                      const G4double& theBR, const G4double& endpointE,
                      const G4double& ex, const G4Ions::G4FloatLevelBase& flb,
                      const G4BetaDecayType& type,
                      const G4double aResonance1,
                      const G4double aWidth1,
                      const G4int aL1,
                      const G4double aResonance2,
                      const G4double aWidth2,
                      const G4int aL2,
                      const G4double aEndpointExcitation);

    virtual ~G4BMDTritonNeutronDecay();

    virtual G4DecayProducts* DecayIt(G4double);

    virtual G4double GetWeight(G4double transferLevel1, G4double transferLevel2);

    virtual std::vector<G4double> GetConfiguration();

    virtual void DumpNuclearInfo();

private:
    void SetUpBetaSpectrumSampler(const G4int& parentZ, const G4int& parentA,
                                  const G4BetaDecayType& type);

    G4double endpointEnergy;
    const G4double resonance1;
    const G4double width1;
    const G4int l1;
    const G4double resonance2;
    const G4double width2;
    const G4int l2;
    const G4double endpointExcitation;
    G4RandGeneral* spectrumSampler;
    G4ParticleDefinition* transferNucleus1;
    G4ParticleDefinition* transferNucleus2;
    G4double transferMass1;
    G4double transferMass2;
    G4double origTransferMass1;
    G4double origTransferMass2;
    G4double bc;
    G4int transferA1;
    G4int transferA2;
    G4int transferZ1;
    G4int transferZ2;
    G4double nucleusMass;
    G4double parentMass;
    G4double tritonMass;
    G4double totalQ;
    G4double eMass;
    G4double neutronMass;
    G4double minLevel1;
    G4double maxLevel1;
    G4double minLevel2;
    G4double maxLevel2;
    G4double maxWeight;
    G4int daughterZ;
    G4int daughterA;
    const G4BetaDecayType betaType;
    G4IonTable* theIonTable;
    CoulombFunctions* cf;
};
#endif

