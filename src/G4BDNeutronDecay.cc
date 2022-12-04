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
//  File:   G4BetaMinusDecay.cc                                               //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   25 October 2014                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4BDNeutronDecay.hh"
#include "G4BetaDecayCorrections.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>
#include "CoulombFunctions.hh"

G4BDNeutronDecay::G4BDNeutronDecay(const G4ParticleDefinition* theParentNucleus,
                                   const G4double& branch, const G4double& e0,
                                   const G4double& excitationE,
                                   const G4Ions::G4FloatLevelBase& flb,
                                   const G4BetaDecayType& aBetaType,
                                   const G4double aResonance,
                                   const G4double aWidth,
                                   const G4int aL,
                                   const G4double aEndpointExcitation)
        : G4NuclearDecay("beta- decay", BetaMinus, excitationE, flb), endpointEnergy(e0), resonance(aResonance),
        width(aWidth), l(aL), endpointExcitation(aEndpointExcitation), betaType(aBetaType)
{
    SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent
    SetBR(branch);

    SetNumberOfDaughters(4);
    theIonTable =
            (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
    daughterZ = theParentNucleus->GetAtomicNumber() + 1;
    daughterA = theParentNucleus->GetAtomicMass() - 1;
    transferNucleus = theIonTable->GetIon(daughterZ, daughterA+1, 0, flb);
    SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, endpointExcitation, flb) );
    SetDaughter(1, "e-");
    SetDaughter(2, "anti_nu_e");
    SetDaughter(3, "neutron");

    parentMass = theParentNucleus -> GetPDGMass();
    transferMass = transferNucleus -> GetPDGMass();
    transferA = daughterA + 1;

    nucleusMass = theIonTable->GetIon(daughterZ, daughterA, endpointExcitation, flb) -> GetPDGMass();
    eMass = 0.510998950;
    neutronMass = 939.56537;

    totalQ = parentMass - ( nucleusMass +  eMass + neutronMass);

    //must find maximal weight. Transferlevel can be minimal nucleusMass + neutronMass - transferMass and maximal parentMass - transferMass + eMass
    minLevel = nucleusMass + neutronMass - transferMass;
    maxLevel = parentMass + eMass - transferMass;

    /*
    G4double resolution = 10000;
    G4double dL = 1.*(maxLevel-minLevel)/resolution;
    G4cout << dL << G4endl;
    //first look with resolution of 1000
    G4double maxWeightLevel = 0;
    maxWeight = 0;
    for(int i = 0; i < resolution; i++){
        G4double w = GetWeight(minLevel + i*dL);
        G4cout << minLevel + i*dL << G4endl;
        G4cout << w << G4endl;
        if(w > maxWeight){
            maxWeightLevel = minLevel + i*dL;
            maxWeight = w;
        }
    }

    G4cout << "first max level " << maxWeightLevel <<  G4endl;

    G4cout << "max Weight " << maxWeight <<  G4endl;

    G4double finerMinLevel = maxWeightLevel - dL;
    G4double finerMaxLevel = maxWeightLevel + dL;
    dL = 1.*(finerMaxLevel-finerMinLevel)/resolution;
    for(int i = 0; i < resolution; i++){
        G4double w = GetWeight(finerMinLevel + i*dL);
        G4cout << finerMinLevel + i*dL << G4endl;
        G4cout << w << G4endl;
        if(w > maxWeight){
            maxWeightLevel = finerMinLevel + i*dL;
            maxWeight = w;
        }
    }*/

    G4double minLevelSearch = minLevel;
    G4double maxLevelSearch = maxLevel;
    G4double resolution = 100;
    G4int noLooks = 5;
    G4double dL = 1.*(maxLevel-minLevel)/resolution;

    G4double maxWeightLevel = 0;
    maxWeight = 0;
    for(int j = 0; j < noLooks; j++){
        for(int i = 0; i < resolution; i++){
            G4double w = GetWeight(minLevelSearch + i*dL);
            if(w > maxWeight){
                maxWeightLevel = minLevelSearch + i*dL;
                maxWeight = w;
            }
        }
        //now search with finer resolution.
        minLevelSearch = maxWeightLevel - dL;
        maxLevelSearch = maxWeightLevel + dL;
        dL = 1.*(maxLevelSearch-minLevelSearch)/resolution;
    }

    SetUpBetaSpectrumSampler(daughterZ, daughterA+1, betaType);
}


G4BDNeutronDecay::~G4BDNeutronDecay()
{
    delete spectrumSampler;
}

G4double G4BDNeutronDecay::GetWeight(G4double transferLevel)
{
    G4double betaQ = parentMass - (transferMass + transferLevel + eMass);
    G4double neutronQ = totalQ - betaQ;
    if(betaQ < 0 || neutronQ < 0) return 0;
    else{
        return GetBetaPhaseSpace(betaQ)*GetNeutronPenetrability(neutronQ,transferA,l)*1/(std::pow(transferLevel-resonance,2) + 1./4.*std::pow(width,2));
    }
}

G4double G4BDNeutronDecay::GetConfiguration()
{
    G4bool found = false;
    G4double transferLevel = 0;
    G4double randomW = 0;
    int i = 0;
    while(!found){
        //find by rejection sampling. First, pick random number between 0 and maximum possible weight
        randomW = maxWeight * G4UniformRand();
        //then pick a transferlevel between the minimum possible value and maximum possible value
        transferLevel = minLevel + (maxLevel-minLevel)*G4UniformRand();
        //if the randomWeight is larger than the weight for this level, reject the sample; else, the sample is good and accept.
        if(GetWeight(transferLevel) > randomW){
            //G4cout << transferLevel << G4endl;
            //G4cout << i << G4endl;
            return transferLevel;
        }
        i++;
    }
}



G4DecayProducts* G4BDNeutronDecay::DecayIt(G4double)
{
    G4double transferLevel = GetConfiguration();
    transferNucleus = theIonTable->GetIon(daughterZ, daughterA+1, transferLevel);
    transferMass = transferNucleus -> GetPDGMass();
    endpointEnergy = parentMass - (transferMass + eMass);
    G4double neutronQ = totalQ - endpointEnergy;

    //must set up new betaSpectrumSampler as betaQ-value has changed.
    SetUpBetaSpectrumSampler(daughterZ, daughterA+1, betaType);

    //first do a beta decay as done in original beta-decay classes with modified Q-value
    // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)
    CheckAndFillParent();

    // Fill G4MT_daughters with e-, nu and residual nucleus (stored by SetDaughter)
    CheckAndFillDaughters();

    // Set up final state
    // parentParticle is set at rest here because boost with correct momentum
    // is done later
    G4DynamicParticle parentParticle(G4MT_parent, G4ThreeVector(0,0,0), 0.0);
    G4DecayProducts* products = new G4DecayProducts(parentParticle);

    // Electron, neutrino and daughter nucleus energies
    G4double eKE = endpointEnergy*spectrumSampler->shoot(G4Random::getTheEngine() );
    G4double eMomentum = std::sqrt(eKE*(eKE + 2.*eMass) );

    G4double cosThetaENu = 2.*G4UniformRand() - 1.;
    G4double eTE = eMass + eKE;
    G4double nuEnergy = ((endpointEnergy - eKE)*(parentMass + transferMass - eTE)
                         - eMomentum*eMomentum)/(parentMass - eTE + eMomentum*cosThetaENu)/2.;

    // Electron 4-vector, isotropic angular distribution
    G4double cosTheta = 2.*G4UniformRand() - 1.0;
    G4double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);

    G4double phi = twopi*G4UniformRand()*rad;
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);

    G4ParticleMomentum eDirection(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
    G4DynamicParticle* dynamicElectron
            = new G4DynamicParticle(G4MT_daughters[1], eDirection*eMomentum);
    products->PushProducts(dynamicElectron);

    // Neutrino 4-vector
    G4double sinThetaENu = std::sqrt(1.0 - cosThetaENu*cosThetaENu);
    phi = twopi*G4UniformRand()*rad;
    G4double sinPhiNu = std::sin(phi);
    G4double cosPhiNu = std::cos(phi);

    G4ParticleMomentum nuDirection;
    nuDirection.setX(sinThetaENu*cosPhiNu*cosTheta*cosPhi -
                     sinThetaENu*sinPhiNu*sinPhi + cosThetaENu*sinTheta*cosPhi);
    nuDirection.setY(sinThetaENu*cosPhiNu*cosTheta*sinPhi +
                     sinThetaENu*sinPhiNu*cosPhi + cosThetaENu*sinTheta*sinPhi);
    nuDirection.setZ(-sinThetaENu*cosPhiNu*sinTheta + cosThetaENu*cosTheta);

    G4DynamicParticle* dynamicNeutrino
            = new G4DynamicParticle(G4MT_daughters[2], nuDirection*nuEnergy);
    products->PushProducts(dynamicNeutrino);

    //we now have a beta decay to the transfer nucleus. I want its velocity so that i can decay the neutron in its rest
    //frame and make a Lorentztransformation back to the rest frame of the parent nucleus.

    G4DynamicParticle dynamicTransfer(transferNucleus, -eDirection*eMomentum - nuDirection*nuEnergy);

    /*G4DynamicParticle* dynamicTransfer =
            new G4DynamicParticle(transferNucleus,
                                  -eDirection*eMomentum - nuDirection*nuEnergy);*/
    G4ThreeVector transferDirection = dynamicTransfer.GetMomentumDirection();

    G4DecayProducts* transferProducts = new G4DecayProducts(dynamicTransfer);

    //do neutron decay as in G4NeutronDecay
    G4double cmMomentum = std::sqrt(neutronQ*(neutronQ + 2.*neutronMass)*
                                    (neutronQ + 2.*nucleusMass)*
                                    (neutronQ + 2.*neutronMass + 2.*nucleusMass) )/
                          (neutronQ + neutronMass + nucleusMass)/2.;

    G4double costheta = 2.*G4UniformRand()-1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    phi  = twopi*G4UniformRand()*rad;
    G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                            costheta);

    G4double KE = std::sqrt(cmMomentum*cmMomentum + neutronMass*neutronMass)
                  - neutronMass;

    G4DynamicParticle* dynamicNeutron =
            new G4DynamicParticle(G4MT_daughters[3], direction, KE, neutronMass);

    transferProducts->PushProducts(dynamicNeutron);

    KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;

    G4DynamicParticle* dynamicDaughter =
            new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);

    transferProducts->PushProducts(dynamicDaughter);

    //boost products to rest frame of parent
    transferProducts->Boost(0,0,0);

    //add the products
    G4DynamicParticle* boostedNeutron = transferProducts -> operator[](0);
    G4DynamicParticle* boostedDaughter = transferProducts -> operator[](1);

    products->PushProducts(boostedNeutron);

    products->PushProducts(boostedDaughter);

    //for checking energy conservation
    /*G4double sumE = 0;
    G4cout << "transferlevel " << transferLevel << G4endl;
    G4cout << "neutronQ " << neutronQ << G4endl;
    G4cout << "betaQ " << endpointEnergy << G4endl;
    G4cout << "sum energy in beta-decay is " << dynamicTransfer.GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy() << G4endl;
    G4cout << "difference is " << endpointEnergy - (dynamicTransfer.GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy()) << G4endl;
    for(int i = 0; i < 4; i++){
        G4cout << "energy of " << products-> operator[](i) -> GetParticleDefinition() -> GetParticleName() << " is " << products -> operator[](i) -> GetKineticEnergy() << G4endl;
        sumE += products -> operator[](i) -> GetKineticEnergy();
    }
    G4cout << "Total energy is " << sumE << G4endl;
    G4cout << "Total energy should be " << totalQ << G4endl;
    G4cout << "Total difference is " << sumE - totalQ << G4endl << G4endl;
*/

    return products;
}


void
G4BDNeutronDecay::SetUpBetaSpectrumSampler(const G4int& daughterZ,
                                           const G4int& daughterA,
                                           const G4BetaDecayType& betaType)
{
    G4double e0 = endpointEnergy/CLHEP::electron_mass_c2;
    G4BetaDecayCorrections corrections(daughterZ, daughterA);
    spectrumSampler = 0;

    if (e0 > 0) {
        // Array to store spectrum pdf
        G4int npti = 100;
        G4double* pdf = new G4double[npti];

        G4double e;  // Total electron energy in units of electron mass
        G4double p;  // Electron momentum in units of electron mass
        G4double f;  // Spectral shape function
        for (G4int ptn = 0; ptn < npti; ptn++) {
            // Calculate simple phase space
            e = 1. + e0*(G4double(ptn) + 0.5)/G4double(npti);
            p = std::sqrt(e*e - 1.);
            f = p*e*(e0 - e + 1.)*(e0 - e + 1.);

            // Apply Fermi factor to get allowed shape
            f *= corrections.FermiFunction(e);

            // Apply shape factor for forbidden transitions
            f *= corrections.ShapeFactor(betaType, p, e0-e+1.);
            pdf[ptn] = f;
        }
        spectrumSampler = new G4RandGeneral(pdf, npti);
        delete[] pdf;
    }
}


void G4BDNeutronDecay::DumpNuclearInfo()
{
    G4cout << " G4BetaMinusDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0) << " , " << GetDaughterName(1)
           << " and " << GetDaughterName(2) << " with branching ratio " << GetBR()
           << "% and endpoint energy " << endpointEnergy/keV << " keV " << G4endl;
}

