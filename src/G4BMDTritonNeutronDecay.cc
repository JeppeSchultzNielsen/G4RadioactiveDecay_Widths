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

#include "G4BMDTritonNeutronDecay.hh"
#include "G4BetaDecayCorrections.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>
#include "CoulombFunctions.hh"

G4BMDTritonNeutronDecay::G4BMDTritonNeutronDecay(const G4ParticleDefinition* theParentNucleus,
                                   const G4double& branch, const G4double& e0,
                                   const G4double& excitationE,
                                   const G4Ions::G4FloatLevelBase& flb,
                                   const G4BetaDecayType& aBetaType,
                                   const G4double aResonance1,
                                   const G4double aWidth1,
                                   const G4int aL1,
                                   const G4double aResonance2,
                                   const G4double aWidth2,
                                   const G4int aL2,
                                   const G4double aEndpointExcitation)
        : G4NuclearDecay("BMDTritonNeutron", BMDTritonNeutron, excitationE, flb), endpointEnergy(e0), resonance1(aResonance1),
          width1(aWidth1), l1(aL1), resonance2(aResonance2),
          width2(aWidth2), l2(aL2), endpointExcitation(aEndpointExcitation), betaType(aBetaType)
{
    SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent
    SetBR(branch);

    SetNumberOfDaughters(5);
    theIonTable =
            (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
    daughterZ = theParentNucleus->GetAtomicNumber();
    daughterA = theParentNucleus->GetAtomicMass() - 4;
    transferNucleus1 = theIonTable->GetIon(daughterZ+1, daughterA+4, 0, flb); //nucleus where 1st resonance exists; after beta- decay
    transferNucleus2 = theIonTable->GetIon(daughterZ, daughterA+1, 0, flb); //nucleus where 2nd resonance exists; after triton decay
    SetDaughter(0, theIonTable->GetIon(daughterZ, daughterA, endpointExcitation, flb) );
    SetDaughter(1, "e-");
    SetDaughter(2, "anti_nu_e");
    SetDaughter(3, "triton");
    SetDaughter(4, "neutron");

    parentMass = theParentNucleus -> GetPDGMass();
    origTransferMass1 = transferNucleus1 -> GetPDGMass();
    origTransferMass2 = transferNucleus2 -> GetPDGMass();
    transferA1 = daughterA + 4;
    transferA2 = daughterA+1;
    transferZ1 = daughterZ +1;
    transferZ2 = daughterZ;

    //initialize Columbfunctions calculator. R0 chosen as 1.4 fm.
    cf = new CoulombFunctions(transferA1, 3, transferZ1, 1, l1, 1.4);
    //bc = cf -> shiftFunction();


    nucleusMass = theIonTable->GetIon(daughterZ, daughterA, endpointExcitation, flb) -> GetPDGMass();
    eMass = CLHEP::electron_mass_c2;
    neutronMass = CLHEP::neutron_mass_c2;
    tritonMass = theIonTable->GetIon(1, 3, 0, flb) -> GetPDGMass();

    totalQ = parentMass - ( nucleusMass +  eMass + neutronMass + tritonMass);

    //the excitation level in the first resonance is minimal when beta-decay takes the entire Q-value. The energy of
    //the transfer nucleus then has energy of triton + neutron + nucleus. The excitation is then
    minLevel1 = nucleusMass + neutronMass + tritonMass - origTransferMass1;
    //the excitation level in the first resonance is minimal when beta-decay takes no Q-value. The energy of
    //the transfer nucleus then has energy of parentMass - eMass. The excitation is then
    maxLevel1 = parentMass - eMass - origTransferMass1;

    //the minimum possible excitation of the second resonance must be when no Q-value is given to the products of the
    //second decay, it then has energy
    minLevel2 = nucleusMass + neutronMass - origTransferMass2;
    //the maximum possible excitation is when beta and triton decays take no kinetic energy.
    maxLevel2 = parentMass - eMass - tritonMass - origTransferMass2;



    G4double minLevelSearch1 = minLevel1;
    G4double maxLevelSearch1 = maxLevel1;
    G4double resolution = 100;
    G4int noLooks = 5;
    G4double dL1 = 1.*(maxLevel1-minLevel1)/resolution;

    G4double currentLevel1 = 0;
    G4double minLevelSearch2 = minLevel2;
    G4double maxLevelSearch2 = maxLevel2;
    G4double currentLevel2 = 0;
    G4double dL2 = 1.*(maxLevelSearch2-minLevelSearch2)/resolution;

    G4double maxWeightLevels[2] = {0,0};
    maxWeight = 0;
    for(int j = 0; j < noLooks; j++){
        for(int i = 0; i < resolution; i++){
            currentLevel1 = minLevelSearch1 + i*dL1;
            for(int k = 0; k < resolution; k++){
                currentLevel2 = minLevelSearch2 + k*dL2;
                G4double w = GetWeight(currentLevel1,currentLevel2);

                /*G4cout << "level 1 " << currentLevel1 << G4endl;
                G4cout << "level 2 " << currentLevel2 << G4endl;
                G4cout << "weight " << w << G4endl;
                G4cout << "max weight " << maxWeight << G4endl << G4endl;*/

                if(w > maxWeight){
                    maxWeightLevels[0] = currentLevel1;
                    maxWeightLevels[1] = currentLevel2;
                    maxWeight = w;
                }
            }
        }
        //now search with finer resolution.
        minLevelSearch1 = maxWeightLevels[0] - dL1;
        maxLevelSearch1 = maxWeightLevels[0] + dL1;
        dL1 = 1.*(maxLevelSearch1-minLevelSearch1)/resolution;
        minLevelSearch2 = maxWeightLevels[1] - dL2;
        maxLevelSearch2 = maxWeightLevels[1] + dL2;
        dL2 = 1.*(maxLevelSearch2-minLevelSearch2)/resolution;
    }

    /*G4cout << "max level 1 " << maxWeightLevels[0] << G4endl;
    G4cout << "max level 2 " << maxWeightLevels[1] << G4endl;
    G4cout << "max weight" << maxWeight << G4endl;

    for(int i = 0; i < 1000; i++){
        G4cout << "energy: " << maxWeightLevels[0]-0.001*i << G4endl;
        G4cout << GetWeight(maxWeightLevels[0]-0.001*i,maxWeightLevels[1]) << G4endl;
    }*/
}


G4BMDTritonNeutronDecay::~G4BMDTritonNeutronDecay()
{
    delete spectrumSampler;
}

G4double G4BMDTritonNeutronDecay::GetWeight(G4double transferLevel1, G4double transferLevel2)
{
    G4double betaQ = parentMass - (origTransferMass1 + transferLevel1 + eMass);
    G4double tritonQ = origTransferMass1 + transferLevel1 - (origTransferMass2 + transferLevel2 + tritonMass);  //decay of transfer1 with excitation transferLevel1 to transfer 2 with excitation transferLevel2;
    G4double neutronQ = totalQ - betaQ - tritonQ;
    if(betaQ < 0.001 || neutronQ < 0.001 || tritonQ < 0.001) return 0; //to avoid numerical errors when calculating penetrations
    else{
        G4double tritonPen = cf -> penetrability(tritonQ);
        //G4double tritonShift = cf -> shiftFunction(tritonQ); GetBetaPhaseSpace(betaQ)*
        return GetBetaPhaseSpace(betaQ)*GetNeutronPenetrability(neutronQ,transferA2,l2)*tritonPen*1/(std::pow(transferLevel1-resonance1,2) + 1./4.*std::pow(width1,2))*1/(std::pow(transferLevel2-resonance2,2) + 1./4.*std::pow(width2,2));
    }
    return 0;
}

std::vector<G4double> G4BMDTritonNeutronDecay::GetConfiguration()
{
    G4bool found = false;
    G4double transferLevel1Conf = 0;
    G4double transferLevel2Conf = 0;
    G4double randomW = 0;
    int i = 0;
    while(!found){
        //find by rejection sampling. First, pick random number between 0 and maximum possible weight
        randomW = maxWeight * G4UniformRand();
        //then pick a transferlevel1 between the minimum possible value and maximum possible value
        transferLevel1Conf = minLevel1 + (maxLevel1-minLevel1)*G4UniformRand();
        //now pick a transferlevel2 between the minimum possible value and maximum possible value
        transferLevel2Conf = minLevel2 + (maxLevel2 - minLevel2)*G4UniformRand();
        //if the randomWeight is larger than the weight for this level, reject the sample; else, the sample is good and accept.
        if(GetWeight(transferLevel1Conf,transferLevel2Conf) > randomW){
            /*G4cout << maxWeight << G4endl;
            G4cout << transferLevel1Conf << " , " << transferLevel2Conf << G4endl;
            G4cout << GetWeight(9.6908,0.0512031) << G4endl;
            G4cout << i << G4endl;
            G4cout << GetWeight(transferLevel1Conf,transferLevel2Conf) << G4endl << G4endl;*/
            return {transferLevel1Conf,transferLevel2Conf};
        }
        i++;
    }
}



G4DecayProducts* G4BMDTritonNeutronDecay::DecayIt(G4double)
{
    std::vector<G4double> transferLevels = GetConfiguration();
    transferNucleus1 = theIonTable->GetIon(transferZ1, transferA1, transferLevels[0]);
    transferNucleus2 = theIonTable->GetIon(transferZ2, transferA2, transferLevels[1]);
    transferMass1 = transferNucleus1 -> GetPDGMass();
    transferMass2 = transferNucleus2 -> GetPDGMass();
    endpointEnergy = parentMass - (transferMass1 + eMass);
    G4double tritonQ = transferMass1 - (transferMass2 + tritonMass);
    G4double neutronQ = totalQ - endpointEnergy - tritonQ;

    //must set up new betaSpectrumSampler as betaQ-value has changed.
    SetUpBetaSpectrumSampler(transferZ1, transferA1+1, betaType);

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
    G4double nuEnergy = ((endpointEnergy - eKE)*(parentMass + transferMass1 - eTE)
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

    //we now have a beta decay to the transfer nucleus. I want its velocity so that i can decay the triton in its rest
    //frame and make a Lorentztransformation back to the rest frame of the parent nucleus.

    G4DynamicParticle dynamicTransfer1(transferNucleus1, -eDirection*eMomentum - nuDirection*nuEnergy);

    G4DecayProducts* transferProducts1 = new G4DecayProducts(dynamicTransfer1);

    //do neutron decay as in G4NeutronDecay
    G4double cmMomentum = cmMomentum = std::sqrt(tritonQ*(tritonQ + 2.*tritonMass)*
                                                 (tritonQ + 2.*transferMass2)*
                                                 (tritonQ + 2.*tritonMass + 2.*transferMass2) )/
                                       (tritonQ + tritonMass + transferMass2)/2.;

    G4double costheta = 2.*G4UniformRand()-1.0;
    G4double sintheta = std::sqrt(1.0 - costheta*costheta);
    phi  = twopi*G4UniformRand()*rad;
    G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),
                            costheta);

    G4double KE = std::sqrt(cmMomentum*cmMomentum + tritonMass*tritonMass)
                  - tritonMass;

    G4DynamicParticle* dynamicTriton =
            new G4DynamicParticle(G4MT_daughters[3], direction, KE, tritonMass);

    transferProducts1->PushProducts(dynamicTriton);

    KE = std::sqrt(cmMomentum*cmMomentum + transferMass2*transferMass2) - transferMass2;

    G4DynamicParticle* dynamicTransfer2 =
            new G4DynamicParticle(transferNucleus2, -1.0*direction, KE, transferMass2);

    transferProducts1->PushProducts(dynamicTransfer2);

    //boost products to rest frame of parent
    transferProducts1->Boost(0,0,0);

    //add the products
    G4DynamicParticle* boostedTriton = transferProducts1 -> operator[](0);
    G4DynamicParticle* boostedTransfer2 = transferProducts1 -> operator[](1);
    products->PushProducts(boostedTriton);


    //now decay of transferNucleus2.
    G4DynamicParticle dynamicTransfer2c(transferNucleus2, boostedTransfer2 -> GetMomentum());
    G4DecayProducts* transferProducts2 = new G4DecayProducts(dynamicTransfer2c);
    //do neutron decay as in G4NeutronDecay
    cmMomentum = std::sqrt(neutronQ*(neutronQ + 2.*neutronMass)*
                                    (neutronQ + 2.*nucleusMass)*
                                    (neutronQ + 2.*neutronMass + 2.*nucleusMass) )/
                          (neutronQ + neutronMass + nucleusMass)/2.;

    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt(1.0 - costheta*costheta);
    phi  = twopi*G4UniformRand()*rad;
    direction = G4ThreeVector(sintheta*std::cos(phi),sintheta*std::sin(phi),
                            costheta);

    KE = std::sqrt(cmMomentum*cmMomentum + neutronMass*neutronMass)
                  - neutronMass;

    G4DynamicParticle* dynamicNeutron =
            new G4DynamicParticle(G4MT_daughters[4], direction, KE, neutronMass);

    transferProducts2->PushProducts(dynamicNeutron);

    KE = std::sqrt(cmMomentum*cmMomentum + nucleusMass*nucleusMass) - nucleusMass;

    G4DynamicParticle* dynamicDaughter =
            new G4DynamicParticle(G4MT_daughters[0], -1.0*direction, KE, nucleusMass);

    transferProducts2->PushProducts(dynamicDaughter);

    //boost products to rest frame of parent
    transferProducts2->Boost(0,0,0);

    //add the products
    G4DynamicParticle* boostedNeutron = transferProducts2 -> operator[](0);
    G4DynamicParticle* boostedDaughter = transferProducts2 -> operator[](1);

    products->PushProducts(boostedNeutron);

    products->PushProducts(boostedDaughter);

    //for checking energy conservation
    /*G4double sumE = 0;
    G4cout << "transferlevels " << transferLevels[0] << " , " << transferLevels[1] << G4endl;
    G4cout << "neutronQ " << neutronQ << G4endl;
    G4cout << "betaQ " << endpointEnergy << G4endl;
    G4cout << "sum energy in beta-decay is " << dynamicTransfer1.GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy() << G4endl;
    G4cout << "difference is " << endpointEnergy - (dynamicTransfer1.GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy()) << G4endl;
    G4cout << "after triton decay sum energy is " << dynamicTransfer2c.GetKineticEnergy() + boostedTriton -> GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy() << G4endl;
    G4cout << "difference is " << endpointEnergy + tritonQ - (dynamicTransfer2c.GetKineticEnergy() + boostedTriton -> GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy() )<< G4endl;
    G4cout << "after neutron decay sum energy is " << boostedDaughter -> GetKineticEnergy() + boostedTriton -> GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy() + boostedNeutron -> GetKineticEnergy()<< G4endl;
    G4cout << "difference is " << endpointEnergy + tritonQ + neutronQ - (boostedDaughter -> GetKineticEnergy() + boostedTriton -> GetKineticEnergy() + dynamicElectron -> GetKineticEnergy() + dynamicNeutrino -> GetKineticEnergy() + boostedNeutron -> GetKineticEnergy()) << G4endl;
    for(int i = 0; i < 5; i++){
        G4cout << "energy of " << products-> operator[](i) -> GetParticleDefinition() -> GetParticleName() << " is " << products -> operator[](i) -> GetKineticEnergy() << G4endl;
        sumE += products -> operator[](i) -> GetKineticEnergy();
    }
    G4cout << "Total energy is " << sumE << G4endl;
    G4cout << "Total energy should be " << totalQ << G4endl;
    G4cout << "Total difference is " << sumE - totalQ << G4endl << G4endl;*/

    return products;
}


void
G4BMDTritonNeutronDecay::SetUpBetaSpectrumSampler(const G4int& daughterZ,
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


void G4BMDTritonNeutronDecay::DumpNuclearInfo()
{
    G4cout << " G4BetaMinusDecay for parent nucleus " << GetParentName() << G4endl;
    G4cout << " decays to " << GetDaughterName(0) << " , " << GetDaughterName(1)
           << " and " << GetDaughterName(2) << " with branching ratio " << GetBR()
           << "% and endpoint energy " << endpointEnergy/keV << " keV " << G4endl;
}

