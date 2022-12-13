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
//  File:   G4NuclearDecay.cc                                                 //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   11 December 2014                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4NuclearDecay.hh"
#include "G4SystemOfUnits.hh"

G4NuclearDecay::G4NuclearDecay(const G4String& channelName,
                               const G4RadioactiveDecayMode& aMode,
                               const G4double& excitationE,
                               const G4Ions::G4FloatLevelBase& flb)
 : G4VDecayChannel(channelName), theMode(aMode), daughterEx(excitationE),
   floatingLevel(flb), halflifeThreshold(nanosecond) 
{}

G4NuclearDecay::~G4NuclearDecay()
{}

G4double G4NuclearDecay::GetNeutronPenetrability(G4double Q, G4double transferA, G4int l){
    //the energy of the emitted neutron in the CM-frame of the decaying particle is then
    // T_n = m_ion / m_n * T_ion, T_n + T_ion = neutronQ => T_n =  m_ion / m_n * (neutronQ - T_n) => T_n = m_ion / m_n * neutronQ / (1 + m_ion / m_n)
    //G4double E = Q * (transferA) / (1 + transferA); actually E is the total CM energy, that is the Q value
    G4double x = 1.4 * std::sqrt(2.*transferA/(transferA+1)*931.502*Q)*(pow(1,1./3.) + pow(transferA,1./3.))/197329.;
    G4cout << x << G4endl;
    if(l == 0){return 1;}
    if(l == 1){return x*x/(1.+x*x);}
    if(l == 2){return std::pow(x,4)/(9.+3.*std::pow(x,3)+std::pow(x,4)); }
    if(l >= 3){return std::pow(x,6)/(225.+45.*std::pow(x,2)+6*std::pow(x,4) + std::pow(x,6)); }
    else return -1;
}

G4double G4NuclearDecay::GetBetaPhaseSpace(G4double Q) {
    //total energy release in decay normalized to electron mass
    G4double e0 = (Q + 0.510998950)/0.510998950;
    //phase space factor can be approximated using
    return std::sqrt(e0*e0-1.)*(std::pow(e0,4)/30.-std::pow(e0,2)/20.-2./15.)+e0/4.*std::log(e0+std::sqrt(e0*e0-1.));
}