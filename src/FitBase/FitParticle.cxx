// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
*    This file is part of NUISANCE.
*
*    NUISANCE is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    NUISANCE is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#include "FitParticle.h"

// NEUT Constructor
#ifdef __NEUT_ENABLED__
FitParticle::FitParticle(NeutPart* part) {

  // Set Momentum
  fP = TLorentzVector(part->fP.X(), part->fP.Y(), part->fP.Z(), part->fP.T());

  fPID = part->fPID;
  fIsAlive = part->fIsAlive;
  fStatus = part->fStatus;
  fMass = part->fMass;
};

// NEUT FSI defined in neutclass/neutfsipart
FitParticle::FitParticle(NeutFsiPart* part) {

  // Set Momentum
  fP = TLorentzVector(part->fDir.X(), part->fDir.Y(), part->fDir.Z(), part->fDir.T());

  fPID = part->fPID;
  // Set these to zero because they don't make sense in NEUT
  fIsAlive = 0;
  fStatus = 0;
  fMass = fP.Mag();
};
#endif

// NUWRO Constructor
#ifdef __NUWRO_ENABLED__
FitParticle::FitParticle(particle* nuwro_particle, Int_t state) {
  // Set Momentum
  this->fP = TLorentzVector(nuwro_particle->p4().x, nuwro_particle->p4().y,
                            nuwro_particle->p4().z, nuwro_particle->p4().t);
  fPID = nuwro_particle->pdg;

  // Set status manually from switch
  switch (state) {
    case 0:
      fIsAlive = 0;
      fStatus = 1;
      break;  // Initial State
    case 1:
      fIsAlive = 1;
      fStatus = 0;
      break;  // Final State
    case 2:
      fIsAlive = 0;
      fStatus = 2;
      break;  // Intermediate State
    default:
      fIsAlive = -1;
      fStatus = 3;
      break;  // Other?
  }

  fMass = nuwro_particle->m();
};
#endif

// GENIE Constructor
#ifdef __GENIE_ENABLED__
FitParticle::FitParticle(genie::GHepParticle* genie_particle) {
  this->fP = TLorentzVector(
      genie_particle->Px() * 1000.0, genie_particle->Py() * 1000.0,
      genie_particle->Pz() * 1000.0, genie_particle->E() * 1000.0);

  fPID = genie_particle->Pdg();

  switch (genie_particle->Status()) {
    case genie::kIStInitialState:
      fIsAlive = 0;
      fStatus = 1;
      break;  // Initial State
    case genie::kIStStableFinalState:
      fIsAlive = 1;
      fStatus = 0;
      break;  // Final State
    case genie::kIStIntermediateState:
      fIsAlive = 0;
      fStatus = 2;
      break;  // Intermediate State
    default:
      fIsAlive = -1;
      fStatus = 3;
      break;
  }

  // Flag to remove nuclear part in genie
  if (fPID > 3000) {
    fIsAlive = -1;
    fStatus = 2;
  }

  fMass = genie_particle->Mass() * 1000.0;

  // Additional flag to remove off shell particles
  if (fabs(fMass - fP.Mag()) > 0.001) {
    fIsAlive = -1;
    fStatus = 2;
  }
};

#endif

#ifdef __GiBUU_ENABLED__
FitParticle::FitParticle(GiBUUStdHepReader* GiRead, Int_t p_it) {
  fPID = GiRead->StdHepPdg[p_it];

  fP = TLorentzVector(
      GiRead->StdHepP4[p_it][0] * 1000.0, GiRead->StdHepP4[p_it][1] * 1000.0,
      GiRead->StdHepP4[p_it][2] * 1000.0, GiRead->StdHepP4[p_it][3] * 1000.0);

  fMass = fP.M();

  // Flag to remove nuclear part
  if (fPID > 100000) {
    fIsAlive = -1;
    fStatus = 2;
  }

  switch (GiRead->StdHepStatus[p_it]) {
    case 0: {  // incoming
      fIsAlive = 0;
      fStatus = 1;
      break;
    }
    case 1: {  // good final state
      fIsAlive = 1;
      fStatus = 0;
      break;
    }
    case 11: {  // struck nucleon
      fIsAlive = 0;
      fStatus = 1;
      break;
    }
  }
}
#endif

FitParticle::FitParticle(UInt_t* i) {
  (void)i;

  // A NULL event has been passed
  //  ERR(FTL)<<"NULL Event Passed to FitEvent.cxx"<<std::endl;

  return;
};

// NUANCE Particle
FitParticle::FitParticle(double x, double y, double z, double t, int pdg, Int_t state){

  // Set Momentum
  this->fP = TLorentzVector(x,
			    y,
			    z,
			    t);
  fPID = pdg;

  // Set status manually from switch
  switch(state){
  case     kInitialState: fIsAlive= 0; fStatus=1; break; // Initial State
  case     kFinalState:   fIsAlive= 1; fStatus=0; break; // Final State
  case     kFSIState:     fIsAlive= 0; fStatus=2; break; // Intermediate State
  default: fIsAlive=-1; fStatus=3; break; // Other?
  }

  fMass = fP.Mag();
};
