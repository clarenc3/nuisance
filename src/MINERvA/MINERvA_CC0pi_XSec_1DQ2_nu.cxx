//Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

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

/*
  Authors: Adrian Orea (v1 2017)
           Clarence Wret (v2 2018)
*/

#include "MINERvA_SignalDef.h"
#include "MINERvA_CC0pi_XSec_1DQ2_nu.h"

//********************************************************************
MINERvA_CC0pi_XSec_1DQ2_nu::MINERvA_CC0pi_XSec_1DQ2_nu(nuiskey samplekey) {
//********************************************************************

  // Sample overview ---------------------------------------------------
  std::string descrip = "MINERvA_CC0pi_XSec_1DQ2_nu sample. \n" \
                        "Target: CH \n" \
                        "Flux: MINERvA Forward Horn Current Low Energy numu\n" \
                        "Signal: Any event with 1 muon, any nucleons, and no other FS particles \n";

  // Setup common settings
  fSettings = LoadSampleSettings(samplekey);
  fSettings.SetDescription(descrip);
  fSettings.SetXTitle("Q^{2}_{QE} (GeV^{2})");
  fSettings.SetYTitle("d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
  fSettings.SetAllowedTypes("FIX,FREE,SHAPE/DIAG,FULL/NORM/MASK", "FIX/FULL");
  fSettings.SetEnuRange(0.0, 100.0);
  fSettings.DefineAllowedTargets("C,H");

  // CCQELike plot information
  fSettings.SetTitle("MINERvA_CC0pi_XSec_1DQ2_nu_2018");
  fSettings.SetDataInput(  FitPar::GetDataBase() + "/MINERvA/CC0pi_ptpz/q2_data.root");
  fSettings.SetCovarInput( FitPar::GetDataBase() + "/MINERvA/CC0pi_ptpz/q2_cov.root");
  fSettings.DefineAllowedSpecies("numu");

  FinaliseSampleSettings();

  fScaleFactor = (GetEventHistogram()->Integral("width") * 1E-38 / (fNEvents + 0.)) / this->TotalIntegratedFlux();

  // Data is __NOT__ bin width normalised, so override the default
  // Set the data file
  SetDataFromFile(fSettings.GetDataInput(), "q2_proj");
  fIsNoWidth = true;

  // Cut out the under and overflow bins
  TMatrixDSym * tempmat = StatUtils::GetCovarFromRootFile(fSettings.GetCovarInput(), "TMatrixDBase");

  // This covariance comes in full units, so need to multiply by scaling factor to ensure good decomp
  double ScalingFactor = 1.0E76;
  //double ScalingFactor = 1;

  fFullCovar = new TMatrixDSym(fDataHist->GetXaxis()->GetNbins());
  int xcnt = 0;
  for (int i = 0; i < tempmat->GetNrows(); ++i) {
    if (i == 0 || i == tempmat->GetNrows()-1) continue;
    int ycnt = 0;
    for (int j = 0; j < tempmat->GetNrows(); ++j) {
      if (j == 0 || j == tempmat->GetNrows()-1) continue;
      (*fFullCovar)(xcnt,ycnt) = (*tempmat)(i,j)*ScalingFactor;
      ycnt++;
    }
    if (ycnt != fDataHist->GetXaxis()->GetNbins()) {
      ERR(FTL) << "ycnt not data dimension" << std::endl;
      throw;
    }
    xcnt++;
  }

  // Now our covariances are in units of 1E38, so inverse covariance is in unit of 1E-38
  covar = StatUtils::GetInvert(fFullCovar);
  fDecomp = StatUtils::GetDecomp(fFullCovar);

  // Final setup  ---------------------------------------------------
  FinaliseMeasurement();
};

//********************************************************************
void MINERvA_CC0pi_XSec_1DQ2_nu::FillEventVariables(FitEvent *event) {
  //********************************************************************
  // Checking to see if there is a Muon
  if (event->NumFSParticle(13) == 0) return;

  // Get the muon kinematics
  TLorentzVector Pmu  = event->GetHMFSParticle(13)->fP;
  TLorentzVector Pnu = event->GetNeutrinoIn()->fP;
  double q2qe = FitUtils::Q2QErec(Pmu, Pnu, 34.0, true);
  fXVar = q2qe;

};

//********************************************************************
bool MINERvA_CC0pi_XSec_1DQ2_nu::isSignal(FitEvent *event) {
  //********************************************************************
  return SignalDef::isCC0pi_MINERvAPTPZ(event, 14, EnuMin, EnuMax);
};

