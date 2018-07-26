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
  fSettings.SetDataInput(  FitPar::GetDataBase() + "/MINERvA/CC0pi_ptpz/data_q2_1d.root");
  fSettings.SetCovarInput( FitPar::GetDataBase() + "/MINERvA/CC0pi_ptpz/q2_Cov.root");
  fSettings.DefineAllowedSpecies("numu");

  FinaliseSampleSettings();

  fScaleFactor = (GetEventHistogram()->Integral("width") * 1E-38 / (fNEvents + 0.)) / this->TotalIntegratedFlux();

  // Data is __NOT__ bin width normalised, so override the default
  // Set the data file
  SetDataFromFile(fSettings.GetDataInput(), "q2_proj");
  fIsNoWidth = true;

  // Set the mapping values and the covariance matrix files
  //SetMapValuesFromText( fSettings.GetMapInput() );
  // Also have to make our own covariance matrix to exclude the under and overflow
  TMatrixDSym * tempmat = StatUtils::GetCovarFromRootFile(fSettings.GetCovarInput(), "TMatrixDBase");

  fFullCovar = new TMatrixDSym(fDataHist->GetXaxis()->GetNbins());
  int xcnt = 0;
  for (int i = 0; i < tempmat->GetNrows(); ++i) {
    if (i == 0 || i == tempmat->GetNrows()-1) continue;
    int ycnt = 0;
    for (int j = 0; j < tempmat->GetNrows(); ++j) {
      if (j == 0 || j == tempmat->GetNrows()-1) continue;
      (*fFullCovar)(xcnt,ycnt) = (*tempmat)(i,j);
      ycnt++;
    }
    if (ycnt != fDataHist->GetXaxis()->GetNbins()) {
      std::cerr << "ycnt not data dimension" << std::endl;
      throw;
    }
    xcnt++;
  }
  //std::cout << "data: " << fDataHist->GetXaxis()->GetNbins() << std::endl;
  //std::cout << "covar: " << xcnt << " x " << xcnt << std::endl;


  /*
     fFullCovar = tempmat;
  // Now we cut out every first and last to exclude under and overflow

  // Count the real covariance matrix x and y
  int xcounter = 0;

  int xbins = fDataHist->GetXaxis()->GetNbins();
  int ybins = fDataHist->GetYaxis()->GetNbins();
  // Loop over the x bins (underflow adds one, overflow adds one)
  for (int i = 0; i < (xbins+2)*(ybins+2); ++i) {
  // Skip the under and overflow
  if (i < (ybins+2) || i % (ybins+1) == 0 || ((i+1)%(ybins+1) == 0) || i > (ybins+2)*(xbins+1)) {
  //std::cout << "Skipping ibin " << i << std::endl;
  continue;
  }

  // The ycounter
  int ycounter = 0;
  // For one bin of pT we have pZ^2 bins
  for (int j = 0; j < (xbins+2)*(ybins+2); ++j) {

  // Skip the under and overflow
  if (j < (ybins+2) || j % (ybins+1) == 0 || ((j+1)%(ybins+1) == 0) || j > (ybins+2)*(xbins+1)) {
  //std::cout << "Skipping jbin " << j << std::endl;
  continue;
  }

  (*fFullCovar)(xcounter, ycounter) = (*tempmat)(i, j);
  //std::cout << xcounter << ", " << ycounter << " === " << i << ", " << j << std::endl;
  ycounter++;
  }
  // Check dimensions
  if (ycounter != xbins*ybins) {
  std::cerr << "Counted " << ycounter << " y bins in cov matrix" << std::endl;
  std::cerr << "Whereas there should be " << xbins*ybins << std::endl;
  }
  xcounter++;
  }
  // Check dimensions
  if (xcounter != xbins*ybins) {
  std::cerr << "Counted " << xcounter << " x bins in cov matrix" << std::endl;
  std::cerr << "Whereas there should be " << xbins*ybins << std::endl;
  }
  // Delete the temporary
  delete tempmat;
  */

  // Now can make the inverted covariance
  covar = StatUtils::GetInvert(fFullCovar);
  for (int i = 0; i < fFullCovar->GetNrows(); ++i) {
    //std::cout << sqrt((*fFullCovar)(i,i)) << " vs " << fDataHist->GetBinError(i+1) << " for " << fDataHist->GetXaxis()->GetBinLowEdge(i+1) << " (" << sqrt((*fFullCovar)(i,i))/fDataHist->GetBinError(i+1) << ")" << std::endl;
  }
  /*
  TFile *tempme = new TFile("temp_invcov.root", "recreate");
  tempme->cd();
  covar->Write("invcov");
  tempme->Close();
  throw;
  */
  fDecomp = StatUtils::GetDecomp(fFullCovar);

  // Use a TH2D version of the covariance to be able to use the global bin numbering scheme
  /*
     covar_th2d = new TH2D((fSettings.Title()+"_th2").c_str(), (fSettings.Title()+"_th2").c_str(), covar->GetNrows(), 0, covar->GetNrows(), covar->GetNcols(), 0, covar->GetNcols());
     for (int i = 0; i < covar_th2d->GetXaxis()->GetNbins(); ++i) {
     for (int j = 0; j < covar_th2d->GetYaxis()->GetNbins(); ++j) {
     covar_th2d->SetBinContent(i+1, j+1, (*covar)(i,j));
     }
     }
     std::cout << "covar is " << covar_th2d->GetXaxis()->GetNbins() << " x " << covar_th2d->GetYaxis()->GetNbins() << " = " << covar_th2d->GetXaxis()->GetNbins()*covar_th2d->GetYaxis()->GetNbins() << std::endl;
     std::cout << "data is " << fDataHist->GetXaxis()->GetNbins() << " x " << fDataHist->GetYaxis()->GetNbins() << " = " << fDataHist->GetXaxis()->GetNbins()*fDataHist->GetYaxis()->GetNbins() << std::endl;
     */

  // Let's make our own mapping histogram
  // The covariance matrix is dominant in Pt and sub-dominant in Pz and includes all bins with under and overflow
  // So we have 13x12 data/MC bins, and including the under/overflow we have 15x14=210 covariance bins
  // i.e. need to cut out the first and last bins of covariance matrix
  // Mapping histogram will have same dimensions as the data
  /*
     fMapHist = (TH2I*)(fDataHist->Clone());
     fMapHist->Reset();
     std::string MapTitle = std::string(fDataHist->GetName())+"_MAP";
     fMapHist->SetNameTitle(MapTitle.c_str(), MapTitle.c_str());
     int counter = 1;
     for (int i = 0; i <= fDataHist->GetXaxis()->GetNbins()+1; ++i) {
     for (int j = 0; j <= fDataHist->GetYaxis()->GetNbins()+1; ++j) {
     if (i == 0 || i == fDataHist->GetXaxis()->GetNbins()+1 || j == 0 || j == fDataHist->GetYaxis()->GetNbins()+1) {
     fMapHist->SetBinContent(i+1, j+1, 0);
     } else {
     fMapHist->SetBinContent(i+1, j+1, counter);
     counter++;
     }
     std::cout << fMapHist->GetBinContent(i+1, j+1) << "    " << fDataHist->GetBinContent(i+1, j+1) << std::endl;
     }
     std::cout << std::endl;
     }
     */

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

/*
//********************************************************************
// Custom likelihood calculator because binning of covariance matrix
double MINERvA_CC0pi_XSec_1DQ2_nu::GetLikelihood() {
//********************************************************************

// The calculated chi2
double chi2 = 0.0;

// Support shape comparisons
double scaleF = fDataHist->Integral() / fMCHist->Integral();
if (fIsShape) {
fMCHist->Scale(scaleF);
fMCFine->Scale(scaleF);
//PlotUtils::ScaleNeutModeArray((TH1**)fMCHist_PDG, scaleF);
}

// Calculate the test-statistic
for (int i = 0; i < covar_th2d->GetXaxis()->GetNbins()+1; ++i) {
// Get the global bin for x
int xbin1, ybin1, zbin1;
fDataHist->GetBinXYZ(i, xbin1, ybin1, zbin1);
double xlo1 = fDataHist->GetXaxis()->GetBinLowEdge(xbin1);
double xhi1 = fDataHist->GetXaxis()->GetBinLowEdge(xbin1+1);
double ylo1 = fDataHist->GetYaxis()->GetBinLowEdge(ybin1);
double yhi1 = fDataHist->GetYaxis()->GetBinLowEdge(ybin1+1);
if (xlo1 < fDataHist->GetXaxis()->GetBinLowEdge(1) ||
ylo1 < fDataHist->GetYaxis()->GetBinLowEdge(1) ||
xhi1 > fDataHist->GetXaxis()->GetBinLowEdge(fDataHist->GetXaxis()->GetNbins()+1) ||
yhi1 > fDataHist->GetYaxis()->GetBinLowEdge(fDataHist->GetYaxis()->GetNbins()+1)) continue;

// Get the data
double data1 = fDataHist->GetBinContent(i);
// Get the MC
double mc1 = fMCHist->GetBinContent(i);

for (int j = 0; j < covar_th2d->GetYaxis()->GetNbins()+1; ++j) {

int xbin2, ybin2, zbin2;
fDataHist->GetBinXYZ(j, xbin2, ybin2, zbin2);
double xlo2 = fDataHist->GetXaxis()->GetBinLowEdge(xbin2);
double xhi2 = fDataHist->GetXaxis()->GetBinLowEdge(xbin2+1);
double ylo2 = fDataHist->GetYaxis()->GetBinLowEdge(ybin2);
double yhi2 = fDataHist->GetYaxis()->GetBinLowEdge(ybin2+1);

if (xlo2 < fDataHist->GetXaxis()->GetBinLowEdge(1) ||
ylo2 < fDataHist->GetYaxis()->GetBinLowEdge(1) ||
xhi2 > fDataHist->GetXaxis()->GetBinLowEdge(fDataHist->GetXaxis()->GetNbins()+1) ||
yhi2 > fDataHist->GetYaxis()->GetBinLowEdge(fDataHist->GetYaxis()->GetNbins()+1)) continue;

//std::cout << "Correlating: (" << xlo1 << "-" << xhi1 << "," << ylo1 << "-" << yhi1 << ") with (" << xlo2 << "-" << xhi2 << "," << ylo2 << "-" << yhi2 << ")" << std::endl;

// Get the data
double data2 = fDataHist->GetBinContent(j);
// Get the MC
double mc2 = fMCHist->GetBinContent(j);
//std::cout << data1 << " " << mc1 << std::endl;
//std::cout << data2 << " " << mc2 << std::endl;
//std::cout << std::endl;

// Get the inverse covariance matrix entry
double coventry = covar_th2d->GetBinContent(i, j);

//std::cout << fDataHist->GetXaxis()->GetBinLowEdge(i+1) << " - " << fDataHist->GetXaxis()->GetBinLowEdge(i+2) << ", " << fDataHist->GetYaxis()->GetBinLowEdge(j+1) << " - " << fDataHist->GetYaxis()->GetBinLowEdge(j+2) << " = " << coventry << " (global = " << global << ")" << std::endl;

chi2 += (data1-mc1)*coventry*(data2-mc2);
}
}

// Normalisation penalty term if included
if (fAddNormPen) {
chi2 +=
(1 - (fCurrentNorm)) * (1 - (fCurrentNorm)) / (fNormError * fNormError);
  LOG(REC) << "Norm penalty = "
  << (1 - (fCurrentNorm)) * (1 - (fCurrentNorm)) /
(fNormError * fNormError)
  << std::endl;
  }

// Adjust the shape back to where it was.
if (fIsShape and !FitPar::Config().GetParB("saveshapescaling")) {
  fMCHist->Scale(1. / scaleF);
  fMCFine->Scale(1. / scaleF);
}

fLikelihood = chi2;

return chi2;
};
*/
