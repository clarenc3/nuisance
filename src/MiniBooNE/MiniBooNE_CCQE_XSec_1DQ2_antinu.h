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

#ifndef MINIBOONE_CCQE_XSEC_1DQ2_ANTINU_H_SEEN
#define MINIBOONE_CCQE_XSEC_1DQ2_ANTINU_H_SEEN

#include "Measurement1D.h"
#include "MeasurementBase.h"
#include "MiniBooNE_Boxes.h"
//********************************************************************
class MiniBooNE_CCQE_XSec_1DQ2_antinu : public Measurement1D {
//********************************************************************

public:

  MiniBooNE_CCQE_XSec_1DQ2_antinu(nuiskey samplekey);
  virtual ~MiniBooNE_CCQE_XSec_1DQ2_antinu() {};

  void FillEventVariables(FitEvent *event);
  // void Write(std::string drawOpt);
  // void FillHistograms();
  bool isSignal(FitEvent *event);
  // void ScaleEvents();
  // void ApplyNormScale(double norm);
  // void ResetAll();
  MeasurementVariableBox* CreateBox();

  void FillExtraHistograms(MeasurementVariableBox* vars,
                              double weight = 1.0);

  TH1D* fMCHist_NONCCPIM[61]; ///< Plots in CCQELike mode to tag PDG of the NONCCPIM background
  TH1D* fMCHist_CCPIM[61]; ///< Plots in CCQELike mode to tag PDG of the CCPIM background
  NuNuBarTrueModeStack* fMCHist_CCQELIKE; ///< Plots in CCQELike mode to tag PDG of the background

private:
  bool fCCQElike; ///< Flag for running in CCQELike mode
  bool fUseCorrectedCTarget; ///< Flag for using corrected `C-Target' data.
  TH1D* fDataHist_CCQELIKE; ///< CCQELike data contribution
  TH1D* fDataHist_CCPIM; ///< CCPIM data contribution
  TH1D* fDataHist_NONCCPIM; ///< NONCCPIM data contribution
};

#endif
