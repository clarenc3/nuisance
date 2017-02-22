#ifndef MINIBOONE_BOXES_H
#define MINIBOONE_BOXES_H
#include "MeasurementBase.h"
#include "MeasurementVariableBox.h"

class MiniBooNE_CCQELike_Box : public MeasurementVariableBox {
public:
	MiniBooNE_CCQELike_Box() { Reset(); };

	void Reset() {
		fTargetTest = -999;
		fNProtons = fNNeutrons = fNIntermediatePions = 0;
		fFSPionMom.clear();
	}

	void FillBoxFromEvent(FitEvent* event) {

		fPDGnu = event->PDGnu();
		fNProtons = event->NumFSParticle(2212);
		fNNeutrons = event->NumFSParticle(2112);
		fTargetTest = event->fTargetPDG;

		for (size_t i = 0; i < event->fNParticles; i++) {
			FitParticle* p = event->PartInfo(i);

			if (p->Status() == kInitialState) {
				continue;
			} else if (p->Status() == kFinalState) {
				fFSPionMom.push_back(p->fP.Vect().Mag());
			} else {
				fNIntermediatePions++;
			}
		}
	};

	int fNProtons;
	int fNNeutrons;
	int fNIntermediatePions;
	std::vector<double> fFSPionMom;
	int fTargetTest;
	int fPDGnu;

};

#endif