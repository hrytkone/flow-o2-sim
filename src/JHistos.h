/*
 *  JHistos.h
 *
 */

#include "TH1D.h"
#include "TH2D.h"

#include "JConst.h"

class JHistos {

public:
	JHistos();
	virtual ~JHistos(){;}

	TH1D *hPt;
    TH1D *hPhi;
    TH1D *hAnisotropicPhi;
    TH1D *hEta;
    TH1D *hMultiplicity;

    TH1D *hSqrtSumWeights[DET_N];

    // Historgrams for resolutions and vobs
    TH1D *hRsub[nCoef][DET_N];
    TH1D *hVnObs[nCoef][DET_N];

    // Histograms for EP-method
    TH1D *hQnQnAEP[nCoef][DET_N];
    TH1D *hQnAQnBEP[nCoef][DET_N];

    // Histograms for SP-method
    TH1D *hQnQnASP[nCoef][DET_N];
    TH1D *hQnAQnBSP[nCoef][DET_N];

};
