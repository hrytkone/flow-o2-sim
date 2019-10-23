#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    int i, j, k;
    double pi = TMath::Pi();

    int NBINS = 150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=100;
    double logBW = (TMath::Log(LimH)-TMath::Log(LimL))/NBINS;
    for (i=0; i<=NBINS; i++) {
        LogBinsX[i] = LimL*exp(i*logBW);
    }

	hPt = new TH1D("hPt","pT - inclusive", NBINS, LogBinsX);
    hPt->Sumw2();
    hPhi = new TH1D("hPhi", "phi - uniform", 129, -pi, pi);
    hPhi->Sumw2();
    hAnisotropicPhi = new TH1D("hAnisotropicPhi", "phi - anisotropic", 129, -2*pi, 2*pi);
    hAnisotropicPhi->Sumw2();
    hCentrality = new TH1D("hCentrality", "centrality", CENTBINS_N, 0.0, 90.0);
    hCentrality->Sumw2();
    hEta = new TH1D("hEta", "pseudorapidity", 401, -6.0, 6.0);

    hMultiplicity = new TH1D("hMultiplicity", "Multiplicity - uniform", 300, 0.0, 30000.);
    hMultiplicity->Sumw2();

    for (i=0; i<DET_N; i++) {
        for (j=0; j<CENTBINS_N; j++){
            hSqrtSumWeights[i][j]= new TH1D(Form("hSqrtSumWeightsD%02iCENT%02i",i,j),"sqrt of sum of weights squares", 480, 0.0, 120.0);
            hSqrtSumWeights[i][j]->Sumw2();
        }
    }

    for (i=0; i<nCoef; i++){
        for (j=0; j<DET_N; j++){
            for (k=0; k<CENTBINS_N; k++){
                hRsub[i][j][k] = new TH1D(Form("hRsubH%02iD%02iCENT%02i",i+1,j,k),Form("hRsubH%02iD%02iCENT%02i",i+1,j,k),404,-1.01,1.01);
                hRsub[i][j][k]->Sumw2();
                hVnObs[i][j][k] = new TH1D(Form("hVnObsH%02iD%02iCENT%02i",i+1,j,k),Form("hVnObsH%02iD%02iCENT%02i",i+1,j,k),401,-1.5,1.5);
                hVnObs[i][j][k]->Sumw2();

                hQnQnAEP[i][j][k] = new TH1D(Form("hQnQnAEPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnQnAEPH%02iD%02iCENT%02i",i+1,j,k),401,-50.0,50.0);
                hQnQnAEP[i][j][k]->Sumw2();
                hQnAQnBEP[i][j][k] = new TH1D(Form("hQnAQnBEPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnAQnBEPH%02iD%02iCENT%02i",i+1,j,k),404,-1.01,1.01);
                hQnAQnBEP[i][j][k]->Sumw2();

                hQnQnASP[i][j][k] = new TH1D(Form("hQnQnASPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnQnASPH%02iD%02iCENT%02i",i+1,j,k),401,-10.0,110.0);
                hQnQnASP[i][j][k]->Sumw2();
                hQnAQnBSP[i][j][k] = new TH1D(Form("hQnAQnBSPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnAQnBSPH%02iD%02iCENT%02i",i+1,j,k),401,-10.0,110.0);
                hQnAQnBSP[i][j][k]->Sumw2();
            }
        }
    }

    for (i=0; i<PTBINS_N; i++) {
        hQnQnAPtBin[i] = new TH1D(Form("hQnQnAPtBin%02i", i+1), Form("hQnQnAPtBin%02i", i+1), 482, -11.0, 11.0);
        hQnQnAPtBin[i]->Sumw2();
        hSqrtSumWeightsPtBins[i] = new TH1D(Form("hSqrtSumWeightsPtBinsH%02i", i+1),Form("sqrt of sum of weights squares for pT bin %02i", i+1), 240, 0.0, 60.0);
        hSqrtSumWeightsPtBins[i]->Sumw2();
    }

    //FOR TESTING
    hV2ComplexPart = new TH1D(Form("hV2ComplexPartH%02i",i+1), Form("hV2ComplexPartH%02i",i+1), 401, -5.0, 5.0);
}
