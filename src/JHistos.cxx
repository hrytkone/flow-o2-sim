#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    int i, j;
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
    hEta = new TH1D("hEta", "pseudorapidity", 401, -6.0, 6.0);

    hMultiplicity = new TH1D("hMultiplicity", "Multiplicity - uniform", 300, 0.0, 30000.);
    hMultiplicity->Sumw2();

    for (i=0; i<DET_N; i++) {
        hSqrtSumWeights[i]= new TH1D(Form("hSqrtSumWeightsD%02i",i),"sqrt of sum of weights squares", 480, 0.0, 120.0);
        hSqrtSumWeights[i]->Sumw2();
    }

    for (i=0; i<nCoef; i++){
        for (j=0; j<DET_N; j++){
            hRsub[i][j] = new TH1D(Form("hRsubH%02iD%02i",i+1,j),Form("hRsubH%02iD%02i",i+1,j),404,-1.01,1.01);
            hRsub[i][j]->Sumw2();
            hVnObs[i][j] = new TH1D(Form("hVnObsH%02iD%02i",i+1,j),Form("hVnObsH%02iD%02i",i+1,j),401,-1.5,1.5);
            hVnObs[i][j]->Sumw2();

            hQnA[i][j] = new TH1D(Form("hQnAH%02iD%02i",i+1,j),Form("hQnAH%02iD%02i",i+1,j),401,-100,100);
            hQnA[i][j]->Sumw2();

            hQnB[i][j] = new TH1D(Form("hQnBH%02iD%02i",i+1,j),Form("hQnBH%02iD%02i",i+1,j),401,-100,100);
            hQnB[i][j]->Sumw2();

            //hQnQnAEP[i][j] = new TH1D(Form("hQnQnAEPH%02iD%02i",i+1,j),Form("hQnQnAEPH%02iD%02i",i+1,j),401,-50.0,50.0);
            //hQnQnAEP[i][j]->Sumw2();
            //hQnAQnBEP[i][j] = new TH1D(Form("hQnAQnBEPH%02iD%02i",i+1,j),Form("hQnAQnBEPH%02iD%02i",i+1,j),404,-1.01,1.01);
            //hQnAQnBEP[i][j]->Sumw2();

            //hQnQnASP[i][j] = new TH1D(Form("hQnQnASPH%02iD%02i",i+1,j),Form("hQnQnASPH%02iD%02i",i+1,j),401,-10.0,110.0);
            //hQnQnASP[i][j]->Sumw2();
            //hQnAQnBSP[i][j] = new TH1D(Form("hQnAQnBSPH%02iD%02i",i+1,j),Form("hQnAQnBSPH%02iD%02i",i+1,j),401,-10.0,110.0);
            //hQnAQnBSP[i][j]->Sumw2();
        }
    }

}
