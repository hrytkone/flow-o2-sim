#include <TMath.h>
#include <TFile.h>

#include "src/ResIter.h"
#include "src/JConst.h"
#include "src/Filipad.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
double CalculateVn(double QnQnA, double QnAQnB, double w);
double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr);
void checkUnderOverFlow( TH1 *h );

void MakeGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root", const int iDet=2) {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    int i, n;
    double pi = TMath::Pi();
    double inputFlow[nCoef] = {0.0};

    TH1D *hInputNumbers = (TH1D*)fIn->Get("hInputNumbers");
    checkUnderOverFlow(hInputNumbers);
    double nEvents = hInputNumbers->GetBinContent(1);
    double nOfFiles = hInputNumbers->GetBinContent(7);
    for(int i = 0; i < nCoef; i++) {
        inputFlow[i] = hInputNumbers->GetBinContent(2+i);
        inputFlow[i] /= nOfFiles;
        cout << inputFlow[i] << endl;
    }

    TH1D *hSqrtSumWeights[DET_N];
    hSqrtSumWeights[iDet] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsD%02i",iDet));
    checkUnderOverFlow(hSqrtSumWeights[iDet]);

    TH1D *hInputFlow = new TH1D("hInputFlow", "hInputFlow", nCoef, 0.5, double(nCoef)+0.5);
    hInputFlow->SetLineStyle(1);
    hInputFlow->SetLineColor(1);
    hInputFlow->SetLineWidth(1);
    for(int i = 0; i < nCoef; i++)
        hInputFlow->Fill(double(i+1), inputFlow[i]);
    hInputFlow->Fill(double(nCoef+1), 0.0);

    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    checkUnderOverFlow(hPhi);

    //=====vn=====
    // Observed vn
    TH1D *hVnObs[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hVnObs[i][iDet] = (TH1D*)fIn->Get(Form("hVnObsH%02iD%02i",i+1,iDet));
        checkUnderOverFlow(hVnObs[i][iDet]);
    }

    // Resolution parameter
    TH1D *hRsub[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hRsub[i][iDet] = (TH1D*)fIn->Get(Form("hRsubH%02iD%02i",i+1,iDet));
        checkUnderOverFlow(hRsub[i][iDet]);
    }

    double vn[nCoef], errorVn[nCoef] = {0};
    double vnRatio[nCoef], errorVnRatio[nCoef] = {0};

    double Rinit;

    double khi0 = 0.5;
    double err = 0.0001;
    double khi;
    double R[nCoef], errorR[nCoef] = {0};
    double Rnonuni[nCoef], errorRnonuni[nCoef] = {0};

    for (i = 0; i < nCoef; i++) {
        n = i+1;
        vn[i] = hVnObs[i][iDet]->GetMean();

        Rinit = TMath::Sqrt(hRsub[i][iDet]->GetMean());
        khi = RIter(khi0, Rinit, err);
        R[i] = R1(TMath::Sqrt(2)*khi); //Because khi sim sqrt(Multi) and full event has twice multiplicity compared to A or B.
        errorR[i] = CalculateRerror(khi, err);

        cout << "R=    " << R[i] <<                       "  err=" << errorR[i] << "\n";

        vn[i] /= R[i];
        errorVn[i] = GetVnError(hVnObs[i][iDet]->GetMean(), hVnObs[i][iDet]->GetMeanError(), R[i], errorR[i]);

        if(inputFlow[i]==0) {
            vnRatio[i] = 0;
            errorVnRatio[i] = 0;

        } else {
            vnRatio[i] = vn[i]/inputFlow[i];
            errorVnRatio[i] = errorVn[i]/inputFlow[i];
        }
    }

    //vn{EP} and vn{SP}
    TH1D *hQnQnAEP[nCoef][DET_N];
    TH1D *hQnAQnBEP[nCoef][DET_N];
    TH1D *hQnQnASP[nCoef][DET_N];
    TH1D *hQnAQnBSP[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hQnQnAEP[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnAEPH%02iD%02i",i+1,iDet));
        checkUnderOverFlow(hQnQnAEP[i][iDet]);
        hQnAQnBEP[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBEPH%02iD%02i",i+1,iDet));
        checkUnderOverFlow(hQnAQnBEP[i][iDet]);
        hQnQnASP[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnASPH%02iD%02i",i+1,iDet));
        checkUnderOverFlow(hQnQnASP[i][iDet]);
        hQnAQnBSP[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBSPH%02iD%02i",i+1,iDet));
        checkUnderOverFlow(hQnAQnBSP[i][iDet]);
    }

    double w = hSqrtSumWeights[iDet]->GetMean();
    double wError = hSqrtSumWeights[iDet]->GetMeanError();

    double vnEP[nCoef], errorVnEP[nCoef] = {0};
    double vnSP[nCoef], errorVnSP[nCoef] = {0};

    double vnEPRatio[nCoef], errorVnEPRatio[nCoef] = {0};
    double vnSPRatio[nCoef], errorVnSPRatio[nCoef] = {0};
    for (i=0; i<nCoef; i++) {
        vnEP[i] = CalculateVn(hQnQnAEP[i][iDet]->GetMean(), hQnAQnBEP[i][iDet]->GetMean(), w);
        errorVnEP[i] = CalculateVnError(hQnQnAEP[i][iDet]->GetMean(), hQnAQnBEP[i][iDet]->GetMean(), hQnQnAEP[i][iDet]->GetMeanError(), hQnAQnBEP[i][iDet]->GetMeanError(),  w, wError);
        vnSP[i] = CalculateVn(hQnQnASP[i][iDet]->GetMean(), hQnAQnBSP[i][iDet]->GetMean(), w);
        errorVnSP[i] = CalculateVnError(hQnQnASP[i][iDet]->GetMean(), hQnAQnBSP[i][iDet]->GetMean(), hQnQnASP[i][iDet]->GetMeanError(), hQnAQnBSP[i][iDet]->GetMeanError(),  w, wError);

        if(inputFlow[i]==0) {
            vnEPRatio[i] = 0;
            errorVnEPRatio[i] = 0;
            vnSPRatio[i] = 0;
            errorVnSPRatio[i] = 0;

        } else {
            vnEPRatio[i] = vnEP[i]/inputFlow[i];
            errorVnEPRatio[i] = errorVnEP[i]/inputFlow[i];
            vnSPRatio[i] = vnSP[i]/inputFlow[i];
            errorVnSPRatio[i] = errorVnSP[i]/inputFlow[i];
        }
    }

    // Make graphs
    TGraphErrors *gVn = new TGraphErrors(nCoef);
    TGraphErrors *gVnEP = new TGraphErrors(nCoef);
    TGraphErrors *gVnSP = new TGraphErrors(nCoef);
    TGraphErrors *gR = new TGraphErrors(nCoef);

    TGraphErrors *gVnRatio = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPRatio = new TGraphErrors(nCoef);
    TGraphErrors *gVnSPRatio = new TGraphErrors(nCoef);

    for(i = 0; i < nCoef; i++){
        gVn->SetPoint(i, double(i+1)-0.1, vn[i]);
        gVn->SetPointError(i, 0.0, errorVn[i]);
        gVnEP->SetPoint(i, double(i+1)+0.1, vnEP[i]);
        gVnEP->SetPointError(i, 0.0, errorVnEP[i]);
        gVnSP->SetPoint(i, double(i+1)+0.3, vnSP[i]);
        gVnSP->SetPointError(i, 0.0, errorVnSP[i]);
        gR->SetPoint(i, double(i+1)+0.05, R[i]);
        gR->SetPointError(i, 0.0, errorR[i]);

        gVnRatio->SetPoint(i, double(i+1)-0.1, vnRatio[i]);
        gVnRatio->SetPointError(i, 0.0, errorVnRatio[i]);
        gVnEPRatio->SetPoint(i, double(i+1)+0.1, vnEPRatio[i]);
        gVnEPRatio->SetPointError(i, 0.0, errorVnEPRatio[i]);
        gVnSPRatio->SetPoint(i, double(i+1)+0.3, vnSPRatio[i]);
        gVnSPRatio->SetPointError(i, 0.0, errorVnSPRatio[i]);
    }

    fOut->cd();
    hPhi->Write("hPhi");
    hInputFlow->Write("hInputFlow");
    gVn->Write("gVn");
    gR->Write("gR");
    gVnEP->Write("gVnEP");
    gVnSP->Write("gVnSP");

    gVnRatio->Write("gVnRatio");
    gVnEPRatio->Write("gVnEPRatio");
    gVnSPRatio->Write("gVnSPRatio");

    fOut->Close();

}

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr) {
    return (vnObs/Rn)*TMath::Sqrt((vnObsErr/vnObs)*(vnObsErr/vnObs) + (RnErr/Rn)*(RnErr/Rn));
}

double CalculateVn(double QnQnA, double QnAQnB, double w) {
    return QnQnA/TMath::Sqrt(QnAQnB)/w;
}

double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr) {
    double vn = CalculateVn(QnQnA, QnAQnB, w);
    return vn*TMath::Sqrt(QnQnAerr*QnQnAerr/QnQnA/QnQnA + 0.5*0.5*QnAQnBerr*QnAQnBerr/QnAQnB/QnAQnB + wErr*wErr/w/w);
}

void checkUnderOverFlow( TH1 *h ){
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}
