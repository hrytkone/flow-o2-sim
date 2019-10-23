#include <TMath.h>
#include <TFile.h>

#include "src/ResIter.h"
#include "src/JConst.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
double CalculateVn(double QnQnA, double QnAQnB, double w);
double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr);
void checkUnderOverFlow( TH1 *h );

void MakeCentralityGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root", const int iDet=3) {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    int nPtBins = 9;

    int i, j, n;
    double pi = TMath::Pi();

    double inputFlow[nCoef][CENTBINS_N-4] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0277401, 0.04488324, 0.06521883, 0.08433443, 0.09597485, 0.10087206, 0.09925828},
        {0.02039728, 0.02369955, 0.02670301, 0.02950095, 0.03118808, 0.03120636, 0.02918556},
        {0.01013229, 0.01171893, 0.0131265, 0.01479335, 0.0159713, 0.01644628, 0.01535014},
        {0.00415816, 0.00467961, 0.00528238, 0.006501, 0.0068885, 0.00690379, 0.00575251}
    };

    TH1D *hInputFlow[nCoef];
    for (i=0; i<nCoef; i++) {
        hInputFlow[i] = new TH1D(Form("hInputFlow%02i",i+1), Form("hInputFlow%02i",i+1), CENTBINS_N-4, centBins);
        hInputFlow[i]->SetLineStyle(1);
        hInputFlow[i]->SetLineColor(1);
        hInputFlow[i]->SetLineWidth(1);
    }
    for(i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hInputFlow[i]->Fill(centBins[j], inputFlow[i][j]);
        }
    }

    //=====vn=====
    // Observed vn
    TH1D *hVnObs[nCoef][DET_N][CENTBINS_N];
    for (i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hVnObs[i][iDet][j] = (TH1D*)fIn->Get(Form("hVnObsH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hVnObs[i][iDet][j]);
        }
    }

    // Resolution parameter
    TH1D *hRsub[nCoef][DET_N][CENTBINS_N];
    TH1D *hRtrue[nCoef][DET_N][CENTBINS_N];
    for (i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hRsub[i][iDet][j] = (TH1D*)fIn->Get(Form("hRsubH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hRsub[i][iDet][j]);
            hRtrue[i][iDet][j] = (TH1D*)fIn->Get(Form("hRtrueH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hRtrue[i][iDet][j]);
        }
    }

    double vn[nCoef][CENTBINS_N], errorVn[nCoef][CENTBINS_N] = {0};
    double vnTrue[nCoef][CENTBINS_N], errorVnTrue[nCoef][CENTBINS_N] = {0};
    double Rinit;

    double khi0 = 0.5;
    double err = 0.0001;
    double khi;
    double R[nCoef][CENTBINS_N], errorR[nCoef][CENTBINS_N] = {0};

    for (i=0; i<nCoef; i++) {
        cout << "n=" << i+1 << "----------------\n";
        for (j=0; j<CENTBINS_N-4; j++) {
            cout << "cent=" << j << "\n";
            n = i+1;
            vn[i][j] = hVnObs[i][iDet][j]->GetMean();
            vnTrue[i][j] = hVnObs[i][iDet][j]->GetMean();

            Rinit = TMath::Sqrt(hRsub[i][iDet][j]->GetMean());
            khi = RIter(khi0, Rinit, err);
            R[i][j] = R1(TMath::Sqrt(2)*khi); //Because khi sim sqrt(Multi) and full event has twice multiplicity compared to A or B.
            errorR[i][j] = CalculateRerror(khi, err);

            cout << "R=    " << R[i][j] <<                       "  err=" << errorR[i][j] << "\n";
            cout << "Rtrue=" << hRtrue[i][iDet][j]->GetMean() << "  err=" << hRtrue[i][iDet][j]->GetMeanError() << "\n";
            cout << "R/Rtrue=" << R[i][j]/hRtrue[i][iDet][j]->GetMean() << "\n\n";

            vn[i][j] /= R[i][j];
            errorVn[i][j] = GetVnError(hVnObs[i][iDet][j]->GetMean(), hVnObs[i][iDet][j]->GetMeanError(), R[i][j], errorR[i][j]);

            vnTrue[i][j] /= hRtrue[i][iDet][j]->GetMean();
            errorVnTrue[i][j] = GetVnError(hVnObs[i][iDet][j]->GetMean(), hVnObs[i][iDet][j]->GetMeanError(), hRtrue[i][iDet][j]->GetMean(), hRtrue[i][iDet][j]->GetMeanError());
        }
    }

    //vn{EP} and vn{SP}
    TH1D *hQnQnAEP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnAQnBEP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnQnASP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnAQnBSP[nCoef][DET_N][CENTBINS_N];
    for (i = 0; i < nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hQnQnAEP[i][iDet][j] = (TH1D*)fIn->Get(Form("hQnQnAEPH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hQnQnAEP[i][iDet][j]);
            hQnAQnBEP[i][iDet][j] = (TH1D*)fIn->Get(Form("hQnAQnBEPH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hQnAQnBEP[i][iDet][j]);
            hQnQnASP[i][iDet][j] = (TH1D*)fIn->Get(Form("hQnQnASPH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hQnQnASP[i][iDet][j]);
            hQnAQnBSP[i][iDet][j] = (TH1D*)fIn->Get(Form("hQnAQnBSPH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hQnAQnBSP[i][iDet][j]);
        }
    }

    TH1D *hSqrtSumWeights[DET_N][CENTBINS_N];
    double w[CENTBINS_N], wError[CENTBINS_N] = {0};
    for (i=0; i<CENTBINS_N-4; i++) {
        hSqrtSumWeights[iDet][i] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsD%02iCENT%02i",iDet,i));
        checkUnderOverFlow(hSqrtSumWeights[iDet][i]);
        w[i] = hSqrtSumWeights[iDet][i]->GetMean();
        wError[i] = hSqrtSumWeights[iDet][i]->GetMeanError();
    }

    double vnEP[nCoef][CENTBINS_N-4], errorVnEP[nCoef][CENTBINS_N-4] = {0};
    double vnSP[nCoef][CENTBINS_N-4], errorVnSP[nCoef][CENTBINS_N-4] = {0};
    for (i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            vnEP[i][j] = CalculateVn(hQnQnAEP[i][iDet][j]->GetMean(), hQnAQnBEP[i][iDet][j]->GetMean(), w[j]);
            errorVnEP[i][j] = CalculateVnError(hQnQnAEP[i][iDet][j]->GetMean(), hQnAQnBEP[i][iDet][j]->GetMean(), hQnQnAEP[i][iDet][j]->GetMeanError(), hQnAQnBEP[i][iDet][j]->GetMeanError(),  w[j], wError[j]);
            vnSP[i][j] = CalculateVn(hQnQnASP[i][iDet][j]->GetMean(), hQnAQnBSP[i][iDet][j]->GetMean(), w[j]);
            errorVnSP[i][j] = CalculateVnError(hQnQnASP[i][iDet][j]->GetMean(), hQnAQnBSP[i][iDet][j]->GetMean(), hQnQnASP[i][iDet][j]->GetMeanError(), hQnAQnBSP[i][iDet][j]->GetMeanError(),  w[j], wError[j]);
        }
    }

    // Make graphs
    TGraphErrors *gVnTrue[nCoef];
    TGraphErrors *gVn[nCoef];
    TGraphErrors *gVnEP[nCoef];
    TGraphErrors *gVnSP[nCoef];
    TGraphErrors *gRtrue[nCoef];
    TGraphErrors *gR[nCoef];

    for(i=0; i<nCoef; i++){
        gVnTrue[i] = new TGraphErrors(CENTBINS_N-4);
        gVn[i] = new TGraphErrors(CENTBINS_N-4);
        gVnEP[i] = new TGraphErrors(CENTBINS_N-4);
        gVnSP[i] = new TGraphErrors(CENTBINS_N-4);
        gRtrue[i] = new TGraphErrors(CENTBINS_N-4);
        gR[i] = new TGraphErrors(CENTBINS_N-4);
        for (j=0; j<CENTBINS_N-4; j++) {
            gVnTrue[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-1.5, vnTrue[i][j]);
            gVnTrue[i]->SetPointError(j, 0.0, errorVnTrue[i][j]);
            gVn[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-0.5, vn[i][j]);
            gVn[i]->SetPointError(j, 0.0, errorVn[i][j]);
            gVnEP[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2+0.5, vnEP[i][j]);
            gVnEP[i]->SetPointError(j, 0.0, errorR[i][j]);
            gVnSP[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2+1.5, vnSP[i][j]);
            gVnSP[i]->SetPointError(j, 0.0, errorR[i][j]);

            gRtrue[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-1.5, hRtrue[i][iDet][j]->GetMean());
            gRtrue[i]->SetPointError(j, 0.0, hRtrue[i][iDet][j]->GetMeanError());
            gR[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-0.5, R[i][j]);
            gR[i]->SetPointError(j, 0.0, errorR[i][j]);
        }
    }

    fOut->cd();
    for (i=0; i<nCoef; i++) {
        hInputFlow[i]->Write(Form("hInputFlow%02i", i+1));

        gVnTrue[i]->Write(Form("gVnTrueH%02i", i+1));
        gVn[i]->Write(Form("gVnH%02i", i+1));
        gVnEP[i]->Write(Form("gVnEPH%02i", i+1));
        gVnSP[i]->Write(Form("gVnSPH%02i", i+1));
        gRtrue[i]->Write(Form("gRtrueH%02i", i+1));
        gR[i]->Write(Form("gRH%02i", i+1));

    }

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
