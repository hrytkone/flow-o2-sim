#include <TMath.h>
#include <vector>

#include "src/ResIter.h"
#include "src/JConst.h"

#define NFILES 6

void checkUnderOverFlow( TH1 *h );

void MakeCentralityGraphs(const char *dirname="/home/heimarry/Desktop/centrality_data/analysis-data-toyflow/", const char *ext=".root") {

    int iFile = 0;
    double pi = TMath::Pi();
    TString datadir(dirname);

    double cent[NFILES] = {0};
    double midCent = 15.0;
    for (int i=0; i<NFILES; i++) {
        cent[i] = midCent;
        midCent += 10.0;
    }

    TH1D *hRsub[nCoef][DET_N];
    double R[DET_N][NFILES], errorR[DET_N][NFILES];

    TGraphErrors *gR[DET_N];
    for (int i=0; i<DET_N; i++)
    gR[i] = new TGraphErrors(NFILES);

    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    TFile *fIn;
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
                cout << fname.Data() << endl;
                fIn = TFile::Open(datadir + fname, "read");

                for (int iDet=0; iDet<DET_N; iDet++) {

                    // Resolution parameter
                    for (int i=0; i<nCoef; i++) {
                        hRsub[i][iDet] = (TH1D*)fIn->Get(Form("hRsubH%02iD%02i",i+1,iDet));
                        checkUnderOverFlow(hRsub[i][iDet]);
                    }

                    double Rinit;
                    double khi0 = 0.5;
                    double err = 0.0001;
                    double khi;

                    Rinit = TMath::Sqrt(hRsub[1][iDet]->GetMean());
                    khi = RIter(khi0, Rinit, err);
                    R[iDet][iFile] = R1(TMath::Sqrt(2)*khi);
                    errorR[iDet][iFile] = CalculateRerror(khi, err);

                    gR[iDet]->SetPoint(iFile, cent[iFile], R[iDet][iFile]);
                    gR[iDet]->SetPointError(iFile, 0.0, errorR[iDet][iFile]);

                }
                iFile++;
            }
        }
    }

    TFile *fOut = TFile::Open("centralityData.root", "recreate");
    fOut->cd();
    for (int i=0; i<DET_N; i++) {
        gR[i]->Write(Form("gRDet%02i", i));
    }

    fOut->Close();
}

void checkUnderOverFlow( TH1 *h ){
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}
