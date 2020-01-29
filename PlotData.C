#include <TMath.h>
#include <TFile.h>

#include "src/JConst.h"

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

TString fileName = "centralityData.root";

TH1D *hInputFlow;

TGraphErrors *gR[DET_N];

double mSize = 1.0;

int mStyle[DET_N] = {24, 25, 27};
int mColor[DET_N] = {1, 2, 4};

char *detNames[DET_N] = {"FT0-A", "FT0-C", "FV0"};

TFile *fIn;

void PlotData() {

    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "c");

    TLegend *leg = new TLegend(0.7, 0.2, 0.85, 0.4,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);

    fIn = TFile::Open(fileName, "read");
    if(fIn==0) ErrorExit(Form("Cannot open file: %s",fileName.Data()));

    for (int i=0; i<DET_N; i++) {

        gR[i] = (TGraphErrors*) fIn->Get(Form("gRDet%02i", i));

        gR[i]->SetMarkerStyle(mStyle[0]);
        gR[i]->SetMarkerColor(mColor[i]);
        gR[i]->SetMarkerSize(mSize);

        if (i==0) {
            gR[i]->Draw("AP");
        } else {
            gR[i]->Draw("P SAME");
        }

        gR[i]->GetXaxis()->SetRangeUser(10.0, 70.0);
        gR[i]->GetYaxis()->SetRangeUser(0.5, 1.01);
        gR[i]->GetXaxis()->CenterTitle(true);
        gR[i]->GetYaxis()->CenterTitle(true);
        gR[i]->SetTitle(Form("Resolution from ToyFlow data, n=1000; centrality; R_{%01i}", 2));
        gR[i]->Draw("P SAME");
        c->Update();

        leg->AddEntry(gR[i], Form("%s", detNames[i]), "p");
    }

    leg->Draw("SAME");
    c->Draw();
}
