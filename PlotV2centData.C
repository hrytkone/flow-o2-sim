#include <TMath.h>
#include <TFile.h>

#include "src/JConst.h"

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

TString fileName = "toyFlowCentralityGraphs.root";

TH1D *hInputFlow;

TGraphErrors *gR;
TGraphErrors *gRtrue;
TGraphErrors *gVn;
TGraphErrors *gVnTrue;
TGraphErrors *gVnEP;
TGraphErrors *gVnSP;

int mMarker = 24;
double mSize = 1.0;

TFile *fIn;

void PlotV2centData() {

    fIn = TFile::Open(fileName, "read");

    int i, j;
    for (i=0; i<nCoef; i++) {
        if(fIn==0) ErrorExit(Form("Cannot open file: %s",fileName.Data()));

        hInputFlow = (TH1D*) fIn->Get(Form("hInputFlow%02i", 2));

        gR = (TGraphErrors*) fIn->Get(Form("gRH%02i", 2));
        gRtrue = (TGraphErrors*) fIn->Get(Form("gRtrueH%02i", 2));

        gVn = (TGraphErrors*) fIn->Get(Form("gVnH%02i", 2));
        gVnTrue = (TGraphErrors*) fIn->Get(Form("gVnTrueH%02i", 2));
        gVnEP = (TGraphErrors*) fIn->Get(Form("gVnEPH%02i", 2));
        gVnSP = (TGraphErrors*) fIn->Get(Form("gVnSPH%02i", 2));

        gR->SetMarkerStyle(mMarker);
        gR->SetMarkerColor(2);
        gR->SetMarkerSize(mSize);

        gRtrue->SetMarkerStyle(mMarker+1);
        gRtrue->SetMarkerColor(4);
        gRtrue->SetMarkerSize(mSize);

        gVn->SetMarkerStyle(mMarker);
        gVn->SetMarkerColor(2);
        gVn->SetMarkerSize(mSize);

        gVnTrue->SetMarkerStyle(mMarker+1);
        gVnTrue->SetMarkerColor(4);
        gVnTrue->SetMarkerSize(mSize);

        gVnEP->SetMarkerStyle(mMarker+3);
        gVnEP->SetMarkerColor(1);
        gVnEP->SetMarkerSize(mSize);

        gVnSP->SetMarkerStyle(mMarker+4);
        gVnSP->SetMarkerColor(8);
        gVnSP->SetMarkerSize(mSize);

    }

    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "c1");

    gVn->Draw("AP");
    gVn->GetYaxis()->SetRangeUser(0.0,0.18);
    gVn->GetXaxis()->CenterTitle(true);
    gVn->GetYaxis()->CenterTitle(true);
    gVn->SetTitle(Form("n=%01i; centrality; v_{n}", 2));
    gVn->Draw("AP");
    gVnTrue->Draw("SAME P");
    gVnEP->Draw("SAME P");
    gVnSP->Draw("SAME P");
    hInputFlow->Draw("SAME HIST");
    c1->Update();

    TLegend *leg1 = new TLegend(0.50,0.6,0.85,0.85,"","brNDC");
    leg1->SetTextSize(0.037);leg1->SetBorderSize(0);

    leg1->AddEntry(gVnTrue, "v_{n}, true RP", "p");
    leg1->AddEntry(gVn, "v_{n}, trad. EP", "p");
    leg1->AddEntry(gVnEP, "v_{n}{EP}", "p");
    leg1->AddEntry(gVnSP, "v_{n}{SP}", "p");
    leg1->Draw("SAME");
    c1->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2");

    gR->Draw("AP");
    gR->GetYaxis()->SetRangeUser(0.0,1.1);
    gR->GetXaxis()->CenterTitle(true);
    gR->GetYaxis()->CenterTitle(true);
    gR->SetTitle(Form("n=%01i; centrality; R_{n}", 2));
    gR->Draw("AP");
    gRtrue->Draw("SAME P");
    c2->Update();

    TLegend *leg2 = new TLegend(0.50,0.6,0.85,0.85,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);

    leg2->AddEntry(gRtrue, "R true, uniform #phi", "p");
    leg2->AddEntry(gR, "R sub event method, uniform #phi", "p");
    leg2->Draw("SAME");
    c2->Draw();

}
