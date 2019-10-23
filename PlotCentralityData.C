#include <TMath.h>
#include <TFile.h>

#include "src/JConst.h"

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

const int nFiles = 3;
const int nRef = 0;
TString fileName[nFiles] = {
                            "toyFlowCentralityGraphs.root"
                           ,"toyFlowCentralityGraphsDoubleGran.root"
                           ,"toyFlowCentralityGraphsGran.root"
                           };
TString sSame[nFiles] = {"AP", "SAME P", "SAME P"};
int gColor[nFiles] = {1,2,4};
bool drawToSameFig = false;


// You can choose how many harmonics are drawn
// by changing this number. Max=4
const int nHarmonics = 4;
int nXH[5] = {0, 1, 2, 2, 2};
int nYH[5] = {0, 1, 1, 2, 2};
int nXP[5] = {0, 500, 1000, 1000, 1000};
int nYP[5] = {0, 500, 500, 1000, 1000};

TH1D *hInputFlow[nCoef][nFiles];

TGraphErrors *gR[nCoef][nFiles];
TGraphErrors *gRtrue[nCoef][nFiles];
TGraphErrors *gVn[nCoef][nFiles];
TGraphErrors *gVnTrue[nCoef][nFiles];
TGraphErrors *gVnEP[nCoef][nFiles];
TGraphErrors *gVnSP[nCoef][nFiles];

int mMarker = 24;
double mSize = 1.0;

TFile *fIn[nFiles];

void PlotCentralityData(int readNFiles = 1) {

    for(int iFil=0; iFil<readNFiles; iFil++)
        fIn[iFil] = TFile::Open(fileName[iFil], "read");

    int i, j;
    for(int iFil=0; iFil<readNFiles; iFil++) {
        for (i=0; i<nCoef; i++) {
            if(fIn[iFil]==0) ErrorExit(Form("Cannot open file: %s",fileName[iFil].Data()));

            hInputFlow[i][iFil] = (TH1D*) fIn[iFil]->Get(Form("hInputFlow%02i", i+1));

            gR[i][iFil] = (TGraphErrors*) fIn[iFil]->Get(Form("gRH%02i", i+1));
            gRtrue[i][iFil] = (TGraphErrors*) fIn[iFil]->Get(Form("gRtrueH%02i", i+1));

            gVn[i][iFil] = (TGraphErrors*) fIn[iFil]->Get(Form("gVnH%02i", i+1));
            gVnTrue[i][iFil] = (TGraphErrors*) fIn[iFil]->Get(Form("gVnTrueH%02i", i+1));
            gVnEP[i][iFil] = (TGraphErrors*) fIn[iFil]->Get(Form("gVnEPH%02i", i+1));
            gVnSP[i][iFil] = (TGraphErrors*) fIn[iFil]->Get(Form("gVnSPH%02i", i+1));

            gR[i][iFil]->SetMarkerStyle(mMarker);
            gR[i][iFil]->SetMarkerColor(gColor[iFil]);
            gR[i][iFil]->SetMarkerSize(mSize);
            gR[i][iFil]->SetFillColor(gColor[iFil]);

            gRtrue[i][iFil]->SetMarkerStyle(mMarker+1);
            gRtrue[i][iFil]->SetMarkerColor(gColor[iFil]);
            gRtrue[i][iFil]->SetMarkerSize(mSize);
            gRtrue[i][iFil]->SetFillColor(gColor[iFil]);

            gVn[i][iFil]->SetMarkerStyle(mMarker);
            gVn[i][iFil]->SetMarkerColor(gColor[iFil]);
            gVn[i][iFil]->SetMarkerSize(mSize);
            gVn[i][iFil]->SetFillColor(gColor[iFil]);

            gVnTrue[i][iFil]->SetMarkerStyle(mMarker+1);
            gVnTrue[i][iFil]->SetMarkerColor(gColor[iFil]);
            gVnTrue[i][iFil]->SetMarkerSize(mSize);
            gVnTrue[i][iFil]->SetFillColor(gColor[iFil]);

            gVnEP[i][iFil]->SetMarkerStyle(mMarker+3);
            gVnEP[i][iFil]->SetMarkerColor(gColor[iFil]);
            gVnEP[i][iFil]->SetMarkerSize(mSize);
            gVnEP[i][iFil]->SetFillColor(gColor[iFil]);

            gVnSP[i][iFil]->SetMarkerStyle(mMarker+4);
            gVnSP[i][iFil]->SetMarkerColor(gColor[iFil]);
            gVnSP[i][iFil]->SetMarkerSize(mSize);
            gVnSP[i][iFil]->SetFillColor(gColor[iFil]);
        }
    }

    gStyle->SetOptStat(0);

    TCanvas *c1;
    if(!drawToSameFig) {
        c1 = new TCanvas("c1", "c1", nXP[nHarmonics], nYP[nHarmonics]);
        c1->Divide(nXH[nHarmonics],nYH[nHarmonics]);
    } else {
        c1 = new TCanvas("c1", "c1", nXP[1], nYP[1]);
    }


    for (i=1; i<nHarmonics+1; i++) {
        for(int iFil=0; iFil<readNFiles; iFil++) {
            if(!drawToSameFig) c1->cd(i);
            gVn[i][iFil]->GetYaxis()->SetRangeUser(0.0,0.11);
            gVn[i][iFil]->GetXaxis()->CenterTitle(true);
            gVn[i][iFil]->GetYaxis()->CenterTitle(true);
            gVn[i][iFil]->SetTitle(Form("n=%01i; centrality; v_{n}", i+1));
            if(drawToSameFig && i!=1) gVn[i][iFil]->Draw("SAME P");
            else     gVn[i][iFil]->Draw(Form("%s",sSame[iFil].Data()));
            gVnTrue[i][iFil]->Draw("SAME P");
            gVnEP[i][iFil]->Draw("SAME P");
            gVnSP[i][iFil]->Draw("SAME P");
            hInputFlow[i][iFil]->Draw("SAME HIST");
            c1->Update();
        }
    }

    TLegend *legVn[nHarmonics+1];
    if(drawToSameFig) {
        double xl[nHarmonics+1] = {0, 0.301205, 0.423695, 0.483936, 0.485944};
        double yl[nHarmonics+1] = {0, 0.67019,  0.350951, 0.230444, 0.150106};
        double xh[nHarmonics+1] = {0, 0.341365, 0.463855, 0.524096, 0.534137};
        double yh[nHarmonics+1] = {0, 0.725159, 0.40592,  0.285412, 0.205074};
        for(i=1; i<nHarmonics+1; i++) {
            legVn[i] = new TLegend(xl[i],yl[i],xh[i],yh[i],"","brNDC");
            legVn[i]->SetTextSize(0.040);legVn[i]->SetBorderSize(0);
            legVn[i]->AddEntry((TH1D*)NULL, Form("v_{%d}",i+1), "");
            legVn[i]->Draw("same");
        }
    }

    TLegend *leg1;
    if(!drawToSameFig) {
        leg1 = new TLegend(0.30,0.5,0.55,0.85,"","brNDC");
    } else {
        leg1 = new TLegend(0.48996,0.348837,0.626506,0.693446,"","brNDC");
    }
    leg1->SetTextSize(0.037);leg1->SetBorderSize(0);

    leg1->AddEntry(gVnTrue[0][nRef], "v_{n}, true RP", "p");
    leg1->AddEntry(gVn[0][nRef], "v_{n}, trad. EP", "p");
    leg1->AddEntry(gVnEP[0][nRef], "v_{n}{EP}", "p");
    leg1->AddEntry(gVnSP[0][nRef], "v_{n}{SP}", "p");
    if(readNFiles>1) leg1->AddEntry(gR[0][0], "Ideal detector", "f");
    if(readNFiles>1) leg1->AddEntry(gR[0][1], "16 sec granular detector", "f");
    if(readNFiles>2) leg1->AddEntry(gR[0][2], "8 sec granular detector", "f");
    leg1->Draw("SAME");
    c1->Draw();

    TCanvas *c2 = new TCanvas("c2", "c2", nXP[nHarmonics], nYP[nHarmonics]);
    c2->Divide(nXH[nHarmonics],nYH[nHarmonics]);

    for (i=1; i<nHarmonics+1; i++) {
        for(int iFil=0; iFil<readNFiles; iFil++) {
            c2->cd(i);
            gR[i][iFil]->GetYaxis()->SetRangeUser(0.35,1.1);
            gR[i][iFil]->GetXaxis()->CenterTitle(true);
            gR[i][iFil]->GetYaxis()->CenterTitle(true);
            gR[i][iFil]->SetTitle(Form("n=%01i; centrality; R_{n}", i+1));
            gR[i][iFil]->Draw(Form("%s",sSame[iFil].Data()));
            gRtrue[i][iFil]->Draw("SAME P");
            c2->Update();
        }
    }

    TLegend *leg2 = new TLegend(0.30,0.7,0.55,0.85,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);

    leg2->AddEntry(gRtrue[0][nRef], "R true", "p");
    leg2->AddEntry(gR[0][nRef], "R sub event method", "p");
    if(readNFiles>1) leg2->AddEntry(gR[0][0], "Ideal detector", "f");
    if(readNFiles>1) leg2->AddEntry(gR[0][1], "16 sec granular detector", "f");
    if(readNFiles>2) leg2->AddEntry(gR[0][2], "8 sec granular detector", "f");
    leg2->Draw("SAME");
    c2->Draw();

}
