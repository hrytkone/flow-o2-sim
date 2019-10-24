#include <TMath.h>
#include <TFile.h>

const int nsets = 1;

void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

TString fileNames[nsets] = {"toyFlowGraphs.root"};

TH1D *hInputFlow[nsets];

TGraphErrors *gR[nsets];
TGraphErrors *gVn[nsets];
TGraphErrors *gVnEP[nsets];
TGraphErrors *gVnSP[nsets];

TGraphErrors *gVnRatio[nsets];
TGraphErrors *gVnEPRatio[nsets];
TGraphErrors *gVnSPRatio[nsets];

TGraphErrors *gPtBin[nsets];

TF1 *fConst;

int mMarker[nsets] = {24};

double mSize = 1.0;

TFile *fIn[nsets];

double VnDist(double *x, double *p);
void hset(TH1& hid, TString xtit="", TString ytit="",
		double titoffx = 1.1, double titoffy = 1.1,
		double titsizex = 0.06, double titsizey = 0.06,
		double labeloffx = 0.01, double labeloffy = 0.001,
		double labelsizex = 0.05, double labelsizey = 0.05,
		int divx = 505, int divy=505);

void PlotVn() {

    TH2F *hfr;

    fConst = new TF1("fConst","1",0.5,5.5);
    fConst->SetLineColor(kBlack);

    int i, j;
    for (i=0; i<nsets; i++) {
        fIn[i] = TFile::Open(fileNames[i], "read");
        if(fIn[i]==0) ErrorExit(Form("Cannot open file: %s",fileNames[i].Data()));

        hInputFlow[i] = (TH1D*) fIn[i]->Get("hInputFlow");

        gR[i] = (TGraphErrors*) fIn[i]->Get("gR");
        gVn[i] = (TGraphErrors*) fIn[i]->Get("gVn");
        gVnEP[i] = (TGraphErrors*) fIn[i]->Get("gVnEP");
        gVnSP[i] = (TGraphErrors*) fIn[i]->Get("gVnSP");

        gVnRatio[i] = (TGraphErrors*) fIn[i]->Get("gVnRatio");
        gVnEPRatio[i] = (TGraphErrors*) fIn[i]->Get("gVnEPRatio");
        gVnSPRatio[i] = (TGraphErrors*) fIn[i]->Get("gVnSPRatio");

        gR[i]->SetMarkerStyle(mMarker[i]+1);
        gR[i]->SetMarkerColor(i+1);
        gR[i]->SetMarkerSize(mSize);

        gVn[i]->SetMarkerStyle(mMarker[i]); gVnRatio[i]->SetMarkerStyle(mMarker[i]);
        gVn[i]->SetMarkerColor(i+1); gVnRatio[i]->SetMarkerColor(i+1);
        gVn[i]->SetMarkerSize(mSize); gVnRatio[i]->SetMarkerSize(mSize);

        gVnEP[i]->SetMarkerStyle(mMarker[i]); gVnEPRatio[i]->SetMarkerStyle(mMarker[i]);
        gVnEP[i]->SetMarkerColor(i+3); gVnEPRatio[i]->SetMarkerColor(i+3);
        gVnEP[i]->SetMarkerSize(mSize); gVnEPRatio[i]->SetMarkerSize(mSize);

        gVnSP[i]->SetMarkerStyle(mMarker[i]); gVnSPRatio[i]->SetMarkerStyle(mMarker[i]);
        gVnSP[i]->SetMarkerColor(i+4); gVnSPRatio[i]->SetMarkerColor(i+4);
        gVnSP[i]->SetMarkerSize(mSize); gVnSPRatio[i]->SetMarkerSize(mSize);

    }

    gStyle->SetOptStat(0);

    //TCanvas *c1 = new TCanvas("c1", "vn values with different methods - uniform phi");
    //c1->cd();
    Filipad *fpad = new Filipad(0, 1.1, 0.3, 100, 100, 0.8, 5);
    fpad->Draw();
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    double c1Min = 0.5, c1Max = 5.5;

    TLegend *leg1 = new TLegend(0.60,0.55,0.80,0.75,"","brNDC");
    leg1->SetTextSize(0.037);leg1->SetBorderSize(0);

    hfr = new TH2F(Form("hfr1%d",0)," ", 1, c1Min, c1Max, 1, -0.019, 0.18);
    hset( *hfr, "n", "v_{n}",1.4,1.1, 0.07,0.06, 0.01,0.001, 0.03,0.03, 5,510);
    hfr->Draw();

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gVn[i]->Draw("SAME P");
            gVn[i]->SetTitle("Vn values with different methods - uniform phi; Vn; n");
            leg1->AddEntry(gVn[i], "v_{n}, trad. EP", "p");
        }


        gVnEP[i]->Draw("SAME P");
        leg1->AddEntry(gVnEP[i], "v_{n}{EP}", "p");

        gVnSP[i]->Draw("SAME P");
        leg1->AddEntry(gVnSP[i], "v_{n}{SP}", "p");

        hInputFlow[i]->Draw("SAME HIST");
        leg1->AddEntry(hInputFlow[i], "Input", "l");
    }

    leg1->Draw("SAME");
    //c1->SaveAs("figures/vn.pdf");

    p = fpad->GetPad(2); //lower pad
    p->SetTickx(); p->SetGridy(0); p->SetLogy(0);p->SetLogx(0); p->cd();

    hfr = new TH2F(Form("hfr1b%d",0)," ", 1, c1Min, c1Max, 1, 0.85, 1.95);
    hset( *hfr, "n", "ratio to input",1.0,0.7, 0.11,0.09, 0.01,0.001, 0.07,0.07, 5,510);
    hfr->Draw();

    fConst->Draw("SAME");
    for (i=0; i<nsets; i++) {
        if (i==0) {
            gVnRatio[i]->Draw("SAME P");
        }

        gVnEPRatio[i]->Draw("SAME P");
        gVnSPRatio[i]->Draw("SAME P");
    }

    fpad->C->SaveAs("figures/vn-with-ratio.pdf");

    TCanvas *c2 = new TCanvas("c2", "R values with different methods");
    c2->cd();

    hfr = new TH2F(Form("hfr%d",1)," ", 1, 0.5, 5.5, 1, -0.021, 1.04);
    hset( *hfr, "n", "R_{n}",1.0,1.0, 0.04,0.04, 0.01,0.01, 0.03,0.03, 510,510);
    hfr->Draw();

    TLegend *leg2 = new TLegend(0.25,0.15,0.55,0.30,"","brNDC");
    leg2->SetTextSize(0.037);leg2->SetBorderSize(0);

    for (i=0; i<nsets; i++) {
        if (i==0) {
            gR[i]->Draw("SAME P");
            gR[i]->SetTitle("Resolution parameter; R; n");
            leg2->AddEntry(gR[i], "R sub event method", "p");
        }

    }

    leg2->Draw("SAME");
    c2->SaveAs("figures/Rn.pdf");

    TH1D *hPhi = (TH1D*)fIn[0]->Get("hPhi");
    TH1D *hPhiNonuni = (TH1D*)fIn[0]->Get("hPhiNonuni");

    TCanvas *c3 = new TCanvas("c3", "Phi distribution");
    hPhi->Draw("HIST");
    c3->SaveAs("figures/phi.pdf");
}

double VnDist(double *x, double *p) {
    double pt = x[0];
    double alpha = p[0];
    double beta = p[1];
    double vnMax = p[2];
    double C = vnMax/(TMath::Power(alpha/beta, alpha)*TMath::Exp(-alpha));
    return C*TMath::Power(pt, alpha)*TMath::Exp(-beta*pt);
}

void hset(TH1& hid, TString xtit="", TString ytit="",
		double titoffx = 1.1, double titoffy = 1.1,
		double titsizex = 0.06, double titsizey = 0.06,
		double labeloffx = 0.01, double labeloffy = 0.001,
		double labelsizex = 0.05, double labelsizey = 0.05,
		int divx = 505, int divy=505)
{
	hid.GetXaxis()->CenterTitle(1);
	hid.GetYaxis()->CenterTitle(1);

	hid.GetXaxis()->SetTitleOffset(titoffx);
	hid.GetYaxis()->SetTitleOffset(titoffy);

	hid.GetXaxis()->SetTitleSize(titsizex);
	hid.GetYaxis()->SetTitleSize(titsizey);

	hid.GetXaxis()->SetLabelOffset(labeloffx);
	hid.GetYaxis()->SetLabelOffset(labeloffy);

	hid.GetXaxis()->SetLabelSize(labelsizex);
	hid.GetYaxis()->SetLabelSize(labelsizey);

	hid.GetXaxis()->SetNdivisions(divx);
	hid.GetYaxis()->SetNdivisions(divy);

	hid.GetXaxis()->SetTitle(xtit);
	hid.GetYaxis()->SetTitle(ytit);

	hid.GetXaxis()->SetLabelFont(42);
	hid.GetYaxis()->SetLabelFont(42);
	hid.GetXaxis()->SetTitleFont(42);
	hid.GetYaxis()->SetTitleFont(42);
}
