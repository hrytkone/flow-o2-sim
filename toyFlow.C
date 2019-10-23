#include <iostream>
#include <stdlib.h>
#include <vector>

// OWN
#include "src/JHistos.h"
#include "src/JEventLists.h"
#include "src/JToyMCTrack.h"
#include "src/JInputs.h"
#include "src/JConst.h"

// ROOT
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TTree.h"

// O2
#include "SimulationDataFormat/MCTrack.h"

// OTHER
#include "TStopwatch.h"

using namespace std;

double GetPhi(double x, double y);
double GetPhi0(double phi, double *vn, double *psi);
double AnisotropicPhiDist(double *x, double *p);
double GetAnisotropicPhi(double x, double y, double phiTemp, double err, double *vn, double *psi, TF1 *fPhiDist);
double GetEta(double px, double py, double pz);

void GetEvent(JHistos *histos, JEventLists *lists, TTree *tree, int iEvent, double *vn, double *psi, TF1 *fPhiDist);
void GetParticleLists(JEventLists *lists);
void AnalyzeEvent(JHistos *histos, JEventLists *lists, bool bUseWeight);

void CalculateQvector(JToyMCTrack *track, TComplex unitVec, TComplex &Qvec, double &norm, int n, double w, bool bUseWeight);
double GetEventPlane(TComplex Qvec, int n);
double GetVnObs(TComplex Qvec, double phi, int n);

double BelongsToA(double phi);

int main(int argc, char **argv) {

    TString inFileName = argc > 1 ? argv[1]:"o2sim.root";
    TString outFileName = argc > 2 ? argv[2]:"toyFlow.root";
    if(inFileName.EqualTo("help",TString::kIgnoreCase)) {
        cout << "Usage: " << argv[0] << " inputFileName.root" << " outputFileName.root" << endl;
        return 0;
    };

    bool bUseWeight = false;

    TStopwatch timer;
    timer.Start();

    TFile *fIn = new TFile(inFileName);
    TTree *tree = (TTree*)fIn->Get("o2sim");
    Int_t nEntries = (Int_t)tree->GetEntries();

    const double scale = 1.0;
    double vn[nCoef] = {scale*0.0, scale*0.15, scale*0.08, scale*0.0, scale*0.0};
    double psi[nCoef] = {0};

    cout << "=========================================== Settings ===========================================" << endl;
    cout << "Input: " << inFileName.Data()
    << ", Events: " << nEntries << endl;
    cout << "Output: " << outFileName.Data() << endl;
    cout << "Vn inputs: ";
    for(int i=0; i<nCoef; i++) cout << vn[i] << ", ";
    cout << endl;
    cout << "================================================================================================" << endl;

    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    gRandom->SetSeed(0);
    TRandom3 *rand = new TRandom3(0);

    JHistos *histos = new JHistos();
    JEventLists *lists = new JEventLists();
    JInputs *inputs = new JInputs();
    inputs->Load();

    // Save input numbers
    TH1D *hInputNumbers = new TH1D("hInputNumbers","hInputNumbers",7, 0.5, 7.5);
    hInputNumbers->Fill(1, double(nEntries));
    hInputNumbers->Fill(2, vn[0]);
    hInputNumbers->Fill(3, vn[1]);
    hInputNumbers->Fill(4, vn[2]);
    hInputNumbers->Fill(5, vn[3]);
    hInputNumbers->Fill(6, vn[4]);
    hInputNumbers->Fill(7, 1.0); // Counting number of files added with hadd.
    fOut->cd();
    hInputNumbers->Write("hInputNumbers");

    TF1 *fPhiDist = new TF1("fPhiDist", AnisotropicPhiDist, -TMath::Pi(), TMath::Pi(), 11);

    int nOutput = nEntries/20;
    if (nOutput<1) nOutput = 1;
    for (int iEvent=0; iEvent<nEntries; iEvent++) {
        if (iEvent % nOutput == 0)
            cout << 100*iEvent/nEntries << " % finished" << endl;

        lists->ClearLists();

        for (int j=0; j<nCoef; j++) {
            psi[j] = rand->Uniform(-TMath::Pi(), TMath::Pi());
        }

        GetEvent(histos, lists, tree, iEvent, vn, psi, fPhiDist);
        GetParticleLists(lists);
        AnalyzeEvent(histos, lists, bUseWeight);
    }

    fOut->Write();
    fOut->Close();
    timer.Print();

    return 0;
}

//======END OF MAIN PROGRAM======
void GetMCEvent(JHistos *histos, JEventLists *lists, TTree *tree, int iEvent, double *vn, double *psi, TF1 *fPhiDist) {
    double pt, phi, eta, energy;
    double x, y, px, py, pz;

    JToyMCTrack track;
    TLorentzVector lVec;

    int nTracks = 0;

    std::vector<o2::MCTrack>* trackArr = nullptr;
    tree->SetBranchAddress("MCTrack", &trackArr);

    tree->GetEntry(iEvent);
    Int_t arrSize = trackArr->size();
    for (Int_t iTrack = 0; iTrack<arrSize; iTrack++) {

        const auto& mcTrack = (*trackArr)[iTrack];

        x = mcTrack.GetStartVertexCoordinatesX();
        y = mcTrack.GetStartVertexCoordinatesY();

        px = mcTrack.GetStartVertexMomentumX();
        py = mcTrack.GetStartVertexMomentumY();
        pz = mcTrack.GetStartVertexMomentumZ();

        energy = mcTrack.GetEnergy();

        pt = mcTrack.GetPt();
        histos->hPt->Fill(pt);

        phi = GetPhi(x, y);
        histos->hPhi->Fill(phi);

        phi = GetAnisotropicPhi(x, y, 0.3, 0.01, vn, psi, fPhiDist);
        histos->hAnisotropicPhi->Fill(phi);

        eta = GetEta(px, py, pz);
        histos->hEta->Fill(eta);

        lVec.SetPxPyPzE(px, py, pz, energy);
        track.SetTrack(lVec);


        lVec.SetPxPyPzE(px, py, pz, energy);
        track.SetTrack(lVec);

        new((*lists->fullEvent)[nTracks]) JToyMCTrack(track);
        nTracks++;
    }

    histos->hMultiplicity->Fill(nTracks);
}

void GetParticleLists(JEventLists *lists) {
    int nMult = lists->fullEvent->GetEntriesFast();

    JToyMCTrack *tempTrack, track, trackA, trackB;
    double eta, phi;
    int detMult[DET_N] = {0};
    int detMultA[DET_N] = {0};
    int detMultB[DET_N] = {0};
    //double detCenter = 0.0;

    int i, j;
    for (i=0; i<nMult; i++) {
        tempTrack = (JToyMCTrack*)lists->fullEvent->At(i);
        eta = tempTrack->GetEta();

        for (j=0; j<DET_N; j++) {
            if (cov[j][0]<eta && eta<cov[j][1]) {
                phi = tempTrack->GetPhi();

                track = *tempTrack;
                new((*lists->GetList(j))[detMult[j]]) JToyMCTrack(track);
                detMult[j]++;

                //detCenter = (cov[j][0]+cov[j][1])/2.0;
                if(/*eta<detCenter*/ i%2==0 /*BelongsToA(phi)*/) { //Later use the function.
                    trackA = *tempTrack;
                    new((*lists->GetList(j,"A"))[detMultA[j]]) JToyMCTrack(trackA);
                    detMultA[j]++;
                } else {
                    trackB = *tempTrack;
                    new((*lists->GetList(j,"B"))[detMultB[j]]) JToyMCTrack(trackB);
                    detMultB[j]++;
                }
            }
        }
    }
}

void AnalyzeEvent(JHistos *histos, JEventLists *lists, bool bUseWeight) {

    int nMult[DET_N], nMultA[DET_N], nMultB[DET_N];

    for(int iDet=0; iDet<DET_N; iDet++) {
        nMult[iDet] = lists->GetList(iDet)->GetEntriesFast();
        nMultA[iDet] = lists->GetList(iDet,"A")->GetEntriesFast();
        nMultB[iDet] = lists->GetList(iDet,"B")->GetEntriesFast();
    }

    int i, j, k, n, jMax;
    double w = 1.0;
    double phi, pt;
    double EventPlaneA, EventPlaneB, Rsub, vobs;

    TComplex Qvec[DET_N], QvecA[DET_N], QvecB[DET_N];
    TComplex unitVec = TComplex(0, 0);
    TComplex autocorr = TComplex(0, 0);

    int centBin = 0;

    double norm[DET_N], normA[DET_N], normB[DET_N];

    double QnQnA, QnAQnB;

    vector<vector<TComplex>> pTBinsQ;
    pTBinsQ.resize(PTBINS_N);

    JToyMCTrack *track;;

    for (i=0; i<nCoef; i++) {

        n = i+1;
        for(int iDet=0; iDet<DET_N; iDet++) {

            if (nMult[iDet] == 0 || nMultA[iDet] == 0 || nMultB[iDet] == 0) continue;

            Qvec[iDet]= TComplex(0, 0);
            QvecA[iDet] = TComplex(0, 0);
            QvecB[iDet] = TComplex(0, 0);

            vobs = 0.0;

            norm[iDet]= 0.0; normA[iDet] = 0.0; normB[iDet] = 0.0;

            // Construct Q-vectors for the detectors
            for (j=0; j<nMult[iDet]; j++) {
                track = (JToyMCTrack*)lists->GetList(iDet)->At(j);
                CalculateQvector(track, unitVec, Qvec[iDet], norm[iDet], n, w, bUseWeight);
            }

            if (nMultA[iDet]<nMultB[iDet]) {
                jMax = nMultA[iDet];
            } else {
                jMax = nMultB[iDet];
            }

            for (j=0; j<jMax; j++) {
            //for (j=0; j<nMult[1]; j++) {
                //track = (JToyMCTrack*)lists->GetList(1)->At(j);
                track = (JToyMCTrack*)lists->GetList(iDet,"A")->At(j);
                CalculateQvector(track, unitVec, QvecA[iDet], normA[iDet], n, w, bUseWeight);
            //}
            //for (j=0; j<nMult[2]; j++) {
                //track = (JToyMCTrack*)lists->GetList(2)->At(j);
                track = (JToyMCTrack*)lists->GetList(iDet,"B")->At(j);
                CalculateQvector(track, unitVec, QvecB[iDet], normB[iDet], n, w, bUseWeight);
            }

            // Calculate vobs from TPC events
            for (j=0; j<nMult[iDet]; j++) {

                track = (JToyMCTrack*)lists->GetList(iDet)->At(j);
                phi = track->GetPhi();
                pt = track->GetPt();

                if (bUseWeight) w = pt;

                autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

                Qvec[iDet] -= autocorr;
                vobs += GetVnObs(Qvec[iDet], phi, n);
                Qvec[iDet] += autocorr;
            }

            vobs /= nMult[iDet];

            // Resolution parameter calculations
            //Rtrue = TMath::Cos(n*(GetEventPlane(Qvec[iDet], n) - Psi[i]));

            EventPlaneA = GetEventPlane(QvecA[iDet], n);
            EventPlaneB = GetEventPlane(QvecB[iDet], n);
            Rsub = TMath::Cos(n*(EventPlaneA - EventPlaneB));

            // EP-method and SP-method
            norm[iDet] = TMath::Sqrt(norm[iDet]);
            normA[iDet] = TMath::Sqrt(normA[iDet]);
            normB[iDet] = TMath::Sqrt(normB[iDet]);

            Qvec[iDet] /= norm[iDet]; QvecA[iDet] /= normA[iDet]; QvecB[iDet] /= normB[iDet];

            QnQnA = Qvec[iDet]*TComplex::Conjugate(QvecA[iDet]);
            QnAQnB = QvecA[iDet]*TComplex::Conjugate(QvecB[iDet]);

            histos->hVnObs[i][iDet][centBin]->Fill(vobs);
            histos->hRsub[i][iDet][centBin]->Fill(Rsub);
            histos->hQnQnAEP[i][iDet][centBin]->Fill(QnQnA/TComplex::Abs(QvecA[iDet]));
            histos->hQnAQnBEP[i][iDet][centBin]->Fill(QnAQnB/(TComplex::Abs(QvecA[iDet])*TComplex::Abs(QvecB[iDet])));
            histos->hQnQnASP[i][iDet][centBin]->Fill(QnQnA);
            histos->hQnAQnBSP[i][iDet][centBin]->Fill(QnAQnB);
        }

        // Divide into pT-bins
        if (n==2) {
            double weight = 0.0;
            double norms[PTBINS_N];
            for (j=0; j<PTBINS_N; j++) norms[j] = 0;

            for (j=0; j<nMult[0]; j++) { //Only do this for TPC.

                track = (JToyMCTrack*)lists->TPClist->At(j);
                phi = track->GetPhi();
                pt = track->GetPt();

                if (bUseWeight) w = pt;

                unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

                for (k=0; k<PTBINS_N; k++) {
                    if ((pTBins[k] <= pt) && (pTBins[k+1] > pt)) {
                        pTBinsQ[k].push_back(unitVec);
                        norms[k] += w*w;
                    }
                }
            }

            double l;
            for (j=0; j<PTBINS_N; j++) {
                Qvec[0] = TComplex(0, 0);
                l = pTBinsQ[j].size();
                for (k=0; k<l; k++) {
                    Qvec[0] += pTBinsQ[j][k];
                }
                weight = TMath::Sqrt(norms[j]);
                if (weight!=0) Qvec[0] /= weight;
                histos->hV2ComplexPart->Fill(Qvec[0].Im()*Qvec[1].Im());
                QnQnA = Qvec[0]*TComplex::Conjugate(Qvec[1]);
                QnQnA /= TComplex::Abs(Qvec[1]);
                histos->hQnQnAPtBin[j]->Fill(QnQnA);
                histos->hSqrtSumWeightsPtBins[j]->Fill(weight);
                pTBinsQ[j].clear();
            }

        }
    }

    for(int iDet=0; iDet<DET_N; iDet++)
        histos->hSqrtSumWeights[iDet][centBin]->Fill(norm[iDet]);
}

double GetPhi(double x, double y) {
    return TMath::ATan2(y, x);
}

double GetPhi0(double phi, double *vn, double *psi) {
    return phi - 2.0*vn[0]*TMath::Sin(phi-psi[0]) + vn[1]*TMath::Sin(2.0*(phi-psi[1])) + (2.0/3.0)*vn[2]*TMath::Sin(3.0*(phi-psi[2])) + (1.0/2.0)*vn[3]*TMath::Sin(4.0*(phi-psi[3])) + (2.0/5.0)*vn[4]*TMath::Sin(5.0*(phi-psi[4]));
}

double AnisotropicPhiDist(double *x, double *p) {
    double phi = x[0];
    double phi0 = p[0];
    double v1 = p[1];
    double v2 = p[2];
    double v3 = p[3];
    double v4 = p[4];
    double v5 = p[5];
    double psi1 = p[6];
    double psi2 = p[7];
    double psi3 = p[8];
    double psi4 = p[9];
    double psi5 = p[10];
    return phi - phi0 + 2.0*v1*TMath::Sin(phi-psi1) + v2*TMath::Sin(2.0*(phi-psi2)) + (2.0/3.0)*v3*TMath::Sin(3.0*(phi-psi3)) + (1.0/2.0)*v4*TMath::Sin(4.0*(phi-psi4)) + (2.0/5.0)*v5*TMath::Sin(5.0*(phi-psi5));
}

double GetAnisotropicPhi(double x, double y, double phiTemp, double err, double *vn, double *psi, TF1 *fPhiDist) {
    double phi0 = GetPhi(x, y);
    double phi = 0;

    fPhiDist->SetParameters(phi0, vn[0], vn[1], vn[2], vn[3], vn[4], psi[0], psi[1], psi[2], psi[3], psi[4]);

    while (TMath::Abs(GetPhi0(phi, vn, psi) - phi0) > err) {
        phi = phiTemp - fPhiDist->Eval(phiTemp)/fPhiDist->Derivative(phiTemp);
        phiTemp = phi;
    }

    if (phi>TMath::Pi())phi -= 2*TMath::Pi();
    if (phi<-TMath::Pi()) phi += 2*TMath::Pi();

    return phi;
}

double GetEta(double px, double py, double pz) {
    double p = TMath::Sqrt(px*px + py*py + pz*pz);
    return TMath::ATanH(pz/p);
}

void CalculateQvector(JToyMCTrack *track, TComplex unitVec, TComplex &Qvec, double &norm, int n, double w, bool bUseWeight) {

    double phi = track->GetPhi();
    double pt = track->GetPt();

    if (bUseWeight) w = pt;
    norm += w*w;

    unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
    //if (bNonuniformPhi) DoCorrections(unitVec, cm, sm, lambdaMinus, lambdaPlus, aMinus, aPlus);
    Qvec += unitVec;
}

double GetEventPlane(TComplex Qvec, int n) {
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/n;
}

double GetVnObs(TComplex Qvec, double phi, int n) {
    return TMath::Cos(n*(phi - GetEventPlane(Qvec, n)));
}

double BelongsToA(double phi) {
    return phi>0;
}
