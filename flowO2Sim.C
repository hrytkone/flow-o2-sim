#include <iostream>
#include <stdlib.h>
#include <vector>

// OWN
#include "src/JHistos.h"
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
#include "DataFormatsFV0/Hit.h"
#include "DataFormatsFT0/HitType.h"

// OTHER
#include "TStopwatch.h"

using namespace std;

double GetPhi(double x, double y);
double GetPhi0(double phi, double *vn, double *psi);
double AnisotropicPhiDist(double *x, double *p);
double GetAnisotropicPhi(double x, double y, double phiTemp, double err, double *vn, double *psi, TF1 *fPhiDist);

void AnalyzeHits(std::vector<o2::ft0::HitType>* hitArrFT0, std::vector<o2::fv0::Hit>* hitArrFV0, double* psi, double* vn, TF1* fPhiDist, JHistos* histos);

void CalculateQvector(double phi, TComplex unitVec, TComplex &Qvec, double &norm, int n, double w);
double GetEventPlane(TComplex Qvec, int n);
double GetVnObs(TComplex Qvec, double phi, int n);

int main(int argc, char **argv) {

    TString inFileName = argc > 1 ? argv[1]:"o2sim.root";
    TString outFileName = argc > 2 ? argv[2]:"toyFlow.root";
    if(inFileName.EqualTo("help",TString::kIgnoreCase)) {
        cout << "Usage: " << argv[0] << " inputFileName.root" << " outputFileName.root" << endl;
        return 0;
    };

    TStopwatch timer;
    timer.Start();

    TFile *fIn = TFile::Open(inFileName);
    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    fIn->cd();
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

    gRandom->SetSeed(0);
    TRandom3 *rand = new TRandom3(0);

    JHistos *histos = new JHistos();

    TF1 *fPhiDist = new TF1("fPhiDist", AnisotropicPhiDist, -TMath::Pi(), TMath::Pi(), 11);

    // Set branch addresses to containers
    std::vector<o2::MCTrack>* trackArr = nullptr;
    std::vector<o2::ft0::HitType>* hitArrFT0 = nullptr;
    std::vector<o2::fv0::Hit>* hitArrFV0 = nullptr;

    tree->SetBranchAddress("MCTrack", &trackArr);
    tree->SetBranchAddress("FT0Hit", &hitArrFT0);
    tree->SetBranchAddress("FV0Hit", &hitArrFV0);

    double x, y, z, pt, phi, eta;

    int nOutput = nEntries/20;
    if (nOutput<1) nOutput = 1;
    for (int iEvent=0; iEvent<nEntries; iEvent++) {
        if (iEvent % nOutput == 0)
            cout << 100*iEvent/nEntries << " % finished" << endl;

        for (int j=0; j<nCoef; j++) {
            psi[j] = rand->Uniform(-TMath::Pi(), TMath::Pi());
        }

        tree->GetEntry(iEvent);

        // MC tracks
        Int_t arrSize = trackArr->size();
        for (Int_t iTrack = 0; iTrack<arrSize; iTrack++) {

            const auto& mcTrack = (*trackArr)[iTrack];

            x = mcTrack.GetStartVertexCoordinatesX();
            y = mcTrack.GetStartVertexCoordinatesY();

            pt = mcTrack.GetPt();
            histos->hPt->Fill(pt);

            phi = GetPhi(x, y);
            histos->hPhi->Fill(phi);

            phi = GetAnisotropicPhi(x, y, 0.3, 0.01, vn, psi, fPhiDist);
            histos->hAnisotropicPhi->Fill(phi);

            eta = mcTrack.GetRapidity();
            histos->hEta->Fill(eta);
        }

        histos->hMultiplicity->Fill(arrSize);

        AnalyzeHits(hitArrFT0, hitArrFV0, vn, psi, fPhiDist, histos);
    }

    fIn->Close();

    fOut->Write();
    fOut->Close();
    timer.Print();

    return 0;
}

//======END OF MAIN PROGRAM======
void AnalyzeHits(std::vector<o2::ft0::HitType>* hitArrFT0, std::vector<o2::fv0::Hit>* hitArrFV0, double* vn, double* psi, TF1* fPhiDist, JHistos* histos) {

    int i, n, iHit, iDet, arrSize;
    double w = 1.0;
    double x, y, z, phi;
    double EventPlaneA, EventPlaneB, Rsub;

    TComplex Qvec[DET_N], QvecA[DET_N], QvecB[DET_N];
    TComplex unitVec = TComplex(0, 0);
    TComplex autocorr = TComplex(0, 0);

    double norm[DET_N], normA[DET_N], normB[DET_N];
    int nMult[DET_N], nMultA[DET_N], nMultB[DET_N];
    double vobs[DET_N];

    double QnQnA, QnAQnB;

    vector<vector<TComplex>> pTBinsQ;
    pTBinsQ.resize(PTBINS_N);

    for (i=0; i<nCoef; i++) {
        n = i+1;

        for (iDet=0; iDet<DET_N; iDet++) {
            vobs[iDet] = 0;
            Qvec[iDet]= TComplex(0, 0);
            QvecA[iDet] = TComplex(0, 0);
            QvecB[iDet] = TComplex(0, 0);
            norm[iDet]= 0.0;
            normA[iDet] = 0.0;
            normB[iDet] = 0.0;
        }

        // Construct Q-vectors for the detectors and their subevents
        // FT0-A & FT0-C
        arrSize = hitArrFT0->size();
        for (iHit=0; iHit<arrSize; iHit++) {
            const auto& ft0Hit = (*hitArrFT0)[iHit];

            x = ft0Hit.GetX();
            y = ft0Hit.GetY();
            z = ft0Hit.GetZ();

            phi = GetPhi(x, y);
            //phi = GetAnisotropicPhi(x, y, 0.3, 0.01, vn, psi, fPhiDist);

            if (z>0) {
                CalculateQvector(phi, unitVec, Qvec[0], norm[0], n, w);
                nMult[0]++;
                if (iHit%2==0) {
                    CalculateQvector(phi, unitVec, QvecA[0], normA[0], n, w);
                    nMultA[0]++;
                } else {
                    CalculateQvector(phi, unitVec, QvecB[0], normB[0], n, w);
                    nMultB[0]++;
                }
            } else {
                CalculateQvector(phi, unitVec, Qvec[1], norm[1], n, w);
                nMult[1]++;
                if (iHit%2==0) {
                    CalculateQvector(phi, unitVec, QvecA[1], normA[1], n, w);
                    nMultA[1]++;
                } else {
                    CalculateQvector(phi, unitVec, QvecB[1], normB[1], n, w);
                    nMultB[1]++;
                }
            }
        }

        // FV0
        arrSize = hitArrFV0->size();
        for (iHit=0; iHit<arrSize; iHit++) {
            const auto& fv0Hit = (*hitArrFV0)[iHit];

            x = fv0Hit.GetX();
            y = fv0Hit.GetY();

            phi = GetPhi(x, y);
            //phi = GetAnisotropicPhi(x, y, 0.3, 0.01, vn, psi, fPhiDist);

            CalculateQvector(phi, unitVec, Qvec[2], norm[2], n, w);
            nMult[2]++;
            if (iHit%2==0) {
                CalculateQvector(phi, unitVec, QvecA[2], normA[2], n, w);
                nMultA[2]++;
            } else {
                CalculateQvector(phi, unitVec, QvecB[2], normB[2], n, w);
                nMultB[2]++;
            }
        }

        // Calculate vobs
        // FT0
        arrSize = hitArrFT0->size();
        for (iHit=0; iHit<arrSize; iHit++) {
            const auto& ft0Hit = (*hitArrFT0)[iHit];

            x = ft0Hit.GetX();
            y = ft0Hit.GetY();
            z = ft0Hit.GetZ();

            phi = GetPhi(x, y);

            autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

            if (z>0) {
                Qvec[0] -= autocorr;
                vobs[0] += GetVnObs(Qvec[0], phi, n);
                Qvec[0] += autocorr;
            } else {
                Qvec[1] -= autocorr;
                vobs[1] += GetVnObs(Qvec[1], phi, n);
                Qvec[1] += autocorr;
            }
        }

        // FV0
        arrSize = hitArrFV0->size();
        for (iHit=0; iHit<arrSize; iHit++) {
            const auto& fv0Hit = (*hitArrFV0)[iHit];

            x = fv0Hit.GetX();
            y = fv0Hit.GetY();

            phi = GetPhi(x, y);

            autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

            Qvec[2] -= autocorr;
            vobs[2] += GetVnObs(Qvec[2], phi, n);
            Qvec[2] += autocorr;
        }

        for (iDet=0; iDet<DET_N; iDet++) {
            vobs[iDet] /= nMult[iDet];

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

            histos->hVnObs[i][iDet]->Fill(vobs[iDet]);
            histos->hRsub[i][iDet]->Fill(Rsub);
            histos->hQnQnAEP[i][iDet]->Fill(QnQnA/TComplex::Abs(QvecA[iDet]));
            histos->hQnAQnBEP[i][iDet]->Fill(QnAQnB/(TComplex::Abs(QvecA[iDet])*TComplex::Abs(QvecB[iDet])));
            histos->hQnQnASP[i][iDet]->Fill(QnQnA);
            histos->hQnAQnBSP[i][iDet]->Fill(QnAQnB);
        }
    }
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

void CalculateQvector(double phi, TComplex unitVec, TComplex &Qvec, double &norm, int n, double w) {
    norm += w*w;
    unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
    Qvec += unitVec;
}

double GetEventPlane(TComplex Qvec, int n) {
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/n;
}

double GetVnObs(TComplex Qvec, double phi, int n) {
    return TMath::Cos(n*(phi - GetEventPlane(Qvec, n)));
}
