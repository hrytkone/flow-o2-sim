
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
#include "TLeaf.h"

// O2
#include "SimulationDataFormat/MCTrack.h"
#include "DataFormatsFV0/Hit.h"
#include "DataFormatsFT0/HitType.h"
#include "FT0Base/Geometry.h"

// OTHER
#include "TStopwatch.h"

using namespace std;

void AnalyzeHits(vector<o2::ft0::HitType>* hitArrFT0, vector<o2::fv0::Hit>* hitArrFV0, TLeaf* ft0ChID, TLeaf* fv0ChID, JHistos* histos, bool bUseGranularity);

double GetPhi(double x, double y);
double getFT0ChannelCenterX(int chID);
double getFT0ChannelCenterY(int chID);
double getFV0Phi(int chID);
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

    bool bUseGranularity = 1;

    TStopwatch timer;
    timer.Start();

    TFile *fIn = TFile::Open(inFileName);
    TFile *fOut = TFile::Open(outFileName, "RECREATE");

    fIn->cd();
    TTree *tree = (TTree*)fIn->Get("o2sim");
    Int_t nEntries = (Int_t)tree->GetEntries();

    cout << "=========================================== Settings ===========================================" << endl;
    cout << "Input: " << inFileName.Data()
    << ", Events: " << nEntries << endl;
    cout << "Output: " << outFileName.Data() << endl;
    cout << endl;
    cout << "================================================================================================" << endl;

    // Save input numbers
    TH1D *hInputNumbers = new TH1D("hInputNumbers", "hInputNumbers", 7, 0.5, 3.5);
    hInputNumbers->Fill(1, double(nEntries));
    hInputNumbers->Fill(2, bUseGranularity);
    hInputNumbers->Fill(3, 1.0); // Counting number of files added with hadd.
    fOut->cd();
    hInputNumbers->Write("hInputNumbers");

    JHistos *histos = new JHistos();

    // Set branch addresses to containers
    vector<o2::MCTrack>* trackArr = nullptr;
    vector<o2::ft0::HitType>* hitArrFT0 = nullptr;
    vector<o2::fv0::Hit>* hitArrFV0 = nullptr;

    tree->SetBranchAddress("MCTrack", &trackArr);
    tree->SetBranchAddress("FT0Hit", &hitArrFT0);
    tree->SetBranchAddress("FV0Hit", &hitArrFV0);

    TH1F* hFV0ChMult = new TH1F("hFV0ChMult", "hFV0ChMult", 41, 0, 41);
    tree->Draw("FV0Hit.mDetectorID>>hFV0ChMult");

    TLeaf* ft0ChID;
    TLeaf* fv0ChID;

    int nOutput = nEntries/20;
    if (nOutput<1) nOutput = 1;
    for (int iEvent=0; iEvent<nEntries; iEvent++) {
        if (iEvent % nOutput == 0)
            cout << 100*iEvent/nEntries << " % finished" << endl;

        tree->GetEntry(iEvent);

        // MC tracks
        Int_t arrSize = trackArr->size();
        histos->hMultiplicity->Fill(arrSize);

        ft0ChID = tree->GetLeaf("FT0Hit.mDetectorID");
        fv0ChID = tree->GetLeaf("FV0Hit.mDetectorID");

        AnalyzeHits(hitArrFT0, hitArrFV0, ft0ChID, fv0ChID, histos, bUseGranularity);
    }

    fIn->Close();

    fOut->Write();
    fOut->Close();
    timer.Print();

    return 0;
}

//======END OF MAIN PROGRAM======
void AnalyzeHits(vector<o2::ft0::HitType>* hitArrFT0, vector<o2::fv0::Hit>* hitArrFV0, TLeaf* ft0ChID, TLeaf* fv0ChID, JHistos* histos, bool bUseGranularity) {

    int i, n, iHit, iDet, arrSize, chID;
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

    for (i=0; i<nCoef; i++) {
        n = i+1;

        for (iDet=0; iDet<DET_N; iDet++) {
            nMult[iDet] = 0;
            vobs[iDet] = 0;
            Qvec[iDet]= TComplex(0, 0);
            QvecA[iDet] = TComplex(0, 0);
            QvecB[iDet] = TComplex(0, 0);
            norm[iDet]= 0.0;
            normA[iDet] = 0.0;
            normB[iDet] = 0.0;
        }

        //______________________________________________________________________
        // Construct Q-vectors for the detectors and their subevents
        // FT0-A & FT0-C
        arrSize = hitArrFT0->size();
        for (iHit=0; iHit<arrSize; iHit++) {
            const auto& ft0Hit = (*hitArrFT0)[iHit];

            if (bUseGranularity) {
                chID = ft0ChID->GetValue(iHit);
                x = getFT0ChannelCenterX(chID);
                y = getFT0ChannelCenterY(chID);
            } else {
                x = ft0Hit.GetX();
                y = ft0Hit.GetY();
            }

            z = ft0Hit.GetZ();
            phi = GetPhi(x, y);

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

            if (bUseGranularity) {
                chID = fv0ChID->GetValue(iHit);
                phi = getFV0Phi(chID);
            } else {
                x = fv0Hit.GetX();
                y = fv0Hit.GetY();
                phi = GetPhi(x, y);
            }

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

        //______________________________________________________________________
        // Calculate vobs
        // FT0
        arrSize = hitArrFT0->size();
        for (iHit=0; iHit<arrSize; iHit++) {
            const auto& ft0Hit = (*hitArrFT0)[iHit];

            if (bUseGranularity) {
                chID = ft0ChID->GetValue(iHit);
                x = getFT0ChannelCenterX(chID);
                y = getFT0ChannelCenterY(chID);
            } else {
                x = ft0Hit.GetX();
                y = ft0Hit.GetY();
            }

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

            if (bUseGranularity) {
                chID = fv0ChID->GetValue(iHit);
                phi = getFV0Phi(chID);
            } else {
                x = fv0Hit.GetX();
                y = fv0Hit.GetY();
                phi = GetPhi(x, y);
            }

            autocorr = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));

            Qvec[2] -= autocorr;
            vobs[2] += GetVnObs(Qvec[2], phi, n);
            Qvec[2] += autocorr;
        }

        //______________________________________________________________________
        for (iDet=0; iDet<DET_N; iDet++) {
            vobs[iDet] /= nMult[iDet];

            EventPlaneA = GetEventPlane(QvecA[iDet], n);
            EventPlaneB = GetEventPlane(QvecB[iDet], n);
            //EventPlaneA = GetEventPlane(Qvec[0], n);
            //EventPlaneB = GetEventPlane(Qvec[1], n);
            Rsub = TMath::Cos(n*(EventPlaneA - EventPlaneB));

            // EP-method and SP-method
            norm[iDet] = TMath::Sqrt(norm[iDet]);
            normA[iDet] = TMath::Sqrt(normA[iDet]);
            normB[iDet] = TMath::Sqrt(normB[iDet]);
            //normA[iDet] = TMath::Sqrt(norm[0]);
            //normB[iDet] = TMath::Sqrt(norm[1]);

            Qvec[iDet] /= norm[iDet]; QvecA[iDet] /= normA[iDet]; QvecB[iDet] /= normB[iDet];

            histos->hVnObs[i][iDet]->Fill(vobs[iDet]);
            histos->hRsub[i][iDet]->Fill(Rsub);
            histos->hQnA[i][iDet]->Fill(QvecA[iDet]);
            histos->hQnB[i][iDet]->Fill(QvecB[iDet]);
        }
    }

    for(int iDet=0; iDet<DET_N; iDet++)
        histos->hSqrtSumWeights[iDet]->Fill(norm[iDet]);
}

//______________________________________________________________________
double GetPhi(double x, double y) {
    return TMath::ATan2(y, x);
}

//______________________________________________________________________
double GetEta(double px, double py, double pz) {
    double p = TMath::Sqrt(px*px + py*py + pz*pz);
    return TMath::ATanH(pz/p);
}

//______________________________________________________________________
void CalculateQvector(double phi, TComplex unitVec, TComplex &Qvec, double &norm, int n, double w) {
    norm += w*w;
    unitVec = TComplex(w*TMath::Cos(n*phi), w*TMath::Sin(n*phi));
    Qvec += unitVec;
}

//______________________________________________________________________
double GetEventPlane(TComplex Qvec, int n) {
    return TMath::ATan2(Qvec.Im(), Qvec.Re())/n;
}

//______________________________________________________________________
double GetVnObs(TComplex Qvec, double phi, int n) {
    return TMath::Cos(n*(phi - GetEventPlane(Qvec, n)));
}

//______________________________________________________________________
// FT0-A: Channels 0-95
// FT0-C: Channels 96-207
double getFT0ChannelCenterX(int chID) {
    if (chID<FT0A_CH_N) {
        return xFT0A[chID];
    } else {
        return xFT0A[chID-FT0A_CH_N-1];
    }
}

//______________________________________________________________________
double getFT0ChannelCenterY(int chID) {
    if (chID<FT0A_CH_N) {
        return yFT0A[chID];
    } else {
        return yFT0A[chID-FT0A_CH_N-1];
    }
}

//______________________________________________________________________
double getFV0Phi(int chID) {
    return TMath::Pi()/8.0 + (chID%8)*TMath::Pi()/4.0;
}
