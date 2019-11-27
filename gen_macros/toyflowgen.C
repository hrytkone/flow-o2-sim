/**
 * Generator based on the older ToyFlow code (includes only flow)
 */

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#endif

#include "JInputs.cpp"

double PtDist(double *x, double *p) {
    return TMath::Exp(-p[0]*x[0]);
};

double PhiDist(double *x, double *p) {
    double phi = x[0];
    double vn[5] = {p[0], p[1], p[2], p[3], p[4]};
    double psi[5] = {p[5], p[6], p[7], p[8], p[9]};
    return 1.0 + 2.0*vn[0]*TMath::Cos(phi - psi[0]) + 2.0*vn[1]*TMath::Cos(2.*(phi - psi[1]))
        + 2.0*vn[2]*TMath::Cos(3.*(phi - psi[2])) + 2.0*vn[3]*TMath::Cos(4.*(phi - psi[3]))
        + 2.0*vn[4]*TMath::Cos(5.*(phi - psi[4]));
};


class ToyFlowGenerator : public FairGenerator
{

public:

    ToyFlowGenerator()
        : FairGenerator("ToyFlowGenerator")
    {
        rand = new TRandom3(0);
        fPhiDist = new TF1("fPhiDist", PhiDist, -TMath::Pi(), TMath::Pi(), 10);
        fPtDist = new TF1("fPtDist", PtDist, 0.0, 10.0, 1);

        double Tdec = 0.12;
        double vr = 0.6;
        double Teff = Tdec * TMath::Sqrt((1.+vr)/(1.-vr));
        fPtDist->SetParameter(0, 1./Teff);

        inputs = new JInputs();
        inputs->Load();
    };

    ~ToyFlowGenerator() override { };

    Bool_t ReadEvent(FairPrimaryGenerator* primGen) override
    {
        double centrality = rand->Uniform(0.0, 60.0);
        int centBin = inputs->GetCentBin(centrality);
        int nParticles = inputs->GetMultiplicity(centBin);

        // Set random flow angles
        for ( int i=0; i<5; i++ )
            psi[i] = rand->Uniform(-TMath::Pi(), TMath::Pi())/(i+1);

        fPhiDist->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4], psi[0], psi[1], psi[2], psi[3], psi[4]);

        for ( int i = 0; i < nParticles; i++ ) {

            double phi = fPhiDist->GetRandom();
            double pT = fPtDist->GetRandom();
            double eta = rand->Uniform(-5.0, 5.0);
            double px = pT*TMath::Cos(phi);
            double py = pT*TMath::Sin(phi);
            double pz = pT*TMath::SinH(eta);

            primGen->AddTrack(211, px, py, pz, 0.0, 0.0, 0.0); // 211 pdg for pion
        }

        return kTRUE;
    };

private:
    TRandom3 *rand;
    TF1 *fPhiDist;
    TF1 *fPtDist;
    JInputs *inputs;

    double vn[5] = {0., 0.15, 0.08, 0.03, 0.01};
    double psi[5] = {TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi()};
};


//______________________________________________________________________________
// CLASSES END HERE
//______________________________________________________________________________

FairGenerator*
    toyflowgen()
{
    auto gen = new ToyFlowGenerator();
    return gen;
}
