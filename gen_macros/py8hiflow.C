/**
 * Pythia8 generator with flow afterburner from macro
 *
 * More about flow afterburner:
 *      DOI: 10.1103/PhysRevC.79.064909
 */

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include "SimConfig/SimConfig.h"
#endif

double GetPhi0(double phi, double *vn, double *psi);
double AnisotropicPhiDist(double *x, double *p);
double GetAnisotropicPhi(double phi0, double phiInit, double err, double *vn, double *psi, TF1 *fPhiDist);

class Pythia8FlowGenerator : public FairGenerator
{

public:

    Pythia8FlowGenerator(TPythia8* pythia8)
        : FairGenerator("Pythia8FlowGenerator"), mPythia(pythia8)
    {
        particles = new TClonesArray("TParticle", 100000);
        fPhiDist = new TF1("fPhiDist", AnisotropicPhiDist, -TMath::Pi(), TMath::Pi(), 11);
        rand = new TRandom3(0);
    };

    ~Pythia8FlowGenerator() override
    {
        //particles->Clear("C");
    };

    Bool_t ReadEvent(FairPrimaryGenerator* primGen) override
    {
        const double mm2cm = 0.1;

        // Set random flow angles
        for (int i=0; i<5; i++) psi[i] = rand->Uniform(-TMath::Pi(), TMath::Pi())/(i+1);

        mPythia->GenerateEvent();
        mPythia->ImportParticles(particles, "All");
        //mPythia->EventListing();

        Int_t nParticles = particles->GetEntriesFast();
        for (Int_t ip = 0; ip < nParticles; ip++) {
            TParticle* part = (TParticle*) particles->At(ip);

            Int_t status = part->GetStatusCode();
            Int_t pdg = part->GetPdgCode();

            if (pdg == 22) {
                cout << "Pb:\n";
                cout << "   px=" << part->Px() << endl;
                cout << "   py=" << part->Py() << endl;
                cout << "   pz=" << part->Pz() << endl;
                cout << "   phi=" << part->Phi() << endl;
                cout << "   phi x=" << part->PhiX() << endl;
                cout << "   phi y=" << part->PhiY() << endl;
                cout << "   phi z=" << part->PhiZ() << endl;
                cout << "   polar phi=" << part->GetPolarPhi() << endl;
                cout << "   polar theta=" << part->GetPolarTheta() << endl;

            }

            if (status <= 0) continue;

            //Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
            //if (charge == 0.) continue;

            double phi0 = part->Phi();
            double phi = GetAnisotropicPhi(phi0, 1., 0.001, vn, psi, fPhiDist);
            double phiDiff = phi - phi0;

            part->ProductionVertex(vec);
            vec.RotateZ(phiDiff);

            Double_t x = vec.X();
            Double_t y = vec.Y();
            Double_t z = vec.Z();
            Double_t px = vec.Px();
            Double_t py = vec.Py();
            Double_t pz = vec.Pz();

            x *= mm2cm;
            y *= mm2cm;
            z *= mm2cm;

            primGen->AddTrack(pdg, px, py, pz, x, y, z);
        }

        particles->Clear("C");
        return kTRUE;
    };

private:
    TPythia8* mPythia;
    TClonesArray* particles;
    TRandom3 *rand;
    TF1 *fPhiDist;

    TLorentzVector vec;
    double vn[5] = {0., 0.15, 0.08, 0.03, 0.01};
    double psi[5] = {TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi(), TMath::Pi()};
};

//______________________________________________________________________________
// CLASS ENDS HERE
//______________________________________________________________________________

FairGenerator*
    py8hiflow()
{
    auto& conf = o2::conf::SimConfig::Instance();

    TPythia8* pythia8 = new TPythia8();
    pythia8->ReadString("Beams:idA 1000822080");                                      // Pb ion
    pythia8->ReadString("Beams:idB 1000822080");                                      // Pb ion
    pythia8->ReadString("Beams:eCM 5520.0");                                          // [GeV]
    pythia8->ReadString("HeavyIon:SigFitNGen 0");                                     // valid for Pb-Pb 5520 only
    pythia8->ReadString("HeavyIon:SigFitDefPar 14.82,1.82,0.25,0.0,0.0,0.0,0.0,0.0"); // valid for Pb-Pb 5520 only
    //pythia8->ReadString(("HeavyIon:bWidth " +  std::to_string(conf.getBMax())).c_str());
    pythia8->ReadString("HeavyIon:bWidth 10");
    pythia8->ReadString("ParticleDecays:tau0Max 0.001");
    pythia8->ReadString("ParticleDecays:limitTau0 on");
    pythia8->Initialize(1000822080, 1000822080, 5520.0);

    auto gen = new Pythia8FlowGenerator(pythia8);

    return gen;
}

//______________________________________________________________________________
// Functions to calculate new phi
//______________________________________________________________________________

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

double GetAnisotropicPhi(double phi0, double phiInit, double err, double *vn, double *psi, TF1 *fPhiDist) {

    double phi = 0;
    fPhiDist->SetParameters(phi0, vn[0], vn[1], vn[2], vn[3], vn[4], psi[0], psi[1], psi[2], psi[3], psi[4]);

    while (TMath::Abs(GetPhi0(phi, vn, psi) - phi0) > err) {
        phi = phiInit - fPhiDist->Eval(phiInit)/fPhiDist->Derivative(phiInit);
        phiInit = phi;
    }

    if (phi>TMath::Pi())phi -= 2*TMath::Pi();
    if (phi<-TMath::Pi()) phi += 2*TMath::Pi();

    return phi;
}
