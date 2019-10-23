#include "JInputs.h"

JInputs::JInputs() {;}

JInputs::~JInputs() {
    for (int i=0; i<CENTBINS_N-1; i++) {
        if(hEtaDist[i]!=0x0) delete hEtaDist[i];
    }
}

// Initializes the eta distribution histogram
// and the total multiplicities for each centrality
// bin.
void JInputs::Load() {
    int i, j;
    for (i=0; i<CENTBINS_N-1; i++) {
        //Initialize histo:
        hEtaDist[i] = new TH1F(Form("hEtaDist%d",i),Form("hEtaDist%d",i),ETADST_N-1,etadst);
        for (j=0; j<ETADST_N; j++) {
            hEtaDist[i]->SetBinContent(j,etanch[i][j]);
        }
        dMulti[i] = hEtaDist[i]->Integral("width");
        cout << "CentBin: " << i << ", Multi: " << dMulti[i] << endl;
    }
}

bool JInputs::CheckCentBin(int centBin) {
    return centBin<CENTBINS_N-1 && centBin>-1;
}

int JInputs::GetMultiplicity(int centrality) {
    double multi;
    if(CheckCentBin(centrality)) {
        multi = dMulti[centrality];
    } else {
        multi = 0;
    }
    return (int)TMath::Floor(multi);
}

double JInputs::GetEta(int centrality) {
    double eta;
    if(CheckCentBin(centrality)) {
        eta = hEtaDist[centrality]->GetRandom();
    } else {
        eta = -999; //ok?
    }
    return eta;
}

int JInputs::GetCentBin(double centrality) {
    for (int i=0; i<CENTBINS_N-1; i++)
        if (centrality>centBins[i] && centrality<centBins[i+1]) return i;
    return -1;
}

double JInputs::GetCentDependVn(int n, double centrality) {
    int centBin = GetCentBin(centrality);
    return centvn[n-1][centBin];
}
