/**
 * ROOT macro for saving ampt.dat to NTuple
 *
 * Needs ampt.dat file
 */

// Function from pythia code
bool isHadron(int particleid) {
    if (particleid <= 100 || (particleid >= 1000000 && particleid <= 9000000) || particleid >= 9900000)
        return false;
    if (particleid == 130 || particleid == 310)
        return true;
    if (particleid%10 == 0 || (particleid/10)%10 == 0 || (particleid/100)%10 == 0)
        return false;
    return true;
}

void AmptToNTuple(int jobNumber, int eventStartId) {

    TString inFileName(Form("ampt%01i.dat", jobNumber));
    FILE *fIn = fopen(inFileName, "r");

    TString outFileName(Form("ampt-output%02i.root", jobNumber));
    TFile *fOut = new TFile(outFileName, "RECREATE");

    //TNtuple *ntuple = new TNtuple("amptEvents", "data from ampt.dat", "eventId:particleId:px:py:pz:x:y:z:isHadron:charge"); // NOT WORKING
    TNtuple *ntuple = new TNtuple("amptEvents", "data from ampt.dat", "eventId:particleId:px:py:pz:x:y:z:isHadron");

    // Event variables
    Int_t ncols = 0, eventid = 0, test = 0, ntracks = 0, npart1 = 0, npart2 = 0,
           npart1_el = 0, npart1_inel = 0, npart2_el = 0, npart2_inel = 0;
    Float_t b = 0.0, time = 0.0;

    // Track variables
    Int_t particleid = 0;
    Float_t px = 0.0, py = 0.0, pz = 0.0, m = 0.0, x = 0.0, y = 0.0, z = 0.0, t = 0.0;
    Double_t charge = 0.0;
    bool ishadron = false;

    while (1) {
        ncols = fscanf(fIn, "%d %d %d %f %d %d %d %d %d %d %f", &eventid, &test,
                        &ntracks, &b, &npart1, &npart2, &npart1_el, &npart2_el,
                        &npart1_inel, &npart2_inel, &time);

        if (ncols<0) break;

        for (Int_t i = 0; i < ntracks; i++) {
            ncols = fscanf(fIn, "%d %f %f %f %f %f %f %f %f",
                            &particleid, &px, &py, &pz, &m, &x, &y, &z, &t);
            //charge = TDatabasePDG::Instance()->GetParticle(particleid)->Charge(); // NOT WORKING
            ishadron = isHadron(particleid);
            //cout << "charge : " << charge << "     ishadron : " << ishadron << endl;
            //ntuple->Fill(eventStartId + eventid, particleid, px, py, pz, x, y, z, ishadron, charge); // NOT WORKING
            ntuple->Fill(eventStartId + eventid, particleid, px, py, pz, x, y, z, ishadron);
        }
    }

   fclose(fIn);

   fOut->Write();
   fOut->Close();
}
