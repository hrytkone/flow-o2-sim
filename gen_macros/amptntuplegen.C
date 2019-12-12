/**
 * AMPT generator from macro
 *
 * "Pythia6Generator.cxx" and "extgen.C" from ALICE O2 framework used as
 * template
 */

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "FairGenerator.h"
#include "FairPrimaryGenerator.h"
#include <iostream>
#include <cstdio>
#endif

class AMPTGenerator : public FairGenerator
{
    public:

        /**
         * Standard constructor
         * @param inputFile The input file name
         */
        AMPTGenerator(const char* inputFile)
            : FairGenerator("AMPTGenerator"), mFileName(inputFile), mInputFile(nullptr)
        {
            std::cout << "Opening input file " << mFileName << std::endl;
            mInputFile = TFile::Open(mFileName);
            if (!mInputFile->IsOpen()) {
                std::cout << "cannot open the input file" << std::endl;
            }
            events = (TNtuple*)mInputFile->Get("amptEvents");

            iStartNewEvent = 0;
            genEventId = 1;
        };

        /**
         * Destructor
         */
        ~AMPTGenerator() override
        {
            CloseInput();
        };

        /**
         * Reads event from the input file and adds the track onto the stack.
         * @param  primGen Pointer to the PrimaryGenerator
         * @return         true if events read, false if no events any more
         */
        Bool_t ReadEvent(FairPrimaryGenerator* primGen) override
        {

            if (!mInputFile->IsOpen()) {
                std::cout << "Input file not opened!" << std::endl;
                return kFALSE;
            }

            // Event variables
            Float_t eventid = 0;

            // Track variables
            Float_t particleid = 0;
            Float_t px = 0.0, py = 0.0, pz = 0.0, x = 0.0, y = 0.0, z = 0.0;

            events->SetBranchAddress("particleId",&particleid);
            events->SetBranchAddress("eventId",&eventid);
            events->SetBranchAddress("px",&px);
            events->SetBranchAddress("py",&py);
            events->SetBranchAddress("pz",&pz);
            events->SetBranchAddress("x",&x);
            events->SetBranchAddress("y",&y);
            events->SetBranchAddress("z",&z);

            Int_t nentries = (Int_t)events->GetEntries();
            for ( Int_t i=iStartNewEvent; i<nentries; i++ ) {
                Int_t entry = events->GetEntry(i);
                //if (i==iStartNewEvent) cout << "entry : " << entry << endl;
                //if (events->GetEntry(i)!=0) {
                //events->GetEntry(i);
                if (entry>0) {
                    if ((Int_t)eventid==genEventId) {
                        primGen->AddTrack(particleid, px, py, pz, x, y, z);
                    } else {
                        genEventId++;
                        iStartNewEvent = i;
                        break;
                    }
                } else if (entry==0) {
                    std::cout << "No entry found" << std::endl;
                    CloseInput();
                    return kFALSE;
                } else {
                    std::cout << "I/O error occured" << std::endl;
                    CloseInput();
                    return kFALSE;
                }
            }

            return kTRUE;
        };

    private:
        const Char_t* mFileName;
        TFile* mInputFile;
        TNtuple* events;
        Int_t genEventId;
        Int_t iStartNewEvent;

        /**
         * Close input file after reading
         */
        void CloseInput()
        {
            if (mInputFile->IsOpen()) {
                std::cout << "Closing input file " << mFileName << std::endl;
                mInputFile->Close();
            }
        };
};

FairGenerator*
    amptntuplegen(const char* inputFile = "ampt.root")
{
    std::cout << "AMPT generator" << std::endl;
    auto gen = new AMPTGenerator(inputFile);
    return gen;
}
