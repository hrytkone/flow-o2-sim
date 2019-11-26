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
            if ((mInputFile = fopen(mFileName, "r")) == nullptr) {
                std::cout << "cannot open the input file" << std::endl;
            }
        };

        /**
         * Destructor
         */
        ~AMPTGenerator() override
        {
            //CloseInput();
        };

        /**
         * Reads event from the input file and adds the track onto the stack.
         * @param  primGen Pointer to the PrimaryGenerator
         * @return         true if events read, false if no events any more
         */
        Bool_t ReadEvent(FairPrimaryGenerator* primGen) override
        {
            if (!mInputFile) {
                std::cout << "Input file not opened!" << std::endl;
                return kFALSE;
            }

            // Event variables
            Int_t ncols = 0, eventid = 0, test = 0, ntracks = 0, npart1 = 0, npart2 = 0,
                    npart1_el = 0, npart1_inel = 0, npart2_el = 0, npart2_inel = 0;
            Float_t b = 0.0, time = 0.0;

            // Track variables
            Int_t particleid = 0;
            Float_t px = 0.0, py = 0.0, pz = 0.0, m = 0.0, x = 0.0, y = 0.0, z = 0.0, t = 0.0;


            ncols = fscanf(mInputFile, "%d %d %d %f %d %d %d %d %d %d %f", &eventid, &test,
                            &ntracks, &b, &npart1, &npart2, &npart1_el, &npart2_el,
                            &npart1_inel, &npart2_inel, &time);

            if (ncols && ntracks > 0) {
                for (Int_t i = 0; i < ntracks; i++) {
                    ncols = fscanf(mInputFile, "%d %f %f %f %f %f %f %f %f",
                                    &particleid, &px, &py, &pz, &m, &x, &y, &z, &t);
                    primGen->AddTrack(particleid, px, py, pz, x, y, z);
                }
            } else {
                std::cout << "End of input file reached" << std::endl;
                //CloseInput();
                return kFALSE;
            }

            if (feof(mInputFile)) {
                std::cout << "End of input file reached" << std::endl;
                //CloseInput();
                return kFALSE;
            }

            return kTRUE;
        };

    private:
        const Char_t* mFileName;
        FILE* mInputFile;

        /**
         * Close input file after reading
         */
        void CloseInput()
        {
            if (mInputFile) {
                std::cout << "Closing input file " << mFileName << std::endl;
                fclose(mInputFile);
            }
            //delete mInputFile;
            //mInputFile = nullptr;
        };
};

FairGenerator*
    amptgen(const char* inputFile = "ampt.dat")
{
    std::cout << "AMPT generator" << std::endl;
    auto gen = new AMPTGenerator(inputFile);
    return gen;
}
