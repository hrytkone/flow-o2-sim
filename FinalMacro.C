#include <TSystem.h>

void FinalMacro(TString sIn="toyFlow.root", TString sOut="toyFlowGraphs.root", int readNFiles=1) {
    gROOT->ProcessLine(".I src/ResIter.h");

    if (gSystem->AccessPathName(sIn)) {
        cout << "File \"" << sIn << "\" not found!\n";
        return 0;
    }

    cout << "Run MakeGraphs.C: " << Form(".x MakeGraphs.C(\"%s\",\"%s\")",sIn.Data(),sOut.Data()) << "\n";
    gROOT->ProcessLine(Form(".x MakeGraphs.C(\"%s\",\"%s\")",sIn.Data(),sOut.Data()));
    cout << "Run PlotVn.C\n";
    gROOT->ProcessLine(".x PlotVn.C");

}
