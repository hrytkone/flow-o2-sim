#include <TSystem.h>

void FinalMacro() {
    gROOT->ProcessLine(".I src/ResIter.h");
    cout << "Run MakeCentralityGraphs.C\n";
    gROOT->ProcessLine(".x MakeCentralityGraphs.C");
    //cout << "Run MakeGraphs.C\n";
    //gROOT->ProcessLine(".x MakeGraphs.C");
    cout << "Run PlotData.C\n";
    gROOT->ProcessLine(".x PlotData.C");
    //cout << "Run PlotVn.C\n";
    //gROOT->ProcessLine(".x PlotVn.C");
}
