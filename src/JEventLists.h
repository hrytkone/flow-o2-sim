#ifndef JEVENTLISTS_H
#define JEVENTLISTS_H

#include "TClonesArray.h"

class JEventLists {

public:
    JEventLists();
    virtual ~JEventLists() {;}

    void ClearLists();
    TClonesArray *GetList(int det_i, TString sAorB="");

    TClonesArray *fullEvent;
    TClonesArray *TPClist;
    TClonesArray *T0PAlist;
    TClonesArray *T0PClist;
    TClonesArray *V0Plist;

    TClonesArray *TPClistA;
    TClonesArray *TPClistB;
    TClonesArray *T0PAlistA;
    TClonesArray *T0PAlistB;
    TClonesArray *T0PClistA;
    TClonesArray *T0PClistB;
    TClonesArray *V0PlistA;
    TClonesArray *V0PlistB;
};

#endif
