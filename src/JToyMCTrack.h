#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TLorentzVector.h"

class JToyMCTrack : public TObject {

public:
    JToyMCTrack();
    JToyMCTrack(TLorentzVector lVec_in);

    virtual ~JToyMCTrack(){;}

    double GetPx(){return lVec.Px();}
	double GetPy(){return lVec.Py();}
	double GetPz(){return lVec.Pz();}
	double GetPt(){return lVec.Perp();}
	double GetPhi(){return lVec.Phi();}
	double GetEta(){return lVec.Eta();}
	double GetMass(){return lVec.Mag();}

    TLorentzVector GetLVector(){return lVec;}

    void SetPx(double px_in){lVec.SetPx(px_in);}
	void SetPy(double py_in){lVec.SetPy(py_in);}
	void SetPz(double pz_in){lVec.SetPz(pz_in);}
	void SetPt(double pt_in){lVec.SetPerp(pt_in);}
	void SetPhi(double phi_in){lVec.SetPhi(phi_in);}
	void SetMass(double m_in);

	void SetLVector(TLorentzVector lVec_in){lVec = lVec_in;}

    void SetTrack(TLorentzVector lVec_in);

private:
    TLorentzVector lVec;

protected:
    ClassDef(JToyMCTrack, 1)

};
