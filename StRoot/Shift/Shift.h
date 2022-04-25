#ifndef EventPlane_h
#define EventPlane_h

#include "StMaker.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/EpProManager/EpProManager.h" // Load STARLibrary header files
#include "TString.h"
#include "TFile.h"
// Configuration file reader
#include "../ConfigReader/ConfigReader.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class CutManager;
class HistManager;
class StRefMultCorr;
class StPicoBTofPidTraits;
class EpProManager;

// run QA
//class TProfile;
//class TH1F;
//class TH2F;

class Shift : public StMaker {
    public:
        Shift(const char *name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName);
        virtual ~Shift();

        virtual Int_t Init();
        virtual Int_t Make();
        virtual void  Clear(Option_t *opt="");
        virtual Int_t Finish();

        //int Centrality(int gRefMult);
        int GetRunIndex(int runId);

        //bool isGoodTrigger(StPicoEvent const*) const;

    private:
	ConfigReader configs;
        static StRefMultCorr *mRefMultCorr;
        StPicoDstMaker *mPicoDstMaker;
        StPicoDst      *mPicoDst;
        StPicoEvent *mPicoEvent;
        CutManager  *mCutManager;
        //HistManager *mHistManager;

	// EPD 
        StEpdGeom *mEpdGeom;
        StEpdEpInfo *mEpdEpInfo;
        EpProManager *mEpProManager;

	Int_t mEnergy;

	// get the recenterpar
        TProfile2D *getep_sub_recen[2][_numSubEvents];
        //  get shift par
        TProfile3D *p_sub_ep_shiftpar_sin[_numSubEvents];
        TProfile3D *p_sub_ep_shiftpar_cos[_numSubEvents];

        TString mOutPut_Shift;

        TFile *mFile_Shift;


        ClassDef(Shift, 1)
};

#endif
