#include "Shift.h"
#include "StRoot/CutManager/CutManager.h"
#include "StRoot/ConstManager/ConstManager.h"
//#include "../HistManager/HistManager.h"
#include "StRoot/HistManager/HistManager.h" // Load STARLibrary header files
#include "StRoot/EpProManager/EpProManager.h" // Load STARLibrary header files
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/Run/run.h"
#include "StThreeVectorF.hh"
//IClass header files
#include "StRoot/IClasses/IEventPlane.h"
#include "StRoot/IClasses/IEvent.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>
//#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//#include "StRoot/StRefMultCorr/CentralityMaker.h"
ClassImp(Shift)

    //StRefMultCorr* Shift::mRefMultCorr = NULL;
    //-----------------------------------------------------------------------------
    Shift::Shift(const char* name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName)
: StMaker(name)
{
    configs.read(configFileName);
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mEnergy = 3;

    mOutPut_Shift=Form("%s_shift.root",jobid);
}

//----------------------------------------------------------------------------- 
Shift::~Shift()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t Shift::Init() 
{
    //if(!mRefMultCorr)
    //{
    //    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    //}
    // qapid
    mCutManager = new CutManager(configs);
    /*mHistManager = new HistManager();
    mHistManager->InitQAPID();*/

    // eventplane 
    // EPD
    mEpdGeom = new StEpdGeom();
    mEpProManager = new EpProManager();
    mEpProManager->InitEP();
    //const int numSubEvents = 3;
    //char subEventModes[_numSubEvents] = {'r','e','e'};

    mFile_Shift = new TFile(mOutPut_Shift.Data(),"RECREATE");
    mFile_Shift->cd();
    
    TFile *mRecenteringInputFile = TFile::Open("/star/data01/pwg/dchen/Ana/3GeV_FXT_2018/evenplane/3gev_2018_recenterpar.root");
    if (mRecenteringInputFile->IsZombie()) {
    	std::cout << "Error opening file with Q vector recentering histograms" << std::endl;
    	std::cout << "No recentering in this analysis." << std::endl;
        for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
        {
            getep_sub_recen[0][sub]     = 0;
	    getep_sub_recen[1][sub]     = 0;
        } // event plane Psi histograms
    }
    else{
      for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
      {
          getep_sub_recen[0][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq1x_Sub_EP_%d",sub + 1));
          getep_sub_recen[1][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq1y_Sub_EP_%d",sub + 1));
      } // event plane Psi histograms
      std::cout << "Recenterpar successfully implemented" << std::endl;
    }
    
    TFile *mShiftInputFile = TFile::Open("/star/u/dchen/ana/3gev_2018/recenter/test_recenter.root");
    if (mShiftInputFile->IsZombie()) {
    	std::cout << "Error opening file with shift correction histograms" << std::endl;
    	std::cout << "No shifting in this analysis." << std::endl;
        for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
	{
          p_sub_ep_shiftpar_sin[sub] = 0;
          p_sub_ep_shiftpar_cos[sub] = 0;
	}
    }
    else{
      for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
      {
        p_sub_ep_shiftpar_sin[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_Sub_%d_sin", sub + 1));
        p_sub_ep_shiftpar_cos[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_Sub_%d_cos", sub + 1));
      }
      std::cout << "Shiftpar successfully implemented" << std::endl;
    }
    return 0;

}

//----------------------------------------------------------------------------- 
Int_t Shift::Finish() 
{


    if(mOutPut_Shift != "")
    {
        mFile_Shift->cd();
        //mHistManager->WriteQAPID();
        mEpProManager->WriteEP();
        mFile_Shift->Close();
    }
    return kStOK;
}

//----------------------------------------------------------------------------- 
void Shift::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t Shift::Make() 
{
    if(!mPicoDstMaker) 
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) 
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
    const Int_t nTracks = mPicoDst->numberOfTracks();
    Int_t TrkMult = 0;
    for(Int_t i = 0; i < nTracks; i++) // TrkMult
    {
        StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
        if(!track->isPrimary()) continue; // Only Primary Tracks
        //mHistManager->FillTrackQA(track,(TVector3)mPicoEvent->primaryVertex());
        //StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(mVertexPos);
	TrkMult ++;
    } // TrkMult

    // RefMult
    Int_t runId = mPicoEvent->runId();

    //cout << "runID = " << runId << endl;
    /*Int_t refMult = mPicoEvent->refMult();
    Float_t vz = mPicoEvent->primaryVertex().Z();
    Float_t vx = mPicoEvent->primaryVertex().X();
    Float_t vy = mPicoEvent->primaryVertex().Y();

    Float_t vzvpd = mPicoEvent->vzVpd();
    Int_t TOF_Mul = mPicoEvent->btofTrayMultiplicity();
    Int_t nMatchedToF = mPicoEvent->nBTOFMatch();
    Float_t zdcX = mPicoEvent->ZDCx();*/
    
    /*mHistManager->FillEventQA(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
    mHistManager->FillEventCut(0);*/

    // runIndex
    const int runIndex = GetRunIndex(runId);
    // event plane IClasses
    //Start with clean events
    IEvent * theEvent = new IEvent;   
    IEvent subEvents[_numSubEvents];
    theEvent->ClearEvent();
    for (int i = 0; i < _numSubEvents; i++)
    {
    	subEvents[i].ClearEvent();
    }	

    // Event Cut
    if(mCutManager->passEventCut(mPicoDst)) // event cut
    {
	
        /*mHistManager->FillEventQaCut(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
        mHistManager->FillEventCut(1);*/
        //mRefMultCorr->init(runId);
        //mRefMultCorr->initEvent(refMult, vz, zdcX);
        //const Int_t cent9 = mRefMultCorr->getCentralityBin9();
        //const Double_t reweight = mRefMultCorr->getWeight();
	//std::cout << "refMult: " << refMult << std::endl;
	//std::cout << "TrkMult: " << TrkMult << std::endl;
        const int cent16 = mCutManager->getCentrality(TrkMult);
	
	//std::cout << "cent16: " << cent16 << std::endl;
        //const double reweight = 1.0;
        if(cent16 >  15 || cent16 < 0) return 0;
        //mHistManager->FillEventCent(cent16);
        //mHistManager->FillEventCut(2);
	
	/// qapid
        //const Int_t nToFMatched = mCutManager->getMatchedToF();
        TVector3 mVertexPos = mPicoDst->event()->primaryVertex();
        //float mField = mPicoEvent->bField();
        //Int_t N_pp = 0, N_pm = 0, N_kp = 0, N_km = 0, N_pr = 0;

	/// eventplane
	// TPC eventplane
        //double qx_TPC_A=0.0, qy_TPC_A=0.0;
        //double qx_TPC_B=0.0, qy_TPC_B=0.0;

        //cout << "nTracks = " << nTracks << endl;
        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
            StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
            if(!track->isPrimary()) continue; // Only Primary Tracks
            //mHistManager->FillTrackCut(0);
            if(!mCutManager->passTrackBasic(track)) continue;
            //mHistManager->FillTrackCut(1);
            //mHistManager->FillTrackPhysics(track );
		
            //StPicoPhysicalHelix helix = track->helix(mField);
            //Float_t dca = helix.geometricSignedDistance(mVertexPos);
            //Short_t  s_charge = track->charge();
            Float_t dca=track->gDCA(mVertexPos).Mag();
	    /*if(mCutManager->isTofTrack(mPicoDst,track)) mHistManager->FillTrackTof(mPicoDst,track);
	    if(mCutManager->isProton(track))
	    {
	   	mHistManager->FillProton(mPicoDst,track,configs.y_mid); 
		N_pr++;
	    }
	    if(mCutManager->isKaon(mPicoDst,track))
	    {
	   	mHistManager->FillKaon(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_kp++;
		else N_km++;
	    }
	    if(mCutManager->isPion(mPicoDst,track))
	    {
	   	mHistManager->FillPion(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_pp++;
		else N_pm++;
	    }*/
	    //TPC EP
            if(!mCutManager->passTrackEP(track,dca)) continue;
	    Double_t eta = track->pMom().Eta();
	    Double_t phi = track->pMom().Phi();
	    Double_t pt = track->pMom().Perp();
	    IEventPlane eventPlane(phi, pt);
	    eventPlane.SetEta(eta);
	    theEvent->AddEPParticle(eventPlane);
            /*if(eta > -2.0 && eta < -1.1)    // sub event A eta(-2,-1.25)   // mid-rapidity is 1.06  -2.0-1.4, -1.05-0.65
            {
	   	mEpProManager->FillTpcAQvec(cent16,runIndex, track); 
            }

            if(eta > -1.0 && eta < 0.)    // sub event B eta(-1.15,0) // -1.3-0.685,  -0.65-0
            {
	   	mEpProManager->FillTpcBQvec(cent16,runIndex,track); 
            }*/
        } // track loop

	/// eventplane
	// EPD information
        Int_t nepdHits = mPicoDst->numberOfEpdHits();
        StPicoEpdHit *epdHit;
        TVector3 StraightLine_center;
        TVector3 StraightLine_random;
        double phi_epd_center = {0.0};
        double phi_epd_random = {0.0};
        double eta_epd_center = {0.0};
        double eta_epd_random = {0.0};
        double mip;
        double TileWeight           = {0};

        /*double Qx_EPD[7]={0};
        double Qy_EPD[7]={0};
        double weight_EPD[7]={0};*/
        for(Int_t iHit=0; iHit<nepdHits; iHit++)
	{ // EPD loop
            epdHit = mPicoDst->epdHit(iHit);
            mip = epdHit->nMIP();
            int iring = epdHit->row() -1;//(1~16)-1 -> 0-15

            if( !epdHit) continue;
            int position = epdHit->position();
            int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner mose
            //std::cout << "ringgroup = " << ringgroup << std::endl;
            if(ringgroup == -1) continue;
            //std::cout << "epdHit id = " << epdHit->id() << std::endl;
            if(epdHit->id() > 0 ) continue; // only East side for FXT 

            StraightLine_center = mEpdGeom->TileCenter(epdHit->id())        - mPicoEvent->primaryVertex();  // the collision is from midille to side of the TPC
            StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - mPicoEvent->primaryVertex();
            //StraightLine_center = mEpdGeom->TileCenter(epdHit->id())       ;
            //StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id());

            phi_epd_center = StraightLine_center.Phi();
            eta_epd_center = StraightLine_center.Eta();
            phi_epd_random = StraightLine_random.Phi();
            eta_epd_random = StraightLine_random.Eta();
	    mEpProManager->FillEpdQa(iring, phi_epd_center, eta_epd_center, phi_epd_random, eta_epd_random);
            if(mip < configs.epd_threshold) continue;
            TileWeight = (mip > configs.epd_max_weight) ? configs.epd_max_weight : mip;
            //std::cout << "tileweight  = " << TileWeight << std::endl;

	    mEpProManager->FillEpdQa(iring, phi_epd_center, eta_epd_center, phi_epd_random, eta_epd_random);
	    mEpProManager->FillEpdMip(eta_epd_center,position,TileWeight);

	    IEventPlane eventPlane(phi_epd_center, TileWeight);
	    eventPlane.SetTileID(epdHit->id());
	    theEvent->AddEPParticle(eventPlane);
            /*for(int i=0; i<4; i++){
                if(i == ringgroup){
                    Qx_EPD[i] += TileWeight * cos(1.0*phi_epd_center);
                    Qy_EPD[i] += TileWeight * sin(1.0*phi_epd_center);
                    weight_EPD[i] += TileWeight;
                }
            }
            if(ringgroup == 0 || ringgroup == 1)
            {
                Qx_EPD[4] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[4] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[4] += TileWeight;
            }
            if(ringgroup == 2 || ringgroup == 3)
            {
                Qx_EPD[5] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[5] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[5] += TileWeight;
            }
            if(ringgroup == 0 || ringgroup == 1 || ringgroup == 2 || ringgroup == 3)
            {
                Qx_EPD[6] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[6] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[6] += TileWeight;
            }*/
        } // EPD loop
    	/*for(int i=0; i<7; i++){
    	    if(weight_EPD[i] == 0. || Qx_EPD[i]== 0.0 || Qy_EPD[i] == 0.0) return 0;
    	    Qx_EPD[i] /= weight_EPD[i];
    	    Qy_EPD[i] /= weight_EPD[i];
    	}
    	for(int i=0; i<7; i++){
	    mEpProManager->FillEpdQvec(i,cent16,runIndex,Qx_EPD[i],Qy_EPD[i]);
	}*/
	//double psi_ievt = theEvent->GetEventPsi(1);
	//std::cout << "ievent psi = " << psi_ievt << std::endl;
	//mHistManager->FillPIDMult(N_pp , N_pm , N_kp , N_km , N_pr ); /// qapid
	for (int sub = 0; sub < _numSubEvents; sub++)
	{
	  //subEvents[sub] = theEvent->GetSubEvent(subEventModes[sub], subEventParams[sub][0] - COMrapidity, subEventParams[sub][1] - COMrapidity);
	  subEvents[sub] = theEvent->GetSubEvent(_subEventModes[sub], _subEventParams[sub][0], _subEventParams[sub][1]);
	}
	//std::cout << "Size of  EPD-C sub EP = " << subEvents[0].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  TPC-A sub EP = " << subEvents[1].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  TPC-B sub EP = " << subEvents[2].GetEPParticles().size() << std::endl << std::endl;
	if (subEvents[0].GetEPParticles().size() < 2
	|| subEvents[1].GetEPParticles().size() < 2
	|| subEvents[2].GetEPParticles().size() < 2)
		{return 0;}	
	
	float subQx[_numSubEvents];
	float subQy[_numSubEvents];
	for (int sub = 0; sub < _numSubEvents; sub++) // raw event plane, store recenter parameter
	{
	  //std::cout<< subEvents[sub].GetEventPsi(1) << std::endl;
	  mEpProManager->FillPsiRaw(sub, subEvents[sub].GetEventPsi(1));
     	  subQx[sub] = subEvents[sub].GetQx(1);
     	  subQy[sub] = subEvents[sub].GetQy(1);
	  mEpProManager->FillSubEpQvec(sub, cent16, runIndex, subQx[sub], subQy[sub]);
	
	} // raw event plane, store recenter parameter
	
	for (int i = 0; i < _numSubEvents; i++)
	{
	  if (getep_sub_recen[0][i] != 0 && getep_sub_recen[1][i] != 0) // recenter event plane, store shift parameter
	     {
	       float subEventQx = (float)getep_sub_recen[0][i]->GetBinContent(cent16+1, runIndex+1);
	       float subEventQy = (float)getep_sub_recen[1][i]->GetBinContent(cent16+1, runIndex+1);
	       subEvents[i].SetQCenter(subEventQx, subEventQy);
	     }
	} // recenter event plane, store shift parameter

	float subQxRec[_numSubEvents];
	float subQyRec[_numSubEvents];
	float psiShift[_numSubEvents];
	for (int i = 0; i < _numSubEvents; i++)
	{ // recenter event plane, store shift parameter
        	if(subEvents[i].GetEventPsi(1) == -99){continue;} //Empty subevent
	        mEpProManager->FillPsiRec(i, subEvents[i].GetEventPsi(1));
     	  	subQxRec[i] = subEvents[i].GetQx(1);
     	  	subQyRec[i] = subEvents[i].GetQy(1);
		// implement shift correction
		float f_Psi_shifted = subEvents[i].GetEventPsi(1);
	        mEpProManager->FillSubEpQvecRec(i, cent16, runIndex, subQxRec[i], subQyRec[i]);
                for(int iorder=0; iorder<20; iorder++)
                {
                    // output first order shift parameter
	            mEpProManager->FillSubEpShiftpar(i, iorder, cent16, runIndex, subEvents[i].GetEventPsi(1));
                    // input first order shift parameter
		    // implement shift correction
	       	    float tmp = (float)(configs.order_m*(iorder+1));
	       	    float sinAve = p_sub_ep_shiftpar_sin[i]->GetBinContent(cent16+1, iorder+1, runIndex+1);
	       	    float cosAve = p_sub_ep_shiftpar_cos[i]->GetBinContent(cent16+1, iorder+1, runIndex+1);
	       	    f_Psi_shifted +=
	       	    2.0*(cosAve*TMath::Sin(tmp* subEvents[i].GetEventPsi(1)) - sinAve*TMath::Cos(tmp* subEvents[i].GetEventPsi(1)))/tmp;
                }
		double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_m;
		double HalfWrapAround = AngleWrapAround/2;
		while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
		while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	        mEpProManager->FillPsiShift(i, f_Psi_shifted);
		psiShift[i] = f_Psi_shifted;
	} // recenter event plane, store shift parameter
	// Fill resolution plots
	mEpProManager->FillPsiResolution(cent16, psiShift[0], psiShift[1], psiShift[2]);

    } // event cut
    return kStOK;
}

int Shift::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<nrun; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!" << endl;
    return runIndex;
}

/*int Shift::Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={4, 9,17,30,50,78, 116,170,205};
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}*/
