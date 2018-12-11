#ifndef AliAnalysisTask_eta_cxx
#define AliAnalysisTask_eta_cxx

class TH1F;
class AliESDEvent;
class AliAODEvent;
class AliAODTrack;
class AliESDtrackCuts;
class TList;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include <vector>
#include "AliAODTrack.h"

class AliAnalysisTask_eta : public AliAnalysisTaskSE {
public:
    AliAnalysisTask_eta();
    AliAnalysisTask_eta(const char *name);
    virtual ~AliAnalysisTask_eta();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    void SetTrackCuts(AliESDtrackCuts* const cuts) { fTrackCuts = cuts; }
    void StoreGlobalTrackReference();
    Double_t Psi_pair(AliAODTrack *neg, AliAODTrack *pos);
    void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
    void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
    void SetMinEEMC(Double_t MinEemcal) { fMinEEMC = MinEemcal; };
    void SetMinEPHS(Double_t MinEphos) { fMinEPHS = MinEphos; };
    void SetMaxE(Double_t MaxEnergy) { fMaxE = MaxEnergy; };
    void SetTrigger(TString Trigger) { fTrigger = Trigger; };
    
    Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };

    
private:
    
    static const int zvtx_bins = 8;
    static const int mult_bins = 7;
    static const unsigned int poolDepth = 80;
    Int_t GetMultBin(Int_t mult);
    Int_t GetZvtxBin(Double_t vertZ);
    
    enum{
        kAODanalysis = BIT(20),
    };
    
    AliVEvent		*fVevent;  		    //! event object
    AliESDEvent		*fESD;    		    //! ESD object
    AliAODEvent		*fAOD;    		    //! AOD object

    TClonesArray	*fTracks;
    TClonesArray	*fCaloClusters;

    Double_t		fMinEEMC; 
    Double_t		fMinEPHS; 
    Double_t		fMaxE;    
    TString		fTrigger;
    
    TList       	*fOutputList; 		    //! Output list
    TH1F		*fHistV0E;		    //! V0's energy distribution
    TH1F        	*fNevents;		    //! no of events
    TH1F        	*fClustStat;		    //! no of clusters
    TH1F		*fV0Stat;		    //! no of V0s
    TH1F        	*fVtxZ;			    //! Vertex z
    TH1F        	*fShapeParam;		    //! Shape Parameter
    TH1F        	*fShapeParam2;		    //! After Cuts  
    TH1F        	*fHistClustE;		    //! cluster energy
    TH1F        	*fHistClustE2;		    //! after cuts
    TH2F        	*fEMCClsEtaPhi;		    //! EMC cluster eta and phi
//    TH1F        	*fHistoNCells;		    //! No of cells per cluster
//    TH1F        	*fHistoNCells2;		    //! after cuts
    TH2F                *ftof;                      //! TOF vs cluster energy
    TH2F                *ftof2;                     //! TOF vs cluster energy after cuts
    	
    TH1I        	*fHistoNtracksMatch;	    //! No of tracks matched to cluster by correction task
    TH2F        	*fHistoTrackMatchedPHOS;    //! delta eta vs delta phi for charged tracks matching in PHOS
    TH2F        	*fHistoTrackMatchedPHOS2;   //! after cuts
    TH2F        	*fHistoTrackMatchedEMC;     //! delta eta vs delta phi for charged tracks matching in EMCal
    TH2F        	*fHistoTrackMatchedEMC2;    //! after cuts

    TH1F		*fHisto_M_V0;		    //! V0 mass after opening angle cut
    TH2F		*fHisto_M_pt_V0;	    //! V0 mass vs pT after opening angle cut
    TH1F		*fHistV0InvMassPi0;	    //! V0 mass after opening angle and pi0 mass cut
    TH2F        	*fHistV0InvMassPtPi0;	    //! V0 mass vs pT after opening angle and pi0 mass cut
    TH2F        	*fHisto_M_pt_EMC;	    //! only EMCal
    TH1F        	*fHisto_M_EMC;		    //!
    TH2F        	*fHisto_M_pt_PHS;	    //! only PHOS
    TH1F        	*fHisto_M_PHS;		    //!    
    TH2F        	*fHisto_M_pt_Pi0;	    //! Calo mass vs. pt with Pi0
    TH1F        	*fHisto_M_Pi0;		    //! Calo mass with Pi0
    TH2F        	*fHisto_M_pt_Eta;	    //! Calo mass vs. pt without Pi0
    TH1F        	*fHisto_M_Eta;		    //! Calo mass without Pi0
    TH2F        	*fHisto_M_pt_All;	    //! mass vs. pt with all pairings allowed
    TH1F        	*fHisto_M_All;		    //! mass with all pairings allowed
    TH1F 		*ftest;			    //!
//    TH2F        	*fHistoE_NCells;	    //!

    AliPIDResponse 	*fPIDResponse;   	    //! PID response object
    


    // persistent members are streamed (copied/stored)
    AliESDtrackCuts *fTrackCuts; // Track cuts
    std::vector<AliAODTrack*>     fGlobalTrackReference;    
    
    std::vector<TLorentzVector> Photons[poolDepth][zvtx_bins][mult_bins]; //!
    
    AliAnalysisTask_eta(const AliAnalysisTask_eta&); // not implemented
    AliAnalysisTask_eta& operator=(const AliAnalysisTask_eta&); // not implemented
    
    ClassDef(AliAnalysisTask_eta, 1); // example of analysis
};

#endif


