/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//////////////////////////////////////////////
//    Service work task for EMCAL //
//////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliExternalTrackParam.h"
#include "AliEMCALRecoUtils.h"
#include "AliTrackerBase.h"

#include <AliESDtrackCuts.h>
#include <AliVVZERO.h>
#include <AliAODv0.h>
#include <AliAODTrack.h>
#include "AliPIDResponse.h"
#include <TLorentzVector.h>

//SM
#include "AliEMCALGeometry.h"

#include "AliAnalysisTask_eta.h"

ClassImp(AliAnalysisTask_eta)
//________________________________________________________________________
AliAnalysisTask_eta::AliAnalysisTask_eta(const char *name)
: AliAnalysisTaskSE(name),
fVevent(0),
fESD(0),
fAOD(0),
fTracks(0),
fHistV0E(0),
fCaloClusters(0),
fOutputList(0),
fNevents(0),
fVtxZ(0),
fShapeParam(0),
fShapeParam2(0),
fHistClustE(0),
fHistClustE2(0),
fEMCClsEtaPhi(0),
//fHistoNCells(0),
//fHistoNCells2(0),
ftof(0),
ftof2(0),
fHistoNtracksMatch(0),
fHistoTrackMatchedPHOS(0),
fHistoTrackMatchedPHOS2(0),
fHistoTrackMatchedEMC(0),
fHistoTrackMatchedEMC2(0),
fHisto_M_V0(0),
fHisto_M_pt_V0(0),
fHistV0InvMassPi0(0),
fHistV0InvMassPtPi0(0),
fHisto_M_pt_EMC(0),
fHisto_M_EMC(0),
fHisto_M_pt_PHS(0),
fHisto_M_PHS(0),
fHisto_M_pt_Pi0(0),
fHisto_M_Pi0(0),
fHisto_M_pt_Eta(0),
fHisto_M_Eta(0),
fHisto_M_pt_All(0),
fHisto_M_All(0),
//fHistoE_NCells(0),
fClustStat(0),
fV0Stat(0),
ftest(0),
fGlobalTrackReference(),
fPIDResponse(0),
fTrackCuts(0)
{
    // Constructor
    
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTask_eta::AliAnalysisTask_eta()
: AliAnalysisTaskSE("DefaultTask_HFEemcQA2"),
fVevent(0),
fESD(0),
fAOD(0),
fTracks(0),
fHistV0E(0),
fCaloClusters(0),
fOutputList(0),
fNevents(0),
fVtxZ(0),
fShapeParam(0),
fShapeParam2(0),
fHistClustE(0),
fHistClustE2(0),
fEMCClsEtaPhi(0),
//fHistoNCells(0),
//fHistoNCells2(0),
ftof(0),
ftof2(0),
fHistoNtracksMatch(0),
fHistoTrackMatchedPHOS(0),
fHistoTrackMatchedPHOS2(0),
fHistoTrackMatchedEMC(0),
fHistoTrackMatchedEMC2(0),
fHisto_M_V0(0),
fHisto_M_pt_V0(0),
fHistV0InvMassPi0(0),
fHistV0InvMassPtPi0(0),
fHisto_M_pt_EMC(0),
fHisto_M_EMC(0),
fHisto_M_pt_PHS(0),
fHisto_M_PHS(0),
fHisto_M_pt_Pi0(0),
fHisto_M_Pi0(0),
fHisto_M_pt_Eta(0),
fHisto_M_Eta(0),
fHisto_M_pt_All(0),
fHisto_M_All(0),
//fHistoE_NCells(0),
fClustStat(0),
fV0Stat(0),
ftest(0),
fGlobalTrackReference(),
fPIDResponse(0),
fTrackCuts(0)
{
    //Default constructor
    // Define input and output slots here
    fGlobalTrackReference.clear();
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTask_eta::~AliAnalysisTask_eta()
{
    //Destructor
    delete fOutputList;
    delete fTracks;
    delete fCaloClusters;
    
}
//________________________________________________________________________
void AliAnalysisTask_eta::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
    AliDebug(3, "Creating Output Objects");
    
    /////////////////////////////////////////////////
    //Automatic determination of the analysis mode//
    ////////////////////////////////////////////////
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");

    fPIDResponse = inputHandler->GetPIDResponse();
    
    
    ////////////////
    //Output list//
    ///////////////
    fOutputList = new TList();
    fOutputList->SetOwner();

    fNevents = new TH1F("fNevents","No of events",4,-0.5,3.5);
    fOutputList->Add(fNevents);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"All events");
    fNevents->GetXaxis()->SetBinLabel(2,"All after trigger selection");
    fNevents->GetXaxis()->SetBinLabel(3,"N of vtx contrib > 1");
    fNevents->GetXaxis()->SetBinLabel(4,"Vtx_{z}<10cm");

    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fOutputList->Add(fVtxZ);

    fClustStat = new TH1F("fClustStat","No. of Clusters",6,0,6);
    fOutputList->Add(fClustStat);
    fClustStat->GetXaxis()->SetBinLabel(1,"All");
    fClustStat->GetXaxis()->SetBinLabel(2,"After NCells cut");
    fClustStat->GetXaxis()->SetBinLabel(3,"After energy cut");
    fClustStat->GetXaxis()->SetBinLabel(4,"After TOF cut");
    fClustStat->GetXaxis()->SetBinLabel(5,"After #lambda0^{2} cut");
    fClustStat->GetXaxis()->SetBinLabel(6,"After trackmatching cuts");

    fV0Stat = new TH1F("fV0Stat", "No. of V0s",12,0,12);
    fOutputList->Add(fV0Stat);
    fV0Stat->GetXaxis()->SetBinLabel(1,"All");
    fV0Stat->GetXaxis()->SetBinLabel(2,"After NProngs cut");
    fV0Stat->GetXaxis()->SetBinLabel(3,"After NDaughters cut");
    fV0Stat->GetXaxis()->SetBinLabel(4,"After charge cut");
    fV0Stat->GetXaxis()->SetBinLabel(5,"After eta cut");
    fV0Stat->GetXaxis()->SetBinLabel(6,"After radius cut");
    fV0Stat->GetXaxis()->SetBinLabel(7,"After line cut");
    fV0Stat->GetXaxis()->SetBinLabel(8,"After Vertex z cut");
    fV0Stat->GetXaxis()->SetBinLabel(9,"After ?");
    fV0Stat->GetXaxis()->SetBinLabel(10,"After A-P cut");
    fV0Stat->GetXaxis()->SetBinLabel(11,"After track cuts");
    fV0Stat->GetXaxis()->SetBinLabel(12,"After PID cuts");

    fHistV0E = new TH1F("fHistV0E", "V0 E distribution", 500, 0., 5.);
    fHistV0E->GetXaxis()->SetTitle("E (GeV)");
    fHistV0E->GetYaxis()->SetTitle("N of V0s");
    fOutputList->Add(fHistV0E);

    fShapeParam = new TH1F("fShapeParam","Shape Parameters of EMCal clusters;#lambda0^{2};counts",1000, 0.0, 20);
    fOutputList->Add(fShapeParam);

    fShapeParam2 = new TH1F("fShapeParam2","Shape Parameters of EMCal clusters after cuts;#lambda0^{2};counts",1000, 0., 20);
    fOutputList->Add(fShapeParam2);

    fHistClustE = new TH1F("fHistClustE", "cluster energy distribution; Cluster E;counts", 5000, 0.0, 100.0);
    fOutputList->Add(fHistClustE);

    fHistClustE2 = new TH1F("fHistClustE2", "cluster energy distribution after cluster cuts; Cluster E;counts", 5000, 0.0, 100.0);
    fOutputList->Add(fHistClustE2);

    fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","Cluster #eta and #phi distribution;#eta;#phi",100,-1,1,200,-4,4);
    fOutputList->Add(fEMCClsEtaPhi);
/*    
    fHistoNCells = new TH1F("fHistoNCells","No of cells in a cluster;N^{EMC}_{cells};counts",30,0,30);
    fOutputList->Add(fHistoNCells);
    
    fHistoNCells2 = new TH1F("fHistoNCells2","No of cells in a cluster after cluster cuts;N^{EMC}_{cells};counts",30,0,30);
    fOutputList->Add(fHistoNCells2);
    
    fHistoE_NCells = new TH2F("fHistoE_NCells","No of cells in a cluster vs Energy;Cluster E;N^{PHI}_{cells}",600,0,100,30,0,30);
    fOutputList->Add(fHistoE_NCells);
*/
    ftof = new TH2F("ftof","Time of Flight vs. Cluster Energy",1000,0,100,1000,-5e-7,5e-7);
    fOutputList->Add(ftof);

    ftof2 = new TH2F("ftof2","Time of Flight vs. Cluster Energy after Cuts",1000,0,100,1000,-5e-7,5e-7);
    fOutputList->Add(ftof2);
    
    fHistoNtracksMatch = new TH1I("fHistoNtracksMatch","No of tracks matched to EMCal clusters by correction task;N_{matched tracks};counts",30,0,30);
    fOutputList->Add(fHistoNtracksMatch);
    
    fHistoTrackMatchedPHOS = new TH2F("fHistoTrackMatchedPHOS","Matching to charged tracks with PHOS clusters;#Delta#eta;#Delta#varphi",500,-55,55,500,-55,55);
    fOutputList->Add(fHistoTrackMatchedPHOS);

    fHistoTrackMatchedPHOS2 = new TH2F("fHistoTrackMatchedPHOS2","Matched Tracks after cluster cuts (PHOS);#Delta#eta;#Delta#varphi", 500,-55,55,500,-55,55);
    fOutputList->Add(fHistoTrackMatchedPHOS2);

    fHistoTrackMatchedEMC = new TH2F("fHistoTrackMatchedEMC","Matching to charged tracks with EMCal clusters;#Delta#eta;#Delta#varphi",500,-0.5,0.5,500,-0.5,0.5);
    fOutputList->Add(fHistoTrackMatchedEMC);

    fHistoTrackMatchedEMC2 = new TH2F("fHistoTrackMatchedEMC2","Matched Tracks after cluster cuts (EMCal);#Delta#eta;#Delta#varphi", 500,-0.5,0.5,500,-0.5,0.5);
    fOutputList->Add(fHistoTrackMatchedEMC2);

    fHisto_M_pt_EMC = new TH2F("fHisto_M_pt_EMC", "Mass vs pT,only EMCal clusters", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_EMC->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_EMC->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_EMC);
    
    fHisto_M_EMC = new TH1F("fHisto_M_EMC", "Mass of only EMCal clusters;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_EMC);

    fHisto_M_pt_PHS = new TH2F("fHisto_M_pt_PHS", "Mass vs pT,only PHOS clusters", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_PHS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_PHS->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_PHS);
    
    fHisto_M_PHS = new TH1F("fHisto_M_PHS", "Mass of only PHOS clusters;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_PHS); 
  
    fHisto_M_pt_Pi0 = new TH2F("fHisto_M_pt_Pi0", "Mass vs pT of calo clusters", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_Pi0->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_Pi0->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_Pi0);
    
    fHisto_M_Pi0 = new TH1F("fHisto_M_Pi0", "Mass of calo clusters;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_Pi0);

    fHisto_M_pt_Eta = new TH2F("fHisto_M_pt_Eta", "Mass vs pT of calo clusters w/o Pi0", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_Eta->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_Eta->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_Eta);
    
    fHisto_M_Eta = new TH1F("fHisto_M_Eta", "Mass of calo clusters w/o Pi0;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_Eta);

    fHisto_M_pt_V0 = new TH2F("fHisto_M_pt_V0", "Mass vs pT of V0s", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_V0->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_V0->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_V0);
    
    fHisto_M_V0 = new TH1F("fHisto_M_V0", "Mass of V0s;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_V0);

    fHistV0InvMassPi0 = new TH1F("fHistV0InvMassPi0", "Mass of V0s w/o Pi0", 1200, 0., 1.4);
    fHistV0InvMassPi0->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fOutputList->Add(fHistV0InvMassPi0);

    fHistV0InvMassPtPi0 = new TH2F("fHistV0InvMassPtPi0", "Inv. mass vs #it{p}_{T} of V0s w/o Pi0", 1200, 0., 1.4, 1000, 0., 100.);
    fHistV0InvMassPtPi0->GetXaxis()->SetTitle("Inv. mass, GeV/#it{c}^{2}");
    fHistV0InvMassPtPi0->GetYaxis()->SetTitle("#it{p}_{T}, GeV/#it{c}");
    fOutputList->Add(fHistV0InvMassPtPi0);

    fHisto_M_pt_All = new TH2F("fHisto_M_pt_All", "Mass vs pT, V0 and calo clusters", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_All->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_All->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_All);
    
    fHisto_M_All = new TH1F("fHisto_M_All", "Mass of V0 and calo clusters;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_All);

    ftest=new TH1F("ftest","Counts PHOS and EMCAL Clusters",3,0,3);
    ftest->GetXaxis()->SetBinLabel(1,"All");
    ftest->GetXaxis()->SetBinLabel(2,"EMCAL");
    ftest->GetXaxis()->SetBinLabel(3,"PHOS");
    fOutputList->Add(ftest);

    PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTask_eta::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    // Post output data.
    
    // Trigger selection according to task->SetCollisionCandidates()
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    //Bool_t isSelected = false;
    //isSelected = ((maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7);
    //if(!isSelected) return;
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
         printf("ERROR: fEvent not available\n");
        return;
    }
    
   if (!IsAODanalysis()) {
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (! fESD) {
         printf("fESD not available\n");
         return;
    }
   }

   if (IsAODanalysis()) {    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) {
         printf("fAOD not available\n");
        return;
    }
   }


        if(IsAODanalysis()) fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
        if(!IsAODanalysis()) fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("ESDFilterTracks"));



    fNevents->Fill(0);

    //////////
    //// Fired online Trigger
    ///////////
    /*
    TString firedTrigger;
    
    TString   MyTrigger1 = "CPHI7";
    TString   MyTrigger2 = "CEMC7";
    TString   MyTrigger3 = "CPHI7"; // CEMC7EG1-B-NOPF for 16l, CEMC7EG2-B-NOPF for 17p/q

    if(fESD){
        firedTrigger = fESD->GetFiredTriggerClasses();
    }
    if(fAOD){
        firedTrigger = fAOD->GetFiredTriggerClasses();
    }
    
    if(firedTrigger.Contains(MyTrigger1)!=1 || firedTrigger.Contains(MyTrigger2)!=1 || firedTrigger.Contains(MyTrigger3)!=1){
        return;
    }*/
    
    fNevents->Fill(1); //events after trigger selection
    
    
    ////////////////
    //Event vertex//
    ///////////////
    Int_t ntracks = -999;
    ntracks = fVevent->GetNumberOfTracks();
    
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<2) {
        PostData(1, fOutputList);
        return;
    }

    fNevents->Fill(2); //vertex with > 2 contributors
    
    Zvertex = pVtx->GetZ();
    if(TMath::Abs(Zvertex) > 10.0){
        PostData(1, fOutputList);
        return;
    }
    fNevents->Fill(3); //events after z vtx cut
    fVtxZ->Fill(Zvertex);
    
    ////////////////////////
    //cluster information//
    ///////////////////////
    Int_t Nclust = -999;
    Int_t NclustEMC;
    Int_t NclustPHS;
    Nclust = fVevent->GetNumberOfCaloClusters();
    
    printf("=============== N of clusters in event: %i =============== \n", Nclust);
    
    /////////////////////////////
    //bins for event mixing/////
    ////////////////////////////
    // may be needed later for event mixing studies
    int izvtx = GetZvtxBin(Zvertex);
    int imult = GetMultBin(Nclust);

    StoreGlobalTrackReference();

    TObjArray* arrayClust = new TObjArray();
    TObjArray* arrayClustEMC = new TObjArray();
    TObjArray* arrayClustPHS = new TObjArray();
    TObjArray* arrayV0 = new TObjArray();
    TObjArray* array = new TObjArray();



  //CLUSTER INITIALIZATION AND CUTS://


    Int_t MinNCells=2;
    Double_t MinLambda=0.1;
    Double_t MaxLambda=0.366;
    Double_t MaxChi2=6.25;
    Double_t MinPi0=0.12;
    Double_t MaxPi0=0.155;
    Double_t TOF=12.5e-9;
    Bool_t IsTrackMatched;




    for(int icl=0; icl<Nclust; icl++)
    {
        AliVCluster *clust = 0x0;
        clust = fVevent->GetCaloCluster(icl);
        if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);
	ftest->Fill(0);
	if(clust && clust->IsEMCAL()) ftest->Fill(1);
        if(clust && clust->IsPHOS())  ftest->Fill(2);
	if(clust)
        {
	    if(clust->IsPHOS() && clust->GetType()!=AliVCluster::kPHOSNeutral) continue;	//reject CPV clusters
            Double_t clustE = clust->E();
            Float_t  posx[3]; // cluster pos
            clust->GetPosition(posx);
            TVector3 clustpos(posx[0],posx[1],posx[2]);
            Double_t emcphi = clustpos.Phi();
            Double_t emceta = clustpos.Eta();
            fHistClustE->Fill(clustE);
            fEMCClsEtaPhi->Fill(emceta,emcphi);
//            fHistoE_NCells->Fill(clustE,clust->GetNCells());
//            fHistoNCells->Fill(clust->GetNCells());
	    ftof->Fill(clustE,clust->GetTOF());
            fClustStat->Fill(0);
            if(clust->GetNCells()<MinNCells) continue;						//NCells cut
	    fClustStat->Fill(1);
            if(clust->IsEMCAL() && clustE<fMinEEMC) continue;					//MinE cut (EMCal/PHOS specific)
	    if(clust->IsPHOS() && clustE<fMinEPHS) continue;
	    //if(clustE>fMaxE) continue;
	    fClustStat->Fill(2);
	    if(TMath::Abs(clust->GetTOF()) > TOF) continue;					//TOF cut
	    fClustStat->Fill(3);

	    if(clust->IsEMCAL()) {
//              fShapeParam->Fill(clust->GetM02());
              if(clust->GetM02() < MinLambda || clust->GetM02() > MaxLambda) continue;		//Shower shape cut (EMCal/PHOS specific)
	    }

	   if(clust->IsPHOS()) {
              fShapeParam->Fill(clust->Chi2());
		if(clust->Chi2() > MaxChi2) continue;
	     }

           fClustStat->Fill(4);


	//////   Track Matching   //////
	    IsTrackMatched=kFALSE;

	    //PHOS track matching
	    if(clust->IsPHOS()) {
              fHistoTrackMatchedPHOS->Fill(clust->GetTrackDz(), clust->GetTrackDx());
              if (sqrt(pow(clust->GetTrackDz(),2)+pow(clust->GetTrackDx(),2))<10) IsTrackMatched = kTRUE;
              if (!IsTrackMatched) fHistoTrackMatchedPHOS2->Fill(clust->GetTrackDz(), clust->GetTrackDx());
	    }

	    //EMCal track matching
	    else {
	    if(IsAODanalysis()) {
             if (clust->GetNTracksMatched() > 0) continue; // track matching tender
             fHistoNtracksMatch->Fill(clust->GetNTracksMatched());}
	    if(!IsAODanalysis()) {
             for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++) {
                
                AliESDtrack* esdTrack = fESD->GetTrack(itrk); //reconstructed track
                if(!esdTrack) {
                    // AliError(Form("ERROR: Could not retrieve any (AOD) track %d",itrk));
                    continue;
                }
                
                Double_t posTrk[3] = {0,0,0};
                esdTrack->GetXYZ(posTrk);
                TVector3 vposTrk(posTrk);
		Double_t pt=esdTrack->Pt();
                
                Double_t fMass          = 0.139;
                Double_t fStepSurface   = 20.;
                Float_t etaproj, phiproj, pttrackproj;
                
                AliExternalTrackParam *trackParam = const_cast<AliExternalTrackParam*>(esdTrack->GetInnerParam());
		//if (trackParam) printf("INFO: Track Param found \n");
                if(!trackParam) {
                //printf("ERROR: Track Param is ZERO \n");
                continue;}
                // magnetic field is loaded by the tender:
                AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(trackParam, 440., fMass, fStepSurface, etaproj, phiproj, pttrackproj);
                fHistoTrackMatchedEMC->Fill(etaproj-emceta, phiproj-emcphi);
                if (sqrt(pow(etaproj-emceta,2)+pow(phiproj-emcphi,2))<0.05) IsTrackMatched = kTRUE;
                if (!IsTrackMatched) fHistoTrackMatchedEMC2->Fill(etaproj-emceta, phiproj-emcphi);
	       }
              }
	    }

	    if (IsTrackMatched) continue;


	    
	    arrayClust->AddLast(clust);
	    if (clust->IsEMCAL()) arrayClustEMC->AddLast(clust);
	    else arrayClustPHS->AddLast(clust);
	    fClustStat->Fill(5);
            fHistClustE2->Fill(clustE);
//            if(clust->IsEMCAL()) fShapeParam2->Fill(clust->GetM02());
            if(clust->IsPHOS()) fShapeParam2->Fill(clust->Chi2());
//            fHistoNCells2->Fill(clust->GetNCells());
	    ftof2->Fill(clustE,clust->GetTOF());
     
        }
    }
    
    Nclust=arrayClust->GetEntries();
    NclustEMC=arrayClustEMC->GetEntries();
    NclustPHS=arrayClustPHS->GetEntries();

    ///////////////////
    //V0 information//
    //////////////////

	Int_t NV0 = fAOD->GetNumberOfV0s();
        printf("=============== N of V0 in event: %i =============== \n", NV0);
    
 
	for (int i = 0; i < NV0; i++) {
		AliAODv0 *v0 = fAOD->GetV0(i);
		if (!v0) continue;
		fV0Stat->Fill(0);
		// Cuts to the V0 selection in order to avoid false pairing
		if (v0->GetNProngs() != 2) continue;
		fV0Stat->Fill(1);
		if (v0->GetNDaughters() != 2) continue;
		fV0Stat->Fill(2);
		if (v0->GetCharge() != 0) continue;
		fV0Stat->Fill(3);
		if (TMath::Abs(v0->Eta()) > 0.9) continue;
		fV0Stat->Fill(4);
		if (v0->RadiusV0() < 5 || v0->RadiusV0() > 180) continue;
		fV0Stat->Fill(5);
		if (v0->RadiusV0() < TMath::Abs(v0->DecayVertexV0Z()) * TMath::Tan(2 * TMath::ATan(TMath::Exp(-0.9))) - 7) continue;            //line cut
		fV0Stat->Fill(6);
		if (TMath::Abs(v0->DecayVertexV0Z()) > 240) continue;
		fV0Stat->Fill(7);
		// daughter tracks
		AliAODTrack *pos =
			static_cast<AliAODTrack *>(fGlobalTrackReference[v0->GetPosID()]);
		AliAODTrack *neg =
			static_cast<AliAODTrack *>(fGlobalTrackReference[v0->GetNegID()]);
		if (!pos || !neg) continue;

		if (TMath::Abs(Psi_pair(neg, pos)) > 0.1) continue;
		fV0Stat->Fill(8);

		// if (v0->GetKFInfo(1,1,2) > 30) continue;		//chi2 cut in ESD data
		// if (v0->Chi2V0() > 30) continue;		        //chi2 cut doesn't work?

		if (v0->GetOnFlyStatus()) continue;  // select only offline v0
		// Get the coordinates of the primary vertex
		Double_t xPV = fAOD->GetPrimaryVertex()->GetX();
		Double_t yPV = fAOD->GetPrimaryVertex()->GetY();
		Double_t zPV = fAOD->GetPrimaryVertex()->GetZ();
		Double_t PV[3] = { xPV, yPV, zPV };
		// Calculate decay vertex variables:
		const float point = v0->CosPointingAngle(PV);
		if (point < 0.99) continue; // in AOD prefilter it's already 0.99 (?)
		//Armenteros Podolanski cuts
		const float armAlpha = v0->AlphaV0();
		const float armQt = v0->PtArmV0();
		if (TMath::Abs(armAlpha) > 0.95) continue;
		if (armQt > 0.05 * TMath::Sqrt(1 - (armAlpha*armAlpha) / (0.95*0.95))) continue;
		fV0Stat->Fill(9);		 //elliptic cut

		// track cuts
		// use AliESDtrackCuts functionality
		if (!(fTrackCuts->IsSelected(pos)) || !(fTrackCuts->IsSelected(neg)))
		{
			continue; // rejected by track cuts
		}

		if (pos->Charge() == neg->Charge()) continue;


		if (pos->Charge() < 0) {
			pos = neg;
			neg = static_cast<AliAODTrack *>(fGlobalTrackReference[v0->GetPosID()]);
		}
		fV0Stat->Fill(10);

		// PID cuts (TPC only)
		Float_t nSigmaTPC_legpos = fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kElectron);
		Float_t nSigmaTPC_legneg = fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kElectron);
		Float_t nSigmaTPCpio_legpos = fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kPion);
		Float_t nSigmaTPCpio_legneg = fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kPion);

		if (pos->P() < 0.4 && nSigmaTPCpio_legpos < 0.5) continue;
		if (pos->P() > 0.4 && nSigmaTPCpio_legpos < 3) continue;
		if (neg->P() < 0.4 && nSigmaTPCpio_legneg < 0.5) continue;
		if (neg->P() > 0.4 && nSigmaTPCpio_legneg < 3) continue;
		if (TMath::Abs(nSigmaTPC_legpos) > 3) continue;
		if (TMath::Abs(nSigmaTPC_legneg) > 3) continue;
		fV0Stat->Fill(11);
/*
		// fill histos for negative V0 leg
		fHistPt_neg->Fill(neg->Pt());
		fHistEta_neg->Fill(neg->Eta());
		fHistPhi_neg->Fill(neg->Phi());
		fHistTPCnSigmaEle_neg->Fill(neg->P(), nSigmaTPC_legneg);

		// fill histos for positive V0 leg
		fHistPt_pos->Fill(pos->Pt());
		fHistEta_pos->Fill(pos->Eta());
		fHistPhi_pos->Fill(pos->Phi());
		fHistTPCnSigmaEle_pos->Fill(pos->P(), nSigmaTPC_legpos);

		// V0 histos
		fHistV0Pt->Fill(v0->Pt());
		fHistV0Eta->Fill(v0->Eta());
		fHistV0Phi->Fill(v0->Phi());
		fHistV0R->Fill(v0->RadiusV0());
		fHistV0CosPA->Fill(point);
		fHistV0Chi2->Fill(v0->Chi2V0()); // doesn't work?
		fHistV0Psi_pair->Fill(Psi_pair(neg, pos));
		fHistV0ArmPod->Fill(armAlpha, armQt);
		fHistV0Z->Fill(v0->DecayVertexV0Z());
		fHistV0RvsZ->Fill(TMath::Abs(v0->DecayVertexV0Z()), v0->RadiusV0());
*/

		TLorentzVector lv;
		lv.SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), 0.0);

		double en = lv.E();
		fHistV0E->Fill(en);

		// fill array with V0 if V0 passes cuts 
		arrayV0->AddLast(v0);

	} // end V0 loop


	NV0=arrayV0->GetEntries();

    TLorentzVector Photon, Photon1, Photon2, Parent;
    Double_t vertex[3];
    Double_t E1=0.0;
    Double_t E2=0.0;
    Double_t E;

    if(fAOD)fAOD->GetVertex()->GetXYZ(vertex);
    if(fESD)fESD->GetVertex()->GetXYZ(vertex);


    // MAIN LOOPS

//V0 LOOPS


	// V0 loop for rejection of V0s with small opening angle
    for (int i=0;i<NV0;i++) {
        AliAODv0 *v0 = static_cast<AliAODv0 *>(arrayV0->At(i));
        if (!v0) continue;
        		   
        TLorentzVector lv1;
        lv1.SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), 0.0);
        
        for (int k=i+1;k<NV0;k++) {
            AliAODv0 *v02 = static_cast<AliAODv0 *>(arrayV0->At(k));
            if (!v02) continue;
            
            TLorentzVector lv2;
            lv2.SetPtEtaPhiM(v02->Pt(), v02->Eta(), v02->Phi(), 0.0);

			double a = lv1.Angle(lv2.Vect());

			if (a < 0.1)
			{
				// reject V0
				arrayV0->Remove(v0);

			}            
        
		} //end second loop
    } //end first loop

	NV0 = arrayV0->GetEntries();

	// V0 loop for rejection of V0s in pi0 mass window
	for (int i = 0; i<NV0; i++) {
		AliAODv0 *v0 = static_cast<AliAODv0 *>(arrayV0->At(i));
		if (!v0) continue;

		TLorentzVector lv1;
		lv1.SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), 0.0);


		for (int k = i + 1; k<NV0; k++) {
			AliAODv0 *v02 = static_cast<AliAODv0 *>(arrayV0->At(k));
			if (!v02) continue;

			TLorentzVector lv2;
			lv2.SetPtEtaPhiM(v02->Pt(), v02->Eta(), v02->Phi(), 0.0);

			double pt = (lv1 + lv2).Pt();
			double m = (lv1 + lv2).M();

			fHisto_M_V0->Fill(m);
			fHisto_M_pt_V0->Fill(m, pt);

			if (m > 0.12 && m < 0.15)
			{
				// reject V0
				arrayV0->Remove(v0);
				arrayV0->Remove(v02);

			}

		} //end second loop
	} //end first loop


	NV0 = arrayV0->GetEntries();

	// fill mass histograms after opening angle and pi0 mass cut
	for (int i = 0; i<NV0; i++) {
		AliAODv0 *v0 = static_cast<AliAODv0 *>(arrayV0->At(i));
		if (!v0) continue;

		TLorentzVector lv1;
		lv1.SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), 0.0);

		for (int k = i + 1; k<NV0; k++) {
			AliAODv0 *v02 = static_cast<AliAODv0 *>(arrayV0->At(k));
			if (!v02) continue;

			TLorentzVector lv2;
			lv2.SetPtEtaPhiM(v02->Pt(), v02->Eta(), v02->Phi(), 0.0);


			double pt = (lv1 + lv2).Pt();
			double m = (lv1 + lv2).M();

			fHistV0InvMassPi0->Fill(m);
			fHistV0InvMassPtPi0->Fill(m, pt);

		} //end second loop
	} //end first loop



    
//CLUSTER LOOPS
    

//Cluster loop to reject cluster pairs in Pi0 mass window:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
            AliVCluster *clust1 = static_cast<AliVCluster *>(arrayClust->At(icl));
            if (!clust1) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //to be used in the event mixing
            Photons[0][izvtx][imult].push_back( TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) );
            //printf("Size of Photons TLorentz Vector = %lu\n", Photons[0][izvtx][imult].size());           


            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2 = static_cast<AliVCluster *>(arrayClust->At(jcl));
                    if (!clust2) continue;
		    if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_Pi0->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_Pi0->Fill(Parent.M());
		    if (Parent.M()>MinPi0 && Parent.M()<MaxPi0) {
			arrayClust->Remove(clust1);
			arrayClust->Remove(clust2);
			arrayClustEMC->Remove(clust1);
			arrayClustEMC->Remove(clust2);
			arrayClustPHS->Remove(clust1);
			arrayClustPHS->Remove(clust2);
		    }
                } 
            }

    Nclust=arrayClust->GetEntries();
    NclustEMC=arrayClustEMC->GetEntries();
    NclustPHS=arrayClustPHS->GetEntries();

//Only EMCal clusters, pi0 cut:

    for(Int_t icl=0; icl<NclustEMC-1; icl++)
    {
            AliVCluster *clust1 = static_cast<AliVCluster *>(arrayClustEMC->At(icl));
            if (!clust1) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < NclustEMC; jcl++) 
	    {
                    AliVCluster *clust2 = static_cast<AliVCluster *>(arrayClustEMC->At(jcl));
                    if (!clust2) continue;
	    	    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_EMC->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_EMC->Fill(Parent.M());
		   
                } //close clust1
            }//close clust2

//Only PHOS clusters, pi0 cut:

    for(Int_t icl=0; icl<NclustPHS-1; icl++)
    {
            AliVCluster *clust1 = static_cast<AliVCluster *>(arrayClustPHS->At(icl));
            if (!clust1) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < NclustPHS; jcl++) 
	    {
            	    AliVCluster *clust2 = static_cast<AliVCluster *>(arrayClustPHS->At(jcl));
            	    if (!clust2) continue;
	    	    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_PHS->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_PHS->Fill(Parent.M());
		   
                } //close clust1
            }//close clust2
/*
//All clusters, pi0 cut, all pairings:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
            AliVCluster *clust1 = static_cast<AliVCluster *>(arrayClust->At(icl));
            if (!clust1) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
            	    AliVCluster *clust2 = static_cast<AliVCluster *>(arrayClust->At(jcl));
            	    if (!clust2) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_All->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_All->Fill(Parent.M());
		   
                } 
            }
*/
//All clusters, pi0 cut, only relevant pairings:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
            AliVCluster *clust1 = static_cast<AliVCluster *>(arrayClust->At(icl));
            if (!clust1) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
            	    AliVCluster *clust2 = static_cast<AliVCluster *>(arrayClust->At(jcl));
            	    if (!clust2) continue;
		    if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_Eta->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_Eta->Fill(Parent.M());
		   
                } 
            }


//Save all Lorentzvectors in array:

    for(int i=0;i<NV0;i++) {
	AliAODv0 *v0 = static_cast<AliAODv0 *>(arrayV0->At(i));
        if (!v0) continue;
	TLorentzVector *lvV0 = new TLorentzVector();
	lvV0->SetPtEtaPhiM(v0->Pt(), v0->Eta(), v0->Phi(), 0.0);
	array->AddLast(lvV0); }
    for(int i=0;i<Nclust;i++) {
	AliVCluster *clust = static_cast<AliVCluster *>(arrayClust->At(i));
        if (!clust) continue;
	E=clust->E();
	TLorentzVector *lvClust = new TLorentzVector();
        clust->GetMomentum(*lvClust,vertex);
        lvClust->SetPxPyPzE(lvClust->Px(), lvClust->Py(), lvClust->Pz(), E);

	array->AddLast(lvClust); }

//All Photons from V0 and Calorimeters
    for(Int_t i=0; i<Nclust+NV0-1; i++)
    {
            TLorentzVector *lv1 = static_cast<TLorentzVector *>(array->At(i));
            if (!lv1) continue;

            //Cluster loop2, for invariant mass.
            for (Int_t j = i+1; j < Nclust+NV0; j++) 
	    {
            	    TLorentzVector *lv2 = static_cast<TLorentzVector *>(array->At(j));
            	    if (!lv2) continue;
		    //if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
                    double pt = (*lv1 + *lv2).Pt();
		    double m = (*lv1 + *lv2).M();
                    fHisto_M_pt_All->Fill(m, pt);
                    fHisto_M_All->Fill(m);
		   
                } 
            }

    PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisTask_eta::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
    
}
//________________________________________________________________________
void AliAnalysisTask_eta::StoreGlobalTrackReference() {
    // This method was inherited form H. Beck & O. Arnold analysis
    // Stores the pointer to the global track
    // Modified to work with vectors
    
    fGlobalTrackReference.clear();
    fGlobalTrackReference.resize(1000);
    AliAODEvent *aodEvent = static_cast<AliAODEvent *>(fInputEvent);
    for (int iTrack = 0; iTrack < aodEvent->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack *>(aodEvent->GetTrack(iTrack));
        if (!track) continue;
        
        // Check that the id is positive
        if (track->GetID() < 0) continue;
        
        // Check id is not too big for buffer
        if (track->GetID() >= static_cast<int>(fGlobalTrackReference.size()))
            fGlobalTrackReference.resize(track->GetID() + 1);
        
        // Warn if we overwrite a track
        auto *trackRef = fGlobalTrackReference[track->GetID()];
        if (trackRef) {
            // Seems like there are FilterMap 0 tracks
            // that have zero TPCNcls, don't store these!
            if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) continue;
            
            // Imagine the other way around, the zero map zero clusters track
            // is stored and the good one wants to be added. We ommit the warning
            // and just overwrite the 'bad' track
            if (fGlobalTrackReference[track->GetID()]->GetFilterMap() ||
                fGlobalTrackReference[track->GetID()]->GetTPCNcls()) {
                // If we come here, there's a problem
                std::cout << "Warning! global track info already there! ";
                std::cout << "TPCNcls track1 "
                << (fGlobalTrackReference[track->GetID()])->GetTPCNcls()
                << " track2 " << track->GetTPCNcls();
                std::cout << " FilterMap track1 "
                << (fGlobalTrackReference[track->GetID()])->GetFilterMap()
                << " track 2" << track->GetFilterMap() << "\n";
                fGlobalTrackReference[track->GetID()] = nullptr;
            }
        }  // Two tracks same id
        
        // Assign the pointer
        (fGlobalTrackReference.at(track->GetID())) = track;
    }
}
//_______________________________________________________________________
Double_t AliAnalysisTask_eta::Psi_pair(AliAODTrack *neg, AliAODTrack *pos)
{
    if (!neg || !pos) return -1;
    Double_t x = TMath::ACos((neg->Px()*pos->Px() + neg->Py()*pos->Py() + neg->Pz()*pos->Pz()) / (neg->P()*pos->P()));
    Double_t t = neg->Eta() - pos->Eta();
    
    return TMath::ASin(t / x);
}
//________________________________________________________________________
Int_t AliAnalysisTask_eta::GetZvtxBin(Double_t vertZ)
{
    
    int izvtx = -1;
    
    if     (vertZ<-3.375)
        izvtx=0;
    else if(vertZ<-1.605)
        izvtx=1;
    else if(vertZ<-0.225)
        izvtx=2;
    else if(vertZ<1.065)
        izvtx=3;
    else if(vertZ<-2.445)
        izvtx=4;
    else if(vertZ<-4.245)
        izvtx=5;
    else
        izvtx=6;
    
    return izvtx;
}
//________________________________________________________________________
Int_t AliAnalysisTask_eta::GetMultBin(Int_t mult){
    
    int imult = -1;
    
    if     (mult<2)
        imult=0;
    else if(mult<3)
        imult=1;
    else if(mult<4)
        imult=2;
    else if(mult<8)
        imult=3;
    else if(mult<15)
        imult=4;
    else
        imult=5;
    
    return imult;  
}
