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
//fHistEAsym(0),
//fHistEAsymPt(0),
//fHisto_M_pt_NCells(0),
//fHisto_M_NCells(0),
//fHisto_M_pt_MinE(0),
//fHisto_M_MinE(0),
//fHisto_M_pt_TOF(0),
//fHisto_M_TOF(0),
//fHisto_M_pt_M02(0),
//fHisto_M_M02(0),
//fHisto_M_pt_Asym(0),
//fHisto_M_Asym(0),
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
ftest(0)
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
//fHistEAsym(0),
//fHistEAsymPt(0),
//fHisto_M_pt_NCells(0),
//fHisto_M_NCells(0),
//fHisto_M_pt_MinE(0),
//fHisto_M_MinE(0),
//fHisto_M_pt_TOF(0),
//fHisto_M_TOF(0),
//fHisto_M_pt_M02(0),
//fHisto_M_M02(0),
//fHisto_M_pt_Asym(0),
//fHisto_M_Asym(0),
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
ftest(0)
{
    //Default constructor
    // Define input and output slots here
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
    
    fClustStat = new TH1F("fClustStat","No. of Clusters",6,0,6);
    fOutputList->Add(fClustStat);
    fClustStat->GetXaxis()->SetBinLabel(1,"All");
    fClustStat->GetXaxis()->SetBinLabel(2,"After NCells cut");
    fClustStat->GetXaxis()->SetBinLabel(3,"After energy cut");
    fClustStat->GetXaxis()->SetBinLabel(4,"After TOF cut");
    fClustStat->GetXaxis()->SetBinLabel(5,"After #lambda0^{2} cut");
    fClustStat->GetXaxis()->SetBinLabel(6,"After trackmatching cuts");
  
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fOutputList->Add(fVtxZ);
    
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
/*   
    fHistEAsym = new TH1F("fHistEAsym","Energy asymmetry of photon pairs",2000,0.,1.1);
    fOutputList->Add(fHistEAsym);

    fHistEAsymPt = new TH2F("fHistEAsymPt","Energy asymmetry vs pT", 1000, 0.,1.1,1000,0.,100.);
    fHistEAsymPt->GetXaxis()->SetTitle("Asymmetry");
    fHistEAsymPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHistEAsymPt);

    fHisto_M_pt_NCells = new TH2F("fHisto_M_pt_NCells", "Mass vs pT after only NCells cut", 1400, 0.0, 1.4, 1200, 0.0, 100.0);
    fHisto_M_pt_NCells->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_NCells->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_NCells);
    
    fHisto_M_NCells = new TH1F("fHisto_M_NCells", "Mass after only NCells cut;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_NCells);

    fHisto_M_pt_MinE = new TH2F("fHisto_M_pt_MinE", "Mass vs pT after NCells and MinE cuts", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_MinE->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_MinE->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_MinE);
    
    fHisto_M_MinE = new TH1F("fHisto_M_MinE", "Mass after NCells and MinE cuts;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_MinE);
   
    fHisto_M_pt_TOF = new TH2F("fHisto_M_pt_TOF", "Mass vs pT after NCells, MinE and TOF cuts", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_TOF->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_TOF->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_TOF);
    
    fHisto_M_TOF = new TH1F("fHisto_M_TOF", "Mass after NCells, MinE and TOF cuts;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_TOF);

    fHisto_M_pt_M02 = new TH2F("fHisto_M_pt_M02", "Mass vs pT after additional #lambda0 cut", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_M02->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_M02->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_M02);
    
    fHisto_M_M02 = new TH1F("fHisto_M_M02", "Mass after additional #lambda0 cut;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_M02); 
*/
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
  
    fHisto_M_pt_Pi0 = new TH2F("fHisto_M_pt_Pi0", "Mass vs pT with Pi0", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_Pi0->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_Pi0->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_Pi0);
    
    fHisto_M_Pi0 = new TH1F("fHisto_M_Pi0", "Mass with Pi0;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_Pi0);

    fHisto_M_pt_All = new TH2F("fHisto_M_pt_All", "Mass vs pT, all pairings", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_All->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_All->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_All);
    
    fHisto_M_All = new TH1F("fHisto_M_All", "Mass with All pairings;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_All);

    fHisto_M_pt_Eta = new TH2F("fHisto_M_pt_Eta", "Mass vs pT", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_Eta->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_Eta->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_Eta);
    
    fHisto_M_Eta = new TH1F("fHisto_M_Eta", "Mass;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_Eta);
/*
    fHisto_M_pt_Asym = new TH2F("fHisto_M_pt_Asym", "Mass vs pT after additional EAsym cut", 1400, 0.0, 1.4, 1400, 0.0, 100.0);
    fHisto_M_pt_Asym->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    fHisto_M_pt_Asym->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutputList->Add(fHisto_M_pt_Asym);
    
    fHisto_M_Asym = new TH1F("fHisto_M_Asym", "Mass after additional EAsym cut;Mass [GeV/c^{2}];#", 1400, 0.0, 1.4);
    fOutputList->Add(fHisto_M_Asym);
*/
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
    
    /////////////////////////////
    //PHOS cluster information//
    ////////////////////////////
    Int_t Nclust = -999;
    Nclust = fVevent->GetNumberOfCaloClusters();
    
    printf("=============== N of clusters in event: %i =============== \n", Nclust);
    
    /////////////////////////////
    //bins for event mixing/////
    ////////////////////////////
    // may be needed later for event mixing studies
    int izvtx = GetZvtxBin(Zvertex);
    int imult = GetMultBin(Nclust);





  //CLUSTER INITIALIZATION AND CUTS://


    Int_t MinNCells=2;
    Double_t MinLambda=0.1;
    Double_t MaxLambda=0.366;
    Double_t MaxChi2=6.25;
    Double_t MinPi0=0.12;
    Double_t MaxPi0=0.155;
    Double_t TOF=12.5e-9;
//    Double_t MaxAsym=0.94;
    Bool_t IsTrackMatched;
    AliVCluster *ClustList[Nclust];
//    AliVCluster *ClustListNCells[Nclust];
//    AliVCluster *ClustListMinE[Nclust];
//    AliVCluster *ClustListM02[Nclust];
//    AliVCluster *ClustListTOF[Nclust];
    AliVCluster *ClustListEMC[Nclust];
    AliVCluster *ClustListPHS[Nclust];


    for(Int_t icl=0; icl<Nclust; icl++)
    {
	ClustList[icl]=0;
//	ClustListNCells[icl]=0;
//	ClustListMinE[icl]=0;
//	ClustListM02[icl]=0;
//	ClustListTOF[icl]=0;
	ClustListEMC[icl]=0;
	ClustListPHS[icl]=0;
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
//	    ClustListNCells[icl]=clust;
	    fClustStat->Fill(1);
            if(clust->IsEMCAL() && clustE<fMinEEMC) continue;					//MinE cut (EMCal/PHOS specific)
	    if(clust->IsPHOS() && clustE<fMinEPHS) continue;
	    //if(clustE>fMaxE) continue;
//	    ClustListMinE[icl]=clust;
	    fClustStat->Fill(2);
	    if(TMath::Abs(clust->GetTOF()) > TOF) continue;					//TOF cut
//	    ClustListTOF[icl]=clust;
	    fClustStat->Fill(3);

	    if(clust->IsEMCAL()) {
//              fShapeParam->Fill(clust->GetM02());
              if(clust->GetM02() < MinLambda || clust->GetM02() > MaxLambda) continue;		//Shower shape cut (EMCal/PHOS specific)
//	      ClustListM02[icl]=clust;
	    }

	   if(clust->IsPHOS()) {
              fShapeParam->Fill(clust->Chi2());
		if(clust->Chi2() > MaxChi2) continue;
//	        ClustListM02[icl]=clust;
	     }

           fClustStat->Fill(4);


	//////   Track Matching   //////
	    IsTrackMatched=kFALSE;

	    //PHOS track matching
	    if(clust->IsPHOS()) {
              fHistoTrackMatchedPHOS->Fill(clust->GetTrackDz(), clust->GetTrackDx());
              if (sqrt(pow(clust->GetTrackDz(),2)+pow(clust->GetTrackDx(),2))<10) IsTrackMatched = kTRUE;
              if (!IsTrackMatched) fHistoTrackMatchedPHOS2->Fill(clust->GetTrackDz(), clust->GetTrackDx());}

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


	    
	    ClustList[icl]=clust;
	    if (clust->IsEMCAL()) ClustListEMC[icl]=clust;
	    else ClustListPHS[icl]=clust;
	    fClustStat->Fill(5);
            fHistClustE2->Fill(clustE);
//            if(clust->IsEMCAL()) fShapeParam2->Fill(clust->GetM02());
            if(clust->IsPHOS()) fShapeParam2->Fill(clust->Chi2());
//            fHistoNCells2->Fill(clust->GetNCells());
	    ftof2->Fill(clustE,clust->GetTOF());
     
        }
    }
    






    // MAIN LOOPS
    
    TLorentzVector Photon1, Photon2, Parent;
    Double_t vertex[3];
    Double_t E1=0.0;
    Double_t E2=0.0;
    Double_t EAsym=0.0;

    if(fAOD)fAOD->GetVertex()->GetXYZ(vertex);
    if(fESD)fESD->GetVertex()->GetXYZ(vertex);
    
/*
//NCells cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustListNCells[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustListNCells[jcl];
	    	    if (clust2==0) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_NCells->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_NCells->Fill(Parent.M());
		   
                } 
            }

//MinE cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustListMinE[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustListMinE[jcl];
	    	    if (clust2==0) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_MinE->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_MinE->Fill(Parent.M());
		   
                } 
            }

//TOF cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustListTOF[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustListTOF[jcl];
	    	    if (clust2==0) continue;
	    	    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_TOF->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_TOF->Fill(Parent.M());
		   
                } //close clust1
            }//close clust2

//Shower shape cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustListM02[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustListM02[jcl];
	    	    if (clust2==0) continue;
	    	    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_M02->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_M02->Fill(Parent.M());
		   
                } //close clust1
            }//close clust2
*/
//Track matching cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustList[icl];
	    if (clust1==0) continue;
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
		    AliVCluster *clust2=ClustList[jcl];
	    	    if (clust2==0) continue;
		    if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
		    E2=clust2->E();
		    //EAsym = TMath::Abs((E1-E2)/(E1+E2));
	            //fHistEAsym->Fill(EAsym);
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
		    //if (Parent.Pt()>3. && Parent.M()<0.2) {
		      //fHistEAsymPt->Fill(EAsym,Parent.Pt());//}
                    fHisto_M_pt_Pi0->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_Pi0->Fill(Parent.M());
		    if (Parent.M()>MinPi0 && Parent.M()<MaxPi0) {
			ClustList[icl]=0;
			ClustList[jcl]=0;
			ClustListEMC[icl]=0;
			ClustListEMC[jcl]=0;
			ClustListPHS[icl]=0;
			ClustListPHS[jcl]=0;
		    }
                } 
            }

//Only EMCal, pi0 cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustListEMC[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustListEMC[jcl];
	    	    if (clust2==0) continue;
		    //if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
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

//Only PHOS, pi0 cut:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustListPHS[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustListPHS[jcl];
	    	    if (clust2==0) continue;
		    //if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
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

//Both, pi0 cut, all pairings:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustList[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustList[jcl];
	    	    if(clust2==0) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_All->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_All->Fill(Parent.M());
		    //EAsym = TMath::Abs((E1-E2)/(E1+E2));
          	    //if(EAsym > MaxAsym) continue;
                    //fHisto_M_pt_Asym->Fill(Parent.M(),Parent.Pt());
                    //fHisto_M_Asym->Fill(Parent.M());
		   
                } 
            }

//Both, pi0 cut, only relevant pairings:

    for(Int_t icl=0; icl<Nclust-1; icl++)
    {
	    AliVCluster *clust1=ClustList[icl];
	    if (clust1==0) continue;
	    E1=clust1->E();
            clust1->GetMomentum(Photon1,vertex);
            Photon1.SetPx(Photon1.Px());
            Photon1.SetPy(Photon1.Py());
            Photon1.SetPz(Photon1.Pz());

            //Cluster loop2, for invariant mass.
            for (Int_t jcl = icl+1; jcl < Nclust; jcl++) 
	    {
		    AliVCluster *clust2=ClustList[jcl];
	    	    if(clust2==0) continue;
		    if((fTrigger=="CPHI7" && clust1->IsEMCAL() && clust2->IsEMCAL()) || (fTrigger=="CEMC7" && clust1->IsPHOS() && clust2->IsPHOS())) continue;
		    E2=clust2->E();
                    clust2->GetMomentum(Photon2,vertex);
                    Photon2.SetPx(Photon2.Px());
                    Photon2.SetPy(Photon2.Py());
                    Photon2.SetPz(Photon2.Pz());

                    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
                    fHisto_M_pt_Eta->Fill(Parent.M(),Parent.Pt());
                    fHisto_M_Eta->Fill(Parent.M());
		    //EAsym = TMath::Abs((E1-E2)/(E1+E2));
          	    //if(EAsym > MaxAsym) continue;
                    //fHisto_M_pt_Asym->Fill(Parent.M(),Parent.Pt());
                    //fHisto_M_Asym->Fill(Parent.M());
		   
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
