#include <string>
#include <stdlib.h>
#include <cmath>
#include "Riostream.h"
#include <cassert>

#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TF1.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TNtuple.h>

int main() {
  gROOT->Reset();
  using namespace std;
  TChain chainData5fb("PATEventTree");
  TChain chainSigMC("PATEventTree");
  const int nRuns = 1;
  const char *FileList[nRuns] = {
    //"2017F_NoTriggerZJUpsSet1_1.txt",
    "Test.txt"
  };
  char filebuff[100];
  for (int iR = 0; iR < nRuns; iR++){
    FILE *RFile;
    RFile = fopen(FileList[iR],"r");
    if (RFile == NULL) std::cout << "Error in File List, No Files Found" << std::endl;
    else{
      while ( ! feof (RFile) ){
	fscanf(RFile,"%s", filebuff);
	
	if(filebuff=="-1") break;
	
	std::cout << "Adding File " << filebuff << " to data chain" << std::endl;
	chainData5fb.Add(filebuff);
	
	
      }
    }
    fclose(RFile);
  }
  
  std::cout<< "\nData files have been loaded" << std::endl;
  
  // create new file
  TFile *fFile = new TFile("SingleMuon2017F_Soft4Mu_Etab_JPsiUpsilonZ_NoTrigger_CompactTest.root","recreate");
  TTree *fTree = new TTree("PATEventTree", "PATEventTree"); 
  // -- general stuff
  unsigned int fRun, fEvent, fLumiBlock;
  int          fBX, fOrbit;
  unsigned int fTimeLo, fTimeHi;
  float        fBz;

  // -- event topology
  float fEvSphericity, fEvAplanarity, fEvLambda[3], fEvThrust, fEvThrust_Major, fEvThrust_Minor, fEvFW[7];

  // -- HLT Path information taken directly from HLT (not PAT)
  // for the array: 0 = path was run or not, 1 = path passed or not, 2 = path encountered exception or not
  bool  fHLTP_DoubleMu3[3], fHLTP_DoubleMu6[3], fHLTP_DoubleMu7[3], fHLTP_Dimuon0_Upsilon_Muon[3],fHLTP_Dimuon0_Jpsi3p5_Muon2[3], fHLTP_QuadMuon0_Dimuon0_Jpsi[3], fHLTP_Dimuon0_Jpsi[3], fHLTP_Dimuon10_Jpsi_Barrel[3], fHLTP_TripleMu5[3], fHLTP_Trimuon5_3p5_2_Upsilon_Muon[3], fHLTP_IsoMu20[3],fHLTP_IsoMu27[3], fHLTP_Mu50[3];//initially ,fHLTP_Dimuon0_Jpsi_Muon[3] changed for 2017 himal
  // HLT prescalersb
  unsigned int fHLTP_DoubleMu3_PS, fHLTP_DoubleMu6_PS, fHLTP_DoubleMu7_PS, fHLTP_Dimuon0_Jpsi3p5_Muon2_PS, fHLTP_Dimuon0_Jpsi_PS, fHLTP_Dimuon10_Jpsi_Barrel_PS, fHLTP_TripleMu5_PS,  fHLTP_Trimuon5_3p5_2_Upsilon_Muon_PS, fHLTP_IsoMu20_PS, fHLTP_Mu50_PS;//initially ,fHLTP_Dimuon0_Jpsi_Muon_PS changed for 2017 himal
  // HLT filters passed
  bool  fHLTP_DoubleMu3_Filters[5], fHLTP_DoubleMu6_Filters[5], fHLTP_DoubleMu7_Filters[5], fHLTP_TripleMu5_Filters[5], fHLTP_Dimuon0_Jpsi_Filters[7], fHLTP_Dimuon0_Jpsi3p5_Muon2_Filters[8], fHLTP_Dimuon10_Jpsi_Barrel_Filters[7],  fHLTP_Trimuon5_3p5_2_Upsilon_Muon_Filters[11], fHLTP_IsoMu20_Filters[5], fHLTP_Mu50_Filters[6];////initially ,fHLTP_Dimuon0_Jpsi_Muon_Filters[8] changed for 2017 himal

  // -- HLT objects from PAT
  const int HLTMAX = 1000;
  int   fHLTN;
  int   fHLT_Index[HLTMAX], fHLT_ToPc[HLTMAX], fHLT_ToJt[HLTMAX], fHLT_PdgId[HLTMAX];
  float fHLT_Mass[HLTMAX], fHLT_Energy[HLTMAX], fHLT_Et[HLTMAX], fHLT_P[HLTMAX], fHLT_Pt[HLTMAX], fHLT_Px[HLTMAX], fHLT_Py[HLTMAX], fHLT_Pz[HLTMAX], fHLT_Theta[HLTMAX], fHLT_Eta[HLTMAX], fHLT_Phi[HLTMAX];


//Himal addition for HLT object
  float	fHLT_Muon_Eta[HLTMAX], fHLT_Muon_Phi[HLTMAX],  fHLT_Muon_Pt[HLTMAX], fHLT_Muon_VertexmumuJpsi[HLTMAX], fHLT_Muon_TripleMuL3[HLTMAX], fHLT_MuonJpsi_L3Matching_Pt[5], fHLT_MuonJpsi_L3Matching_Eta[5], fHLT_MuonJpsi_L3Matching_Phi[5], fHLT_MuonJpsi_VtxMatching_Pt[5], fHLT_MuonJpsi_VtxMatching_Eta[5], fHLT_MuonJpsi_VtxMatching_Phi[5], fHLT_MuonUpsilon_L3Matching_Pt[5], fHLT_MuonUpsilon_L3Matching_Eta[5], fHLT_MuonUpsilon_L3Matching_Phi[5], fHLT_MuonUpsilon_VtxMatching_Pt[5], fHLT_MuonUpsilon_VtxMatching_Eta[5], fHLT_MuonUpsilon_VtxMatching_Phi[5], fHLT_MuonZ_L3Matching_Pt[5],  fHLT_MuonZ_L3Matching_Eta[5],  fHLT_MuonZ_L3Matching_Phi[5],  fHLT_MuonZ_VtxMatching_Pt[5],  fHLT_MuonZ_VtxMatching_Eta[5],  fHLT_MuonZ_VtxMatching_Phi[5];
//end




  bool  fHLT_Mu[HLTMAX][2], fHLT_Mu12[HLTMAX][2], fHLT_Mu15[HLTMAX][2], fHLT_Mu20[HLTMAX][2], fHLT_Mu24[HLTMAX][2], fHLT_Mu30[HLTMAX][2], fHLT_IsoMu12[HLTMAX][2], fHLT_IsoMu15[HLTMAX][2], fHLT_IsoMu17[HLTMAX][2], fHLT_IsoMu24[HLTMAX][2], fHLT_IsoMu30[HLTMAX][2], fHLT_DoubleMu3[HLTMAX][2], fHLT_DoubleMu6[HLTMAX][2], fHLT_DoubleMu7[HLTMAX][2], fHLT_Dimuon0_Upsilon_Muon[HLTMAX][2],fHLT_Dimuon0_Jpsi3p5_Muon2[HLTMAX][2], fHLT_Dimuon0_Jpsi[HLTMAX][2], fHLT_Dimuon7_Jpsi_Displaced[HLTMAX][2], fHLT_Dimuon7_Jpsi_X_Barrel[HLTMAX][2], fHLT_Dimuon10_Jpsi_Barrel[HLTMAX][2], fHLT_TripleMu5[HLTMAX][2], fHLT_Jet[HLTMAX][2];//initially fHLT_Dimuon0_Jpsi_Muon[HLTMAX][2] changed for 2017 by himal

  // -- Particles
  const int PARTMAX = 10000;
  int   fPcN, fTkN, fMuN, fElecN, fMiscTkN, fPhotN; //fHadrN, 
  int   fPcIndex[PARTMAX], fMuIndex[PARTMAX], fElecIndex[PARTMAX], fMiscTkIndex[PARTMAX], fPhotIndex[PARTMAX], fPcToGn[PARTMAX], fPcToTk[PARTMAX], fTkToPc[PARTMAX], fPcToPv[PARTMAX], fPcTkQuality[PARTMAX], fPcJtN[PARTMAX], fPcPdgId[PARTMAX], fPcPixHitN[PARTMAX], fPcPixLayN[PARTMAX], fPcStripHitN[PARTMAX], fPcStripLayN[PARTMAX]; //fHadrIndex[PARTMAX]
  float fPcCharge[PARTMAX], fPcChi2[PARTMAX], fPcNdof[PARTMAX], fPcEnergy[PARTMAX], fPcEt[PARTMAX], fPcP[PARTMAX], fPcPt[PARTMAX], fPcPx[PARTMAX], fPcPy[PARTMAX], fPcPz[PARTMAX], fPcTheta[PARTMAX], fPcEta[PARTMAX], fPcPhi[PARTMAX], fPcD0[PARTMAX], fPcDz[PARTMAX], fPcEtaErr[PARTMAX], fPcPhiErr[PARTMAX], fPcD0Err[PARTMAX], fPcDzErr[PARTMAX], fPcVx[PARTMAX], fPcVy[PARTMAX], fPcVz[PARTMAX], fPcEcalIso[PARTMAX], fPcHcalIso[PARTMAX], fPcTrackIso[PARTMAX], fPcIP[PARTMAX], fPcIPxy[PARTMAX];
  // -- muons
  int   fMuHitN[PARTMAX], fMuMatchedN[PARTMAX], fMuMatchedNSegArb[PARTMAX], fMuMatchedNSegTrkArb[PARTMAX], fMuMatchedNSegTrkArbClean[PARTMAX], fMuHLTN[PARTMAX], fMuToHLT[PARTMAX], fMuBestProbI[4], fTkBestProbI, fMuBestProbByPtI[4];
  float fMuChi2[PARTMAX], fMuNdof[PARTMAX], fMuTkKink[PARTMAX], fMuGlbKink[PARTMAX], fMuGlbProb[PARTMAX], fMuTkSADist[PARTMAX], fMuTkSAdR[PARTMAX], fMuECALEnergy[PARTMAX], fMuHCALEnergy[PARTMAX], fMuCalCompat[PARTMAX];
  bool  fMuIsGlobal[PARTMAX], fMuIsSoft[PARTMAX], fMuIsLoose[PARTMAX], fMuIsTracker[PARTMAX], fMuIsStandalone[PARTMAX], fMuIsCalo[PARTMAX], fMuArbitrated[PARTMAX], fMuLastStationLoose[PARTMAX], fMuLastStationTight[PARTMAX], fMu2DCompatibilityLoose[PARTMAX], fMu2DCompatibilityTight[PARTMAX], fMuOneStationLoose[PARTMAX], fMuOneStationTight[PARTMAX], fMuHLTMatch[PARTMAX][2], fMuL3Match[PARTMAX], fMuTightMatch[PARTMAX], fPcTPFilter[PARTMAX], fPcBasicFilter[PARTMAX];
  int ntkplusloose, ntkminusloose, nplussoft, nminussoft;

  // -- Primary Vertices
  const int PVMAX = 300;
  int   fPvN, fRePvN, fAllPvN;
  int   fPvIndex[PVMAX], fPvTkN[PVMAX], fPvClosestI[PVMAX];
  float fPvX[PVMAX], fPvY[PVMAX],  fPvZ[PVMAX], fPvXe[PVMAX], fPvYe[PVMAX],  fPvZe[PVMAX], fPvPx[PVMAX], fPvPy[PVMAX], fPvPz[PVMAX], fPvPt[PVMAX], fPvEta[PVMAX],  fPvChi2[PVMAX], fPvNdof[PVMAX], fPvMass[PVMAX];
  bool  fPvIsFake[PVMAX], fPvIsRefit[PVMAX];

  // -- Jets
  const int JETMAX = 2000;
  const int PARTINJETMAX = 15;
  const int PARTTOJETMAX = 5;

  int JtShift; //used to differentiate jet collections
  int fJtN, fJtStandN, fJtFatN, fJtSubN, fJtFiltN, fJtBdRN, fJtBFlavN;
  int fJtIndex[JETMAX], fJtStandIndex[JETMAX], fJtFatIndex[JETMAX], fJtSubIndex[JETMAX], fJtFiltIndex[JETMAX], fJtToPv[JETMAX], fJtHLTN[JETMAX], fJtToHLT[JETMAX], fJtTkN[JETMAX], fJtSsvN[JETMAX], fJtGtvN[JETMAX], fJtnConstituents[JETMAX], fJtFlavour[JETMAX], fJtn60[JETMAX], fJtn90[JETMAX], fJtnChargedParticles[JETMAX], fJtnNeutralParticles[JETMAX], fJtnChargedHadrons[JETMAX], fJtnNeutralHadrons[JETMAX], fJtnPhotons[JETMAX], fJtnElectrons[JETMAX], fJtnMuons[JETMAX], fJtnHFHadrons[JETMAX], fJtnHFEMParticles[JETMAX], fJtRankTCHE[JETMAX], fJtRankTCHP[JETMAX], fJtRankP[JETMAX], fJtRankBP[JETMAX], fJtRankSSVHE[JETMAX], fJtRankSSVHP[JETMAX], fJtRankCSV[JETMAX], fJtRankCSVMVA[JETMAX], fJtRankGT[JETMAX], fJtRankSE[JETMAX], fJtRankSM[JETMAX];
  float fJtCharge[JETMAX], fJtDiscTCHE[JETMAX], fJtDiscTCHP[JETMAX], fJtDiscP[JETMAX], fJtDiscBP[JETMAX], fJtDiscSSVHE[JETMAX], fJtDiscSSVHP[JETMAX], fJtDiscCSV[JETMAX], fJtDiscCSVMVA[JETMAX], fJtDiscGT[JETMAX], fJtDiscSE[JETMAX], fJtDiscSM[JETMAX], fJtMaxDist[JETMAX], fJtPhi[JETMAX], fJtTheta[JETMAX], fJtEta[JETMAX], fJtRapidity[JETMAX], fJtP[JETMAX], fJtPt[JETMAX], fJtPx[JETMAX], fJtPy[JETMAX], fJtPz[JETMAX], fJtEnergy[JETMAX], fJtEt[JETMAX], fJtMass[JETMAX], fJtMt[JETMAX], fJtVx[JETMAX], fJtVy[JETMAX], fJtVz[JETMAX], fJtChargedEmEnergy[JETMAX], fJtNeutralEmEnergy[JETMAX], fJtChargedHadronEnergy[JETMAX], fJtNeutralHadronEnergy[JETMAX], fJtPhotonEnergy[JETMAX], fJtElectronEnergy[JETMAX], fJtMuonEnergy[JETMAX], fJtHFHadronEnergy[JETMAX], fJtHFEMEnergy[JETMAX], fJtdRMean[JETMAX], fJtdRMax[JETMAX], fJtPtRelMean[JETMAX], fJtPtRelMax[JETMAX], fJtPtRelSum[JETMAX], fJtPullPx[JETMAX], fJtPullPy[JETMAX], fJtPullPz[JETMAX];
  bool fJtIsStandard[JETMAX], fJtIsFat[JETMAX], fJtIsSub[JETMAX], fJtIsFilt[JETMAX], fJtBFromH[JETMAX], fJtBdRMatch[JETMAX], fJtHLTMatch[JETMAX][2], fJtVeto[JETMAX];

  // -- References between particle and jet
  int   fPcToJt[PARTMAX][PARTTOJETMAX], fJtToPc[JETMAX][PARTINJETMAX];

  // -- Fat and subjet specific variables
  int   fFatSubN[JETMAX], fFatFiltN[JETMAX], fSubToFat[JETMAX], fFiltToFat[JETMAX];

  // for jet matching to B
  int fJtBdRIndex[JETMAX], fJtToBdRIndex[JETMAX], fJtBFlavIndex[JETMAX];
  float fJtBdR[JETMAX];  // dR to B hadron

  // -- Secondary Vertices
  const int SVMAX = 20;
  const int TKINSVMAX = 20;
  const int SVINJETMAX = 4;

  int   fSvN, fSsvN, fGtvN;
  int   fSvIndex[SVMAX], fSvTkN[SVMAX], fSvSeqInJt[SVMAX];
  float fSvX[SVMAX], fSvY[SVMAX], fSvZ[SVMAX], fSvXe[SVMAX], fSvYe[SVMAX],  fSvZe[SVMAX], fSvPx[SVMAX], fSvPy[SVMAX], fSvPz[SVMAX], fSvPt[SVMAX], fSvEta[SVMAX], fSvChi2[SVMAX], fSvNdof[SVMAX], fSvDist[SVMAX], fSvDistCM[SVMAX], fSvMass[SVMAX], fSvTau[SVMAX], fSvTauCM[SVMAX];
  bool  fSvIsGTV[SVMAX];

  // -- References between vertex and tracks
  int   fSvToPc[SVMAX][TKINSVMAX], fPcToSsv[PARTMAX], fPcToGtv[PARTMAX];
  // -- References between jet and SV
  int   fSvToJt[SVMAX], fJtToSsv[JETMAX][SVINJETMAX], fJtToGtv[JETMAX][SVINJETMAX];

  // -- MET information
  const int METMAX = 10;

  int fMETN;
  int fMETIndex[METMAX];
  float fMETPhi[METMAX], fMETTheta[METMAX], fMETEta[METMAX], fMETRapidity[METMAX], fMETCharge[METMAX], fMETP[METMAX], fMETPt[METMAX], fMETPx[METMAX], fMETPy[METMAX], fMETPz[METMAX], fMETEnergy[METMAX], fMETEt[METMAX], fMETMass[METMAX], fMETMt[METMAX], fMETVx[METMAX], fMETVy[METMAX], fMETVz[METMAX];

  // -- JPsi->MuMu candidates
  const int JPsiMAX = 1000;
  int fJPsiN, fJPsiMuMuN, fJPsiMuTkN, fBaseJPsiI[2], fJPsiBestProbI[2], fJPsiVtxBestProbI[2];
  int fJPsiIndex[JPsiMAX], fJPsiClosestPVinZ[JPsiMAX], fJPsiMuI[JPsiMAX][2], fJPsiMuCategory[JPsiMAX][2];
  float fJPsiCharge[JPsiMAX], fJPsiPhi[JPsiMAX], fJPsiTheta[JPsiMAX], fJPsiEta[JPsiMAX], fJPsiRapidity[JPsiMAX], fJPsiP[JPsiMAX], fJPsiPt[JPsiMAX], fJPsiPx[JPsiMAX], fJPsiPy[JPsiMAX], fJPsiPz[JPsiMAX], fJPsiEnergy[JPsiMAX], fJPsiEt[JPsiMAX], fJPsiMass[JPsiMAX], fJPsiMt[JPsiMAX], fJPsiChi2[JPsiMAX], fJPsiNdof[JPsiMAX], fJPsiVx[JPsiMAX], fJPsiVy[JPsiMAX], fJPsiVz[JPsiMAX], fJPsiVxE[JPsiMAX], fJPsiVyE[JPsiMAX], fJPsiVzE[JPsiMAX], fJPsiVtxPhi[JPsiMAX], fJPsiVtxTheta[JPsiMAX], fJPsiVtxEta[JPsiMAX], fJPsiVtxRapidity[JPsiMAX], fJPsiVtxP[JPsiMAX], fJPsiVtxPt[JPsiMAX], fJPsiVtxPx[JPsiMAX], fJPsiVtxPy[JPsiMAX], fJPsiVtxPz[JPsiMAX], fJPsiVtxEnergy[JPsiMAX], fJPsiVtxEt[JPsiMAX], fJPsiVtxMass[JPsiMAX], fJPsiVtxMt[JPsiMAX];
  bool fJPsiMuCutKin[JPsiMAX][2], fJPsiMuCutHLT[JPsiMAX][2], fJPsiMuCutIso[JPsiMAX][2], fJPsiMuCutSA[JPsiMAX][2], fJPsiMuCutTrk[JPsiMAX][2], fJPsiMuType[JPsiMAX][2][5], fJPsiBasicFilter[JPsiMAX], fJPsiVtxBasicFilter[JPsiMAX];
  // -- isolation information
  int fJPsiFromClosestI[JPsiMAX][5], fNJPsiSharedTk2[JPsiMAX], fNJPsiWithin05[JPsiMAX], fNJPsiWithin03[JPsiMAX];
  float fJPsidRToClosest[JPsiMAX][5];

  // -- Eta_b->2 J/Psi candidates
  const int ETABMAX = 500;

  int fEtabN, fBaseEtabI, fEtabBestMassI, fEtabVtxBestMassI, fEtabBestProbI, fEtabBest4MuProbI, fEtabVtxBestProbI, fEtabVtxBestPtI;
  int fEtabIndex[ETABMAX], fEtabDuplicatesI[ETABMAX], fEtabJPsiI[ETABMAX][2], fEtabMuI[ETABMAX][4], fEtabMuN[ETABMAX], fEtabL3MatchMuN[ETABMAX], fEtabToRePvI[ETABMAX];
  float fEtabCharge[ETABMAX], fEtabPhi[ETABMAX], fEtabTheta[ETABMAX], fEtabEta[ETABMAX], fEtabRapidity[ETABMAX], fEtabP[ETABMAX], fEtabPt[ETABMAX], fEtabPx[ETABMAX], fEtabPy[ETABMAX], fEtabPz[ETABMAX], fEtabEnergy[ETABMAX], fEtabEt[ETABMAX], fEtabMass[ETABMAX], fEtabMt[ETABMAX], fEtabChi2[ETABMAX], fEtabNdof[ETABMAX], fEtabVx[ETABMAX], fEtabVy[ETABMAX], fEtabVz[ETABMAX], fEtabVxE[ETABMAX], fEtabVyE[ETABMAX], fEtabVzE[ETABMAX], fEtabVtxPhi[ETABMAX], fEtabVtxTheta[ETABMAX], fEtabVtxEta[ETABMAX], fEtabVtxRapidity[ETABMAX], fEtabVtxP[ETABMAX], fEtabVtxPt[ETABMAX], fEtabVtxPx[ETABMAX], fEtabVtxPy[ETABMAX], fEtabVtxPz[ETABMAX], fEtabVtxEnergy[ETABMAX], fEtabVtxEt[ETABMAX], fEtabVtxMass[ETABMAX], fEtabVtxMt[ETABMAX], fEtabCT[ETABMAX], fEtabCTxy[ETABMAX], fEtabVtxCT[ETABMAX], fEtabVtxCTxy[ETABMAX], fEtabJPsiDeltaL[ETABMAX], fEtabJPsiDeltaT[ETABMAX], fEtabJPsiVtxErr[ETABMAX], fEtabJPsiVtxErrxy[ETABMAX], fEtabJPsiProjX[ETABMAX][2], fEtabJPsiProjY[ETABMAX][2], fEtabJPsiProjZ[ETABMAX][2], fEtabJPsiCT[ETABMAX][2], fEtabJPsiCTxy[ETABMAX][2], fEtabJPsiToPVVtxErr[ETABMAX][2], fEtabJPsiToPVVtxErrxy[ETABMAX][2], fEtabJPsiVtxCT[ETABMAX][2], fEtabJPsiVtxCTxy[ETABMAX][2];
  bool fEtabBasicFilter[ETABMAX], fEtabVtxBasicFilter[ETABMAX];
//Himal additions
  bool  fEtabBasicFilter_Muon[ETABMAX];


  // isolation information
  int fEtabJPsiIsoTkN[ETABMAX][2];
  float fEtabJPsiIso7PV[ETABMAX][2], fEtabJPsiIsoTkCA[ETABMAX][2];

  // -- H->bb candidates
  const int HMAX = 100;

  int fHN;
  int fHIndex[HMAX], fHJtI[HMAX][2];
  float fHCharge[HMAX], fHPhi[HMAX], fHTheta[HMAX], fHEta[HMAX], fHRapidity[HMAX], fHP[HMAX], fHPt[HMAX], fHPx[HMAX], fHPy[HMAX], fHPz[HMAX], fHEnergy[HMAX], fHEt[HMAX], fHMass[HMAX], fHMt[HMAX], fHVx[HMAX], fHVy[HMAX], fHVz[HMAX];

  // -- Gen Particles
  const int GENMAX = 3000;
  const int MOTHERMAX = 5;
  const int DAUGHTERMAX = 20;

  int fGnN, fGnBN;
  int fGnIndex[GENMAX], fGnBIndex[GENMAX], fGnIsJet[GENMAX], fGnLongLived[GENMAX], fGnPdgId[GENMAX], fGnNMothers[GENMAX], fGnNDaughters[GENMAX];
  float fGnCharge[GENMAX], fGnPhi[GENMAX], fGnTheta[GENMAX], fGnEta[GENMAX], fGnRapidity[GENMAX], fGnP[GENMAX], fGnPt[GENMAX], fGnPx[GENMAX], fGnPy[GENMAX], fGnPz[GENMAX], fGnEnergy[GENMAX], fGnEt[GENMAX], fGnMass[GENMAX], fGnMt[GENMAX], fGnVx[GENMAX], fGnVy[GENMAX], fGnVz[GENMAX], fGnDist[GENMAX], fGnDistCM[GENMAX], fGnTau[GENMAX], fGnTauCM[GENMAX];

  // -- References to daughter candidates
  int   fGnMotherIndex[GENMAX][MOTHERMAX], fGnDaughterIndex[GENMAX][DAUGHTERMAX];

  // -- Signal quantities (Etab and J/Psi information)
  int   fEtabI, fEtabdRMatchI, fJPsiI[2], fJPsidRMatchI[2], fMuI[2][2], fMudRMatchI[2][2], fMuByPtI[4], fMuByPtMatchI[4];
  float fMu11_4MomJPsi1CM[4], fMu12_4MomJPsi1CM[4], fMu21_4MomJPsi2CM[4], fMu22_4MomJPsi2CM[4], fJPsi1_4MomEtabCM[4], fJPsi2_4MomEtabCM[4], fEtab_4MomJPsi1CM[4], fEtab_4MomJPsi2CM[4];

  // create the branches in the new tree
  fTree->Branch("run",          &fRun,          "run/i");
  fTree->Branch("lumiblock",    &fLumiBlock,    "lumiblock/i");
  fTree->Branch("event",        &fEvent,        "event/i");
  fTree->Branch("bx",           &fBX,           "bx/I");
  fTree->Branch("orbit",        &fOrbit,        "orbit/I");
  fTree->Branch("bz",           &fBz,           "bz/F");
  fTree->Branch("tlo",          &fTimeLo,       "tlo/i");
  fTree->Branch("thi",          &fTimeHi,       "thi/i");

  fTree->Branch("EvSphericity",   &fEvSphericity,   "EvSphericity/F");
  fTree->Branch("EvAplanarity",   &fEvAplanarity,   "EvAplanarity/F");
  fTree->Branch("EvLambda",       fEvLambda,        "EvLambda[3]/F");
  fTree->Branch("EvThrust",       &fEvThrust,       "EvThrust/F");
  fTree->Branch("EvThrust_Major", &fEvThrust_Major, "EvThrust_Major/F");
  fTree->Branch("EvThrust_Minor", &fEvThrust_Minor, "EvThrust_Minor/F");
  fTree->Branch("EvFW",           fEvFW,            "EvFW[7]/F");

  fTree->Branch("HLTP_DoubleMu3",         fHLTP_DoubleMu3,       "HLTP_DoubleMu3[3]/O");
  fTree->Branch("HLTP_DoubleMu6",         fHLTP_DoubleMu6,       "HLTP_DoubleMu6[3]/O");
  fTree->Branch("HLTP_DoubleMu7",         fHLTP_DoubleMu7,       "HLTP_DoubleMu7[3]/O");
  fTree->Branch("HLTP_Dimuon0_Jpsi",      fHLTP_Dimuon0_Jpsi,    "HLTP_Dimuon0_Jpsi[3]/O");
  //fTree->Branch("HLTP_Dimuon0_Jpsi_Muon", fHLTP_Dimuon0_Jpsi_Muon, "HLTP_Dimuon0_Jpsi_Muon[3]/O");
  fTree->Branch("HLTP_Dimuon0_Jpsi3p5_Muon2", fHLTP_Dimuon0_Jpsi3p5_Muon2, "HLTP_Dimuon0_Jpsi3p5_Muon2[3]/O");
  fTree->Branch("HLTP_Trimuon5_3p5_2_Upsilon_Muon", fHLTP_Trimuon5_3p5_2_Upsilon_Muon, "HLTP_Trimuon5_3p5_2_Upsilon_Muon[3]/O");
  fTree->Branch("HLTP_Trimuon5_3p5_2_Upsilon_Muon_Filters", fHLTP_Trimuon5_3p5_2_Upsilon_Muon_Filters, "HLTP_Trimuon5_3p5_2_Upsilon_Muon_Filters[11]/O");
  fTree->Branch("HLTP_Trimuon5_3p5_2_Upsilon_Muon_PS", &fHLTP_Trimuon5_3p5_2_Upsilon_Muon_PS, "HLTP_Trimuon5_3p5_2_Upsilon_Muon_PS/i");
  fTree->Branch("HLTP_QuadMuon0_Dimuon0_Jpsi", fHLTP_QuadMuon0_Dimuon0_Jpsi, "HLTP_QuadMuon0_Dimuon0_Jpsi[3]/O");
  fTree->Branch("HLTP_Dimuon0_Upsilon_Muon", fHLTP_Dimuon0_Upsilon_Muon, "HLTP_Dimuon0_Upsilon_Muon[3]/O");
  fTree->Branch("HLTP_Dimuon10_Jpsi_Barrel", fHLTP_Dimuon10_Jpsi_Barrel, "HLTP_Dimuon10_Jpsi_Barrel[3]/O");
  fTree->Branch("HLTP_TripleMu5",         fHLTP_TripleMu5,       "HLTP_TripleMu5[3]/O");
  fTree->Branch("HLTP_DoubleMu3_PS",         &fHLTP_DoubleMu3_PS,       "HLTP_DoubleMu3_PS/i");
  fTree->Branch("HLTP_DoubleMu6_PS",         &fHLTP_DoubleMu6_PS,       "HLTP_DoubleMu6_PS/i");
  fTree->Branch("HLTP_DoubleMu7_PS",         &fHLTP_DoubleMu7_PS,       "HLTP_DoubleMu7_PS/i");
  //fTree->Branch("HLTP_Dimuon0_Jpsi_Muon_PS", &fHLTP_Dimuon0_Jpsi_Muon_PS, "HLTP_Dimuon0_Jpsi_Muon_PS/i");
  fTree->Branch("HLTP_Dimuon0_Jpsi3p5_Muon2_PS", &fHLTP_Dimuon0_Jpsi3p5_Muon2_PS, "HLTP_Dimuon0_Jpsi3p5_Muon2_PS/i");
  fTree->Branch("HLTP_Dimuon0_Jpsi_PS",      &fHLTP_Dimuon0_Jpsi_PS,    "HLTP_Dimuon0_Jpsi_PS/i");
  fTree->Branch("HLTP_Dimuon10_Jpsi_Barrel_PS", &fHLTP_Dimuon10_Jpsi_Barrel_PS, "HLTP_Dimuon10_Jpsi_Barrel_PS/i");
  fTree->Branch("HLTP_TripleMu5_PS",         &fHLTP_TripleMu5_PS,       "HLTP_TripleMu5_PS/i");
  fTree->Branch("HLTP_DoubleMu3_Filters",    fHLTP_DoubleMu3_Filters,  "HLTP_DoubleMu3_Filters[5]/O");
  fTree->Branch("HLTP_DoubleMu6_Filters",    fHLTP_DoubleMu6_Filters,  "HLTP_DoubleMu6_Filters[5]/O");
  fTree->Branch("HLTP_DoubleMu7_Filters",    fHLTP_DoubleMu7_Filters,  "HLTP_DoubleMu7_Filters[5]/O");
  fTree->Branch("HLTP_TripleMu5_Filters",    fHLTP_TripleMu5_Filters,  "HLTP_TripleMu5_Filters[5]/O");
  fTree->Branch("HLTP_Dimuon0_Jpsi_Filters", fHLTP_Dimuon0_Jpsi_Filters, "HLTP_Dimuon0_Jpsi_Filters[7]/O");
  //fTree->Branch("HLTP_Dimuon0_Jpsi_Muon_Filters", fHLTP_Dimuon0_Jpsi_Muon_Filters, "HLTP_Dimuon0_Jpsi_Muon_Filters[8]/O");
  fTree->Branch("HLTP_Dimuon0_Jpsi3p5_Muon2_Filters", fHLTP_Dimuon0_Jpsi3p5_Muon2_Filters, "HLTP_Dimuon0_Jpsi3p5_Muon2_Filters[8]/O");
  fTree->Branch("HLTP_Dimuon10_Jpsi_Barrel_Filters", fHLTP_Dimuon10_Jpsi_Barrel_Filters, "HLTP_Dimuon10_Jpsi_Barrel_Filters[7]/O");
  fTree->Branch("HLTP_IsoMu20", fHLTP_IsoMu20, "HLTP_IsoMu20[3]/O");
  fTree->Branch("HLTP_IsoMu27", fHLTP_IsoMu27, "HLTP_IsoMu27[3]/O");
  fTree->Branch("HLTP_Mu50", fHLTP_Mu50, "HLTP_Mu50[3]/O");
  fTree->Branch("HLTP_IsoMu20_Filters", fHLTP_IsoMu20_Filters, "HLTP_IsoMu20_Filters[7]/O");
  fTree->Branch("HLTP_Mu50_Filters", fHLTP_Mu50_Filters, "HLTP_Mu50_Filters[6]/O");
  fTree->Branch("HLTP_IsoMu20_PS", &fHLTP_IsoMu20_PS, "HLTP_IsoMu20_PS/i");
  fTree->Branch("HLTP_Mu50_PS", &fHLTP_Mu50_PS, "HLTP_Mu50_PS/i");

  fTree->Branch("HLTN",           &fHLTN,           "HLTN/I");
//himal addition of branch
  fTree->Branch("HLT_Muon_Eta",        fHLT_Muon_Eta,         "HLT_Muon_Eta[HLTN]/F");
  fTree->Branch("HLT_Muon_Phi",        fHLT_Muon_Phi,         "HLT_Muon_Phi[HLTN]/F");
  fTree->Branch("HLT_Muon_Pt",         fHLT_Muon_Pt,          "HLT_Muon_Pt[HLTN]/F");
  fTree->Branch("HLT_Muon_VertexmumuJpsi", fHLT_Muon_VertexmumuJpsi, "HLT_Muon_VertexmumuJpsi[HLTN]/F");
  fTree->Branch("HLT_Muon_TripleMuL3",  fHLT_Muon_TripleMuL3,          "HLT_Muon_TripleMuL3[HLTN]/F");
//end
  //trigger matching stuff
  //---Jpsi Trigger---------------
  fTree->Branch("HLT_MuonJpsi_L3Matching_Pt",        fHLT_MuonJpsi_L3Matching_Pt,         "HLT_MuonJpsi_L3Matching_Pt[5]/F");
  fTree->Branch("HLT_MuonJpsi_L3Matching_Eta",        fHLT_MuonJpsi_L3Matching_Eta,         "HLT_MuonJpsi_L3Matching_Eta[5]/F");
  fTree->Branch("HLT_MuonJpsi_L3Matching_Phi",        fHLT_MuonJpsi_L3Matching_Phi,         "HLT_MuonJpsi_L3Matching_Phi[5]/F");
  fTree->Branch("HLT_MuonJpsi_VtxMatching_Pt",        fHLT_MuonJpsi_VtxMatching_Pt,         "HLT_MuonJpsi_VtxMatching_Pt[5]/F");
  fTree->Branch("HLT_MuonJpsi_VtxMatching_Eta",        fHLT_MuonJpsi_VtxMatching_Eta,         "HLT_MuonJpsi_VtxMatching_Eta[5]/F");
  fTree->Branch("HLT_MuonJpsi_VtxMatching_Phi",        fHLT_MuonJpsi_VtxMatching_Phi,         "HLT_MuonJpsi_VtxMatching_Phi[5]/F");
  //------Upsilon Trigger-------------
  fTree->Branch("HLT_MuonUpsilon_L3Matching_Pt",        fHLT_MuonUpsilon_L3Matching_Pt,         "HLT_MuonUpsilon_L3Matching_Pt[5]/F");
  fTree->Branch("HLT_MuonUpsilon_L3Matching_Eta",        fHLT_MuonUpsilon_L3Matching_Eta,         "HLT_MuonUpsilon_L3Matching_Eta[5]/F");
  fTree->Branch("HLT_MuonUpsioln_L3Matching_Phi",        fHLT_MuonUpsilon_L3Matching_Phi,         "HLT_MuonUpsilon_L3Matching_Phi[5]/F");
  fTree->Branch("HLT_MuonUpsilon_VtxMatching_Pt",        fHLT_MuonUpsilon_VtxMatching_Pt,         "HLT_MuonUpsilon_VtxMatching_Pt[5]/F");
  fTree->Branch("HLT_MuonUpsilon_VtxMatching_Eta",        fHLT_MuonUpsilon_VtxMatching_Eta,         "HLT_MuonUpsilon_VtxMatching_Eta[5]/F");
  fTree->Branch("HLT_MuonUpsilon_VtxMatching_Phi",        fHLT_MuonUpsilon_VtxMatching_Phi,         "HLT_MuonUpsilon_VtxMatching_Phi[5]/F");
  //-----Z Trigger-------------------
  fTree->Branch("HLT_MuonZ_L3Matching_Pt",        fHLT_MuonZ_L3Matching_Pt,         "HLT_MuonZ_L3Matching_Pt[5]/F");
  fTree->Branch("HLT_MuonZ_L3Matching_Eta",        fHLT_MuonZ_L3Matching_Eta,         "HLT_MuonZ_L3Matching_Eta[5]/F");
  fTree->Branch("HLT_MuonZ_L3Matching_Phi",        fHLT_MuonZ_L3Matching_Phi,         "HLT_MuonZ_L3Matching_Phi[5]/F");
  fTree->Branch("HLT_MuonZ_VtxMatching_Pt",        fHLT_MuonZ_VtxMatching_Pt,         "HLT_MuonZ_VtxMatching_Pt[5]/F");
  fTree->Branch("HLT_MuonZ_VtxMatching_Eta",        fHLT_MuonZ_VtxMatching_Eta,         "HLT_MuonZ_VtxMatching_Eta[5]/F");
  fTree->Branch("HLT_MuonZ_VtxMatching_Phi",        fHLT_MuonZ_VtxMatching_Phi,         "HLT_MuonZ_VtxMatching_Phi[5]/F");
  
  fTree->Branch("HLT_Index",      fHLT_Index,       "HLT_Index[HLTN]/I");
  fTree->Branch("HLT_ToPc",       fHLT_ToPc,        "HLT_ToPc[HLTN]/I");
  fTree->Branch("HLT_ToJt",       fHLT_ToJt,        "HLT_ToJt[HLTN]/I");
  fTree->Branch("HLT_PdgId",      fHLT_PdgId,       "HLT_PdgId[HLTN]/I");
  fTree->Branch("HLT_Mass",       fHLT_Mass,        "HLT_Mass[HLTN]/F");
  fTree->Branch("HLT_Energy",     fHLT_Energy,      "HLT_Energy[HLTN]/F");
  fTree->Branch("HLT_Et",         fHLT_Et,          "HLT_Et[HLTN]/F");
  fTree->Branch("HLT_P",          fHLT_P,           "HLT_P[HLTN]/F");
  fTree->Branch("HLT_Pt",         fHLT_Pt,          "HLT_Pt[HLTN]/F");
  fTree->Branch("HLT_Px",         fHLT_Px,          "HLT_Px[HLTN]/F");
  fTree->Branch("HLT_Py",         fHLT_Py,          "HLT_Py[HLTN]/F");
  fTree->Branch("HLT_Pz",         fHLT_Pz,          "HLT_Pz[HLTN]/F");
  fTree->Branch("HLT_Theta",      fHLT_Theta,       "HLT_Theta[HLTN]/F");
  fTree->Branch("HLT_Eta",        fHLT_Eta,         "HLT_Eta[HLTN]/F");
  fTree->Branch("HLT_Phi",        fHLT_Phi,         "HLT_Phi[HLTN]/F");

  fTree->Branch("HLT_Mu",         fHLT_Mu,          "HLT_Mu[HLTN][2]/O");
  fTree->Branch("HLT_Mu12",       fHLT_Mu12,        "HLT_Mu12[HLTN][2]/O");
  fTree->Branch("HLT_Mu15",       fHLT_Mu15,        "HLT_Mu15[HLTN][2]/O");
  fTree->Branch("HLT_Mu20",       fHLT_Mu20,        "HLT_Mu20[HLTN][2]/O");
  fTree->Branch("HLT_Mu24",       fHLT_Mu24,        "HLT_Mu24[HLTN][2]/O");
  fTree->Branch("HLT_Mu30",       fHLT_Mu30,        "HLT_Mu30[HLTN][2]/O");
  fTree->Branch("HLT_IsoMu12",    fHLT_IsoMu12,     "HLT_IsoMu12[HLTN][2]/O");
  fTree->Branch("HLT_IsoMu15",    fHLT_IsoMu15,     "HLT_IsoMu15[HLTN][2]/O");
  fTree->Branch("HLT_IsoMu17",    fHLT_IsoMu17,     "HLT_IsoMu17[HLTN][2]/O");
  fTree->Branch("HLT_IsoMu24",    fHLT_IsoMu24,     "HLT_IsoMu24[HLTN][2]/O");
  fTree->Branch("HLT_IsoMu30",    fHLT_IsoMu30,     "HLT_IsoMu30[HLTN][2]/O");
  fTree->Branch("HLT_DoubleMu3",  fHLT_DoubleMu3,   "HLT_DoubleMu3[HLTN][2]/O");
  fTree->Branch("HLT_DoubleMu6",  fHLT_DoubleMu6,   "HLT_DoubleMu6[HLTN][2]/O");
  fTree->Branch("HLT_DoubleMu7",  fHLT_DoubleMu7,   "HLT_DoubleMu7[HLTN][2]/O");
  //fTree->Branch("HLT_Dimuon0_Jpsi_Muon", fHLT_Dimuon0_Jpsi_Muon, "HLT_Dimuon0_Jpsi_Muon[HLTN][2]/O");
  fTree->Branch("HLT_Dimuon0_Jpsi3p5_Muon2", fHLT_Dimuon0_Jpsi3p5_Muon2, "HLT_Dimuon0_Jpsi3p5_Muon2[HLTN][2]/O");//changed foe 2017 himal
  fTree->Branch("HLT_Dimuon0_Upsilon_Muon", fHLT_Dimuon0_Upsilon_Muon, "HLT_Dimuon0_Upsilon_Muon[HLTN][2]/O");
  fTree->Branch("HLT_Dimuon0_Jpsi", fHLT_Dimuon0_Jpsi, "HLT_Dimuon0_Jpsi[HLTN][2]/O");
  fTree->Branch("HLT_Dimuon7_Jpsi_Displaced", fHLT_Dimuon7_Jpsi_Displaced, "HLT_Dimuon7_Jpsi_Displaced[HLTN][2]/O");
  fTree->Branch("HLT_Dimuon7_Jpsi_X_Barrel", fHLT_Dimuon7_Jpsi_X_Barrel, "HLT_Dimuon7_Jpsi_X_Barrel[HLTN][2]/O");
  fTree->Branch("HLT_Dimuon10_Jpsi_Barrel", fHLT_Dimuon10_Jpsi_Barrel, "HLT_Dimuon10_Jpsi_Barrel[HLTN][2]/O");
  fTree->Branch("HLT_TripleMu5", fHLT_TripleMu5, "HLT_TripleMu5[HLTN][2]/O");
  fTree->Branch("HLT_Jet",        fHLT_Jet,         "HLT_Jet[HLTN][2]/O");

  fTree->Branch("PcN",          &fPcN,          "PcN/I");
  fTree->Branch("TkN",          &fTkN,          "TkN/I");
  fTree->Branch("MuN",          &fMuN,          "MuN/I");
  fTree->Branch("ElecN",        &fElecN,        "ElecN/I");
  //fTree->Branch("HadrN",        &fHadrN,        "HadrN/I");
  fTree->Branch("MiscTkN",      &fMiscTkN,      "MiscTkN/I");
  fTree->Branch("PhotN",        &fPhotN,        "PhotN/I");
  fTree->Branch("PcPlusLoose",  ntkplusloose,   "PcPlusLoose/I");
  fTree->Branch("PcMinusLoose", ntkminusloose,  "PcMinusLoose/I");
  fTree->Branch("PcPlusSoft",   nplussoft,      "PcPlusSoft/I");
  fTree->Branch("PcMinusSoft",  nminussoft,     "PcMinusSoft/I");
  fTree->Branch("PcIndex",      fPcIndex,       "PcIndex[PcN]/I");
  fTree->Branch("MuIndex",      fMuIndex,       "MuIndex[MuN]/I");
  fTree->Branch("ElecIndex",    fElecIndex,     "ElecIndex[ElecN]/I");
  //fTree->Branch("HadrIndex",    fHadrIndex,     "HadrIndex[HadrN]/I");
  fTree->Branch("MiscTkIndex",  fMiscTkIndex,   "MiscTkIndex[MiscTkN]/I");
  fTree->Branch("PhotIndex",    fPhotIndex,     "PhotIndex[PhotN]/I");
  fTree->Branch("PcToTk",       fPcToTk,        "PcToTk[PcN]/I");
  fTree->Branch("TkToPc",       fTkToPc,        "TkToPc[TkN]/I");
  fTree->Branch("PcToPv",       fPcToPv,        "PcToPv[PcN]/I");
  fTree->Branch("PcToSsv",      fPcToSsv,       "PcToSsv[PcN]/I");
  fTree->Branch("PcToGtv",      fPcToGtv,       "PcToGtv[PcN]/I");
  fTree->Branch("PcTkQuality",  fPcTkQuality,   "PcTkQuality[PcN]/I");
  fTree->Branch("PcJtN",        fPcJtN,         "PcJtN[PcN]/I");
  fTree->Branch("PcPdgId",      fPcPdgId,       "PcPdgId[PcN]/I");
  fTree->Branch("PcPixHitN",    fPcPixHitN,     "PcPixHitN[PcN]/I");
  fTree->Branch("PcPixLayN",    fPcPixLayN,     "PcPixLayN[PcN]/I");
  fTree->Branch("PcStripHitN",  fPcStripHitN,   "PcStripHitN[PcN]/I");
  fTree->Branch("PcStripLayN",  fPcStripLayN,   "PcStripLayN[PcN]/I");
  fTree->Branch("PcCharge",     fPcCharge,      "PcCharge[PcN]/F");
  fTree->Branch("PcChi2",       fPcChi2,        "PcChi2[PcN]/F");
  fTree->Branch("PcNdof",       fPcNdof,        "PcNdof[PcN]/F");
  fTree->Branch("PcEnergy",     fPcEnergy,      "PcEnergy[PcN]/F");
  fTree->Branch("PcEt",         fPcEt,          "PcEt[PcN]/F");
  fTree->Branch("PcP",          fPcP,           "PcP[PcN]/F");
  fTree->Branch("PcPt",         fPcPt,          "PcPt[PcN]/F");
  fTree->Branch("PcPx",         fPcPx,          "PcPx[PcN]/F");
  fTree->Branch("PcPy",         fPcPy,          "PcPy[PcN]/F");
  fTree->Branch("PcPz",         fPcPz,          "PcPz[PcN]/F");
  fTree->Branch("PcTheta",      fPcTheta,       "PcTheta[PcN]/F");
  fTree->Branch("PcEta",        fPcEta,         "PcEta[PcN]/F");
  fTree->Branch("PcPhi",        fPcPhi,         "PcPhi[PcN]/F");
  fTree->Branch("PcD0",         fPcD0,          "PcD0[PcN]/F");
  fTree->Branch("PcDz",         fPcDz,          "PcDz[PcN]/F");
  fTree->Branch("PcEtaErr",     fPcEtaErr,      "PcEtaErr[PcN]/F");
  fTree->Branch("PcPhiErr",     fPcPhiErr,      "PcPhiErr[PcN]/F");
  fTree->Branch("PcD0Err",      fPcD0Err,       "PcD0Err[PcN]/F");
  fTree->Branch("PcDzErr",      fPcDzErr,       "PcDzErr[PcN]/F");
  fTree->Branch("PcVx",         fPcVx,          "PcVx[PcN]/F");
  fTree->Branch("PcVy",         fPcVy,          "PcVy[PcN]/F");
  fTree->Branch("PcVz",         fPcVz,          "PcVz[PcN]/F");
  fTree->Branch("PcEcalIso",    fPcEcalIso,     "PcEcalIso[PcN]/F");
  fTree->Branch("PcHcalIso",    fPcHcalIso,     "PcHcalIso[PcN]/F");
  fTree->Branch("PcTrackIso",   fPcTrackIso,    "PcTrackIso[PcN]/F");
  fTree->Branch("PcIP",         fPcIP,          "PcIP[PcN]/F");
  fTree->Branch("PcIPxy",       fPcIPxy,        "PcIPxy[PcN]/F");
  fTree->Branch("MuHitN",       fMuHitN,        "MuHitN[MuN]/I");
  fTree->Branch("MuMatchedN",   fMuMatchedN,    "MuMatchedN[MuN]/I");
  fTree->Branch("MuMatchedNSegArb",   fMuMatchedNSegArb,    "MuMatchedNSegArb[MuN]/I");
  fTree->Branch("MuMatchedNSegTrkArb",   fMuMatchedNSegTrkArb,    "MuMatchedNSegTrkArb[MuN]/I");
  fTree->Branch("MuMatchedNSegTrkArbClean",   fMuMatchedNSegTrkArbClean,    "MuMatchedNSegTrkArbClean[MuN]/I");

  fTree->Branch("MuHLTN",       fMuHLTN,        "MuHLTN[MuN]/I");
  fTree->Branch("MuToHLT",      fMuToHLT,       "MuToHLT[MuN]/I");
  fTree->Branch("MuChi2",       fMuChi2,        "MuChi2[MuN]/F");
  fTree->Branch("MuNdof",       fMuNdof,        "MuNdof[MuN]/F");
  fTree->Branch("MuTkKink",     fMuTkKink,      "MuTkKink[MuN]/F");
  fTree->Branch("MuGlbKink",    fMuGlbKink,     "MuGlbKink[MuN]/F");
  fTree->Branch("MuGlbProb",    fMuGlbProb,     "MuGlbProb[MuN]/F");
  fTree->Branch("MuTkSADist",   fMuTkSADist,    "MuTkSADist[MuN]/F");
  fTree->Branch("MuTkSAdR",     fMuTkSAdR,      "MuTkSAdR[MuN]/F");
  fTree->Branch("MuECALEnergy", fMuECALEnergy,  "MuECALEnergy[MuN]/F");
  fTree->Branch("MuHCALEnergy", fMuHCALEnergy,  "MuHCALEnergy[MuN]/F");
  fTree->Branch("MuCalCompat",  fMuCalCompat,   "MuCalCompat[MuN]/F");
  fTree->Branch("MuIsSoft",     fMuIsSoft,      "MuIsSoft[MuN]/O");
  fTree->Branch("MuIsLoose",     fMuIsLoose,      "MuIsLoose[MuN]/O");
  fTree->Branch("MuIsGlobal",   fMuIsGlobal,    "MuIsGlobal[MuN]/O");
  fTree->Branch("MuIsTracker",  fMuIsTracker,   "MuIsTracker[MuN]/O");
  fTree->Branch("MuIsStandalone",fMuIsStandalone,"MuIsStandalone[MuN]/O");
  fTree->Branch("MuIsCalo",     fMuIsCalo,      "MuIsCalo[MuN]/O");
  fTree->Branch("MuArbitrated", fMuArbitrated,  "MuArbitrated[MuN]/O");
  fTree->Branch("MuLastStationLoose", fMuLastStationLoose, "MuLastStationLoose[MuN]/O");
  fTree->Branch("MuLastStationTight", fMuLastStationTight, "MuLastStationTight[MuN]/O");
  fTree->Branch("Mu2DCompatibilityLoose", fMu2DCompatibilityLoose, "Mu2DCompatibilityLoose[MuN]/O");
  fTree->Branch("Mu2DCompatibilityTight", fMu2DCompatibilityTight, "Mu2DCompatibilityTight[MuN]/O");
  fTree->Branch("MuOneStationLoose", fMuOneStationLoose, "MuOneStationLoose[MuN]/O");
  fTree->Branch("MuOneStationTight", fMuOneStationTight, "MuOneStationTight[MuN]/O");
  fTree->Branch("MuHLTMatch",   fMuHLTMatch,    "MuHLTMatch[MuN][2]/O");
  fTree->Branch("MuL3Match",    fMuL3Match,     "MuL3Match[MuN]/O");
  fTree->Branch("MuTightMatch", fMuTightMatch,  "MuTightMatch[MuN]/O");
  fTree->Branch("PcTPFilter",   fPcTPFilter,    "PcTPFilter[PcN]/O");
  fTree->Branch("PcBasicFilter",fPcBasicFilter, "PcBasicFilter[PcN]/O");
  fTree->Branch("MuBestProbI",  fMuBestProbI,   "MuBestProbI[4]/I");
  fTree->Branch("MuBestProbByPtI",fMuBestProbByPtI,"MuBestProbByPtI[4]/I");
  fTree->Branch("TkBestProbI",  &fTkBestProbI,  "TkBestProbI/I");

  fTree->Branch("PcToJt",       fPcToJt,        "PcToJt[PcN][5]/I");

  fTree->Branch("PvN",          &fPvN,          "PvN/I");
  fTree->Branch("RePvN",        &fRePvN,        "RePvN/I");
  fTree->Branch("AllPvN",       &fAllPvN,       "AllPvN/I");
  fTree->Branch("PvIndex",      fPvIndex,       "PvIndex[AllPvN]/I");
  fTree->Branch("PvTkN",        fPvTkN,         "PvTkN[AllPvN]/I");
  fTree->Branch("PvClosestI",   fPvClosestI,    "PvClosestI[AllPvN]/I");
  fTree->Branch("PvX",          fPvX,           "PvX[AllPvN]/F");
  fTree->Branch("PvY",          fPvY,           "PvY[AllPvN]/F");
  fTree->Branch("PvZ",          fPvZ,           "PvZ[AllPvN]/F");
  fTree->Branch("PvXe",         fPvXe,          "PvXe[AllPvN]/F");
  fTree->Branch("PvYe",         fPvYe,          "PvYe[AllPvN]/F");
  fTree->Branch("PvZe",         fPvZe,          "PvZe[AllPvN]/F");
  fTree->Branch("PvPx",         fPvPx,          "PvPx[AllPvN]/F");
  fTree->Branch("PvPy",         fPvPy,          "PvPy[AllPvN]/F");
  fTree->Branch("PvPz",         fPvPz,          "PvPz[AllPvN]/F");
  fTree->Branch("PvPt",         fPvPt,          "PvPt[AllPvN]/F");
  fTree->Branch("PvEta",        fPvEta,         "PvEta[AllPvN]/F");
  fTree->Branch("PvChi2",       fPvChi2,        "PvChi2[AllPvN]/F");
  fTree->Branch("PvNdof",       fPvNdof,        "PvNdof[AllPvN]/F");
  fTree->Branch("PvMass",       fPvMass,        "PvMass[AllPvN]/F");
  fTree->Branch("PvIsFake",     fPvIsFake,      "PvIsFake[AllPvN]/O");
  fTree->Branch("PvIsRefit",    fPvIsRefit,     "PvIsRefit[AllPvN]/O");

  fTree->Branch("JtN",          &fJtN,          "JtN/I");
  fTree->Branch("JtStandN",     &fJtStandN,     "JtStandN/I");
  fTree->Branch("JtTkN",        fJtTkN,         "JtTkN[JtN]/I");
  fTree->Branch("JtSsvN",       fJtSsvN,        "JtSsvN[JtN]/I");
  fTree->Branch("JtGtvN",       fJtGtvN,        "JtGtvN[JtN]/I");
  fTree->Branch("JtIndex",      fJtIndex,       "JtIndex[JtN]/I");
  fTree->Branch("JtStandIndex", fJtStandIndex,  "JtStandIndex[JtStandN]/I");
  fTree->Branch("JtToPv",       fJtToPv,        "JtToPv[JtN]/I");
  fTree->Branch("JtnConstituents", fJtnConstituents, "JtnConstituents[JtN]/I");
  fTree->Branch("Jtn60",        fJtn60,         "Jtn60[JtN]/I");
  fTree->Branch("Jtn90",        fJtn90,         "Jtn90[JtN]/I");
  fTree->Branch("JtnChargedParticles", fJtnChargedParticles, "JtnChargedParticles[JtN]/I");
  fTree->Branch("JtnNeutralParticles", fJtnNeutralParticles, "JtnNeutralParticles[JtN]/I");
  fTree->Branch("JtnChargedHadrons", fJtnChargedHadrons, "JtnChargedHadrons[JtN]/I");
  fTree->Branch("JtnNeutralHadrons", fJtnNeutralHadrons, "JtnNeutralHadrons[JtN]/I");
  fTree->Branch("JtnPhotons",   fJtnPhotons,    "JtnPhotons[JtN]/I");
  fTree->Branch("JtnElectrons", fJtnElectrons,  "JtnElectrons[JtN]/I");
  fTree->Branch("JtnMuons",     fJtnMuons,      "JtnMuons[JtN]/I");
  fTree->Branch("JtnHFHadrons", fJtnHFHadrons,  "JtnHFHadrons[JtN]/I");
  fTree->Branch("JtnHFEMParticles", fJtnHFEMParticles, "JtnHFEMParticles[JtN]/I");
  fTree->Branch("JtRankTCHE",   fJtRankTCHE,    "JtRankTCHE[JtN]/I");
  fTree->Branch("JtRankTCHP",   fJtRankTCHP,    "JtRankTCHP[JtN]/I");
  fTree->Branch("JtRankP",      fJtRankP,       "JtRankP[JtN]/I");
  fTree->Branch("JtRankBP",     fJtRankBP,      "JtRankBP[JtN]/I");
  fTree->Branch("JtRankSSVHE",  fJtRankSSVHE,   "JtRankSSVHE[JtN]/I");
  fTree->Branch("JtRankSSVHP",  fJtRankSSVHP,   "JtRankSSVHP[JtN]/I");
  fTree->Branch("JtRankCSV",    fJtRankCSV,     "JtRankCSV[JtN]/I");
  fTree->Branch("JtRankCSVMVA", fJtRankCSVMVA,  "JtRankCSVMVA[JtN]/I");
  fTree->Branch("JtRankGT",     fJtRankGT,      "JtRankGT[JtN]/I");
  fTree->Branch("JtRankSE",     fJtRankSE,      "JtRankSE[JtN]/I");
  fTree->Branch("JtRankSM",     fJtRankSM,      "JtRankSM[JtN]/I");
  fTree->Branch("JtCharge",     fJtCharge,      "JtCharge[JtN]/F");
  fTree->Branch("JtDiscTCHE",   fJtDiscTCHE,    "JtDiscTCHE[JtN]/F");
  fTree->Branch("JtDiscTCHP",   fJtDiscTCHP,    "JtDiscTCHP[JtN]/F");
  fTree->Branch("JtDiscP",      fJtDiscP,       "JtDiscP[JtN]/F");
  fTree->Branch("JtDiscBP",     fJtDiscBP,      "JtDiscBP[JtN]/F");
  fTree->Branch("JtDiscSSVHE",  fJtDiscSSVHE,   "JtDiscSSVHE[JtN]/F");
  fTree->Branch("JtDiscSSVHP",  fJtDiscSSVHP,   "JtDiscSSVHP[JtN]/F");
  fTree->Branch("JtDiscCSV",    fJtDiscCSV,     "JtDiscCSV[JtN]/F");
  fTree->Branch("JtDiscCSVMVA", fJtDiscCSVMVA,  "JtDiscCSVMVA[JtN]/F");
  fTree->Branch("JtDiscGT",     fJtDiscGT,      "JtDiscGT[JtN]/F");
  fTree->Branch("JtDiscSE",     fJtDiscSE,      "JtDiscSE[JtN]/F");
  fTree->Branch("JtDiscSM",     fJtDiscSM,      "JtDiscSM[JtN]/F");
  fTree->Branch("JtMaxDist",    fJtMaxDist,     "JtMaxDist[JtN]/F");
  fTree->Branch("JtPhi",        fJtPhi,         "JtPhi[JtN]/F");
  fTree->Branch("JtTheta",      fJtTheta,       "JtTheta[JtN]/F");
  fTree->Branch("JtEta",        fJtEta,         "JtEta[JtN]/F");
  fTree->Branch("JtRapidity",   fJtRapidity,    "JtRapidity[JtN]/F");
  fTree->Branch("JtP",          fJtP,           "JtP[JtN]/F");
  fTree->Branch("JtPt",         fJtPt,          "JtPt[JtN]/F");
  fTree->Branch("JtPx",         fJtPx,          "JtPx[JtN]/F");
  fTree->Branch("JtPy",         fJtPy,          "JtPy[JtN]/F");
  fTree->Branch("JtPz",         fJtPz,          "JtPz[JtN]/F");
  fTree->Branch("JtEnergy",     fJtEnergy,      "JtEnergy[JtN]/F");
  fTree->Branch("JtEt",         fJtEt,          "JtEt[JtN]/F");
  fTree->Branch("JtMass",       fJtMass,        "JtMass[JtN]/F");
  fTree->Branch("JtMt",         fJtMt,          "JtMt[JtN]/F");
  fTree->Branch("JtVx",         fJtVx,          "JtVx[JtN]/F");
  fTree->Branch("JtVy",         fJtVy,          "JtVy[JtN]/F");
  fTree->Branch("JtVz",         fJtVz,          "JtVz[JtN]/F");
  fTree->Branch("JtChargedEmEnergy", fJtChargedEmEnergy, "JtChargedEmEnergy[JtN]/F");
  fTree->Branch("JtNeutralEmEnergy", fJtNeutralEmEnergy, "JtNeutralEmEnergy[JtN]/F");
  fTree->Branch("JtChargedHadronEnergy", fJtChargedHadronEnergy, "JtChargedHadronEnergy[JtN]/F");
  fTree->Branch("JtNeutralHadronEnergy", fJtNeutralHadronEnergy, "JtNeutralHadronEnergy[JtN]/F");
  fTree->Branch("JtPhotonEnergy", fJtPhotonEnergy, "JtPhotonEnergy[JtN]/F");
  fTree->Branch("JtElectronEnergy", fJtElectronEnergy, "JtElectronEnergy[JtN]/F");
  fTree->Branch("JtMuonEnergy", fJtMuonEnergy,  "JtMuonEnergy[JtN]/F");
  fTree->Branch("JtHFHadronEnergy", fJtHFHadronEnergy, "JtHFHadronEnergy[JtN]/F");
  fTree->Branch("JtHFEMEnergy", fJtHFEMEnergy,  "JtHFEMEnergy[JtN]/F");
  fTree->Branch("JtdRMean",     fJtdRMean,      "JtdRMean[JtN]/F");
  fTree->Branch("JtdRMax",      fJtdRMax,       "JtdRMax[JtN]/F");
  fTree->Branch("JtPtRelMean",  fJtPtRelMean,   "JtPtRelMean[JtN]/F");
  fTree->Branch("JtPtRelMax",   fJtPtRelMax,    "JtPtRelMax[JtN]/F");
  fTree->Branch("JtPtRelSum",   fJtPtRelSum,    "JtPtRelSum[JtN]/F");
  fTree->Branch("JtPullPx",     fJtPullPx,      "JtPullPx[JtN]/F");
  fTree->Branch("JtPullPy",     fJtPullPy,      "JtPullPy[JtN]/F");
  fTree->Branch("JtPullPz",     fJtPullPz,      "JtPullPz[JtN]/F");
  fTree->Branch("JtIsStandard", fJtIsStandard,  "JtIsStandard[JtN]/O");
  fTree->Branch("JtIsFat",      fJtIsFat,       "JtIsFat[JtN]/O");
  fTree->Branch("JtIsSub",      fJtIsSub,       "JtIsSub[JtN]/O");
  fTree->Branch("JtIsFilt",     fJtIsFilt,      "JtIsFilt[JtN]/O");
  fTree->Branch("JtHLTMatch",   fJtHLTMatch,    "JtHLTMatch[JtN][2]/O");
  fTree->Branch("JtVeto",       fJtVeto,        "JtVeto[JtN]/O");

  fTree->Branch("JtHLTN",       fJtHLTN,        "JtHLTN[JtN]/I");
  fTree->Branch("JtToHLT",      fJtToHLT,       "JtToHLT[JtN]/I");
  fTree->Branch("JtToPc",       fJtToPc,        "JtToPc[JtN][20]/I");
  fTree->Branch("JtToSsv",      fJtToSsv,       "JtToSsv[JtN][2]/I");
  fTree->Branch("JtToGtv",      fJtToGtv,       "JtToGtv[JtN][4]/I");

  fTree->Branch("SvN",          &fSvN,          "SvN/I");
  fTree->Branch("SsvN",         &fSsvN,         "SsvN/I");
  fTree->Branch("GtvN",         &fGtvN,         "GtvN/I");
  fTree->Branch("SvIndex",      fSvIndex,       "SvIndex[SvN]/I");
  fTree->Branch("SvTkN",        fSvTkN,         "SvTkN[SvN]/I");
  fTree->Branch("SvToJt",       fSvToJt,        "SvToJt[SvN]/I");
  fTree->Branch("SvSeqInJt",    fSvSeqInJt,     "SvSeqInJt[SvN]/I");
  fTree->Branch("SvX",          fSvX,           "SvX[SvN]/F");
  fTree->Branch("SvY",          fSvY,           "SvY[SvN]/F");
  fTree->Branch("SvZ",          fSvZ,           "SvZ[SvN]/F");
  fTree->Branch("SvXe",         fSvXe,          "SvXe[SvN]/F");
  fTree->Branch("SvYe",         fSvYe,          "SvYe[SvN]/F");
  fTree->Branch("SvZe",         fSvZe,          "SvZe[SvN]/F");
  fTree->Branch("SvPx",         fSvPx,          "SvPx[SvN]/F");
  fTree->Branch("SvPy",         fSvPy,          "SvPy[SvN]/F");
  fTree->Branch("SvPz",         fSvPz,          "SvPz[SvN]/F");
  fTree->Branch("SvPt",         fSvPt,          "SvPt[SvN]/F");
  fTree->Branch("SvEta",        fSvEta,         "SvEta[SvN]/F");
  fTree->Branch("SvChi2",       fSvChi2,        "SvChi2[SvN]/F");
  fTree->Branch("SvNdof",       fSvNdof,        "SvNdof[SvN]/F");
  fTree->Branch("SvDist",       fSvDist,        "SvDist[SvN]/F");
  fTree->Branch("SvDistCM",     fSvDistCM,      "SvDistCM[SvN]/F");
  fTree->Branch("SvMass",       fSvMass,        "SvMass[SvN]/F");
  fTree->Branch("SvTau",        fSvTau,         "SvTau[SvN]/F");
  fTree->Branch("SvTauCM",      fSvTauCM,       "SvTauCM[SvN]/F");
  fTree->Branch("SvIsGTV",      fSvIsGTV,       "SvIsGTV[SvN]/O");

  fTree->Branch("SvToPc",       fSvToPc,        "SvToPc[SvN][20]/I");

  fTree->Branch("METN",         &fMETN,         "METN/I");
  fTree->Branch("METIndex",     fMETIndex,      "METIndex[METN]/I");
  fTree->Branch("METCharge",    fMETCharge,     "METCharge[METN]/F");
  fTree->Branch("METPhi",       fMETPhi,        "METPhi[METN]/F");
  fTree->Branch("METTheta",     fMETTheta,      "METTheta[METN]/F");
  fTree->Branch("METEta",       fMETEta,        "METEta[METN]/F");
  fTree->Branch("METRapidity",  fMETRapidity,   "METRapidity[METN]/F");
  fTree->Branch("METP",         fMETP,          "METP[METN]/F");
  fTree->Branch("METPt",        fMETPt,         "METPt[METN]/F");
  fTree->Branch("METPx",        fMETPx,         "METPx[METN]/F");
  fTree->Branch("METPy",        fMETPy,         "METPy[METN]/F");
  fTree->Branch("METPz",        fMETPz,         "METPz[METN]/F");
  fTree->Branch("METEnergy",    fMETEnergy,     "METEnergy[METN]/F");
  fTree->Branch("METEt",        fMETEt,         "METEt[METN]/F");
  fTree->Branch("METMass",      fMETMass,       "METMass[METN]/F");
  fTree->Branch("METMt",        fMETMt,         "METMt[METN]/F");
  fTree->Branch("METVx",        fMETVx,         "METVx[METN]/F");
  fTree->Branch("METVy",        fMETVy,         "METVy[METN]/F");
  fTree->Branch("METVz",        fMETVz,         "METVz[METN]/F");

  fTree->Branch("JPsiN",           &fJPsiN,           "JPsiN/I");
  fTree->Branch("JPsiMuMuN",       &fJPsiMuMuN,       "JPsiMuMuN/I");
  fTree->Branch("JPsiMuTkN",       &fJPsiMuTkN,       "JPsiMuTkN/I");
  fTree->Branch("BaseJPsiI",       fBaseJPsiI,        "BaseJPsiI[2]/I");
  fTree->Branch("JPsiBestProbI",   fJPsiBestProbI,    "JPsiBestProbI[2]/I");
  fTree->Branch("JPsiVtxBestProbI",fJPsiVtxBestProbI, "JPsiVtxBestProbI[2]/I");
  fTree->Branch("JPsiIndex",       fJPsiIndex,        "JPsiIndex[JPsiN]/I");
  fTree->Branch("JPsiClosestPVinZ", fJPsiClosestPVinZ, "JPsiClosestPVinZ[JPsiN]/I");
  fTree->Branch("JPsiCharge",      fJPsiCharge,       "JPsiCharge[JPsiN]/F");
  fTree->Branch("JPsiPhi",         fJPsiPhi,          "JPsiPhi[JPsiN]/F");
  fTree->Branch("JPsiTheta",       fJPsiTheta,        "JPsiTheta[JPsiN]/F");
  fTree->Branch("JPsiEta",         fJPsiEta,          "JPsiEta[JPsiN]/F");
  fTree->Branch("JPsiRapidity",    fJPsiRapidity,     "JPsiRapidity[JPsiN]/F");
  fTree->Branch("JPsiP",           fJPsiP,            "JPsiP[JPsiN]/F");
  fTree->Branch("JPsiPt",          fJPsiPt,           "JPsiPt[JPsiN]/F");
  fTree->Branch("JPsiPx",          fJPsiPx,           "JPsiPx[JPsiN]/F");
  fTree->Branch("JPsiPy",          fJPsiPy,           "JPsiPy[JPsiN]/F");
  fTree->Branch("JPsiPz",          fJPsiPz,           "JPsiPz[JPsiN]/F");
  fTree->Branch("JPsiEnergy",      fJPsiEnergy,       "JPsiEnergy[JPsiN]/F");
  fTree->Branch("JPsiEt",          fJPsiEt,           "JPsiEt[JPsiN]/F");
  fTree->Branch("JPsiMass",        fJPsiMass,         "JPsiMass[JPsiN]/F");
  fTree->Branch("JPsiMt",          fJPsiMt,           "JPsiMt[JPsiN]/F");
  fTree->Branch("JPsiChi2",        fJPsiChi2,         "JPsiChi2[JPsiN]/F");
  fTree->Branch("JPsiNdof",        fJPsiNdof,         "JPsiNdof[JPsiN]/F");
  fTree->Branch("JPsiVx",          fJPsiVx,           "JPsiVx[JPsiN]/F");
  fTree->Branch("JPsiVy",          fJPsiVy,           "JPsiVy[JPsiN]/F");
  fTree->Branch("JPsiVz",          fJPsiVz,           "JPsiVz[JPsiN]/F");
  fTree->Branch("JPsiVxE",         fJPsiVxE,          "JPsiVxE[JPsiN]/F");
  fTree->Branch("JPsiVyE",         fJPsiVyE,          "JPsiVyE[JPsiN]/F");
  fTree->Branch("JPsiVzE",         fJPsiVzE,          "JPsiVzE[JPsiN]/F");
  fTree->Branch("JPsiVtxPhi",      fJPsiVtxPhi,       "JPsiVtxPhi[JPsiN]/F");
  fTree->Branch("JPsiVtxTheta",    fJPsiVtxTheta,     "JPsiVtxTheta[JPsiN]/F");
  fTree->Branch("JPsiVtxEta",      fJPsiVtxEta,       "JPsiVtxEta[JPsiN]/F");
  fTree->Branch("JPsiVtxRapidity", fJPsiVtxRapidity,  "JPsiVtxRapidity[JPsiN]/F");
  fTree->Branch("JPsiVtxP",        fJPsiVtxP,         "JPsiVtxP[JPsiN]/F");
  fTree->Branch("JPsiVtxPt",       fJPsiVtxPt,        "JPsiVtxPt[JPsiN]/F");
  fTree->Branch("JPsiVtxPx",       fJPsiVtxPx,        "JPsiVtxPx[JPsiN]/F");
  fTree->Branch("JPsiVtxPy",       fJPsiVtxPy,        "JPsiVtxPy[JPsiN]/F");
  fTree->Branch("JPsiVtxPz",       fJPsiVtxPz,        "JPsiVtxPz[JPsiN]/F");
  fTree->Branch("JPsiVtxEnergy",   fJPsiVtxEnergy,    "JPsiVtxEnergy[JPsiN]/F");
  fTree->Branch("JPsiVtxEt",       fJPsiVtxEt,        "JPsiVtxEt[JPsiN]/F");
  fTree->Branch("JPsiVtxMass",     fJPsiVtxMass,      "JPsiVtxMass[JPsiN]/F");
  fTree->Branch("JPsiVtxMt",       fJPsiVtxMt,        "JPsiVtxMt[JPsiN]/F");

  fTree->Branch("JPsiMuI",         fJPsiMuI,          "JPsiMuI[JPsiN][2]/I");
  fTree->Branch("JPsiMuCategory",  fJPsiMuCategory,   "JPsiMuCategory[JPsiN][2]/I");
  fTree->Branch("JPsiMuCutKin",    fJPsiMuCutKin,     "JPsiMuCutKin[JPsiN][2]/O");
  fTree->Branch("JPsiMuCutHLT",    fJPsiMuCutHLT,     "JPsiMuCutHLT[JPsiN][2]/O");
  fTree->Branch("JPsiMuCutIso",    fJPsiMuCutIso,     "JPsiMuCutIso[JPsiN][2]/O");
  fTree->Branch("JPsiMuCutSA",     fJPsiMuCutSA,      "JPsiMuCutSA[JPsiN][2]/O");
  fTree->Branch("JPsiMuCutTrk",    fJPsiMuCutTrk,     "JPsiMuCutTrk[JPsiN][2]/O");
  fTree->Branch("JPsiMuType",      fJPsiMuType,       "JPsiMuType[JPsiN][2][5]/O");
  fTree->Branch("JPsiBasicFilter", fJPsiBasicFilter,  "JPsiBasicFilter[JPsiN]/O");
  fTree->Branch("JPsiVtxBasicFilter",fJPsiVtxBasicFilter,"JPsiVtxBasicFilter[JPsiN]/O");

  // isolation information
  fTree->Branch("JPsiFromClosestI",fJPsiFromClosestI, "JPsiFromClosestI[JPsiN][5]/I");
  fTree->Branch("NJPsiSharedTk2",  fNJPsiSharedTk2,   "NJPsiSharedTk2[JPsiN]/I");
  fTree->Branch("NJPsiWithin05",   fNJPsiWithin05,    "NJPsiWithin05[JPsiN]/I");
  fTree->Branch("NJPsiWithin03",   fNJPsiWithin03,    "NJPsiWithin03[JPsiN]/I");
  fTree->Branch("JPsidRToClosest", fJPsidRToClosest,  "JPsidRToClosest[JPsiN][5]/F");

  fTree->Branch("EtabN",        &fEtabN,        "EtabN/I");
  fTree->Branch("EtabIndex",    fEtabIndex,     "EtabIndex[EtabN]/I");
  fTree->Branch("EtabDuplicatesI", fEtabDuplicatesI, "EtabDuplicatesI[EtabN]/I");
  fTree->Branch("EtabCharge",   fEtabCharge,    "EtabCharge[EtabN]/F");
  fTree->Branch("EtabPhi",      fEtabPhi,       "EtabPhi[EtabN]/F");
  fTree->Branch("EtabTheta",    fEtabTheta,     "EtabTheta[EtabN]/F");
  fTree->Branch("EtabEta",      fEtabEta,       "EtabEta[EtabN]/F");
  fTree->Branch("EtabRapidity", fEtabRapidity,  "EtabRapidity[EtabN]/F");
  fTree->Branch("EtabP",        fEtabP,         "EtabP[EtabN]/F");
  fTree->Branch("EtabPt",       fEtabPt,        "EtabPt[EtabN]/F");
  fTree->Branch("EtabPx",       fEtabPx,        "EtabPx[EtabN]/F");
  fTree->Branch("EtabPy",       fEtabPy,        "EtabPy[EtabN]/F");
  fTree->Branch("EtabPz",       fEtabPz,        "EtabPz[EtabN]/F");
  fTree->Branch("EtabEnergy",   fEtabEnergy,    "EtabEnergy[EtabN]/F");
  fTree->Branch("EtabEt",       fEtabEt,        "EtabEt[EtabN]/F");
  fTree->Branch("EtabMass",     fEtabMass,      "EtabMass[EtabN]/F");
  fTree->Branch("EtabMt",       fEtabMt,        "EtabMt[EtabN]/F");
  fTree->Branch("EtabChi2",     fEtabChi2,      "EtabChi2[EtabN]/F");
  fTree->Branch("EtabNdof",     fEtabNdof,      "EtabNdof[EtabN]/F");
  fTree->Branch("EtabVx",       fEtabVx,        "EtabVx[EtabN]/F");
  fTree->Branch("EtabVy",       fEtabVy,        "EtabVy[EtabN]/F");
  fTree->Branch("EtabVz",       fEtabVz,        "EtabVz[EtabN]/F");
  fTree->Branch("EtabVxE",      fEtabVxE,       "EtabVxE[EtabN]/F");
  fTree->Branch("EtabVyE",      fEtabVyE,       "EtabVyE[EtabN]/F");
  fTree->Branch("EtabVzE",      fEtabVzE,       "EtabVzE[EtabN]/F");
  fTree->Branch("EtabVtxPhi",   fEtabVtxPhi,    "EtabVtxPhi[EtabN]/F");
  fTree->Branch("EtabVtxTheta", fEtabVtxTheta,  "EtabVtxTheta[EtabN]/F");
  fTree->Branch("EtabVtxEta",   fEtabVtxEta,    "EtabVtxEta[EtabN]/F");
  fTree->Branch("EtabVtxRapidity",fEtabVtxRapidity,"EtabVtxRapidity[EtabN]/F");
  fTree->Branch("EtabVtxP",     fEtabVtxP,      "EtabVtxP[EtabN]/F");
  fTree->Branch("EtabVtxPt",    fEtabVtxPt,     "EtabVtxPt[EtabN]/F");
  fTree->Branch("EtabVtxPx",    fEtabVtxPx,     "EtabVtxPx[EtabN]/F");
  fTree->Branch("EtabVtxPy",    fEtabVtxPy,     "EtabVtxPy[EtabN]/F");
  fTree->Branch("EtabVtxPz",    fEtabVtxPz,     "EtabVtxPz[EtabN]/F");
  fTree->Branch("EtabVtxEnergy",fEtabVtxEnergy, "EtabVtxEnergy[EtabN]/F");
  fTree->Branch("EtabVtxEt",    fEtabVtxEt,     "EtabVtxEt[EtabN]/F");
  fTree->Branch("EtabVtxMass",  fEtabVtxMass,   "EtabVtxMass[EtabN]/F");
  fTree->Branch("EtabVtxMt",    fEtabVtxMt,     "EtabVtxMt[EtabN]/F");

  fTree->Branch("BaseEtabI",    &fBaseEtabI,    "BaseEtabI/I");
  fTree->Branch("EtabBestMassI", &fEtabBestMassI, "EtabBestMassI/I");
  fTree->Branch("EtabVtxBestMassI", &fEtabVtxBestMassI, "EtabVtxBestMassI/I");
  fTree->Branch("EtabBestProbI", &fEtabBestProbI, "EtabBestProbI/I");
  fTree->Branch("EtabBest4MuProbI", &fEtabBest4MuProbI, "EtabBest4MuProbI/I");
  fTree->Branch("EtabVtxBestProbI", &fEtabVtxBestProbI, "EtabVtxBestProbI/I");
  fTree->Branch("EtabVtxBestPtI", &fEtabVtxBestPtI, "EtabVtxBestPtI/I");
  fTree->Branch("EtabJPsiI",    fEtabJPsiI,     "EtabJPsiI[EtabN][2]/I");
  fTree->Branch("EtabMuI",      fEtabMuI,       "EtabMuI[EtabN][4]/I");
  fTree->Branch("EtabMuN",      fEtabMuN,       "EtabMuN[EtabN]/I");
  fTree->Branch("EtabL3MatchMuN",fEtabL3MatchMuN,"EtabL3MatchMuN[EtabN]/I");
  fTree->Branch("EtabToRePvI",  fEtabToRePvI,   "EtabToRePvI[EtabN]/I");
  fTree->Branch("EtabCT",       fEtabCT,        "EtabCT[EtabN]/F");
  fTree->Branch("EtabCTxy",     fEtabCTxy,      "EtabCTxy[EtabN]/F");
  fTree->Branch("EtabVtxCT",    fEtabVtxCT,     "EtabVtxCT[EtabN]/F");
  fTree->Branch("EtabVtxCTxy",  fEtabVtxCTxy,   "EtabVtxCTxy[EtabN]/F");
  fTree->Branch("EtabJPsiDeltaL", fEtabJPsiDeltaL, "EtabJPsiDeltaL[EtabN]/F");
  fTree->Branch("EtabJPsiDeltaT", fEtabJPsiDeltaT, "EtabJPsiDeltaT[EtabN]/F");
  fTree->Branch("EtabJPsiVtxErr", fEtabJPsiVtxErr, "EtabJPsiVtxErr[EtabN]/F");
  fTree->Branch("EtabJPsiVtxErrxy",fEtabJPsiVtxErrxy,"EtabJPsiVtxErrxy[EtabN]/F");
  fTree->Branch("EtabJPsiProjX", fEtabJPsiProjX, "EtabJPsiProjX[EtabN][2]/F");
  fTree->Branch("EtabJPsiProjY", fEtabJPsiProjY, "EtabJPsiProjY[EtabN][2]/F");
  fTree->Branch("EtabJPsiProjZ", fEtabJPsiProjZ, "EtabJPsiProjZ[EtabN][2]/F");
  fTree->Branch("EtabJPsiCT",   fEtabJPsiCT,    "EtabJPsiCT[EtabN][2]/F");
  fTree->Branch("EtabJPsiCTxy", fEtabJPsiCTxy,  "EtabJPsiCTxy[EtabN][2]/F");
  fTree->Branch("EtabJPsiVtxCT",fEtabJPsiVtxCT, "EtabJPsiVtxCT[EtabN][2]/F");
  fTree->Branch("EtabJPsiVtxCTxy",fEtabJPsiVtxCTxy,"EtabJPsiVtxCTxy[EtabN][2]/F");
  fTree->Branch("EtabJPsiToPVVtxErr",fEtabJPsiToPVVtxErr,"EtabJPsiToPVVtxErr[EtabN][2]/F");
  fTree->Branch("EtabJPsiToPVVtxErrxy",fEtabJPsiToPVVtxErrxy,"EtabJPsiToPVVtxErrxy[EtabN][2]/F");
  fTree->Branch("EtabBasicFilter", fEtabBasicFilter, "EtabBasicFilter[EtabN]/O");

//himal additions
//  fTree->Branch(" fEtabBasicFilter", fEtabBasicFilter,  "fEtabBasicFilter[ETABMAX]/0")
//end
//
  fTree->Branch("EtabVtxBasicFilter", fEtabVtxBasicFilter, "EtabVtxBasicFilter[EtabN]/O");
  // isolation information
  fTree->Branch("EtabJPsiIsoTkN", fEtabJPsiIsoTkN, "EtabJPsiIsoTkN[EtabN][2]/I");
  fTree->Branch("EtabJPsiIso7PV", fEtabJPsiIso7PV, "EtabJPsiIso7PV[EtabN][2]/F");
  fTree->Branch("EtabJPsiIsoTkCA",fEtabJPsiIsoTkCA,"EtabJPsiIsoTkCA[EtabN][2]/F");

  fTree->Branch("HN",           &fHN,           "HN/I");
  fTree->Branch("HIndex",       fHIndex,        "HIndex[HN]/I");
  fTree->Branch("HCharge",      fHCharge,       "HCharge[HN]/F");
  fTree->Branch("HPhi",         fHPhi,          "HPhi[HN]/F");
  fTree->Branch("HTheta",       fHTheta,        "HTheta[HN]/F");
  fTree->Branch("HEta",         fHEta,          "HEta[HN]/F");
  fTree->Branch("HRapidity",    fHRapidity,     "HRapidity[HN]/F");
  fTree->Branch("HP",           fHP,            "HP[HN]/F");
  fTree->Branch("HPt",          fHPt,           "HPt[HN]/F");
  fTree->Branch("HPx",          fHPx,           "HPx[HN]/F");
  fTree->Branch("HPy",          fHPy,           "HPy[HN]/F");
  fTree->Branch("HPz",          fHPz,           "HPz[HN]/F");
  fTree->Branch("HEnergy",      fHEnergy,       "HEnergy[HN]/F");
  fTree->Branch("HEt",          fHEt,           "HEt[HN]/F");
  fTree->Branch("HMass",        fHMass,         "HMass[HN]/F");
  fTree->Branch("HMt",          fHMt,           "HMt[HN]/F");
  fTree->Branch("HVx",          fHVx,           "HVx[HN]/F");
  fTree->Branch("HVy",          fHVy,           "HVy[HN]/F");
  fTree->Branch("HVz",          fHVz,           "HVz[HN]/F");

  fTree->Branch("HJtI",         fHJtI,          "HJtI[HN][2]/I");

  // read in the following variables
  chainData5fb.SetBranchAddress("run",          &fRun);
  chainData5fb.SetBranchAddress("lumiblock",    &fLumiBlock);
  chainData5fb.SetBranchAddress("event",        &fEvent);
  chainData5fb.SetBranchAddress("bx",           &fBX);
  chainData5fb.SetBranchAddress("orbit",        &fOrbit);
  chainData5fb.SetBranchAddress("bz",           &fBz);
  chainData5fb.SetBranchAddress("tlo",          &fTimeLo);
  chainData5fb.SetBranchAddress("thi",          &fTimeHi);

  chainData5fb.SetBranchAddress("EvSphericity",   &fEvSphericity);
  chainData5fb.SetBranchAddress("EvAplanarity",   &fEvAplanarity);
  chainData5fb.SetBranchAddress("EvLambda",       fEvLambda);
  chainData5fb.SetBranchAddress("EvThrust",       &fEvThrust);
  chainData5fb.SetBranchAddress("EvThrust_Major", &fEvThrust_Major);
  chainData5fb.SetBranchAddress("EvThrust_Minor", &fEvThrust_Minor);
  chainData5fb.SetBranchAddress("EvFW",           fEvFW);

  chainData5fb.SetBranchAddress("HLTP_DoubleMu3",         fHLTP_DoubleMu3);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu6",         fHLTP_DoubleMu6);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu7",         fHLTP_DoubleMu7);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi",      fHLTP_Dimuon0_Jpsi);
  //chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_Muon", fHLTP_Dimuon0_Jpsi_Muon);
  //--------addition to add JPsi, Upsilon Z trigger-------------------------------------------
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi3p5_Muon2", fHLTP_Dimuon0_Jpsi3p5_Muon2); //himal addition for 2017 data
  chainData5fb.SetBranchAddress("HLTP_Trimuon5_3p5_2_Upsilon_Muon", fHLTP_Trimuon5_3p5_2_Upsilon_Muon); 
  chainData5fb.SetBranchAddress("HLTP_Trimuon5_3p5_2_Upsilon_Muon_PS", &fHLTP_Trimuon5_3p5_2_Upsilon_Muon_PS); 
  chainData5fb.SetBranchAddress("HLTP_Trimuon5_3p5_2_Upsilon_Muon_Filters", fHLTP_Trimuon5_3p5_2_Upsilon_Muon_Filters);
  chainData5fb.SetBranchAddress("HLTP_IsoMu20", fHLTP_IsoMu20);
  chainData5fb.SetBranchAddress("HLTP_IsoMu27", fHLTP_IsoMu27);
  chainData5fb.SetBranchAddress("HLTP_IsoMu20_PS", &fHLTP_IsoMu20_PS);
  chainData5fb.SetBranchAddress("HLTP_IsoMu20_Filters", fHLTP_IsoMu20_Filters);
  chainData5fb.SetBranchAddress("HLTP_Mu50", fHLTP_Mu50);
  chainData5fb.SetBranchAddress("HLTP_Mu50_PS", &fHLTP_Mu50_PS);
  chainData5fb.SetBranchAddress("HLTP_Mu50_Filters", fHLTP_Mu50_Filters);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_L3Matching_Pt", fHLT_MuonJpsi_L3Matching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_L3Matching_Eta", fHLT_MuonJpsi_L3Matching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_L3Matching_Phi", fHLT_MuonJpsi_L3Matching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_VtxMatching_Pt", fHLT_MuonJpsi_VtxMatching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_VtxMatching_Eta", fHLT_MuonJpsi_VtxMatching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_VtxMatching_Phi", fHLT_MuonJpsi_VtxMatching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonUpsilon_L3Matching_Pt", fHLT_MuonUpsilon_L3Matching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonUpsilon_L3Matching_Eta", fHLT_MuonUpsilon_L3Matching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonUpsilon_L3Matching_Phi", fHLT_MuonUpsilon_L3Matching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonUpsilon_VtxMatching_Pt", fHLT_MuonUpsilon_VtxMatching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonUpsilon_VtxMatching_Eta", fHLT_MuonUpsilon_VtxMatching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonUpsilon_VtxMatching_Phi", fHLT_MuonUpsilon_VtxMatching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonZ_L3Matching_Pt", fHLT_MuonZ_L3Matching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonZ_L3Matching_Eta", fHLT_MuonZ_L3Matching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonZ_L3Matching_Phi", fHLT_MuonZ_L3Matching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonZ_VtxMatching_Pt", fHLT_MuonZ_VtxMatching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonZ_VtxMatching_Eta", fHLT_MuonZ_VtxMatching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonZ_VtxMatching_Phi", fHLT_MuonZ_VtxMatching_Phi);
  //------------end----------------------------
  chainData5fb.SetBranchAddress("HLTP_QuadMuon0_Dimuon0_Jpsi", fHLTP_QuadMuon0_Dimuon0_Jpsi);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Upsilon_Muon", fHLTP_Dimuon0_Upsilon_Muon);
  chainData5fb.SetBranchAddress("HLTP_Dimuon10_Jpsi_Barrel", fHLTP_Dimuon10_Jpsi_Barrel);
  chainData5fb.SetBranchAddress("HLTP_TripleMu5",         fHLTP_TripleMu5);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu3_PS",         &fHLTP_DoubleMu3_PS);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu6_PS",         &fHLTP_DoubleMu6_PS);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu7_PS",         &fHLTP_DoubleMu7_PS);
  //chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_Muon_PS", &fHLTP_Dimuon0_Jpsi_Muon_PS);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi3p5_Muon2_PS", &fHLTP_Dimuon0_Jpsi3p5_Muon2_PS);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_PS",      &fHLTP_Dimuon0_Jpsi_PS);
  chainData5fb.SetBranchAddress("HLTP_Dimuon10_Jpsi_Barrel_PS", &fHLTP_Dimuon10_Jpsi_Barrel_PS);
  chainData5fb.SetBranchAddress("HLTP_TripleMu5_PS",         &fHLTP_TripleMu5_PS);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu3_Filters",    fHLTP_DoubleMu3_Filters);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu6_Filters",    fHLTP_DoubleMu6_Filters);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu7_Filters",    fHLTP_DoubleMu7_Filters);
  chainData5fb.SetBranchAddress("HLTP_TripleMu5_Filters",    fHLTP_TripleMu5_Filters);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_Filters", fHLTP_Dimuon0_Jpsi_Filters);
  //chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_Muon_Filters", fHLTP_Dimuon0_Jpsi_Muon_Filters);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi3p5_Muon2_Filters", fHLTP_Dimuon0_Jpsi3p5_Muon2_Filters);
  chainData5fb.SetBranchAddress("HLTP_Dimuon10_Jpsi_Barrel_Filters", fHLTP_Dimuon10_Jpsi_Barrel_Filters);
//Himal addition for  setbranchaddress
  chainData5fb.SetBranchAddress("HLT_Muon_Eta",        fHLT_Muon_Eta);
  chainData5fb.SetBranchAddress("HLT_Muon_Phi",        fHLT_Muon_Phi);
  chainData5fb.SetBranchAddress("HLT_Muon_Pt",         fHLT_Muon_Pt);
  chainData5fb.SetBranchAddress("HLT_Muon_VertexmumuJpsi",       fHLT_Muon_VertexmumuJpsi);
  chainData5fb.SetBranchAddress("HLT_Muon_TripleMuL3",         fHLT_Muon_TripleMuL3);
 //end


  chainData5fb.SetBranchAddress("HLTN",           &fHLTN);
  chainData5fb.SetBranchAddress("HLT_Index",      fHLT_Index);
  chainData5fb.SetBranchAddress("HLT_ToPc",       fHLT_ToPc);
  chainData5fb.SetBranchAddress("HLT_ToJt",       fHLT_ToJt);
  chainData5fb.SetBranchAddress("HLT_PdgId",      fHLT_PdgId);
  chainData5fb.SetBranchAddress("HLT_Mass",       fHLT_Mass);
  chainData5fb.SetBranchAddress("HLT_Energy",     fHLT_Energy);
  chainData5fb.SetBranchAddress("HLT_Et",         fHLT_Et);
  chainData5fb.SetBranchAddress("HLT_P",          fHLT_P);
  chainData5fb.SetBranchAddress("HLT_Pt",         fHLT_Pt);
  chainData5fb.SetBranchAddress("HLT_Px",         fHLT_Px);
  chainData5fb.SetBranchAddress("HLT_Py",         fHLT_Py);
  chainData5fb.SetBranchAddress("HLT_Pz",         fHLT_Pz);
  chainData5fb.SetBranchAddress("HLT_Theta",      fHLT_Theta);
  chainData5fb.SetBranchAddress("HLT_Eta",        fHLT_Eta);
  chainData5fb.SetBranchAddress("HLT_Phi",        fHLT_Phi);

  chainData5fb.SetBranchAddress("HLT_Mu",         fHLT_Mu);
  chainData5fb.SetBranchAddress("HLT_Mu12",       fHLT_Mu12);
  chainData5fb.SetBranchAddress("HLT_Mu15",       fHLT_Mu15);
  chainData5fb.SetBranchAddress("HLT_Mu20",       fHLT_Mu20);
  chainData5fb.SetBranchAddress("HLT_Mu24",       fHLT_Mu24);
  chainData5fb.SetBranchAddress("HLT_Mu30",       fHLT_Mu30);
  chainData5fb.SetBranchAddress("HLT_IsoMu12",    fHLT_IsoMu12);
  chainData5fb.SetBranchAddress("HLT_IsoMu15",    fHLT_IsoMu15);
  chainData5fb.SetBranchAddress("HLT_IsoMu17",    fHLT_IsoMu17);
  chainData5fb.SetBranchAddress("HLT_IsoMu24",    fHLT_IsoMu24);
  chainData5fb.SetBranchAddress("HLT_IsoMu30",    fHLT_IsoMu30);
  chainData5fb.SetBranchAddress("HLT_DoubleMu3",  fHLT_DoubleMu3);
  chainData5fb.SetBranchAddress("HLT_DoubleMu6",  fHLT_DoubleMu6);
  chainData5fb.SetBranchAddress("HLT_DoubleMu7",  fHLT_DoubleMu7);
  //chainData5fb.SetBranchAddress("HLT_Dimuon0_Jpsi_Muon", fHLT_Dimuon0_Jpsi_Muon);
  chainData5fb.SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", fHLT_Dimuon0_Jpsi3p5_Muon2);//changed for 2017 by himal
  chainData5fb.SetBranchAddress("HLT_Dimuon0_Upsilon_Muon", fHLT_Dimuon0_Upsilon_Muon);
  chainData5fb.SetBranchAddress("HLT_Dimuon0_Jpsi", fHLT_Dimuon0_Jpsi);
  chainData5fb.SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced", fHLT_Dimuon7_Jpsi_Displaced);
  chainData5fb.SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel", fHLT_Dimuon7_Jpsi_X_Barrel);
  chainData5fb.SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel", fHLT_Dimuon10_Jpsi_Barrel);
  chainData5fb.SetBranchAddress("HLT_TripleMu5", fHLT_TripleMu5);
  chainData5fb.SetBranchAddress("HLT_Jet",        fHLT_Jet);

  chainData5fb.SetBranchAddress("PcN",          &fPcN);
  chainData5fb.SetBranchAddress("TkN",          &fTkN);
  chainData5fb.SetBranchAddress("MuN",          &fMuN);
  chainData5fb.SetBranchAddress("ElecN",        &fElecN);
  chainData5fb.SetBranchAddress("MiscTkN",      &fMiscTkN);
  chainData5fb.SetBranchAddress("PhotN",        &fPhotN);
  chainData5fb.SetBranchAddress("PcIndex",      fPcIndex);
  chainData5fb.SetBranchAddress("MuIndex",      fMuIndex);
  chainData5fb.SetBranchAddress("ElecIndex",    fElecIndex);
  chainData5fb.SetBranchAddress("MiscTkIndex",  fMiscTkIndex);
  chainData5fb.SetBranchAddress("PhotIndex",    fPhotIndex);
  chainData5fb.SetBranchAddress("PcToTk",       fPcToTk);
  chainData5fb.SetBranchAddress("TkToPc",       fTkToPc);
  chainData5fb.SetBranchAddress("PcToPv",       fPcToPv);
  chainData5fb.SetBranchAddress("PcToSsv",      fPcToSsv);
  chainData5fb.SetBranchAddress("PcToGtv",      fPcToGtv);
  chainData5fb.SetBranchAddress("PcTkQuality",  fPcTkQuality);
  chainData5fb.SetBranchAddress("PcJtN",        fPcJtN);
  chainData5fb.SetBranchAddress("PcPdgId",      fPcPdgId);
  chainData5fb.SetBranchAddress("PcPixHitN",    fPcPixHitN);
  chainData5fb.SetBranchAddress("PcPixLayN",    fPcPixLayN);
  chainData5fb.SetBranchAddress("PcStripHitN",  fPcStripHitN);
  chainData5fb.SetBranchAddress("PcStripLayN",  fPcStripLayN);
  chainData5fb.SetBranchAddress("PcCharge",     fPcCharge);
  chainData5fb.SetBranchAddress("PcChi2",       fPcChi2);
  chainData5fb.SetBranchAddress("PcNdof",       fPcNdof);
  chainData5fb.SetBranchAddress("PcEnergy",     fPcEnergy);
  chainData5fb.SetBranchAddress("PcEt",         fPcEt);
  chainData5fb.SetBranchAddress("PcP",          fPcP);
  chainData5fb.SetBranchAddress("PcPt",         fPcPt);
  chainData5fb.SetBranchAddress("PcPx",         fPcPx);
  chainData5fb.SetBranchAddress("PcPy",         fPcPy);
  chainData5fb.SetBranchAddress("PcPz",         fPcPz);
  chainData5fb.SetBranchAddress("PcTheta",      fPcTheta);
  chainData5fb.SetBranchAddress("PcEta",        fPcEta);
  chainData5fb.SetBranchAddress("PcPhi",        fPcPhi);
  chainData5fb.SetBranchAddress("PcD0",         fPcD0);
  chainData5fb.SetBranchAddress("PcDz",         fPcDz);
  chainData5fb.SetBranchAddress("PcEtaErr",     fPcEtaErr);
  chainData5fb.SetBranchAddress("PcPhiErr",     fPcPhiErr);
  chainData5fb.SetBranchAddress("PcD0Err",      fPcD0Err);
  chainData5fb.SetBranchAddress("PcDzErr",      fPcDzErr);
  chainData5fb.SetBranchAddress("PcVx",         fPcVx);
  chainData5fb.SetBranchAddress("PcVy",         fPcVy);
  chainData5fb.SetBranchAddress("PcVz",         fPcVz);
  chainData5fb.SetBranchAddress("PcEcalIso",    fPcEcalIso);
  chainData5fb.SetBranchAddress("PcHcalIso",    fPcHcalIso);
  chainData5fb.SetBranchAddress("PcTrackIso",   fPcTrackIso);
  chainData5fb.SetBranchAddress("PcIP",         fPcIP);
  chainData5fb.SetBranchAddress("PcIPxy",       fPcIPxy);
  chainData5fb.SetBranchAddress("MuHitN",       fMuHitN);
  chainData5fb.SetBranchAddress("MuMatchedN",   fMuMatchedN);
  chainData5fb.SetBranchAddress("MuMatchedNSegArb",   fMuMatchedNSegArb);
  chainData5fb.SetBranchAddress("MuMatchedNSegTrkArb",   fMuMatchedNSegTrkArb);
  chainData5fb.SetBranchAddress("MuMatchedNSegTrkArbClean",   fMuMatchedNSegTrkArbClean);
  chainData5fb.SetBranchAddress("MuHLTN",       fMuHLTN);
  chainData5fb.SetBranchAddress("MuToHLT",      fMuToHLT);
  chainData5fb.SetBranchAddress("MuChi2",       fMuChi2);
  chainData5fb.SetBranchAddress("MuNdof",       fMuNdof);
  chainData5fb.SetBranchAddress("MuTkKink",     fMuTkKink);
  chainData5fb.SetBranchAddress("MuGlbKink",    fMuGlbKink);
  chainData5fb.SetBranchAddress("MuGlbProb",    fMuGlbProb);
  chainData5fb.SetBranchAddress("MuTkSADist",   fMuTkSADist);
  chainData5fb.SetBranchAddress("MuTkSAdR",     fMuTkSAdR);
  chainData5fb.SetBranchAddress("MuECALEnergy", fMuECALEnergy);
  chainData5fb.SetBranchAddress("MuHCALEnergy", fMuHCALEnergy);
  chainData5fb.SetBranchAddress("MuCalCompat",  fMuCalCompat);
  chainData5fb.SetBranchAddress("MuIsSoft",     fMuIsSoft);
  chainData5fb.SetBranchAddress("MuIsGlobal",   fMuIsGlobal);
  chainData5fb.SetBranchAddress("MuIsTracker",  fMuIsTracker);
  chainData5fb.SetBranchAddress("MuIsStandalone", fMuIsStandalone);
  chainData5fb.SetBranchAddress("MuIsCalo",     fMuIsCalo);
  chainData5fb.SetBranchAddress("MuArbitrated", fMuArbitrated);
  chainData5fb.SetBranchAddress("MuLastStationLoose", fMuLastStationLoose);
  chainData5fb.SetBranchAddress("MuLastStationTight", fMuLastStationTight);
  chainData5fb.SetBranchAddress("Mu2DCompatibilityLoose", fMu2DCompatibilityLoose);
  chainData5fb.SetBranchAddress("Mu2DCompatibilityTight", fMu2DCompatibilityTight);
  chainData5fb.SetBranchAddress("MuOneStationLoose", fMuOneStationLoose);
  chainData5fb.SetBranchAddress("MuOneStationTight", fMuOneStationTight);
  chainData5fb.SetBranchAddress("MuHLTMatch",   fMuHLTMatch);
  chainData5fb.SetBranchAddress("MuL3Match",    fMuL3Match);
  chainData5fb.SetBranchAddress("MuTightMatch", fMuTightMatch);

  chainData5fb.SetBranchAddress("PcToJt",       fPcToJt);

  chainData5fb.SetBranchAddress("PvN",          &fPvN);
  chainData5fb.SetBranchAddress("RePvN",        &fRePvN);
  chainData5fb.SetBranchAddress("AllPvN",       &fAllPvN);
  chainData5fb.SetBranchAddress("PvIndex",      fPvIndex);
  chainData5fb.SetBranchAddress("PvTkN",        fPvTkN);
  chainData5fb.SetBranchAddress("PvX",          fPvX);
  chainData5fb.SetBranchAddress("PvY",          fPvY);
  chainData5fb.SetBranchAddress("PvZ",          fPvZ);
  chainData5fb.SetBranchAddress("PvXe",         fPvXe);
  chainData5fb.SetBranchAddress("PvYe",         fPvYe);
  chainData5fb.SetBranchAddress("PvZe",         fPvZe);
  chainData5fb.SetBranchAddress("PvPx",         fPvPx);
  chainData5fb.SetBranchAddress("PvPy",         fPvPy);
  chainData5fb.SetBranchAddress("PvPz",         fPvPz);
  chainData5fb.SetBranchAddress("PvPt",         fPvPt);
  chainData5fb.SetBranchAddress("PvEta",        fPvEta);
  chainData5fb.SetBranchAddress("PvChi2",       fPvChi2);
  chainData5fb.SetBranchAddress("PvNdof",       fPvNdof);
  chainData5fb.SetBranchAddress("PvMass",       fPvMass);
  chainData5fb.SetBranchAddress("PvIsFake",     fPvIsFake);
  chainData5fb.SetBranchAddress("PvIsRefit",    fPvIsRefit);

  chainData5fb.SetBranchAddress("JtN",          &fJtN);
  chainData5fb.SetBranchAddress("JtStandN",     &fJtStandN);
  chainData5fb.SetBranchAddress("JtTkN",        fJtTkN);
  chainData5fb.SetBranchAddress("JtSsvN",       fJtSsvN);
  chainData5fb.SetBranchAddress("JtGtvN",       fJtGtvN);
  chainData5fb.SetBranchAddress("JtIndex",      fJtIndex);
  chainData5fb.SetBranchAddress("JtStandIndex", fJtStandIndex);
  chainData5fb.SetBranchAddress("JtToPv",       fJtToPv);
  chainData5fb.SetBranchAddress("JtnConstituents", fJtnConstituents);
  chainData5fb.SetBranchAddress("Jtn60",        fJtn60);
  chainData5fb.SetBranchAddress("Jtn90",        fJtn90);
  chainData5fb.SetBranchAddress("JtnChargedParticles", fJtnChargedParticles);
  chainData5fb.SetBranchAddress("JtnNeutralParticles", fJtnNeutralParticles);
  chainData5fb.SetBranchAddress("JtnChargedHadrons", fJtnChargedHadrons);
  chainData5fb.SetBranchAddress("JtnNeutralHadrons", fJtnNeutralHadrons);
  chainData5fb.SetBranchAddress("JtnPhotons",   fJtnPhotons);
  chainData5fb.SetBranchAddress("JtnElectrons", fJtnElectrons);
  chainData5fb.SetBranchAddress("JtnMuons",     fJtnMuons);
  chainData5fb.SetBranchAddress("JtnHFHadrons", fJtnHFHadrons);
  chainData5fb.SetBranchAddress("JtnHFEMParticles", fJtnHFEMParticles);
  chainData5fb.SetBranchAddress("JtRankTCHE",   fJtRankTCHE);
  chainData5fb.SetBranchAddress("JtRankTCHP",   fJtRankTCHP);
  chainData5fb.SetBranchAddress("JtRankP",      fJtRankP);
  chainData5fb.SetBranchAddress("JtRankBP",     fJtRankBP);
  chainData5fb.SetBranchAddress("JtRankSSVHE",  fJtRankSSVHE);
  chainData5fb.SetBranchAddress("JtRankSSVHP",  fJtRankSSVHP);
  chainData5fb.SetBranchAddress("JtRankCSV",    fJtRankCSV);
  chainData5fb.SetBranchAddress("JtRankCSVMVA", fJtRankCSVMVA);
  chainData5fb.SetBranchAddress("JtRankGT",     fJtRankGT);
  chainData5fb.SetBranchAddress("JtRankSE",     fJtRankSE);
  chainData5fb.SetBranchAddress("JtRankSM",     fJtRankSM);
  chainData5fb.SetBranchAddress("JtCharge",     fJtCharge);
  chainData5fb.SetBranchAddress("JtDiscTCHE",   fJtDiscTCHE);
  chainData5fb.SetBranchAddress("JtDiscTCHP",   fJtDiscTCHP);
  chainData5fb.SetBranchAddress("JtDiscP",      fJtDiscP);
  chainData5fb.SetBranchAddress("JtDiscBP",     fJtDiscBP);
  chainData5fb.SetBranchAddress("JtDiscSSVHE",  fJtDiscSSVHE);
  chainData5fb.SetBranchAddress("JtDiscSSVHP",  fJtDiscSSVHP);
  chainData5fb.SetBranchAddress("JtDiscCSV",    fJtDiscCSV);
  chainData5fb.SetBranchAddress("JtDiscCSVMVA", fJtDiscCSVMVA);
  chainData5fb.SetBranchAddress("JtDiscGT",     fJtDiscGT);
  chainData5fb.SetBranchAddress("JtDiscSE",     fJtDiscSE);
  chainData5fb.SetBranchAddress("JtDiscSM",     fJtDiscSM);
  chainData5fb.SetBranchAddress("JtMaxDist",    fJtMaxDist);
  chainData5fb.SetBranchAddress("JtPhi",        fJtPhi);
  chainData5fb.SetBranchAddress("JtTheta",      fJtTheta);
  chainData5fb.SetBranchAddress("JtEta",        fJtEta);
  chainData5fb.SetBranchAddress("JtRapidity",   fJtRapidity);
  chainData5fb.SetBranchAddress("JtP",          fJtP);
  chainData5fb.SetBranchAddress("JtPt",         fJtPt);
  chainData5fb.SetBranchAddress("JtPx",         fJtPx);
  chainData5fb.SetBranchAddress("JtPy",         fJtPy);
  chainData5fb.SetBranchAddress("JtPz",         fJtPz);
  chainData5fb.SetBranchAddress("JtEnergy",     fJtEnergy);
  chainData5fb.SetBranchAddress("JtEt",         fJtEt);
  chainData5fb.SetBranchAddress("JtMass",       fJtMass);
  chainData5fb.SetBranchAddress("JtMt",         fJtMt);
  chainData5fb.SetBranchAddress("JtVx",         fJtVx);
  chainData5fb.SetBranchAddress("JtVy",         fJtVy);
  chainData5fb.SetBranchAddress("JtVz",         fJtVz);
  chainData5fb.SetBranchAddress("JtChargedEmEnergy", fJtChargedEmEnergy);
  chainData5fb.SetBranchAddress("JtNeutralEmEnergy", fJtNeutralEmEnergy);
  chainData5fb.SetBranchAddress("JtChargedHadronEnergy", fJtChargedHadronEnergy);
  chainData5fb.SetBranchAddress("JtNeutralHadronEnergy", fJtNeutralHadronEnergy);
  chainData5fb.SetBranchAddress("JtPhotonEnergy", fJtPhotonEnergy);
  chainData5fb.SetBranchAddress("JtElectronEnergy", fJtElectronEnergy);
  chainData5fb.SetBranchAddress("JtMuonEnergy", fJtMuonEnergy);
  chainData5fb.SetBranchAddress("JtHFHadronEnergy", fJtHFHadronEnergy);
  chainData5fb.SetBranchAddress("JtHFEMEnergy", fJtHFEMEnergy);
  chainData5fb.SetBranchAddress("JtdRMean",     fJtdRMean);
  chainData5fb.SetBranchAddress("JtdRMax",      fJtdRMax);
  chainData5fb.SetBranchAddress("JtPtRelMean",  fJtPtRelMean);
  chainData5fb.SetBranchAddress("JtPtRelMax",   fJtPtRelMax);
  chainData5fb.SetBranchAddress("JtPtRelSum",   fJtPtRelSum);
  chainData5fb.SetBranchAddress("JtPullPx",     fJtPullPx);
  chainData5fb.SetBranchAddress("JtPullPy",     fJtPullPy);
  chainData5fb.SetBranchAddress("JtPullPz",     fJtPullPz);
  chainData5fb.SetBranchAddress("JtIsStandard", fJtIsStandard);
  chainData5fb.SetBranchAddress("JtIsFat",      fJtIsFat);
  chainData5fb.SetBranchAddress("JtIsSub",      fJtIsSub);
  chainData5fb.SetBranchAddress("JtIsFilt",     fJtIsFilt);
  chainData5fb.SetBranchAddress("JtHLTMatch",   fJtHLTMatch);
  chainData5fb.SetBranchAddress("JtVeto",       fJtVeto);

  chainData5fb.SetBranchAddress("JtHLTN",       fJtHLTN);
  chainData5fb.SetBranchAddress("JtToHLT",      fJtToHLT);
  chainData5fb.SetBranchAddress("JtToPc",       fJtToPc);
  chainData5fb.SetBranchAddress("JtToSsv",      fJtToSsv);
  chainData5fb.SetBranchAddress("JtToGtv",      fJtToGtv);

  chainData5fb.SetBranchAddress("SvN",          &fSvN);
  chainData5fb.SetBranchAddress("SsvN",         &fSsvN);
  chainData5fb.SetBranchAddress("GtvN",         &fGtvN);
  chainData5fb.SetBranchAddress("SvIndex",      fSvIndex);
  chainData5fb.SetBranchAddress("SvTkN",        fSvTkN);
  chainData5fb.SetBranchAddress("SvToJt",       fSvToJt);
  chainData5fb.SetBranchAddress("SvSeqInJt",    fSvSeqInJt);
  chainData5fb.SetBranchAddress("SvX",          fSvX);
  chainData5fb.SetBranchAddress("SvY",          fSvY);
  chainData5fb.SetBranchAddress("SvZ",          fSvZ);
  chainData5fb.SetBranchAddress("SvXe",         fSvXe);
  chainData5fb.SetBranchAddress("SvYe",         fSvYe);
  chainData5fb.SetBranchAddress("SvZe",         fSvZe);
  chainData5fb.SetBranchAddress("SvPx",         fSvPx);
  chainData5fb.SetBranchAddress("SvPy",         fSvPy);
  chainData5fb.SetBranchAddress("SvPz",         fSvPz);
  chainData5fb.SetBranchAddress("SvPt",         fSvPt);
  chainData5fb.SetBranchAddress("SvEta",        fSvEta);
  chainData5fb.SetBranchAddress("SvChi2",       fSvChi2);
  chainData5fb.SetBranchAddress("SvNdof",       fSvNdof);
  chainData5fb.SetBranchAddress("SvDist",       fSvDist);
  chainData5fb.SetBranchAddress("SvDistCM",     fSvDistCM);
  chainData5fb.SetBranchAddress("SvMass",       fSvMass);
  chainData5fb.SetBranchAddress("SvTau",        fSvTau);
  chainData5fb.SetBranchAddress("SvTauCM",      fSvTauCM);
  chainData5fb.SetBranchAddress("SvIsGTV",      fSvIsGTV);

  chainData5fb.SetBranchAddress("SvToPc",       fSvToPc);

  chainData5fb.SetBranchAddress("METN",         &fMETN);
  chainData5fb.SetBranchAddress("METIndex",     fMETIndex);
  chainData5fb.SetBranchAddress("METCharge",    fMETCharge);
  chainData5fb.SetBranchAddress("METPhi",       fMETPhi);
  chainData5fb.SetBranchAddress("METTheta",     fMETTheta);
  chainData5fb.SetBranchAddress("METEta",       fMETEta);
  chainData5fb.SetBranchAddress("METRapidity",  fMETRapidity);
  chainData5fb.SetBranchAddress("METP",         fMETP);
  chainData5fb.SetBranchAddress("METPt",        fMETPt);
  chainData5fb.SetBranchAddress("METPx",        fMETPx);
  chainData5fb.SetBranchAddress("METPy",        fMETPy);
  chainData5fb.SetBranchAddress("METPz",        fMETPz);
  chainData5fb.SetBranchAddress("METEnergy",    fMETEnergy);
  chainData5fb.SetBranchAddress("METEt",        fMETEt);
  chainData5fb.SetBranchAddress("METMass",      fMETMass);
  chainData5fb.SetBranchAddress("METMt",        fMETMt);
  chainData5fb.SetBranchAddress("METVx",        fMETVx);
  chainData5fb.SetBranchAddress("METVy",        fMETVy);
  chainData5fb.SetBranchAddress("METVz",        fMETVz);

  chainData5fb.SetBranchAddress("JPsiN",           &fJPsiN);
  chainData5fb.SetBranchAddress("JPsiMuMuN",       &fJPsiMuMuN);
  chainData5fb.SetBranchAddress("JPsiMuTkN",       &fJPsiMuTkN);
  chainData5fb.SetBranchAddress("BaseJPsiI",       fBaseJPsiI);
  chainData5fb.SetBranchAddress("JPsiIndex",       fJPsiIndex);
  chainData5fb.SetBranchAddress("JPsiClosestPVinZ",fJPsiClosestPVinZ);
  chainData5fb.SetBranchAddress("JPsiCharge",      fJPsiCharge);
  chainData5fb.SetBranchAddress("JPsiPhi",         fJPsiPhi);
  chainData5fb.SetBranchAddress("JPsiTheta",       fJPsiTheta);
  chainData5fb.SetBranchAddress("JPsiEta",         fJPsiEta);
  chainData5fb.SetBranchAddress("JPsiRapidity",    fJPsiRapidity);
  chainData5fb.SetBranchAddress("JPsiP",           fJPsiP);
  chainData5fb.SetBranchAddress("JPsiPt",          fJPsiPt);
  chainData5fb.SetBranchAddress("JPsiPx",          fJPsiPx);
  chainData5fb.SetBranchAddress("JPsiPy",          fJPsiPy);
  chainData5fb.SetBranchAddress("JPsiPz",          fJPsiPz);
  chainData5fb.SetBranchAddress("JPsiEnergy",      fJPsiEnergy);
  chainData5fb.SetBranchAddress("JPsiEt",          fJPsiEt);
  chainData5fb.SetBranchAddress("JPsiMass",        fJPsiMass);
  chainData5fb.SetBranchAddress("JPsiMt",          fJPsiMt);
  chainData5fb.SetBranchAddress("JPsiChi2",        fJPsiChi2);
  chainData5fb.SetBranchAddress("JPsiNdof",        fJPsiNdof);
  chainData5fb.SetBranchAddress("JPsiVx",          fJPsiVx);
  chainData5fb.SetBranchAddress("JPsiVy",          fJPsiVy);
  chainData5fb.SetBranchAddress("JPsiVz",          fJPsiVz);
  chainData5fb.SetBranchAddress("JPsiVxE",         fJPsiVxE);
  chainData5fb.SetBranchAddress("JPsiVyE",         fJPsiVyE);
  chainData5fb.SetBranchAddress("JPsiVzE",         fJPsiVzE);
  chainData5fb.SetBranchAddress("JPsiVtxPhi",      fJPsiVtxPhi);
  chainData5fb.SetBranchAddress("JPsiVtxTheta",    fJPsiVtxTheta);
  chainData5fb.SetBranchAddress("JPsiVtxEta",      fJPsiVtxEta);
  chainData5fb.SetBranchAddress("JPsiVtxRapidity", fJPsiVtxRapidity);
  chainData5fb.SetBranchAddress("JPsiVtxP",        fJPsiVtxP);
  chainData5fb.SetBranchAddress("JPsiVtxPt",       fJPsiVtxPt);
  chainData5fb.SetBranchAddress("JPsiVtxPx",       fJPsiVtxPx);
  chainData5fb.SetBranchAddress("JPsiVtxPy",       fJPsiVtxPy);
  chainData5fb.SetBranchAddress("JPsiVtxPz",       fJPsiVtxPz);
  chainData5fb.SetBranchAddress("JPsiVtxEnergy",   fJPsiVtxEnergy);
  chainData5fb.SetBranchAddress("JPsiVtxEt",       fJPsiVtxEt);
  chainData5fb.SetBranchAddress("JPsiVtxMass",     fJPsiVtxMass);
  chainData5fb.SetBranchAddress("JPsiVtxMt",       fJPsiVtxMt);
  chainData5fb.SetBranchAddress("JPsiMuI",         fJPsiMuI);
  chainData5fb.SetBranchAddress("JPsiMuCategory",  fJPsiMuCategory);
  chainData5fb.SetBranchAddress("JPsiMuCutKin",    fJPsiMuCutKin);
  chainData5fb.SetBranchAddress("JPsiMuCutHLT",    fJPsiMuCutHLT);
  chainData5fb.SetBranchAddress("JPsiMuCutIso",    fJPsiMuCutIso);
  chainData5fb.SetBranchAddress("JPsiMuCutSA",     fJPsiMuCutSA);
  chainData5fb.SetBranchAddress("JPsiMuCutTrk",    fJPsiMuCutTrk);
  chainData5fb.SetBranchAddress("JPsiMuType",      fJPsiMuType);

  chainData5fb.SetBranchAddress("EtabN",        &fEtabN);
  chainData5fb.SetBranchAddress("EtabIndex",    fEtabIndex);
  chainData5fb.SetBranchAddress("EtabDuplicatesI", fEtabDuplicatesI);
  chainData5fb.SetBranchAddress("EtabCharge",   fEtabCharge);
  chainData5fb.SetBranchAddress("EtabPhi",      fEtabPhi);
  chainData5fb.SetBranchAddress("EtabTheta",    fEtabTheta);
  chainData5fb.SetBranchAddress("EtabEta",      fEtabEta);
  chainData5fb.SetBranchAddress("EtabRapidity", fEtabRapidity);
  chainData5fb.SetBranchAddress("EtabP",        fEtabP);
  chainData5fb.SetBranchAddress("EtabPt",       fEtabPt);
  chainData5fb.SetBranchAddress("EtabPx",       fEtabPx);
  chainData5fb.SetBranchAddress("EtabPy",       fEtabPy);
  chainData5fb.SetBranchAddress("EtabPz",       fEtabPz);
  chainData5fb.SetBranchAddress("EtabEnergy",   fEtabEnergy);
  chainData5fb.SetBranchAddress("EtabEt",       fEtabEt);
  chainData5fb.SetBranchAddress("EtabMass",     fEtabMass);
  chainData5fb.SetBranchAddress("EtabMt",       fEtabMt);
  chainData5fb.SetBranchAddress("EtabChi2",     fEtabChi2);
  chainData5fb.SetBranchAddress("EtabNdof",     fEtabNdof);
  chainData5fb.SetBranchAddress("EtabVx",       fEtabVx);
  chainData5fb.SetBranchAddress("EtabVy",       fEtabVy);
  chainData5fb.SetBranchAddress("EtabVz",       fEtabVz);
  chainData5fb.SetBranchAddress("EtabVxE",      fEtabVxE);
  chainData5fb.SetBranchAddress("EtabVyE",      fEtabVyE);
  chainData5fb.SetBranchAddress("EtabVzE",      fEtabVzE);
  chainData5fb.SetBranchAddress("EtabVtxPhi",   fEtabVtxPhi);
  chainData5fb.SetBranchAddress("EtabVtxTheta", fEtabVtxTheta);
  chainData5fb.SetBranchAddress("EtabVtxEta",   fEtabVtxEta);
  chainData5fb.SetBranchAddress("EtabVtxRapidity",fEtabVtxRapidity);
  chainData5fb.SetBranchAddress("EtabVtxP",     fEtabVtxP);
  chainData5fb.SetBranchAddress("EtabVtxPt",    fEtabVtxPt);
  chainData5fb.SetBranchAddress("EtabVtxPx",    fEtabVtxPx);
  chainData5fb.SetBranchAddress("EtabVtxPy",    fEtabVtxPy);
  chainData5fb.SetBranchAddress("EtabVtxPz",    fEtabVtxPz);
  chainData5fb.SetBranchAddress("EtabVtxEnergy",fEtabVtxEnergy);
  chainData5fb.SetBranchAddress("EtabVtxEt",    fEtabVtxEt);
  chainData5fb.SetBranchAddress("EtabVtxMass",  fEtabVtxMass);
  chainData5fb.SetBranchAddress("EtabVtxMt",    fEtabVtxMt);

  chainData5fb.SetBranchAddress("BaseEtabI",    &fBaseEtabI);
  chainData5fb.SetBranchAddress("EtabJPsiI",    fEtabJPsiI);
  chainData5fb.SetBranchAddress("EtabMuI",      fEtabMuI);
  chainData5fb.SetBranchAddress("EtabToRePvI",  fEtabToRePvI);
  chainData5fb.SetBranchAddress("EtabMuN",      fEtabMuN);
  chainData5fb.SetBranchAddress("EtabCT",       fEtabCT);
  chainData5fb.SetBranchAddress("EtabCTxy",     fEtabCTxy);
  chainData5fb.SetBranchAddress("EtabVtxCT",    fEtabVtxCT);
  chainData5fb.SetBranchAddress("EtabVtxCTxy",  fEtabVtxCTxy);
  chainData5fb.SetBranchAddress("EtabJPsiDeltaL", fEtabJPsiDeltaL);
  chainData5fb.SetBranchAddress("EtabJPsiDeltaT", fEtabJPsiDeltaT);
  chainData5fb.SetBranchAddress("EtabJPsiVtxErr", fEtabJPsiVtxErr);
  chainData5fb.SetBranchAddress("EtabJPsiVtxErrxy", fEtabJPsiVtxErrxy);
  chainData5fb.SetBranchAddress("EtabJPsiProjX", fEtabJPsiProjX);
  chainData5fb.SetBranchAddress("EtabJPsiProjY", fEtabJPsiProjY);
  chainData5fb.SetBranchAddress("EtabJPsiProjZ", fEtabJPsiProjZ);
  chainData5fb.SetBranchAddress("EtabJPsiCT",   fEtabJPsiCT);
  chainData5fb.SetBranchAddress("EtabJPsiCTxy", fEtabJPsiCTxy);
  chainData5fb.SetBranchAddress("EtabJPsiVtxCT",   fEtabJPsiVtxCT);
  chainData5fb.SetBranchAddress("EtabJPsiVtxCTxy", fEtabJPsiVtxCTxy);
  chainData5fb.SetBranchAddress("EtabJPsiToPVVtxErr",fEtabJPsiToPVVtxErr);
  chainData5fb.SetBranchAddress("EtabJPsiToPVVtxErrxy",fEtabJPsiToPVVtxErrxy);
  chainData5fb.SetBranchAddress("EtabJPsiIsoTkN",fEtabJPsiIsoTkN);
  chainData5fb.SetBranchAddress("EtabJPsiIso7PV",fEtabJPsiIso7PV);
  chainData5fb.SetBranchAddress("EtabJPsiIsoTkCA",fEtabJPsiIsoTkCA);

  chainData5fb.SetBranchAddress("HN",           &fHN);
  chainData5fb.SetBranchAddress("HIndex",       fHIndex);
  chainData5fb.SetBranchAddress("HCharge",      fHCharge);
  chainData5fb.SetBranchAddress("HPhi",         fHPhi);
  chainData5fb.SetBranchAddress("HTheta",       fHTheta);
  chainData5fb.SetBranchAddress("HEta",         fHEta);
  chainData5fb.SetBranchAddress("HRapidity",    fHRapidity);
  chainData5fb.SetBranchAddress("HP",           fHP);
  chainData5fb.SetBranchAddress("HPt",          fHPt);
  chainData5fb.SetBranchAddress("HPx",          fHPx);
  chainData5fb.SetBranchAddress("HPy",          fHPy);
  chainData5fb.SetBranchAddress("HPz",          fHPz);
  chainData5fb.SetBranchAddress("HEnergy",      fHEnergy);
  chainData5fb.SetBranchAddress("HEt",          fHEt);
  chainData5fb.SetBranchAddress("HMass",        fHMass);
  chainData5fb.SetBranchAddress("HMt",          fHMt);
  chainData5fb.SetBranchAddress("HVx",          fHVx);
  chainData5fb.SetBranchAddress("HVy",          fHVy);
  chainData5fb.SetBranchAddress("HVz",          fHVz);

  chainData5fb.SetBranchAddress("HJtI",         fHJtI);

std::cout << "Branches have been set" << std::endl;

  // efficiency cut variables
 Int_t cutHLT, cutHLTUps,cutHLTZ, cut4MuN, cutSoftMuN, softDiMu, cutJPsiVtxMass, cutJPsiVtxProb, cutDiJPsi, cutEtabN, cutEtabHLTN, cutEtabProb, nFill;
  cutHLT = cutHLTUps = cut4MuN = cutSoftMuN = softDiMu = cutJPsiVtxMass = cutJPsiVtxProb = cutDiJPsi = cutEtabN = cutEtabHLTN = cutEtabProb = nFill = 0;

  TH1F *DimuCandidateMass= new TH1F("DimuCandidateMass","DimuCandidateMass",10,0,10);
  
  Float_t aa,bb,cc,a1,b1,c1;//defining Etabmass,JPsiMassI,JPsiMassII,EtabP,JPsiMassIP,JPsiMassIIP    
  TNtuple *nhimal = new TNtuple("ntuple","Ntuple analysis","etab:j1:j2:etaP:j1P:j2P");
  
  //-------------------------------------------------------------
  //loop over events and get entries
  Int_t nevent = chainData5fb.GetEntries();
  
  for (Int_t i=0; i<nevent; i++) { 
    if( i % 1000 == 0 ) cout<<i+1<<" of "<<nevent<<endl;
    //    if (i == 100000) break;
    chainData5fb.GetEvent(i);    //read complete accepted event in memory

    fJPsiBestProbI[0]=fJPsiBestProbI[1]=fJPsiVtxBestProbI[0]=fJPsiVtxBestProbI[1]=fMuBestProbI[3]=fMuBestProbI[2]=fMuBestProbI[1]=fMuBestProbI[0]=fTkBestProbI=fMuBestProbByPtI[3]=fMuBestProbByPtI[2]=fMuBestProbByPtI[1]=fMuBestProbByPtI[0]=-9999;

    //    int DiMuCandCounter = 0;
    //    int UpsCandCounter = 0;
    //   int EtabCandCounter = 0;
    /////// HLT
    if(fHLTP_Dimuon0_Jpsi3p5_Muon2[1]) cutHLT++;
    else if (fHLTP_Trimuon5_3p5_2_Upsilon_Muon[1]) cutHLTUps++;
    else if (fHLTP_IsoMu27[1]) cutHLTZ++;
    //else continue;//Himal change for removing trigger cut


    ////// muon cuts
    Int_t nmusoft, ntkplussoft, ntkminussoft;
    nmusoft=0; ntkplussoft=0; ntkminussoft=0;
    if (fMuN >= 4) cut4MuN++;
    else continue;
    for (Int_t imu=0; imu<fMuN; imu++) {
      fPcTPFilter[imu] = fPcBasicFilter[imu] = false;
      /// new soft cut
      if (fMuIsSoft[imu]){
        nmusoft++;
       if( fPcCharge[imu]>0 ) ntkplussoft++;
       else if( fPcCharge[imu]<0 ) ntkminussoft++;
       fPcBasicFilter[imu] = true;
        }
    } // end loop over muons

    if( nmusoft>=4 ) {
      cutSoftMuN++;
      if( ntkplussoft>=2 && ntkminussoft>=2 ) softDiMu++;
      else continue;
    }
    else continue;
    


    // J/Psi cuts
    Int_t njpsimuptetamumu=0, njpsimassmumu=0, njpsiptmumu=0, njpsiprobmumu=0, njpsidispmumu=0, njpsikinmumu=0, njpsimassmutk=0, njpsiptmutk=0, njpsiprobmutk=0, njpsivtxmass=0;
    float bestMuProb1 = 0.0, bestMuProb2 = 0.0, bestMuProb3 = 0.0, bestMuProb4 = 0.0, bestTkProb = 0.0;
    float bestPt = 0.0, secondPt = 0.0;	
    // J/Psi cuts
    Int_t njpsivtxptmumu=0, njpsivtxprobmumu=0, njpsivtxdispmumu=0, njpsivtxkinmumu=0, njpsimuvtxptetamumu=0, njpsivtxptmutk=0, njpsivtxprobmutk=0;


    // loop over jpsi candidates - find JPsiN (does it count Jpsi without Vtx fit?)
    // plot Fill(fJPsiMuMuN) == multiplicity   THF1 histogram ....

    // loop over J/Psi using vtx kinematics
    for (Int_t ijp=0; ijp<fJPsiMuMuN; ijp++) {


      // plot JPsi Mass
      DimuCandidateMass->Fill(fJPsiVtxMass[ijp]);

      fJPsiBasicFilter[ijp] = fJPsiVtxBasicFilter[ijp] = false;

      // make sure J/Psi tracks exist
      if( fJPsiMuI[ijp][0]==-9999 || fJPsiMuI[ijp][1]==-9999 )  continue;


      // if both tracks pass looser mu criteria
     if( fPcBasicFilter[fJPsiMuI[ijp][0]] && fPcBasicFilter[fJPsiMuI[ijp][1]] ) {
       // J/Psi cuts
       if( fJPsiChi2[ijp]!=-9999 && TMath::Prob(fJPsiChi2[ijp],fJPsiNdof[ijp])>0.005 ) {
         njpsivtxprobmumu++;
         if( fPcPt[fJPsiMuI[ijp][0]]>=2.0&&fPcPt[fJPsiMuI[ijp][1]]>=2.0&&abs(fPcEta[fJPsiMuI[ijp][0]])<=2.4&&abs(fPcEta[fJPsiMuI[ijp][1]])<=2.4){
           njpsimuptetamumu++;
           njpsidispmumu++;
           fJPsiBasicFilter[ijp] = true;}
	 if((fJPsiVtxMass[ijp]>2.8&&fJPsiVtxMass[ijp]<4.0)||(fJPsiVtxMass[ijp]>8.5&&fJPsiVtxMass[ijp]<12)||(fJPsiVtxMass[ijp]>70&&fJPsiVtxMass[ijp]<110) ){
           ++njpsivtxmass; // add valid JPsi candidate
           }
         }
       }
     } // end loop over J/Psi

    if( njpsivtxprobmumu > 0 ) {
      cutJPsiVtxProb++;
      if( njpsivtxmass > 0 ) {
        cutJPsiVtxMass++;
        if( njpsivtxprobmumu >= 2 && njpsivtxmass >= 2 ) cutDiJPsi++;
        else continue;
      }
      else continue;
    }
    else continue;

    // Eta_b cuts 
    Int_t netab=0, netabhlt=0, netabsignif=0, netabprob=0, netabmass=0, netabmult=0, netabvtx=0, netabvtxhlt=0, netabvtxsignif=0;
    fEtabBestMassI = fEtabVtxBestMassI = fEtabBest4MuProbI = fEtabBestProbI = fEtabVtxBestProbI = fEtabVtxBestPtI = -9999;
    float bestMass = 9999.0, bestVtxMass = 9999.0, bestProb = 0.0, best4muProb = 0.0, bestJP1Prob = 0.0, bestJP2Prob = 0.0, bestJP1VtxProb = 0.0, bestJP2VtxProb = 0.0, biggestPt = 0.0;

    /////////////////////////////////////////////////////////
    //loop for matching of Etab muon and HLT muon         //
    //Himal Acharya                                       //
    //--To determine the appropriate position for the cut //                              
    ////////////////////////////////////////////////////////

    /* 
   //end 
   */


   ///////////////////////////////////////
   //Himal code for the selected Muon////
   //////////////////////////////////////
   
   /*
   //end
    */


    // loop over Eta_b
    for (Int_t ie=0; ie<fEtabN; ie++) {

      fEtabBasicFilter[ie]=false;
      fEtabVtxBasicFilter[ie]=false;

      fEtabL3MatchMuN[ie] = 0;
      for (Int_t imu=0; imu<fEtabMuN[ie]; imu++) {
        if( fEtabMuI[ie][imu]!=-9999 && fPcBasicFilter[fEtabMuI[ie][imu]] && fMuL3Match[fEtabMuI[ie][imu]] ) fEtabL3MatchMuN[ie] = fEtabL3MatchMuN[ie]+1;
      }

      if( fEtabMuN[ie]<4 ) continue;
      // make sure Eta_b tracks match up

      if( fEtabMuI[ie][0]==-9999 || fEtabMuI[ie][1]==-9999 || fEtabMuI[ie][2]==-9999 || fEtabMuI[ie][3]==-9999 ) continue;
      // make sure Eta_b J/Psi match up

      if( fEtabJPsiI[ie][0]==-9999 || fEtabJPsiI[ie][1]==-9999 ) continue;
      // Eta_b cuts
      // can we make 4-muon from looser cuts?

      if( fJPsiBasicFilter[fEtabJPsiI[ie][0]] && fJPsiBasicFilter[fEtabJPsiI[ie][1]] ) {
        netab++;

	//	if( fEtabL3MatchMuN[ie] >=3 ){
	  netabhlt++;
	
	  if( fEtabChi2[ie]!=-9999 && TMath::Prob(fEtabChi2[ie],fEtabNdof[ie])>0.005 ){
	    netabprob++;


	    fEtabBasicFilter[ie]=true;
         

         
            // record best jpsi and eta_b candidates
            if( TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][0]],fJPsiNdof[fEtabJPsiI[ie][0]])>bestJP1Prob ) { fJPsiBestProbI[0] = fEtabJPsiI[ie][0]; bestJP1Prob = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][0]],fJPsiNdof[fEtabJPsiI[ie][0]]); }
            if( TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][1]],fJPsiNdof[fEtabJPsiI[ie][1]])>bestJP1Prob ) { fJPsiBestProbI[0] = fEtabJPsiI[ie][1]; bestJP1Prob = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][1]],fJPsiNdof[fEtabJPsiI[ie][1]]); }
          }
	  //        }
      } 



    } // end loop over Eta_b
//bracket for himal loop
	
    // loop again for best prob
    for (Int_t ie=0; ie<fEtabN; ie++) {

      if( fEtabBasicFilter[ie] && fJPsiBestProbI[0]==fEtabJPsiI[ie][0] && TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][1]],fJPsiNdof[fEtabJPsiI[ie][1]])>bestJP2Prob ) { fEtabBestProbI = ie; fJPsiBestProbI[1] = fEtabJPsiI[ie][1]; bestJP2Prob = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][1]],fJPsiNdof[fEtabJPsiI[ie][1]]); }
      else if( fEtabBasicFilter[ie] && fJPsiBestProbI[0]==fEtabJPsiI[ie][1] && TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][0]],fJPsiNdof[fEtabJPsiI[ie][0]])>bestJP2Prob ) { fEtabBestProbI = ie; fJPsiBestProbI[1] = fEtabJPsiI[ie][0]; bestJP2Prob = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][0]],fJPsiNdof[fEtabJPsiI[ie][0]]); }
      if( fEtabVtxBasicFilter[ie] && fJPsiVtxBestProbI[0]==fEtabJPsiI[ie][0] && TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][1]],fJPsiNdof[fEtabJPsiI[ie][1]])>bestJP2VtxProb ) { fEtabVtxBestProbI = ie; fJPsiVtxBestProbI[1] = fEtabJPsiI[ie][1]; bestJP2VtxProb = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][1]],fJPsiNdof[fEtabJPsiI[ie][1]]); }
      else if( fEtabVtxBasicFilter[ie] && fJPsiVtxBestProbI[0]==fEtabJPsiI[ie][1] && TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][0]],fJPsiNdof[fEtabJPsiI[ie][0]])>bestJP2VtxProb ) { fEtabVtxBestProbI = ie; fJPsiVtxBestProbI[1] = fEtabJPsiI[ie][0]; bestJP2VtxProb = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][0]],fJPsiNdof[fEtabJPsiI[ie][0]]); }
      if( fEtabBasicFilter[ie] && sqrt(pow(fJPsiVtxMass[fEtabJPsiI[ie][0]]-3.096,2)+pow(fJPsiVtxMass[fEtabJPsiI[ie][1]]-3.096,2))<bestMass ) { fEtabBestMassI = ie; bestMass = sqrt(pow(fJPsiVtxMass[fEtabJPsiI[ie][0]]-3.096,2)+pow(fJPsiVtxMass[fEtabJPsiI[ie][1]]-3.096,2)); }
      if( fEtabVtxBasicFilter[ie] && sqrt(pow(fJPsiVtxMass[fEtabJPsiI[ie][0]]-3.096,2)+pow(fJPsiVtxMass[fEtabJPsiI[ie][1]]-3.096,2))<bestVtxMass ) { fEtabVtxBestMassI = ie; bestVtxMass = sqrt(pow(fJPsiVtxMass[fEtabJPsiI[ie][0]]-3.096,2)+pow(fJPsiVtxMass[fEtabJPsiI[ie][1]]-3.096,2)); }
      //Grant additions to pick etab on pT
      if( fEtabVtxBasicFilter[ie] && (fJPsiVtxPt[fEtabJPsiI[ie][0]])+(fJPsiVtxPt[fEtabJPsiI[ie][1]])>biggestPt) { fEtabVtxBestPtI = ie;biggestPt = (fJPsiVtxPt[fEtabJPsiI[ie][0]])+(fJPsiVtxPt[fEtabJPsiI[ie][1]]); /*cout<<(fJPsiVtxPt[fEtabJPsiI[ie][0]])<<" + "<<(fJPsiVtxPt[fEtabJPsiI[ie][1]])<<" = "<<biggestPt<<endl;*/}
      if( fEtabBasicFilter[ie] && (TMath::Prob(fEtabChi2[ie],fEtabNdof[ie])>best4muProb)) {fEtabBest4MuProbI = ie; best4muProb = TMath::Prob(fEtabChi2[ie],fEtabNdof[ie]);}

    } // end loop over Eta_b




    if(fEtabBest4MuProbI!=-9999){
      int ie123 = fEtabBest4MuProbI;
      
      aa = fEtabVtxMass[ie123];//EtabMass
      bb = fJPsiVtxMass[fEtabJPsiI[ie123][0]];//JPsiIMass
      cc = fJPsiVtxMass[fEtabJPsiI[ie123][1]];//JpsiIIMass
      a1 = best4muProb;//EtabP
      b1 = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie123][0]],fJPsiNdof[fEtabJPsiI[ie123][0]]);//JPsiIP
      c1 = TMath::Prob(fJPsiChi2[fEtabJPsiI[ie123][1]],fJPsiNdof[fEtabJPsiI[ie123][1]]);//JPsiIIP


      nhimal->Fill(aa,bb,cc,a1,b1,c1);

    }



    // sort best muons by pT
    if( fEtabBestProbI!=-9999 ) {

      int tempI;
      fMuBestProbByPtI[0]=fEtabMuI[fEtabVtxBestPtI][0];
      fMuBestProbByPtI[1]=fEtabMuI[fEtabVtxBestPtI][1];
      fMuBestProbByPtI[2]=fEtabMuI[fEtabVtxBestPtI][2];
      fMuBestProbByPtI[3]=fEtabMuI[fEtabVtxBestPtI][3];
      for( int im = 0; im < 4; im++ ) {
        if( fMuBestProbByPtI[im]==-9999 ) continue;
        for( int im2 = im+1; im2 < 4; im2++ ) {
          if( fMuBestProbByPtI[im2]==-9999 ) continue;
          if( fPcPt[fMuBestProbByPtI[im]]<fPcPt[fMuBestProbByPtI[im2]] ) {
            tempI=fMuBestProbByPtI[im]; fMuBestProbByPtI[im]=fMuBestProbByPtI[im2]; fMuBestProbByPtI[im2]=tempI;
          }
        }
      }
    }


    if( netab>0 ) {
      cutEtabN++;
      if( netabhlt>0 ) {
        cutEtabHLTN++;
        if( netabprob>0) cutEtabProb++;
         else continue;
      }
     else continue;
    }
    else continue;
    fTree->Fill();
    nFill++;
    
    
  } // end loop over events

  cout << "Final Counter Output" << endl; 
  cout << " HLT Muon with fHLTP_QuadMuon0_Dimuon0_Jpsi or HLTP_Dimuon0_Jpsi3p5_Muon2 Trigger fired are : " << cutHLT <<endl;
  cout<< " Event with MuN greater than or equal to 4 are :  " << cut4MuN <<endl;
  cout<< " Event with nmusoft greater than or equal to 4 are : " << cutSoftMuN  <<endl;
  cout<< " Event with ntkplussoft and ntkinusoft muons both greater than or equal to 2 are :" << softDiMu << endl;
  cout<< " Event with JPsi Vertex Probability greater than 0.005  are : " << cutJPsiVtxProb << endl;
  cout<<" Event with JPsi Mass greater than 2.8 and less than 3.35 are :  "  << cutJPsiVtxMass << endl;
  cout<<" Event having JPsi Vertex Probability and JPsi Vertex Mass both  for 2 or more than 2 muon : " << cutDiJPsi << endl;
  cout<<" Etab candidate (Event in which both EtabJPsiCandidates passes Basic Filter are : " << cutEtabN << endl;
  cout<<" Etab candidate in which  both EtabJPsiCandidates passes JPsiBasic Filter are : " << cutEtabHLTN << endl;
  cout<<" Etab candidate with vertex probability greater than 0.005 are : "<< cutEtabProb << endl;
  cout<<" Number of Etab event filled in Tree are : " << nFill << endl;

  fFile->Write();
  fFile->Close();
  delete fFile;
} // end macro
