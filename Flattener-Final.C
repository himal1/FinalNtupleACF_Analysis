#include <string>
#include <stdlib.h>
#include <cmath>
#include "Riostream.h"

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

//Defining the helicity function
Float_t coshel(TLorentzVector particle, TLorentzVector parent,
	       TLorentzVector grandparent) {

  TVector3 boosttoparent = -(parent.BoostVector());

  particle.Boost(boosttoparent);
  grandparent.Boost(boosttoparent);

  TVector3 particle3 = particle.Vect();
  TVector3 grandparent3 = grandparent.Vect();
  Float_t numerator = particle3.Dot(grandparent3);
  Float_t denominator = (particle3.Mag())*(grandparent3.Mag());
  Float_t temp = numerator/denominator;

  return temp;

}



int main() {
  gROOT->Reset();
  using namespace std;
  TChain chainData5fb("PATEventTree");
  TChain chainSigMC("PATEventTree");
  chainData5fb.Add("SingleMuon2017F_Soft4Mu_Etab_JPsiUpsilonZ_NoTrigger_CompactTest.root");

  // make some histogram
  TH1F *h1 = new TH1F("h1", "h1", 11, 0, 10);
  //end

  // create new file
  TFile *fFile = new TFile("SingleMuon2017F_NoTrigger_FlattenH125ZJpsi.root","recreate");
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
  bool  fHLTP_DoubleMu3[3], fHLTP_DoubleMu6[3], fHLTP_DoubleMu7[3], fHLTP_Dimuon0_Jpsi3p5_Muon2[3], fHLTP_Dimuon0_Jpsi[3], fHLTP_Dimuon10_Jpsi_Barrel[3], fHLTP_TripleMu5[3], fHLTP_QuadMuon0_Dimuon0_Jpsi[3],fHLTP_IsoMu20[3],fHLTP_IsoMu27[3],fHLTP_Mu50[3]; //initially HLTP_Dimuon0_Jpsi_Muon changed for 2017 by himal
	//fHLTP_QuadMuon0_Dimuon0_Upsilon[3];
  // HLT prescalersb
  unsigned int fHLTP_DoubleMu3_PS, fHLTP_DoubleMu6_PS, fHLTP_DoubleMu7_PS, fHLTP_Dimuon0_Jpsi_PS, fHLTP_Dimuon10_Jpsi_Barrel_PS, fHLTP_TripleMu5_PS;
  // HLT filters passed
  bool  fHLTP_DoubleMu3_Filters[5], fHLTP_DoubleMu6_Filters[5], fHLTP_DoubleMu7_Filters[5], fHLTP_TripleMu5_Filters[5], fHLTP_Dimuon0_Jpsi_Filters[7], fHLTP_Dimuon10_Jpsi_Barrel_Filters[7];
  //himal add
  int iMu0, iMu1;
  //end
  // -- HLT objects from PAT
  const int HLTMAX = 1000;
  int   fHLTN;
  int   fHLT_Index[HLTMAX], fHLT_ToPc[HLTMAX], fHLT_ToJt[HLTMAX], fHLT_PdgId[HLTMAX];
  float fHLT_Mass[HLTMAX], fHLT_Energy[HLTMAX], fHLT_Et[HLTMAX], fHLT_P[HLTMAX], fHLT_Pt[HLTMAX], fHLT_Px[HLTMAX], fHLT_Py[HLTMAX], fHLT_Pz[HLTMAX], fHLT_Theta[HLTMAX], fHLT_Eta[HLTMAX], fHLT_Phi[HLTMAX];
  bool  fHLT_Mu[HLTMAX][2], fHLT_Mu12[HLTMAX][2], fHLT_Mu15[HLTMAX][2], fHLT_Mu20[HLTMAX][2], fHLT_Mu24[HLTMAX][2], fHLT_Mu30[HLTMAX][2], fHLT_IsoMu12[HLTMAX][2], fHLT_IsoMu15[HLTMAX][2], fHLT_IsoMu17[HLTMAX][2], fHLT_IsoMu24[HLTMAX][2], fHLT_IsoMu30[HLTMAX][2], fHLT_DoubleMu3[HLTMAX][2], fHLT_DoubleMu6[HLTMAX][2], fHLT_DoubleMu7[HLTMAX][2], fHLT_Dimuon0_Jpsi3p5_Muon2[HLTMAX][2], fHLT_Dimuon0_Jpsi[HLTMAX][2], fHLT_Dimuon7_Jpsi_Displaced[HLTMAX][2], fHLT_Dimuon7_Jpsi_X_Barrel[HLTMAX][2], fHLT_Dimuon10_Jpsi_Barrel[HLTMAX][2], fHLT_TripleMu5[HLTMAX][2], fHLT_Jet[HLTMAX][2]; //initially fHLT_Dimuon0_Jpsi_Muon changed for 2017 by himal

  // -- Particles
  const int PARTMAX = 10000;
  int   fPcN, fTkN, fMuN, fElecN, fMiscTkN, fPhotN; //fHadrN, 
  int   fPcIndex[PARTMAX], fMuIndex[PARTMAX], fElecIndex[PARTMAX], fMiscTkIndex[PARTMAX], fPhotIndex[PARTMAX], fPcToGn[PARTMAX], fPcToTk[PARTMAX], fTkToPc[PARTMAX], fPcToPv[PARTMAX], fPcTkQuality[PARTMAX], fPcJtN[PARTMAX], fPcPdgId[PARTMAX], fPcPixHitN[PARTMAX], fPcPixLayN[PARTMAX], fPcStripHitN[PARTMAX], fPcStripLayN[PARTMAX]; //fHadrIndex[PARTMAX]
  float fPcCharge[PARTMAX], fPcChi2[PARTMAX], fPcNdof[PARTMAX], fPcEnergy[PARTMAX], fPcEt[PARTMAX], fPcP[PARTMAX], fPcPt[PARTMAX], fPcPx[PARTMAX], fPcPy[PARTMAX], fPcPz[PARTMAX], fPcTheta[PARTMAX], fPcEta[PARTMAX], fPcPhi[PARTMAX], fPcD0[PARTMAX], fPcDz[PARTMAX], fPcEtaErr[PARTMAX], fPcPhiErr[PARTMAX], fPcD0Err[PARTMAX], fPcDzErr[PARTMAX], fPcVx[PARTMAX], fPcVy[PARTMAX], fPcVz[PARTMAX], fPcEcalIso[PARTMAX], fPcHcalIso[PARTMAX], fPcTrackIso[PARTMAX], fPcIP[PARTMAX], fPcIPxy[PARTMAX];
  // -- muons
  int   fMuHitN[PARTMAX], fMuMatchedN[PARTMAX], fMuHLTN[PARTMAX], fMuToHLT[PARTMAX], fMuBestProbI[4], fTkBestProbI, fMuBestProbByPtI[4];
  float fMuChi2[PARTMAX], fMuNdof[PARTMAX], fMuTkKink[PARTMAX], fMuGlbKink[PARTMAX], fMuGlbProb[PARTMAX], fMuTkSADist[PARTMAX], fMuTkSAdR[PARTMAX], fMuECALEnergy[PARTMAX], fMuHCALEnergy[PARTMAX], fMuCalCompat[PARTMAX];
  bool  fMuIsGlobal[PARTMAX], fMuIsTracker[PARTMAX], fMuIsStandalone[PARTMAX], fMuIsCalo[PARTMAX], fMuArbitrated[PARTMAX], fMuLastStationLoose[PARTMAX], fMuLastStationTight[PARTMAX], fMu2DCompatibilityLoose[PARTMAX], fMu2DCompatibilityTight[PARTMAX], fMuOneStationLoose[PARTMAX], fMuOneStationTight[PARTMAX], fMuHLTMatch[PARTMAX][2], fMuL3Match[PARTMAX], fMuTightMatch[PARTMAX], fPcTPFilter[PARTMAX], fPcBasicFilter[PARTMAX];

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
  const int JPsiMAX = 500;

  int fJPsiN, fJPsiMuMuN, fJPsiMuTkN, fBaseJPsiI[2], fJPsiBestProbI[2], fJPsiVtxBestProbI[2];
  int fJPsiIndex[JPsiMAX], fJPsiClosestPVinZ[JPsiMAX], fJPsiMuI[JPsiMAX][2], fJPsiMuCategory[JPsiMAX][2];
  float fJPsiCharge[JPsiMAX], fJPsiPhi[JPsiMAX], fJPsiTheta[JPsiMAX], fJPsiEta[JPsiMAX], fJPsiRapidity[JPsiMAX], fJPsiP[JPsiMAX], fJPsiPt[JPsiMAX], fJPsiPx[JPsiMAX], fJPsiPy[JPsiMAX], fJPsiPz[JPsiMAX], fJPsiEnergy[JPsiMAX], fJPsiEt[JPsiMAX], fJPsiMass[JPsiMAX], fJPsiMt[JPsiMAX], fJPsiChi2[JPsiMAX], fJPsiNdof[JPsiMAX], fJPsiVx[JPsiMAX], fJPsiVy[JPsiMAX], fJPsiVz[JPsiMAX], fJPsiVxE[JPsiMAX], fJPsiVyE[JPsiMAX], fJPsiVzE[JPsiMAX], fJPsiVtxPhi[JPsiMAX], fJPsiVtxTheta[JPsiMAX], fJPsiVtxEta[JPsiMAX], fJPsiVtxRapidity[JPsiMAX], fJPsiVtxP[JPsiMAX], fJPsiVtxPt[JPsiMAX], fJPsiVtxPx[JPsiMAX], fJPsiVtxPy[JPsiMAX], fJPsiVtxPz[JPsiMAX], fJPsiVtxEnergy[JPsiMAX], fJPsiVtxEt[JPsiMAX], fJPsiVtxMass[JPsiMAX], fJPsiVtxMt[JPsiMAX];
  bool fJPsiMuCutKin[JPsiMAX][2], fJPsiMuCutHLT[JPsiMAX][2], fJPsiMuCutIso[JPsiMAX][2], fJPsiMuCutSA[JPsiMAX][2], fJPsiMuCutTrk[JPsiMAX][2], fJPsiMuType[JPsiMAX][2][5], fJPsiBasicFilter[JPsiMAX], fJPsiVtxBasicFilter[JPsiMAX];
  // -- isolation information
  int fJPsiFromClosestI[JPsiMAX][5], fNJPsiSharedTk2[JPsiMAX], fNJPsiWithin05[JPsiMAX], fNJPsiWithin03[JPsiMAX];
  float fJPsidRToClosest[JPsiMAX][5];

  // -- Eta_b->2 J/Psi candidates
  const int ETABMAX = 200;

  int fEtabN, fBaseEtabI, fEtabBestMassI, fEtabBestProbI, fEtabVtxBestMassI, fEtabVtxBestProbI, fEtabBest4MuProbI;
  int fEtabIndex[ETABMAX], fEtabDuplicatesI[ETABMAX], fEtabJPsiI[ETABMAX][2], fEtabMuI[ETABMAX][4], fEtabMuN[ETABMAX], fEtabL3MatchMuN[ETABMAX], fEtabToRePvI[ETABMAX];
  float fEtabCharge[ETABMAX], fEtabPhi[ETABMAX], fEtabTheta[ETABMAX], fEtabEta[ETABMAX], fEtabRapidity[ETABMAX], fEtabP[ETABMAX], fEtabPt[ETABMAX], fEtabPx[ETABMAX], fEtabPy[ETABMAX], fEtabPz[ETABMAX], fEtabEnergy[ETABMAX], fEtabEt[ETABMAX], fEtabMass[ETABMAX], fEtabMt[ETABMAX], fEtabChi2[ETABMAX], fEtabNdof[ETABMAX], fEtabVx[ETABMAX], fEtabVy[ETABMAX], fEtabVz[ETABMAX], fEtabVxE[ETABMAX], fEtabVyE[ETABMAX], fEtabVzE[ETABMAX], fEtabVtxPhi[ETABMAX], fEtabVtxTheta[ETABMAX], fEtabVtxEta[ETABMAX], fEtabVtxRapidity[ETABMAX], fEtabVtxP[ETABMAX], fEtabVtxPt[ETABMAX], fEtabVtxPx[ETABMAX], fEtabVtxPy[ETABMAX], fEtabVtxPz[ETABMAX], fEtabVtxEnergy[ETABMAX], fEtabVtxEt[ETABMAX], fEtabVtxMass[ETABMAX], fEtabVtxMt[ETABMAX], fEtabCT[ETABMAX], fEtabCTxy[ETABMAX], fEtabVtxCT[ETABMAX], fEtabVtxCTxy[ETABMAX], fEtabJPsiDeltaL[ETABMAX], fEtabJPsiDeltaT[ETABMAX], fEtabJPsiVtxErr[ETABMAX], fEtabJPsiVtxErrxy[ETABMAX], fEtabJPsiProjX[ETABMAX][2], fEtabJPsiProjY[ETABMAX][2], fEtabJPsiProjZ[ETABMAX][2], fEtabJPsiCT[ETABMAX][2], fEtabJPsiCTxy[ETABMAX][2], fEtabJPsiToPVVtxErr[ETABMAX][2], fEtabJPsiToPVVtxErrxy[ETABMAX][2], fEtabJPsiVtxCT[ETABMAX][2], fEtabJPsiVtxCTxy[ETABMAX][2];
  bool fEtabBasicFilter[ETABMAX], fEtabVtxBasicFilter[ETABMAX];
  // isolation information
  int fEtabJPsiIsoTkN[ETABMAX][2];
  float fEtabJPsiIso7PV[ETABMAX][2], fEtabJPsiIsoTkCA[ETABMAX][2];

  ////Himal addition for HLT object                                                                                                                              
  float fHLT_Muon_Eta[HLTMAX], fHLT_Muon_Phi[HLTMAX],  fHLT_Muon_Pt[HLTMAX], fHLT_Muon_VertexmumuJpsi[HLTMAX], fHLT_Muon_TripleMuL3[HLTMAX];
  //end                                                                                                                                                        
  float fGenMuonPt[15], fGenMuonEta[15],  fGenMuonPhi[15];




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

  // -- Variables to output
  float FourMu_CT, FourMu_CTxy, FourMu_Rapidity, FourMu_Mass, FourMu_pT, FourMu_VtxProb, Psi1To2DistTot, Psi1To2Significance, Psi1To2DistT, Psi1To2DistL, Psi1To2_S, Psi1To2_dY, Psi_Mass[2], PsiW_Mass[4], Psi_y[2], Psi1_eta, Psi1_phi, Psi_pT[2], Psi1_p, Psi1_px, Psi1_py, Psi1_pz, Psi1_CT, Psi_CTxy[2], Psi_VtxProb[2], Psi1_CTErr, Psi1_CTErrxy, Psi1_CTSig, Psi1_CTxySig, Psi2_y, Psi2_eta, Psi2_phi, Psi2_p, Psi2_px, Psi2_py, Psi2_pz, Psi2_CT, Psi2_CTErr, Psi2_CTErrxy, Psi2_CTSig, Psi2_CTxySig, Mu11_eta, Mu11_phi, Mu_px[4], Mu_py[4], Mu_pz[4], Mu_pT[4], Mu_Eta[4], Mu12_phi, Mu12_pT, Mu21_eta, Mu21_phi, Mu21_pT, Mu22_eta, Mu22_phi, Mu22_pT, Mu11_charge, Mu12_charge, Mu21_charge, Mu22_charge, Mu_SAdist[4], Mu_SAdR[4], Mu_d0[4], Mu_dz[4],Mu_Iso[4],Mu_EcalIso[4], Mu_HcalIso[4],Mu_TrakIso[4],Mu_TIcorr,Mu_helicity[2], fHLT_MuonJpsi_L3Matching_Pt[5],fHLT_MuonJpsi_L3Matching_Eta[5],fHLT_MuonJpsi_L3Matching_Phi[5],fHLT_MuonJpsi_VtxMatching_Eta[5],fHLT_MuonJpsi_VtxMatching_Pt[5],fHLT_MuonJpsi_VtxMatching_Phi[5], fHLT_MuonZ_L3Matching_Pt[5],fHLT_MuonZ_L3Matching_Eta[5],fHLT_MuonZ_L3Matching_Phi[5],fHLT_MuonZ_VtxMatching_Eta[5],fHLT_MuonZ_VtxMatching_Pt[5],fHLT_MuonZ_VtxMatching_Phi[5];
  
  int   Mu_charge[4], muonI[4];
  bool  PsiNearSamePV, Psi_is_Seagull[2], Mu_Arbitrated[4], Mu_LastTight[4], Mu_OneTight[4], Mu_2DCompatTight[4], Mu_Tight[4];


  // create the branches in the new tree
  fTree->Branch("run",          &fRun,          "run/i");
  fTree->Branch("lumiblock",    &fLumiBlock,    "lumiblock/i");
  fTree->Branch("event",        &fEvent,        "event/i");
  fTree->Branch("bx",           &fBX,           "bx/I");
  fTree->Branch("orbit",        &fOrbit,        "orbit/I");

  fTree->Branch("FourMu_CT",      &FourMu_CT,           "FourMu_CT/F");
  fTree->Branch("FourMu_CTxy",    &FourMu_CTxy,         "FourMu_CTxy/F");
  fTree->Branch("FourMu_Mass",    &FourMu_Mass,         "FourMu_Mass/F");
  fTree->Branch("FourMu_Rapidity",    &FourMu_Rapidity,         "FourMu_Rapidity/F");
  fTree->Branch("FourMu_pT",      &FourMu_pT,           "FourMu_pT/F");
  fTree->Branch("FourMu_VtxProb", &FourMu_VtxProb,      "FourMu_VtxProb/F");
  fTree->Branch("Psi1To2DistTot", &Psi1To2DistTot,      "Psi1To2DistTot/F");
  fTree->Branch("Psi1To2Significance",&Psi1To2Significance,"Psi1To2Significance/F");
  fTree->Branch("Psi1To2DistT",   &Psi1To2DistT,        "Psi1To2DistT/F");
  fTree->Branch("Psi1To2DistL",   &Psi1To2DistL,        "Psi1To2DistL/F");
  fTree->Branch("Psi1To2_S",      &Psi1To2_S,           "Psi1To2_S/F");
  fTree->Branch("Psi1To2_dY",     &Psi1To2_dY,          "Psi1To2_dY/F");
  fTree->Branch("Psi_Mass",       Psi_Mass,             "Psi_Mass[2]/F");
  fTree->Branch("PsiW_Mass",       PsiW_Mass,             "PsiW_Mass[4]/F");
  fTree->Branch("Psi_y",          Psi_y,                "Psi_y[2]/F");
  fTree->Branch("Psi1_eta",       &Psi1_eta,            "Psi1_eta/F");
  fTree->Branch("Psi1_phi",       &Psi1_phi,            "Psi1_phi/F");
  fTree->Branch("Psi_pT",         Psi_pT,               "Psi_pT[2]/F");
  fTree->Branch("Psi1_p",         &Psi1_p,              "Psi1_p/F");
  fTree->Branch("Psi1_px",        &Psi1_px,             "Psi1_px/F");
  fTree->Branch("Psi1_py",        &Psi1_py,             "Psi1_py/F");
  fTree->Branch("Psi1_pz",        &Psi1_pz,             "Psi1_pz/F");
  fTree->Branch("Psi1_CT",        &Psi1_CT,             "Psi1_CT/F");
  fTree->Branch("Psi_CTxy",       Psi_CTxy,             "Psi_CTxy[2]/F");
  fTree->Branch("Psi_VtxProb",    Psi_VtxProb,          "Psi_VtxProb[2]/F");
  fTree->Branch("Psi1_CTErr",     &Psi1_CTErr,          "Psi1_CTErr/F");
  fTree->Branch("Psi1_CTErrxy",   &Psi1_CTErrxy,        "Psi1_CTErrxy/F");
  fTree->Branch("Psi1_CTSig",     &Psi1_CTSig,          "Psi1_CTSig/F");
  fTree->Branch("Psi1_CTxySig",   &Psi1_CTxySig,        "Psi1_CTxySig/F");
  fTree->Branch("Psi2_eta",       &Psi2_eta,            "Psi2_eta/F");
  fTree->Branch("Psi2_phi",       &Psi2_phi,            "Psi2_phi/F");
  fTree->Branch("Psi2_p",         &Psi2_p,              "Psi2_p/F");
  fTree->Branch("Psi2_px",        &Psi2_px,             "Psi2_px/F");
  fTree->Branch("Psi2_py",        &Psi2_py,             "Psi2_py/F");
  fTree->Branch("Psi2_pz",        &Psi2_pz,             "Psi2_pz/F");
  fTree->Branch("Psi2_CT",        &Psi2_CT,             "Psi2_CT/F");
  fTree->Branch("Psi2_CTErr",     &Psi2_CTErr,          "Psi2_CTErr/F");
  fTree->Branch("Psi2_CTErrxy",   &Psi2_CTErrxy,        "Psi2_CTErrxy/F");
  fTree->Branch("Psi2_CTSig",     &Psi2_CTSig,          "Psi2_CTSig/F");
  fTree->Branch("Psi2_CTxySig",   &Psi2_CTxySig,        "Psi2_CTxySig/F");
  fTree->Branch("Mu_charge",      Mu_charge,            "Mu_charge[4]/I");
  fTree->Branch("Mu_Eta",         Mu_Eta,               "Mu_Eta[4]/F");
  fTree->Branch("Mu_px",          Mu_px,                "Mu_px[4]/F");
  fTree->Branch("Mu_py",          Mu_py,                "Mu_py[4]/F");
  fTree->Branch("Mu_pz",          Mu_pz,                "Mu_pz[4]/F");
  fTree->Branch("Mu_pT",          Mu_pT,                "Mu_pT[4]/F");
  fTree->Branch("Mu_SAdist",      Mu_SAdist,            "Mu_SAdist[4]/F");
  fTree->Branch("Mu_SAdR",        Mu_SAdR,              "Mu_SAdR[4]/F");
  fTree->Branch("Mu_d0",          Mu_d0,                "Mu_d0[4]/F");
  fTree->Branch("Mu_dz",          Mu_dz,                "Mu_dz[4]/F");
  fTree->Branch("Mu11_phi",       &Mu11_phi,            "Mu11_phi/F");
  fTree->Branch("Mu12_phi",       &Mu12_phi,            "Mu12_phi/F");
  fTree->Branch("Mu21_phi",       &Mu21_phi,            "Mu21_phi/F");
  fTree->Branch("Mu22_phi",       &Mu22_phi,            "Mu22_phi/F");
  fTree->Branch("Mu_Arbitrated",  Mu_Arbitrated,        "Mu_Arbitrated[4]/O");
  fTree->Branch("Mu_LastTight",   Mu_LastTight,         "Mu_LastTight[4]/O");
  fTree->Branch("Mu_OneTight",    Mu_OneTight,          "Mu_OneTight[4]/O");
  fTree->Branch("Mu_2DCompatTight",Mu_2DCompatTight,    "Mu_2DCompatTight[4]/O");
  fTree->Branch("Mu_Tight",       Mu_Tight,             "Mu_Tight[4]/O");
  fTree->Branch("PsiNearSamePV",  &PsiNearSamePV,       "PsiNearSamePV/O");
  fTree->Branch("Psi_is_Seagull", Psi_is_Seagull,       "Psi_is_Seagull[2]/O");
  fTree->Branch("PvN",            &fPvN,                "PvN/I");
  fTree->Branch("Mu_Iso",         &Mu_Iso,              "Mu_Iso[4]/F");
  fTree->Branch("Mu_EcalIso",     &Mu_EcalIso,          "Mu_EcalIso[4]/F");
  fTree->Branch("Mu_HcalIso",     &Mu_HcalIso,          "Mu_HcalIso[4]/F");
  fTree->Branch("Mu_TrakIso",    &Mu_TrakIso,         "Mu_TrakIso[4]/F");
  fTree->Branch("Mu_TIcorr",    &Mu_TIcorr,         "Mu_TIcorr/F");
  fTree->Branch("Mu_helicity",    &Mu_helicity,         "Mu_helicity[2]/F");
  //Himal addition for trigger matching
  
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
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi3p5_Muon2", fHLTP_Dimuon0_Jpsi3p5_Muon2);
  //hima add
  //JPsi
  chainData5fb.SetBranchAddress("HLTP_IsoMu20", fHLTP_IsoMu20);
  chainData5fb.SetBranchAddress("HLTP_IsoMu27", fHLTP_IsoMu27);
  chainData5fb.SetBranchAddress("HLTP_Mu50", fHLTP_Mu50);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_L3Matching_Pt",        fHLT_MuonJpsi_L3Matching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_L3Matching_Eta",        fHLT_MuonJpsi_L3Matching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_L3Matching_Phi",        fHLT_MuonJpsi_L3Matching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_VtxMatching_Pt",        fHLT_MuonJpsi_VtxMatching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_VtxMatching_Eta",        fHLT_MuonJpsi_VtxMatching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonJpsi_VtxMatching_Phi",        fHLT_MuonJpsi_VtxMatching_Phi);
  //Z
  chainData5fb.SetBranchAddress("HLT_MuonZ_L3Matching_Pt",        fHLT_MuonZ_L3Matching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonZ_L3Matching_Eta",        fHLT_MuonZ_L3Matching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonZ_L3Matching_Phi",        fHLT_MuonZ_L3Matching_Phi);
  chainData5fb.SetBranchAddress("HLT_MuonZ_VtxMatching_Pt",        fHLT_MuonZ_VtxMatching_Pt);
  chainData5fb.SetBranchAddress("HLT_MuonZ_VtxMatching_Eta",        fHLT_MuonZ_VtxMatching_Eta);
  chainData5fb.SetBranchAddress("HLT_MuonZ_VtxMatching_Phi",        fHLT_MuonZ_VtxMatching_Phi);
  //end                                                                                                                
  chainData5fb.SetBranchAddress("HLTP_QuadMuon0_Dimuon0_Jpsi", fHLTP_QuadMuon0_Dimuon0_Jpsi);
  //  chainData5fb.SetBranchAddress("HLTP_QuadMuon0_Dimuon0_Upsilon", fHLTP_QuadMuon0_Dimuon0_Upsilon);
  chainData5fb.SetBranchAddress("HLTP_Dimuon10_Jpsi_Barrel", fHLTP_Dimuon10_Jpsi_Barrel);
  chainData5fb.SetBranchAddress("HLTP_TripleMu5",         fHLTP_TripleMu5);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu3_PS",         &fHLTP_DoubleMu3_PS);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu6_PS",         &fHLTP_DoubleMu6_PS);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu7_PS",         &fHLTP_DoubleMu7_PS);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_PS",      &fHLTP_Dimuon0_Jpsi_PS);
  chainData5fb.SetBranchAddress("HLTP_Dimuon10_Jpsi_Barrel_PS", &fHLTP_Dimuon10_Jpsi_Barrel_PS);
  chainData5fb.SetBranchAddress("HLTP_TripleMu5_PS",         &fHLTP_TripleMu5_PS);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu3_Filters",    fHLTP_DoubleMu3_Filters);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu6_Filters",    fHLTP_DoubleMu6_Filters);
  chainData5fb.SetBranchAddress("HLTP_DoubleMu7_Filters",    fHLTP_DoubleMu7_Filters);
  chainData5fb.SetBranchAddress("HLTP_TripleMu5_Filters",    fHLTP_TripleMu5_Filters);
  chainData5fb.SetBranchAddress("HLTP_Dimuon0_Jpsi_Filters", fHLTP_Dimuon0_Jpsi_Filters);
  chainData5fb.SetBranchAddress("HLTP_Dimuon10_Jpsi_Barrel_Filters", fHLTP_Dimuon10_Jpsi_Barrel_Filters);

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
  chainData5fb.SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", fHLT_Dimuon0_Jpsi3p5_Muon2);
  chainData5fb.SetBranchAddress("HLT_Dimuon0_Jpsi", fHLT_Dimuon0_Jpsi);
  chainData5fb.SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced", fHLT_Dimuon7_Jpsi_Displaced);
  chainData5fb.SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel", fHLT_Dimuon7_Jpsi_X_Barrel);
  chainData5fb.SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel", fHLT_Dimuon10_Jpsi_Barrel);
  chainData5fb.SetBranchAddress("HLT_TripleMu5", fHLT_TripleMu5);
  chainData5fb.SetBranchAddress("HLT_Jet",        fHLT_Jet);
  //Gen Matching branch  
  //chainData5fb.SetBranchAddress("GenMuonPt",           fGenMuonPt);
  //chainData5fb.SetBranchAddress("GenMuonEta",           fGenMuonEta);
  //chainData5fb.SetBranchAddress("GenMuonPhi",           fGenMuonPhi);
  //end                                                       
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
  chainData5fb.SetBranchAddress("JPsiClosestPVinZ", fJPsiClosestPVinZ);
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
  //added by Himal                                                                                                                       
  chainData5fb.SetBranchAddress("HLT_Muon_Eta",        fHLT_Muon_Eta);
  chainData5fb.SetBranchAddress("HLT_Muon_Phi",        fHLT_Muon_Phi);
  chainData5fb.SetBranchAddress("HLT_Muon_Pt",         fHLT_Muon_Pt);
  chainData5fb.SetBranchAddress("HLT_Muon_VertexmumuJpsi",       fHLT_Muon_VertexmumuJpsi);
  chainData5fb.SetBranchAddress("HLT_Muon_TripleMuL3",         fHLT_Muon_TripleMuL3);
  //end   


  chainData5fb.SetBranchAddress("HJtI",         fHJtI);

  chainData5fb.SetBranchAddress("JPsiBestProbI",    fJPsiBestProbI);
  chainData5fb.SetBranchAddress("EtabBestMassI",    &fEtabBestMassI);
  chainData5fb.SetBranchAddress("EtabBestProbI",    &fEtabBestProbI);
  chainData5fb.SetBranchAddress("EtabBest4MuProbI",    &fEtabBest4MuProbI);
  chainData5fb.SetBranchAddress("JPsiVtxBestProbI", fJPsiVtxBestProbI);
  chainData5fb.SetBranchAddress("EtabVtxBestMassI", &fEtabVtxBestMassI);
  chainData5fb.SetBranchAddress("EtabVtxBestProbI", &fEtabVtxBestProbI);


  //loop over events and get entries
  Int_t nevent = chainData5fb.GetEntries();
  //initialiaziation
  float mu_mass=0.1056583745;
  float Mu_mass=0.1056583745;

  Int_t ncand(0);
  Int_t nev(0),oev(0);
  Int_t nerr(0);
  Int_t nerr1(0);
  Int_t nerr2(0);
  Int_t nerr3(0);
  Int_t nerr4(0);
  Int_t nerr5(0);
  int muacc=0;
  int totaletab=0;
  for (int i=0;i<4;i++){ 
    Mu_Iso[i] = -9999;
    Mu_EcalIso[i]=-9999;
    Mu_HcalIso[i]=-9999;
    Mu_TrakIso[i]=-9999;
  }


  for (Int_t i=0; i<nevent; i++) { 
    if( i % 10 == 0 ) cout<<i+1<<" of "<<nevent<<endl;
    
    chainData5fb.GetEvent(i);    //read complete accepted event in memory
    
    totaletab++;
    if (fHLTP_Dimuon0_Jpsi3p5_Muon2[1]==false) continue;
    //if (fHLTP_Trimuon5_3p5_2_Upsilon_Muon[1]==false) continue;
    //if (fHLTP_IsoMu20[1]==false) continue;
    // if (fHLTP_IsoMu27[1]==false) continue;
    nerr1++;
    // <sms> 
    //    if(fEtabN>1) continue;
    for (Int_t ie=0; ie<fEtabN; ie++) {
      // Best Candidate decision      if( fEtabJPsiVtxErr[ie]>=0 && ie==fEtabBest4MuProbI ) {
      //<sms>
      if( fEtabJPsiVtxErr[ie]<0) nerr++;
      if( fEtabJPsiVtxErr[ie]>=0 ) {
	nerr2++;
        Int_t ijp1, ijp2, ejp1, ejp2;
        if( fJPsiVtxPt[fEtabJPsiI[ie][0]]>fJPsiVtxPt[fEtabJPsiI[ie][1]] ) { ijp1 = fEtabJPsiI[ie][0]; ijp2 = fEtabJPsiI[ie][1]; ejp1 = 0; ejp2 = 1;}
	else { ijp1 = fEtabJPsiI[ie][1]; ijp2 = fEtabJPsiI[ie][0]; ejp1 = 1; ejp2 = 0; }
        
	FourMu_CT    = fEtabVtxCT[ie];
        FourMu_CTxy  = fEtabVtxCTxy[ie];
	
	//addition of Isolation calculation 
	for (int imu=0;imu<4;imu++){
	  Mu_Iso[imu] = fPcTrackIso[fJPsiMuI[fEtabJPsiI[ie][ejp1]][imu]]/fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][imu]];
	  Mu_EcalIso[imu]= fPcEcalIso[fJPsiMuI[fEtabJPsiI[ie][ejp1]][imu]];
	  Mu_HcalIso[imu]= fPcHcalIso[fJPsiMuI[fEtabJPsiI[ie][ejp1]][imu]];
	  Mu_TrakIso[imu]= fPcTrackIso[fJPsiMuI[fEtabJPsiI[ie][ejp1]][imu]];
	}
	

        FourMu_Mass  = fEtabVtxMass[ie];
        FourMu_Rapidity = fEtabVtxRapidity[ie];
        FourMu_pT    = fEtabVtxPt[ie];
        FourMu_VtxProb=TMath::Prob(fEtabChi2[ie],fEtabNdof[ie]);
        Psi1To2DistTot = sqrt(pow(fJPsiVx[ijp1]-fJPsiVx[ijp2],2)+pow(fJPsiVy[ijp1]-fJPsiVy[ijp2],2)+pow(fJPsiVz[ijp1]-fJPsiVz[ijp2],2));
        Psi1To2Significance = Psi1To2DistTot/fEtabJPsiVtxErr[ie];
        Psi1To2DistT = fEtabJPsiDeltaT[ie];
        Psi1To2DistL = fEtabJPsiDeltaL[ie];
        Psi1To2_S    = 0.5*( sqrt( pow(fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]]+fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]], 2)+pow(fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]]+fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]], 2) )/(fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]]+fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]]) + sqrt( pow(fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]]+fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]], 2)+pow(fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]]+fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]], 2) )/(fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]]+fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]]) );
        Psi1To2_dY   = fabs(fJPsiVtxRapidity[fEtabJPsiI[ie][ejp1]]-fJPsiVtxRapidity[fEtabJPsiI[ie][ejp2]]);
        Psi_Mass[0]  = fJPsiVtxMass[fEtabJPsiI[ie][ejp1]];
        Psi_y[0]     = fJPsiVtxRapidity[fEtabJPsiI[ie][ejp1]];
        Psi1_eta     = fJPsiVtxEta[fEtabJPsiI[ie][ejp1]];
        Psi1_phi     = fJPsiVtxPhi[fEtabJPsiI[ie][ejp1]];
        Psi_pT[0]    = fJPsiVtxPt[fEtabJPsiI[ie][ejp1]];
        Psi1_p       = fJPsiVtxP[fEtabJPsiI[ie][ejp1]];
        Psi1_px      = fJPsiVtxPx[fEtabJPsiI[ie][ejp1]];
        Psi1_py      = fJPsiVtxPy[fEtabJPsiI[ie][ejp1]];
        Psi1_pz      = fJPsiVtxPz[fEtabJPsiI[ie][ejp1]];
        Psi_VtxProb[0]=TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][ejp1]], fJPsiNdof[fEtabJPsiI[ie][ejp1]]);
        Psi1_CT      = fEtabJPsiVtxCT[ie][ejp1];
        Psi_CTxy[0]  = fEtabJPsiVtxCTxy[ie][ejp1];
        Psi1_CTErr   = (fJPsiVtxMass[fEtabJPsiI[ie][ejp1]]/fJPsiVtxP[fEtabJPsiI[ie][ejp1]])*fEtabJPsiToPVVtxErr[ie][ejp1];
        Psi1_CTErrxy = (fJPsiVtxMass[fEtabJPsiI[ie][ejp1]]/fJPsiVtxPt[fEtabJPsiI[ie][ejp1]])*fEtabJPsiToPVVtxErrxy[ie][ejp1];
        Psi_Mass[1]  = fJPsiVtxMass[fEtabJPsiI[ie][ejp2]];
        Psi_y[1]     = fJPsiVtxRapidity[fEtabJPsiI[ie][ejp2]];
        Psi2_eta     = fJPsiVtxEta[fEtabJPsiI[ie][ejp2]];
        Psi2_phi     = fJPsiVtxPhi[fEtabJPsiI[ie][ejp2]];
        Psi_pT[1]    = fJPsiVtxPt[fEtabJPsiI[ie][ejp2]];
        Psi2_p       = fJPsiVtxP[fEtabJPsiI[ie][ejp2]];
        Psi2_px      = fJPsiVtxPx[fEtabJPsiI[ie][ejp2]];
        Psi2_py      = fJPsiVtxPy[fEtabJPsiI[ie][ejp2]];
        Psi2_pz      = fJPsiVtxPz[fEtabJPsiI[ie][ejp2]];
        Psi_VtxProb[1]=TMath::Prob(fJPsiChi2[fEtabJPsiI[ie][ejp2]], fJPsiNdof[fEtabJPsiI[ie][ejp2]]);
        Psi2_CT      = fEtabJPsiVtxCT[ie][ejp2];
        Psi_CTxy[1]  = fEtabJPsiVtxCTxy[ie][ejp2];
        Psi2_CTErr   = (fJPsiVtxMass[fEtabJPsiI[ie][ejp2]]/fJPsiVtxP[fEtabJPsiI[ie][ejp2]])*fEtabJPsiToPVVtxErr[ie][ejp2];
        Psi2_CTErrxy = (fJPsiVtxMass[fEtabJPsiI[ie][ejp2]]/fJPsiVtxPt[fEtabJPsiI[ie][ejp2]])*fEtabJPsiToPVVtxErrxy[ie][ejp2];
	
        if (fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]]>=fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]]) {iMu0 = 0; iMu1 = 1;}
        else{iMu0 = 1; iMu1 = 0;}
        Mu_Eta[0]    = fPcEta[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];
        Mu11_phi     = fPcPhi[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];
        Mu_pT[0]     = fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];
        Mu_Eta[1]    = fPcEta[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];
        Mu12_phi     = fPcPhi[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];
        Mu_pT[1]     = fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];
        Mu_Eta[2]    = fPcEta[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];
        Mu21_phi     = fPcPhi[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];
        Mu_pT[2]     = fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];
        Mu_Eta[3]    = fPcEta[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];
        Mu22_phi     = fPcPhi[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];
        Mu_pT[3]     = fPcPt[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];

        Mu11_charge  = fPcCharge[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];

        Mu_px[0]      = fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];
        Mu_py[0]      = fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];
        Mu_pz[0]      = fPcPz[fJPsiMuI[fEtabJPsiI[ie][ejp1]][0]];

        Mu12_charge  = fPcCharge[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];

        Mu_px[1]      = fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];
        Mu_py[1]      = fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];
        Mu_pz[1]      = fPcPz[fJPsiMuI[fEtabJPsiI[ie][ejp1]][1]];

        Mu21_charge  = fPcCharge[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];

        Mu_px[2]      = fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];
        Mu_py[2]      = fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];
        Mu_pz[2]      = fPcPz[fJPsiMuI[fEtabJPsiI[ie][ejp2]][0]];

        Mu22_charge   = fPcCharge[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];

        Mu_px[3]      = fPcPx[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];
        Mu_py[3]      = fPcPy[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];
        Mu_pz[3]      = fPcPz[fJPsiMuI[fEtabJPsiI[ie][ejp2]][1]];

        muonI[0]      = fJPsiMuI[fEtabJPsiI[ie][ejp1]][0];
        muonI[1]      = fJPsiMuI[fEtabJPsiI[ie][ejp1]][1];

        muonI[2]      = fJPsiMuI[fEtabJPsiI[ie][ejp2]][0];
        muonI[3]      = fJPsiMuI[fEtabJPsiI[ie][ejp2]][1];

        Psi1_CTSig   = Psi1_CT/Psi1_CTErr;
        Psi1_CTxySig = Psi_CTxy[0]/Psi1_CTErrxy;
        Psi2_CTSig   = Psi2_CT/Psi2_CTErr;
        Psi2_CTxySig = Psi_CTxy[1]/Psi2_CTErrxy;
	
	
	float Mu_phi[4] = {Mu11_phi, Mu12_phi, Mu21_phi,Mu22_phi};

        if( fJPsiClosestPVinZ[ijp1]==fJPsiClosestPVinZ[ijp2] ) PsiNearSamePV=true;
        else PsiNearSamePV=false;

        TVector3 mu1(0,0,0);
        TVector3 mu2(0,0,0);
        TVector3 mu3(0,0,0);
        TVector3 mu4(0,0,0);
        TVector3 xprod1(0,0,0);
        TVector3 xprod2(0,0,0);

        Psi_is_Seagull[0]=false;
        Psi_is_Seagull[1]=false;

        if( Mu11_charge > 0 && Mu12_charge < 0 ) {
          mu1.SetX(Mu_px[0]);
          mu1.SetY(Mu_py[0]);
          mu1.SetZ(Mu_pz[0]);


          mu2.SetX(Mu_px[1]);
          mu2.SetY(Mu_py[1]);
          mu2.SetZ(Mu_pz[1]);
        } else if( Mu11_charge < 0 && Mu12_charge > 0 ) {
          mu1.SetX(Mu_px[1]);
          mu1.SetY(Mu_py[1]);
          mu1.SetZ(Mu_pz[1]);

          mu2.SetX(Mu_px[0]);
          mu2.SetY(Mu_py[0]);
          mu2.SetZ(Mu_pz[0]);
        }
        xprod1 = mu1.Cross(mu2);
        if( xprod1.z() > 0 ) Psi_is_Seagull[0]=true;
        if( Mu21_charge > 0 && Mu22_charge < 0 ) {
          mu3.SetX(Mu_px[2]);
          mu3.SetY(Mu_py[2]);
          mu3.SetZ(Mu_pz[2]);

          mu4.SetX(Mu_px[3]);
          mu4.SetY(Mu_py[3]);
          mu4.SetZ(Mu_pz[3]);
        } else if( Mu21_charge < 0 && Mu22_charge > 0 ) {
          mu3.SetX(Mu_px[3]);
          mu3.SetY(Mu_py[3]);
          mu3.SetZ(Mu_pz[3]);

          mu4.SetX(Mu_px[2]);
          mu4.SetY(Mu_py[2]);
          mu4.SetZ(Mu_pz[2]);
        }
	//adding helicity and isolation Himal
	//helicity
	//defining TLorentz Vector


	TLorentzVector mu4Vector1;
	TLorentzVector mu4Vector2;
	TLorentzVector psi4Vector1;


	mu4Vector1.SetXYZM(mu1.X(),mu1.Y(),mu1.Z(),mu_mass);
	mu4Vector2.SetXYZM(mu2.X(),mu2.Y(),mu2.Z(),mu_mass);


	psi4Vector1 = mu4Vector1 + mu4Vector2;     


	//-----------

	TLorentzVector mu4Vector3;
	TLorentzVector mu4Vector4;
	mu4Vector3.SetXYZM(mu3.X(),mu3.Y(),mu3.Z(),mu_mass);
	mu4Vector4.SetXYZM(mu4.X(),mu4.Y(),mu4.Z(),mu_mass);

	TLorentzVector psi4Vector2;
	psi4Vector2 = mu4Vector3 + mu4Vector4;     


	TLorentzVector FourMu4Vector;
	FourMu4Vector = psi4Vector1 + psi4Vector2;


	// find other wrong neutral combination

	// 1-4   2-3 
        TLorentzVector psi4Vector14;
        TLorentzVector psi4Vector23;
        TLorentzVector psi4Vector13;
        TLorentzVector psi4Vector24;
	psi4Vector14 = mu4Vector1 + mu4Vector4;
	psi4Vector23 = mu4Vector2 + mu4Vector3;

	PsiW_Mass[0] = psi4Vector14.M();
	PsiW_Mass[1] = psi4Vector23.M();

	psi4Vector13 = mu4Vector1 + mu4Vector3;
	psi4Vector24 = mu4Vector2 + mu4Vector4;


	PsiW_Mass[2] = psi4Vector13.M();
	PsiW_Mass[3] = psi4Vector24.M();



	// ----


        Mu_helicity[0] = coshel(mu4Vector1,  psi4Vector1,  FourMu4Vector);
        Mu_helicity[1] = coshel(mu4Vector3,  psi4Vector2,  FourMu4Vector);
	
	//---------------------------------------------
	//.....
	//isolation - incorrect muon association
	Mu_TIcorr = 0.;
	if(Mu_TrakIso[0]>0){

	  // make vector corresponding to
	  TLorentzVector P0,P1,P2,P3;	  
	  P0.SetXYZM(Mu_px[0],Mu_py[0], Mu_pz[0],mu_mass); 
	  P1.SetXYZM(Mu_px[1],Mu_py[1], Mu_pz[1],mu_mass); 
	  P2.SetXYZM(Mu_px[2],Mu_py[2], Mu_pz[2],mu_mass); 
	  P3.SetXYZM(Mu_px[3],Mu_py[3], Mu_pz[3],mu_mass); 
	  
	  float DeltaR12 = P0.DeltaR(P1);
	  float DeltaR13 = P0.DeltaR(P2);
	  float DeltaR14 = P0.DeltaR(P3);
	  
	  if (DeltaR12<0.3) {
	    //	    cout<<" Himal deltaR12 "<< DeltaR12<< "  " << Mu_TrakIso[0] << " - " << P1.Pt() << endl;

	    // really subtract?
	    Mu_TIcorr = Mu_TrakIso[0] - P1.Pt(); 

	    if(P1.Pt()/Mu_TrakIso[0] >10.) Mu_TIcorr = Mu_TrakIso[0] ;

	  }
	  //	  if (DeltaR13<0.3) {cout<<" Himal deltaR13 "<< DeltaR13<< "  " << P2.Pt() << " - " << Mu_TrakIso[0] << endl;}
	  //	  if (DeltaR14<0.3) {cout<<" Himal deltaR14 "<< DeltaR14<< "  " << P3.Pt() << " - " << Mu_TrakIso[0] << endl;}
	}

	//---------------------------------

        xprod2 = mu3.Cross(mu4);
        if( xprod2.z() > 0 ) Psi_is_Seagull[1]=true;

        for(Int_t imu=0; imu<4; imu++ ) {
          Mu_charge[imu] = (int)fPcCharge[muonI[imu]];
          Mu_SAdist[imu] = fMuTkSADist[muonI[imu]];
          Mu_SAdR[imu] = fMuTkSAdR[muonI[imu]];
          Mu_d0[imu] = fPcD0[muonI[imu]];
          Mu_dz[imu] = fPcDz[muonI[imu]];
          Mu_Arbitrated[imu] = fMuArbitrated[muonI[imu]];
          Mu_LastTight[imu] = fMuLastStationTight[muonI[imu]];
          Mu_OneTight[imu] = fMuOneStationTight[muonI[imu]];
          Mu_2DCompatTight[imu] = fMu2DCompatibilityTight[muonI[imu]];
          Mu_Tight[imu] = fMuTightMatch[muonI[imu]];
	}



	nerr3++;

	// some final criteria
	//if( FourMu_VtxProb <= 0.05 ) continue;
	nerr5++;
	

	//Trigger and Gen Matching Stuff
	bool Mu_Matched_Himal[4];
	int matchedmuon=0;
        int  genmatchedMuon=0;
        float Mu_mass=0.1056583745;
	TLorentzVector MuReco,MuTrig,MuGen;
	bool Mu_MatchedZ_Himal[4];
	int matchedmuonZ=0;
        int  genmatchedMuonZ=0;
        
	TLorentzVector MuRecoZ,MuTrigZ,MuGenZ;
	Mu_Matched_Himal[0]=Mu_Matched_Himal[1]=Mu_Matched_Himal[2]=Mu_Matched_Himal[3]=false;
        for (int ihe =0;ihe<4;ihe++){
          for (int ihh =0;ihh<5;ihh++){
            if (fHLT_MuonJpsi_L3Matching_Pt[ihh]>=-10){
	      if (abs(fHLT_MuonJpsi_L3Matching_Pt[ihh] - Mu_pT[ihe])<10){
                MuReco.SetPtEtaPhiM(Mu_pT[ihe],Mu_Eta[ihe],Mu_phi[ihe],Mu_mass);
                MuTrig.SetPtEtaPhiM(fHLT_MuonJpsi_L3Matching_Pt[ihh],fHLT_MuonJpsi_L3Matching_Eta[ihh],fHLT_MuonJpsi_L3Matching_Phi[ihh],Mu_mass);
		if (MuReco.DeltaR(MuTrig)<0.01){
		  Mu_Matched_Himal[ihe]=true;
                  //cout<<"Found the matched muon"<<"The diff is "<<abs(fHLT_MuonJpsi_L3Matching_Pt[ihh]-Mu_pT[ihe])<<endl;
                  matchedmuon++;
                }
              }
            }
          }
	}
	Mu_MatchedZ_Himal[0]=Mu_MatchedZ_Himal[1]=Mu_MatchedZ_Himal[2]=Mu_MatchedZ_Himal[3]=false;
        for (int ihe =0;ihe<4;ihe++){
          for (int ihh =0;ihh<5;ihh++){
            if (fHLT_MuonZ_L3Matching_Pt[ihh]>=-10){
	      if (abs(fHLT_MuonZ_L3Matching_Pt[ihh] - Mu_pT[ihe])<10){
                MuReco.SetPtEtaPhiM(Mu_pT[ihe],Mu_Eta[ihe],Mu_phi[ihe],Mu_mass);
                MuTrig.SetPtEtaPhiM(fHLT_MuonZ_L3Matching_Pt[ihh],fHLT_MuonZ_L3Matching_Eta[ihh],fHLT_MuonZ_L3Matching_Phi[ihh],Mu_mass);
		if (MuReco.DeltaR(MuTrig)<0.01){
		  Mu_MatchedZ_Himal[ihe]=true;
                  //cout<<"Found the matched muon"<<"The diff is "<<abs(fHLT_MuonJpsi_L3Matching_Pt[ihh]-Mu_pT[ihe])<<endl;
                  matchedmuonZ++;
                }
              }
            }
          }
	}
	//only open when MC
	//for (int ihe =0;ihe<4;ihe++){
	//  for (int ihi =0;ihi<4;ihi++){
        //    if (fGenMuonPt[ihi]>-100 && fGenMuonEta[ihi]>-100 && fGenMuonPhi[ihi]>-100){
	//      MuReco.SetPtEtaPhiM(Mu_pT[ihe],Mu_Eta[ihe],Mu_phi[ihe],Mu_mass);
        //      MuGen.SetPtEtaPhiM(fGenMuonPt[ihi],fGenMuonEta[ihi],fGenMuonPhi[ihi],Mu_mass);
        //      if (MuReco.DeltaR(MuGen)<0.01){
	//       genmatchedMuon++;
	//      cout<<"DeltaR is given by ---> "<<MuReco.DeltaR(MuGen)<<endl;
	//    }
	//  }
	//}
        //}
        nerr3++;
        //cout <<"Number of gen matched muon = "<<genmatchedMuon<<endl;
	//end
	
	
	//if (Mu_pT[0]>2 && Mu_pT[1]>2 && Mu_pT[2]>2 && Mu_pT[3]>2 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 &&Psi_VtxProb[0]>0.005&&Psi_VtxProb[1]>0.005 && FourMu_VtxProb>0.005){ //Plot for cut for JPsi Pt(to address comments)
	if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 &&Psi_VtxProb[0]>0.005&&Psi_VtxProb[1]>0.005 && FourMu_VtxProb>0.005){ //Plot for cut for JPsi Pt(to address comments)
	//if (Mu_pT[0]>4 && Mu_pT[1]>4 && Mu_pT[2]>4 && Mu_pT[3]>4 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 &&Psi_VtxProb[0]>0.005&&Psi_VtxProb[1]>0.005 && FourMu_VtxProb>0.005){ //Plot for cut for JPsi Pt(to address comments)
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5){ //baseline cuts                                                                                           
	//addind dist1to2 cut <0.05
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && Psi_Mass[0]<3.2 && Psi_Mass[0] > 3.0 && Psi_Mass[1]<3.25 && Psi_Mass[1] > 2.95 && Psi_VtxProb[0]>0.005 && Psi_VtxProb[1]>0.005 && Psi1To2Significance<10){
	//<himsms>  
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && Psi_Mass[0]<3.2 && Psi_Mass[0] > 3.0 && Psi_Mass[1]<3.2 && Psi_Mass[1] > 3.0 && Psi_VtxProb[0]>0.01 && Psi_VtxProb[1]>0.01 && Psi1To2Significance<3 && FourMu_VtxProb>0.05){
	//end of distance cut
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4 && abs(Mu_Eta[1]) < 2.4 && abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && FourMu_VtxProb>0.05){                                                                                                  
	  
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4&& Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && FourMu_Mass>30 && FourMu_VtxProb>0.05){                                                                                
	  //blinding
	  //if (FourMu_Mass>83.5 && FourMu_Mass<98.5)continue;
	  //if (FourMu_Mass>120.0 && FourMu_Mass<130.0)continue;
	
	// + 7                                                                                                                                    
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && FourMu_VtxProb>0.05 && FourMu_Mass>30 && Psi_Mass[0]<3.2 && Psi_Mass[0] > 3.0 && Psi_Mass[1]< 3.25 && Psi_Mass[1] > 2.95){  
	// + trackiso(new +8)                                                                    
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && FourMu_VtxProb>0.05 && FourMu_Mass>30 && Psi_Mass[0]<3.2 && Psi_Mass[0] > 3.0 && Psi_Mass[1]< 3.25 && Psi_Mass[1] > 2.95 && Mu_TIcorr/Mu_pT[0] < 0.5){                                                                                                   
	// + deltaY                                                                                                             
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && FourMu_VtxProb>0.05 && Psi_Mass[0]<3.2 && Psi_Mass[0] > 3.0 && Psi_Mass[1]< 3.25 && Psi_Mass[1] > 2.95 &&Psi1To2_dY<3.0 && Mu_TIcorr/Mu_pT[0] < 0.5 && FourMu_Mass>30){                                                                                                   
	
	//+ All cut)    
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && FourMu_VtxProb>0.05 &&Psi1To2_dY<3.0 && Mu_TIcorr/Mu_pT[0] < 0.5){ 
	  //if  (Mu_TIcorr/Mu_pT[0] > 0.1) continue;
	  if (((Psi_Mass[0] < 2.8 || Psi_Mass[0] > 4.0)&&(Psi_Mass[0] < 70 || Psi_Mass[0] > 110)) || ((Psi_Mass[1] < 2.8 || Psi_Mass[1] > 4.0)&&(Psi_Mass[1] < 70 || Psi_Mass[1] > 110))) continue; //Himal change for ZJ
	  //if ((Psi_Mass[0] < 8.0 || Psi_Mass[0] > 12.0) && (Psi_Mass[1] < 8.0 || Psi_Mass[1] > 12.0)) continue; //Himal change for YY
	  //if (Psi_Mass[0] < 8. || Psi_Mass[0] > 12.) continue;//for YY
	  //if (Psi_Mass[1] < 8. || Psi_Mass[1] > 12.) continue;
	  
	  if (abs(Psi_Mass[0]-Psi_Mass[1])<50) continue;
	  //if (Mu_TIcorr/Mu_pT[0] > 0.5) continue;
	//if (matchedmuon<3) continue;
	//if ((fJPsiMass[ijp1]&&fJPsiMass[ijp2])>50.) continue;
	//TRIGGER MATCHING
	  //if (Psi_Mass[0]>3.4 && Psi_Mass[0]<80.) continue;
	  //if (Psi_Mass[1]>3.4 && Psi_Mass[1]<80.) continue;
	
	
	//new criteria for the efficiency study                                                                                            
	//if (Mu_pT[0]>3 && Mu_pT[1]>3 && Mu_pT[2]>3 && Mu_pT[3]>3 && abs(Mu_Eta[0])<2.4&& abs(Mu_Eta[1])<2.4&&abs(Mu_Eta[2])<2.4&& abs(Mu_Eta[3])<2.4 && Psi_pT[0]>3.5 && Psi_pT[1]>3.5 && Psi_Mass[0]<3.2 && Psi_Mass[0] > 3.0 && Psi_Mass[1]< 3.25 && Psi_Mass[1] > 2.95 ){                 
	  
	  //GENMatching
	  //if (genmatchedMuon<4) continue;

	  nerr4++;
	  
	  if(i != oev){
	    nev++;
	    oev = i;
	    if(ncand>0) h1->Fill((float) ncand); 
            ncand=1;
	  }else{
	    ncand++;
	  }
	  fTree->Fill();
	  //if( fRun>=160431 && fRun<=167913 ) cout<<fRun<<":"<<fLumiBlock<<":"<<fEvent<<" ="<<endl;
	  //if( fRun>=170826 && fRun<=172619 ) cout<<fRun<<":"<<fLumiBlock<<":"<<fEvent<<" ="<<endl;
	  //if( fRun>=172620 && fRun<=173692 ) cout<<fRun<<":"<<fLumiBlock<<":"<<fEvent<<" ="<<endl;
	  //if( fRun>=175860 && fRun<=180252 ) cout<<fRun<<":"<<fLumiBlock<<":"<<fEvent<<" ="<<endl;
	  //}
	  //		}
	}
      }
    }
  }
  
  // end loop over Eta_b 
  
  // goto next event

  
  // end loop over events
  
  h1->Fill((float) ncand); 
  h1->Write();

  fFile->Write();
  fFile->Close();

  delete fFile;
  
  
  cout << "-------------------------" << endl;
  cout << "Error     : " << nerr << endl;
  cout << "Trigger   : " << nerr1 << endl;
  cout << "All eta   : " << nerr3 << endl;
  cout << "Pvts>0.05 : " << nerr5 << endl;
  

  cout << "# Cands : " << nerr4 << endl;
  cout << "# Evts  : " << nev << endl;

  cout<<"himaltotal etab latest"<<totaletab<<endl;
  cout<<"himal total muon acceptance"<<muacc<<endl;

  
} // end macro
