#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "Pythia8/Pythia.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>

using namespace Pythia8;
using namespace fastjet;
#include <iostream>

double Absdphi(double phi1, double phi2) {
  // ***** This is the function to calcuate |phi_l-phi_j-pi| ***** 
  double Pi=3.14159265;
  double temp1,temp2,temp3;
  temp1=TMath::Abs(phi1-phi2);
  temp2=TMath::Abs(phi1-phi2-2*Pi);
  temp3=TMath::Abs(phi1-phi2+2*Pi);
  temp1=TMath::Min(temp1,temp2);
  temp1=TMath::Min(temp1,temp3);
  return temp1;
}

int main() {
  // ***** This is the main function *****
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  //PYTHIA SETTINGS
  Pythia pythia;
  Event& event = pythia.event;

  //energy to energy frame
  pythia.readString("Beams:frameType = 2");
  
  //proton
  pythia.readString("Beams:idA = 2212");

  //muon
  pythia.readString("Beams:idB = 13");

  //specify the beam energies
  //muon ion collider: 960 GeV muon, 275 GeV proton
  pythia.readString("Beams:eA = 275");
  pythia.readString("Beams:eB = 960");
  
  //neutral current
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  //charged current
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  
  //minimum Q2 cut
  pythia.settings.parm("PhaseSpace:Q2Min", 1);

  //needed setting for DIS + showering
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");
  
  //Some effects' switch
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  
  pythia.init();

  // ***** ROOT SETUP  *****
  
  //pt of jet cuts, a jet with pt lower than this num[GeV/c] will not be recorded as a jet
  const float jetptcut=5;
  
  //eta cut from etacutn to etacutp
  const float etacutn=-5;
  const float etacutp=2.4;
  
  //only a particle has a pt bigger than Parptcut can be recorded
  const float Parptcut=0.2;
  
  //y cut due to the uncertainty results of DIS Q2 x y
  const float ycutmin=0.001;
  const float ycutmax=0.95;

  float pi=3.14159266;
  
  std::vector< float > h2dformuouteta2e_eta;
  std::vector< float > h2dformuouteta2e_pt;

  float Q2,W2,x,y;
  
  float tempjetnum;
  float tempdeltaR;
  float tempparticle;
  float tempphi;
  float tempphi1,tempeta1,tempdeltaR1;
  float tempphi2,tempeta2,tempdeltaR2;
  
  std::vector< float > genJetPt;
  std::vector< float > genJetEta;
  std::vector< float > genJetPhi;
  std::vector< int >   genJetChargedMultiplicity;
  std::vector< float > gendPhiej;
  std::vector< float > genoneJetPt;
  
  std::vector< int > genJetNum; //**************
  std::vector< int > genParticleNum;
  std::vector< float > gendiJetdR;
  std::vector< float > gendijetdPhi;
  
  std::vector< float > F2x;
  std::vector< float > F2y;
  std::vector< float > F2Q2;

  TFile * f = TFile::Open("MuIC.root","recreate");
  TTree * trackTree = new TTree("trackTree","v1");

  trackTree->Branch("Q2",&Q2);
  trackTree->Branch("W2",&W2);
  trackTree->Branch("x",&x);
  trackTree->Branch("y",&y);
  
  trackTree->Branch("F2Q2",&F2Q2);
  trackTree->Branch("F2x",&F2x);
  trackTree->Branch("F2y",&F2y);
  
  trackTree->Branch("genJetEta",&genJetEta);
  trackTree->Branch("genJetPt",&genJetPt);
  trackTree->Branch("genJetPhi",&genJetPhi);
  trackTree->Branch("genJetChargedMultiplicity",&genJetChargedMultiplicity);
  
  trackTree->Branch("genoneJetPt",&genoneJetPt);
  
  trackTree->Branch("genJetNum",&genJetNum);
  trackTree->Branch("genParticleNum",&genParticleNum); 
  trackTree->Branch("gendiJetdR",&gendiJetdR);
  trackTree->Branch("gendijetdPhi",&gendijetdPhi);
  trackTree->Branch("gendPhiej",&gendPhiej);
  
  // ***** ANALYZER *****

  //jet clustering
  // choose a jet definition and R
  double R = 1;
  JetDefinition jet_def(antikt_algorithm, R);

  // Begin event loop.
  int zoombeta=100;
  int nEvent = 500000*zoombeta;
  
  for (int iQ=0;iQ<=0;iQ++) {
    //This process allow you to combine differnet Q2 events
    pythia.settings.parm("PhaseSpace:Q2Min",TMath::Power(10,iQ));
    pythia.init();
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;
      // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
      Vec4 pProton = event[1].p();
      Vec4 peIn    = event[4].p();
      Vec4 peOut   = event[6].p();
      Vec4 pPhoton = peIn - peOut;
    
      Vec4 tempp1;

      Q2    = - pPhoton.m2Calc();
      W2    = (pProton + pPhoton).m2Calc();
      x     = Q2 / (2. * pProton * pPhoton);
      y     = (pProton * pPhoton) / (pProton * peIn);
      
      if(y<ycutmin || y>ycutmax) continue;
 
      //getting pythia particle info
      vector<PseudoJet> particles;
      vector<PseudoJet> antiparticles;
      vector<PseudoJet> chargedparticles;
      int multiplicity = 0;
      //start loop at i=7 because we want to exclude the scattered lepton (which is i=6 index)
      for (int i = 7; i < pythia.event.size(); ++i){
        if (pythia.event[i].isFinal()&&((pythia.event[i].id()!=13)||(pythia.event[i].id()==13&&pythia.event[pythia.event[i].mother1()].id()!=13))) {
          if( pythia.event[i].eta() < etacutn || pythia.event[i].eta() > etacutp || pythia.event[i].pT() < Parptcut){
            continue;
          }
        
          PseudoJet pj = PseudoJet(   pythia.event[i].px(),  pythia.event[i].py(),  pythia.event[i].pz(), pythia.event[i].e() );
          pj.set_user_index(i);
          particles.push_back( pj );
        
          if( pythia.event[i].isCharged() ) {
            multiplicity++;
            chargedparticles.push_back( pj );
          }
        }
      }
    
      // run the clustering, extract the jets
      ClusterSequence cs(particles, jet_def);
      vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

      //loop over jets
      tempjetnum=0;
      for (unsigned i = 0; i < jets.size(); i++) {
        if(jets[i].pt() < jetptcut ) continue;
        tempjetnum++;
        genJetPt.push_back(jets[i].pt());
        genJetEta.push_back(jets[i].eta());
        genJetPhi.push_back(jets[i].phi());

        //loop over constituents
        vector<PseudoJet> constituents = jets[i].constituents();
        int chMult = 0;
        for (unsigned j = 0; j < constituents.size(); j++) {
          if( ( pythia.event[ constituents[j].user_index() ] ).isCharged()){
            chMult++;
          }
        }
        genJetChargedMultiplicity.push_back(chMult);
        genParticleNum.push_back(constituents.size());
      
        if (i==0&&jets[1].pt() < jetptcut) {
        
          genoneJetPt.push_back(jets[i].pt());
        
          tempeta1=TMath::Abs(jets[0].eta()-pythia.event[5].eta());
          tempphi1=Absdphi(jets[0].phi(),pythia.event[5].phi());
          tempdeltaR1=sqrt(tempeta1*tempeta1+tempphi1*tempphi1);
          tempdeltaR=tempdeltaR1;
        
          gendiJetdR.push_back(tempdeltaR);
          gendPhiej.push_back(TMath::Abs(Absdphi(jets[0].phi(),pythia.event[6].phi())-pi));
        }
      
        if (i==1&&jets[2].pt() < jetptcut) {
          for (unsigned j=0;j<chargedparticles.size();j++) {
            tempphi=jets[0].phi();
            gendijetdPhi.push_back(Absdphi(tempphi,chargedparticles[j].phi()));
          }
        }
      }
      genJetNum.push_back(tempjetnum);

      trackTree->Fill();
      x=0;
      y=0;
      Q2=0;
      W2=0;
    
      genJetPt.clear();
      genJetEta.clear();
      genJetPhi.clear();
      genJetChargedMultiplicity.clear();
      genoneJetPt.clear();
    
      F2Q2.clear();
      F2x.clear();
      F2y.clear();
    
      genJetNum.clear();
      gendiJetdR.clear();
      genParticleNum.clear();
      gendijetdPhi.clear();
      gendPhiej.clear();
    }
  }
  trackTree->Write();
  f->Close();
  return 0;
}
    