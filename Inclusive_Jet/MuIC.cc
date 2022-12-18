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
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

  //*********************************************** PYTHIA SETTINGS **************************************************
  Pythia pythia;
  Event& event = pythia.event;

  pythia.readString("Beams:frameType = 2"); 
  //proton
  pythia.readString("Beams:idA = 2212");
  //electron
  //pythia.readString("Beams:idB = 11");
  //pythia.readString("Beams:idB = 2212");
  //muon
  pythia.readString("Beams:idB = 13");

  //specify the beam energies
  //electron: 5-30 GeV for EIC
  //proton:   50-250 GeV for EIC
  //muon ion collider: 960 GeV muon, 275 GeV proton
  pythia.readString("Beams:eA = 275");//275
  //pythia.readString("Beams:eB = 18");//18
  //pythia.readString("Beams:eB = 250");//18
  pythia.readString("Beams:eB = 960");
  
  //neutral current
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  //charged current
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  
  //minimum Q2 cut
  pythia.settings.parm("PhaseSpace:Q2Min", 10000);
  const double Q2Max=10000;

  //needed setting for DIS + showering
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");
  
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  
  
  //cout<<"sadkljkfghodsalkfhdlk"<<endl;
  
  pythia.init();

  //cout<<"1132145678985764357687980"<<endl;
  //*********************************************** ROOT SETUP  **************************************************

  const float jetptcut=5;
  const float etacut=9;
  const float Parptcut=0.2;
  const float ycutmin=0.001;
  const float ycutmax=0.95;

  float pi;
  
  pi=3.14159266;
  
  TH1 *h1d = new TH1D("Jetnum","Jetnum",10,0,10);
  TH1 *h1dParnum = new TH1D("Parnum","Parnum",30,0,30);
  TH1 *h1dR = new TH1D("JetdR","JetdR",30,0,5);
  TH1 *h1dPhi = new TH1D("JetdPhi","JetdPhi",30,0,3.14159265);
  TH1 *h1dPhiej = new TH1D("dPhiej","dPhiej",30,0,0.5);

  TH2 *h2d = new TH2D("Qre2Qp","Qre2Qp",50,0,800,50,0,800);
  TH2 *h2dx = new TH2D("Qre2Qp","Qre2Qp",50,0,20,50,0,20);
  TH2 *h2dy = new TH2D("Qre2Qp","Qre2Qp",50,0,20,50,0,20);
  TH2 *h2dz = new TH2D("Qre2Qp","Qre2Qp",50,0,800,50,0,800);
  TH2 *h2d2 = new TH2D("eta2x","eta2x",100,-7,0,100,-10,10);
  TH2 *h2d3 = new TH2D("eta2eta","eta2eta",100,-10,10,100,-10,10);
  TH2 *h2d4 = new TH2D("pt2phi","pt2phi",100,0,pi,100,0,3);
  
  const float etarangemin=-8;
  const float etarangemax=2;
  
  const float binmax=500;
  
  TH2 *h2dformuouteta2e = new TH2D("mueta2E","mueta2E",100,etarangemin,etarangemax,100,0,binmax);
  std::vector< float > h2dformuouteta2e_eta;
  std::vector< float > h2dformuouteta2e_pt;
  
  const float jetarangemin=-7;
  const float jetarangemax=4;
  
  TH2 *jh2dformuouteta2e = new TH2D("mueta2E","mueta2E",100,jetarangemin,jetarangemax,100,0,binmax);
  std::vector< float > jh2dformuouteta2e_eta;
  std::vector< float > jh2dformuouteta2e_pt;
  

  float Q2,W2,x,y;
  
  float Qre,Qp,tempjetnum;
  
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
  

  
  
  //std::vector< std::vector<int> > gendau_chg;
  //std::vector< std::vector<int> > gendau_pid;
  //std::vector< std::vector<float> > gendau_pt;
  //std::vector< std::vector<float> > gendau_eta;
  //std::vector< std::vector<float> > gendau_phi;

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
  
  trackTree->Branch("genJetNum",&genJetNum);  //**************
  trackTree->Branch("genParticleNum",&genParticleNum); 
  trackTree->Branch("gendiJetdR",&gendiJetdR);
  trackTree->Branch("gendijetdPhi",&gendijetdPhi);
  trackTree->Branch("gendPhiej",&gendPhiej);
  
  
  trackTree->Branch("Qre",&Qre);
  trackTree->Branch("Qp",&Qp);
  
  trackTree->Branch("h2dformuouteta2e_eta",&h2dformuouteta2e_eta);
  trackTree->Branch("h2dformuouteta2e_pt",&h2dformuouteta2e_pt);


  trackTree->Branch("jh2dformuouteta2e_eta",&jh2dformuouteta2e_eta);
  trackTree->Branch("jh2dformuouteta2e_pt",&jh2dformuouteta2e_pt);

  
  //******************************************** ANALYZER ********************************************************

  //jet clustering
  // choose a jet definition
  double R = 1;
  JetDefinition jet_def(antikt_algorithm, R);

  // Begin event loop.
  int zoombeta=1;
  int nEvent = 1000000*zoombeta;
  
for (int iQ=4;iQ<=4;iQ++) {
  pythia.settings.parm("PhaseSpace:Q2Min",TMath::Power(10,iQ));
  pythia.init();
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
    //cout<<iEvent<<"!!!!!!!!!!!!!!!!!"<<endl;
    if (!pythia.next()) continue;
    //pythia.event.list();
    // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
    Vec4 pProton = event[1].p();
    Vec4 peIn    = event[4].p();
    Vec4 peOut   = event[6].p();
    Vec4 pPhoton = peIn - peOut;
    
    Vec4 tempp1;

    // Q2, W2, Bjorken x, y.
    Q2    = - pPhoton.m2Calc();
    W2    = (pProton + pPhoton).m2Calc();
    x     = Q2 / (2. * pProton * pPhoton);
    y     = (pProton * pPhoton) / (pProton * peIn);
    if(y<ycutmin || y>ycutmax) continue;
    //if(Q2>Q2Max) continue;
 
    //getting pythia particle info
    vector<PseudoJet> particles;
    vector<PseudoJet> antiparticles;
    vector<PseudoJet> chargedparticles;
    int multiplicity = 0;
    //start loop at i=7 because we want to exclude the scattered lepton (which is i=6 index)
    //TMath::Abs(pythia.event[i].pT()-pythia.event[6].pT())!=0
    for (int i = 7; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal()&&((pythia.event[i].id()!=13)||(pythia.event[i].id()==13&&pythia.event[pythia.event[i].mother1()].id()!=13))){
      //if (pythia.event[i].isFinal()&&pythia.event[i].id()!=11){
        //basic kinematic cuts
        if( TMath::Abs( pythia.event[i].eta() ) > etacut || pythia.event[i].pT() < Parptcut){
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
    
    F2x.push_back(x);
    F2y.push_back(y);
    F2Q2.push_back(Q2);
      
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
      
      genParticleNum.push_back(constituents.size()); //**************
      h1dParnum->Fill(constituents.size());
      
      if (i==0&&jets[1].pt() < jetptcut) {
        
        genoneJetPt.push_back(jets[i].pt());
        
        tempp1=pythia.event[5].p();
        Qp=tempp1.e();
        PseudoJet tempp2=jets[0];
        Qre=tempp2.e();
        h2d->Fill(Qp,Qre);
        
        Qp=tempp1.px();
        Qre=tempp2.px();
        h2dx->Fill(Qp,Qre);
        
        Qp=tempp1.py();
        Qre=tempp2.py();
        h2dy->Fill(Qp,Qre);
        
        Qp=-tempp1.pz();
        Qre=-tempp2.pz();
        h2dz->Fill(Qp,Qre);
        
        tempeta1=TMath::Abs(jets[0].eta()-pythia.event[5].eta());
        tempphi1=Absdphi(jets[0].phi(),pythia.event[5].phi());
        tempdeltaR1=sqrt(tempeta1*tempeta1+tempphi1*tempphi1);
        tempdeltaR=tempdeltaR1;
        
        gendiJetdR.push_back(tempdeltaR);
        h1dR->Fill(tempdeltaR);
        
        if (9<pythia.event[6].pT()&&pythia.event[6].pT()<11) {
          gendPhiej.push_back(TMath::Abs(Absdphi(jets[0].phi(),pythia.event[6].phi())-pi));
          h1dPhiej->Fill(TMath::Abs(Absdphi(jets[0].phi(),pythia.event[6].phi())-pi));
        }
        
        
        //$$$$$$$$$$$$$$$$$$$$$$QXcut2squre$$$$$$$$$$$$$$$$$$$$$$$
        float eventeta=pythia.event[6].eta();
        float eventnum2=pythia.event[6].pT();
        
        float jeteta=jets[0].eta();
        float jetnum2=jets[0].pt();
        
        h2dformuouteta2e->Fill(eventeta,eventnum2);
        jh2dformuouteta2e->Fill(jeteta,jetnum2);
        
        h2dformuouteta2e_eta.push_back(eventeta);
        jh2dformuouteta2e_eta.push_back(jeteta);
        h2dformuouteta2e_pt.push_back(eventnum2);
        jh2dformuouteta2e_pt.push_back(jetnum2);
        
        
      }
      
      if (i==1&&jets[2].pt() < jetptcut) {
        //pythia.event.list();
        h2d3->Fill(jets[0].eta(),jets[1].eta());
        h2d2->Fill(TMath::Log10(x),jets[0].eta());
        h2d2->Fill(TMath::Log10(x),jets[1].eta());
        for (unsigned j=0;j<chargedparticles.size();j++) {
          tempphi=jets[0].phi();
          gendijetdPhi.push_back(Absdphi(tempphi,chargedparticles[j].phi()));
          h1dPhi->Fill(Absdphi(tempphi,chargedparticles[j].phi()));
          h2d4->Fill(Absdphi(tempphi,chargedparticles[j].phi()),chargedparticles[j].pt());
        }
      }
      
    }
    genJetNum.push_back(tempjetnum);
    h1d->Fill(tempjetnum);

    //fill and clear branches
    
    trackTree->Fill();
    x=0;
    y=0;
    Q2=0;
    W2=0;
    Qre=0;
    Qp=0;
    
    genJetPt.clear();
    genJetEta.clear();
    genJetPhi.clear();
    genJetChargedMultiplicity.clear();
    genoneJetPt.clear();
    
    F2Q2.clear();
    F2x.clear();
    F2y.clear();
    
    genJetNum.clear();  //**************
    gendiJetdR.clear();
    genParticleNum.clear();
    gendijetdPhi.clear();
    gendPhiej.clear();
    
        h2dformuouteta2e_eta.clear();
    h2dformuouteta2e_pt.clear();
    
    jh2dformuouteta2e_eta.clear();
    jh2dformuouteta2e_pt.clear();

  }
}
  
    

    
  trackTree->Write();
  
  f->Close();
  return 0;
}
    