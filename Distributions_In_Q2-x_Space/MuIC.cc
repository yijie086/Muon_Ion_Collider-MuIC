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
  pythia.settings.parm("PhaseSpace:Q2Min", 1);

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
  const float ycutmin=0.00;
  const float ycutmax=1.00;

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
  std::vector< float > genJetQ2;
  std::vector< float > genJetx;
  std::vector< int >   genJetChargedMultiplicity;
  std::vector< float > gendPhiej;
  
  std::vector< float > genoneJetPt;
  std::vector< float > genoneJetEta;
  std::vector< float > genoneJetx;
  std::vector< float > genoneJetQ2;
  
  std::vector< int > genJetNum; //**************
  std::vector< int > genParticleNum;
  std::vector< float > gendiJetdR;
  std::vector< float > gendijetdPhi;
  
  std::vector< float > F2x;
  std::vector< float > F2y;
  std::vector< float > F2Q2;
  std::vector< float > F2W2;
  
  std::vector< float > proton_eta;
  std::vector< float > proton_pt;
  std::vector< float > proton_Q2;
  std::vector< float > proton_x;
  std::vector< float > antiproton_eta;
  std::vector< float > antiproton_pt;
  std::vector< float > antiproton_Q2;
  std::vector< float > antiproton_x;
  
  std::vector< float > n0_eta;
  std::vector< float > n0_pt;
  std::vector< float > n0_Q2;
  std::vector< float > n0_x;
  std::vector< float > antin0_eta;
  std::vector< float > antin0_pt;
  std::vector< float > antin0_Q2;
  std::vector< float > antin0_x;
  
  std::vector< float > pionp_eta;
  std::vector< float > pionp_pt;
  std::vector< float > pionp_Q2;
  std::vector< float > pionp_x;
  std::vector< float > pionn_eta;
  std::vector< float > pionn_pt;
  std::vector< float > pionn_Q2;
  std::vector< float > pionn_x;
  
  std::vector< float > Kp_eta;
  std::vector< float > Kp_pt;
  std::vector< float > Kp_Q2;
  std::vector< float > Kp_x;
  std::vector< float > Kn_eta;
  std::vector< float > Kn_pt;
  std::vector< float > Kn_Q2;
  std::vector< float > Kn_x;
  std::vector< float > KL_eta;
  std::vector< float > KL_pt;
  std::vector< float > KL_Q2;
  std::vector< float > KL_x;
  
  std::vector< float > gamma_eta;
  std::vector< float > gamma_pt;
  std::vector< float > gamma_Q2;
  std::vector< float > gamma_x;

  
  
  //std::vector< std::vector<int> > gendau_chg;
  //std::vector< std::vector<int> > gendau_pid;
  //std::vector< std::vector<float> > gendau_pt;
  //std::vector< std::vector<float> > gendau_eta;
  //std::vector< std::vector<float> > gendau_phi;

  TFile * f = TFile::Open("MuIC111.root","recreate");
  TTree * trackTree = new TTree("trackTree","v1");

  trackTree->Branch("proton_eta",&proton_eta);
  trackTree->Branch("proton_pt",&proton_pt);
  trackTree->Branch("proton_Q2",&proton_Q2);
  trackTree->Branch("proton_x",&proton_x);
  trackTree->Branch("antiproton_eta",&antiproton_eta);
  trackTree->Branch("antiproton_pt",&antiproton_pt);
  trackTree->Branch("antiproton_Q2",&antiproton_Q2);
  trackTree->Branch("antiproton_x",&antiproton_x);
  
  trackTree->Branch("n0_eta",&n0_eta);
  trackTree->Branch("n0_pt",&n0_pt);
  trackTree->Branch("n0_Q2",&n0_Q2);
  trackTree->Branch("n0_x",&n0_x);
  trackTree->Branch("antin0_eta",&antin0_eta);
  trackTree->Branch("antin0_pt",&antin0_pt);
  trackTree->Branch("antin0_Q2",&antin0_Q2);
  trackTree->Branch("antin0_x",&antin0_x);
  
  trackTree->Branch("pionp_eta",&pionp_eta);
  trackTree->Branch("pionp_pt",&pionp_pt);
  trackTree->Branch("pionp_Q2",&pionp_Q2);
  trackTree->Branch("pionp_x",&pionp_x);
  trackTree->Branch("pionn_eta",&pionn_eta);
  trackTree->Branch("pionn_pt",&pionn_pt);
  trackTree->Branch("pionn_Q2",&pionn_Q2);
  trackTree->Branch("pionn_x",&pionn_x);
  
  trackTree->Branch("Kp_eta",&Kp_eta);
  trackTree->Branch("Kp_pt",&Kp_pt);
  trackTree->Branch("Kp_Q2",&Kp_Q2);
  trackTree->Branch("Kp_x",&Kp_x);
  trackTree->Branch("Kn_eta",&Kn_eta);
  trackTree->Branch("Kn_pt",&Kn_pt);
  trackTree->Branch("Kn_Q2",&Kn_Q2);
  trackTree->Branch("Kn_x",&Kn_x);
  trackTree->Branch("KL_eta",&KL_eta);
  trackTree->Branch("KL_pt",&KL_pt);
  trackTree->Branch("KL_Q2",&KL_Q2);
  trackTree->Branch("KL_x",&KL_x);
  
  trackTree->Branch("gamma_eta",&gamma_eta);
  trackTree->Branch("gamma_pt",&gamma_pt);
  trackTree->Branch("gamma_Q2",&gamma_Q2);
  trackTree->Branch("gamma_x",&gamma_x);

  
  trackTree->Branch("Q2",&Q2);
  trackTree->Branch("W2",&W2);
  trackTree->Branch("x",&x);
  trackTree->Branch("y",&y);
  
  trackTree->Branch("F2Q2",&F2Q2);
  trackTree->Branch("F2x",&F2x);
  trackTree->Branch("F2y",&F2y);
  trackTree->Branch("F2W2",&F2W2);
  
  trackTree->Branch("genJetEta",&genJetEta);
  trackTree->Branch("genJetPt",&genJetPt);
  trackTree->Branch("genJetPhi",&genJetPhi);
  trackTree->Branch("genJetQ2",&genJetQ2);
  trackTree->Branch("genJetx",&genJetx);
  trackTree->Branch("genJetChargedMultiplicity",&genJetChargedMultiplicity);
  
  trackTree->Branch("genoneJetPt",&genoneJetPt);
  trackTree->Branch("genoneJetEta",&genoneJetEta);
  trackTree->Branch("genoneJetx",&genoneJetx);
  trackTree->Branch("genoneJetQ2",&genoneJetQ2);
  
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
  int nEvent = 2000*zoombeta;
  
for (int iQ=0;iQ<=0;iQ++) {
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
        
        //proton
        if (pythia.event[i].id()==2212) {
          proton_eta.push_back(pythia.event[i].eta());
          proton_pt.push_back(pythia.event[i].pT());
          proton_Q2.push_back(Q2);
          proton_x.push_back(x);
        }
        //antiproton
        else if (pythia.event[i].id()==-2212) {
          antiproton_eta.push_back(pythia.event[i].eta());
          antiproton_pt.push_back(pythia.event[i].pT());
          antiproton_Q2.push_back(Q2);
          antiproton_x.push_back(x);
        }
        //n0
        else if (pythia.event[i].id()==2112) {
          n0_eta.push_back(pythia.event[i].eta());
          n0_pt.push_back(pythia.event[i].pT());
          n0_Q2.push_back(Q2);
          n0_x.push_back(x);
        }
        //antin0
        else if (pythia.event[i].id()==-2112) {
          antin0_eta.push_back(pythia.event[i].eta());
          antin0_pt.push_back(pythia.event[i].pT());
          antin0_Q2.push_back(Q2);
          antin0_x.push_back(x);
        }
        //pionp
        else if (pythia.event[i].id()==211) {
          pionp_eta.push_back(pythia.event[i].eta());
          pionp_pt.push_back(pythia.event[i].pT());
          pionp_Q2.push_back(Q2);
          pionp_x.push_back(x);
        }
        //pionn
        else if (pythia.event[i].id()==-211) {
          pionn_eta.push_back(pythia.event[i].eta());
          pionn_pt.push_back(pythia.event[i].pT());
          pionn_Q2.push_back(Q2);
          pionn_x.push_back(x);
        }
        //Kp
        else if (pythia.event[i].id()==321) {
          Kp_eta.push_back(pythia.event[i].eta());
          Kp_pt.push_back(pythia.event[i].pT());
          Kp_Q2.push_back(Q2);
          Kp_x.push_back(x);
        }
        //Kn
        else if (pythia.event[i].id()==-321) {
          Kn_eta.push_back(pythia.event[i].eta());
          Kn_pt.push_back(pythia.event[i].pT());
          Kn_Q2.push_back(Q2);
          Kn_x.push_back(x);
        }
        //KL
        else if (pythia.event[i].id()==130) {
          KL_eta.push_back(pythia.event[i].eta());
          KL_pt.push_back(pythia.event[i].pT());
          KL_Q2.push_back(Q2);
          KL_x.push_back(x);
        }
        //gamma
        else if (pythia.event[i].id()==22) {
          gamma_eta.push_back(pythia.event[i].eta());
          gamma_pt.push_back(pythia.event[i].pT());
          gamma_Q2.push_back(Q2);
          gamma_x.push_back(x);
        }
        
        if( pythia.event[i].isCharged() ) {
          multiplicity++;
          chargedparticles.push_back( pj );
        }
      }
    }
    
    F2x.push_back(x);
    F2y.push_back(y);
    F2Q2.push_back(Q2);
    F2W2.push_back(W2);
    float eventeta=pythia.event[6].eta();
    float eventnum2=pythia.event[6].pT();
    h2dformuouteta2e->Fill(eventeta,eventnum2);
    h2dformuouteta2e_eta.push_back(eventeta);
    h2dformuouteta2e_pt.push_back(eventnum2);
      
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
      genJetQ2.push_back(Q2);
      genJetx.push_back(x);

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
        genoneJetEta.push_back(jets[i].eta());
        genoneJetx.push_back(x);
        genoneJetQ2.push_back(Q2);
        
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
        
        
        
        float jeteta=jets[0].eta();
        float jetnum2=jets[0].pt();
        jh2dformuouteta2e->Fill(jeteta,jetnum2);
        jh2dformuouteta2e_eta.push_back(jeteta);
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
    genJetx.clear();
    genJetQ2.clear();
    genJetChargedMultiplicity.clear();
    
    genoneJetPt.clear();
    genoneJetEta.clear();
    genoneJetQ2.clear();
    genoneJetx.clear();
    
    F2Q2.clear();
    F2W2.clear();
    F2x.clear();
    F2y.clear();
    
    genJetNum.clear();  //**************
    gendiJetdR.clear();
    genParticleNum.clear();
    gendijetdPhi.clear();
    gendPhiej.clear();
    
    proton_eta.clear();
    proton_pt.clear();
    proton_Q2.clear();
    proton_x.clear();
    antiproton_eta.clear();
    antiproton_pt.clear();
    antiproton_Q2.clear();
    antiproton_x.clear();
    
    n0_eta.clear();
    n0_pt.clear();
    n0_Q2.clear();
    n0_x.clear();
    antin0_eta.clear();
    antin0_pt.clear();
    antin0_Q2.clear();
    antin0_x.clear();
    
    pionp_eta.clear();
    pionp_pt.clear();
    pionp_Q2.clear();
    pionp_x.clear();
    pionn_eta.clear();
    pionn_pt.clear();
    pionn_Q2.clear();
    pionn_x.clear();
    
    Kp_eta.clear();
    Kp_pt.clear();
    Kp_Q2.clear();
    Kp_x.clear();
    Kn_eta.clear();
    Kn_pt.clear();
    Kn_Q2.clear();
    Kn_x.clear();
    KL_eta.clear();
    KL_pt.clear();
    KL_Q2.clear();
    KL_x.clear();
    
    gamma_eta.clear();
    gamma_pt.clear();
    gamma_Q2.clear();
    gamma_x.clear();
    
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
    