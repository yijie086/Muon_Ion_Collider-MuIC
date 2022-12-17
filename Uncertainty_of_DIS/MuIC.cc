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
  const float etacut1=5;
  const float etacut2=2.4;
  const float Parptcut=0.2;
  const float ycutmin=0.00;
  const float ycutmax=1.0;

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
  std::vector< float > jh2dformuouteta2e_e;
  
  std::vector< float > ph2dformuouteta2e_eta;
  std::vector< float > ph2dformuouteta2e_pt;
  std::vector< float > ph2dformuouteta2e_e;

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
  std::vector< float > F2W2;
  
  std::vector< float > jF2x;
  std::vector< float > jF2y;
  std::vector< float > jF2Q2;
  
  std::vector< float > proton_eta;
  std::vector< float > proton_phi;
  std::vector< float > proton_pt;
  std::vector< float > proton_e;
  std::vector< float > proton_Q2;
  std::vector< float > proton_x;
  std::vector< float > proton_num;
  double proton_numt=0;
  std::vector< float > antiproton_eta;
  std::vector< float > antiproton_phi;
  std::vector< float > antiproton_pt;
  std::vector< float > antiproton_e;
  std::vector< float > antiproton_Q2;
  std::vector< float > antiproton_x;
  std::vector< float > antiproton_num;
  double antiproton_numt=0;
  
  std::vector< float > n0_eta;
  std::vector< float > n0_phi;
  std::vector< float > n0_pt;
  std::vector< float > n0_e;
  std::vector< float > n0_Q2;
  std::vector< float > n0_x;
  std::vector< float > n0_num;
  double n0_numt=0;
  std::vector< float > antin0_eta;
  std::vector< float > antin0_phi;
  std::vector< float > antin0_pt;
  std::vector< float > antin0_e;
  std::vector< float > antin0_Q2;
  std::vector< float > antin0_x;
  std::vector< float > antin0_num;
  double antin0_numt=0;
  
  std::vector< float > pionp_eta;
  std::vector< float > pionp_phi;
  std::vector< float > pionp_pt;
  std::vector< float > pionp_e;
  std::vector< float > pionp_Q2;
  std::vector< float > pionp_x;
  std::vector< float > pionp_num;
  double pionp_numt=0;
  std::vector< float > pionn_eta;
  std::vector< float > pionn_phi;
  std::vector< float > pionn_pt;
  std::vector< float > pionn_e;
  std::vector< float > pionn_Q2;
  std::vector< float > pionn_x;
  std::vector< float > pionn_num;
  double pionn_numt=0;
  
  std::vector< float > Kp_eta;
  std::vector< float > Kp_phi;
  std::vector< float > Kp_pt;
  std::vector< float > Kp_e;
  std::vector< float > Kp_Q2;
  std::vector< float > Kp_x;
  std::vector< float > Kp_num;
  double Kp_numt=0;
  std::vector< float > Kn_eta;
  std::vector< float > Kn_phi;
  std::vector< float > Kn_pt;
  std::vector< float > Kn_e;
  std::vector< float > Kn_Q2;
  std::vector< float > Kn_x;
  std::vector< float > Kn_num;
  double Kn_numt=0;
  std::vector< float > KL_eta;
  std::vector< float > KL_phi;
  std::vector< float > KL_pt;
  std::vector< float > KL_e;
  std::vector< float > KL_Q2;
  std::vector< float > KL_x;
  std::vector< float > KL_num;
  double KL_numt=0;
  
  std::vector< float > gamma_eta;
  std::vector< float > gamma_phi;
  std::vector< float > gamma_pt;
  std::vector< float > gamma_e;
  std::vector< float > gamma_Q2;
  std::vector< float > gamma_x;
  std::vector< float > gamma_num;
  double gamma_numt=0;

  std::vector< float > muon_eta;
  std::vector< float > muon_phi;
  std::vector< float > muon_pt;
  std::vector< float > muon_e;
  std::vector< float > muon_Q2;
  std::vector< float > muon_x;
  
  std::vector< float > ep_eta;
  std::vector< float > ep_phi;
  std::vector< float > ep_pt;
  std::vector< float > ep_e;
  std::vector< float > ep_Q2;
  std::vector< float > ep_x;
  std::vector< float > ep_num;
  double ep_numt=0;
  std::vector< float > en_eta;
  std::vector< float > en_phi;
  std::vector< float > en_pt;
  std::vector< float > en_e;
  std::vector< float > en_Q2;
  std::vector< float > en_x;
  std::vector< float > en_num;
  double en_numt=0;
  
  //std::vector< std::vector<int> > gendau_chg;
  //std::vector< std::vector<int> > gendau_pid;
  //std::vector< std::vector<float> > gendau_pt;
  //std::vector< std::vector<float> > gendau_eta;
  //std::vector< std::vector<float> > gendau_phi;

  TFile * f = TFile::Open("MuIC.root","recreate");
  TTree * trackTree = new TTree("trackTree","v1");

  trackTree->Branch("muon_eta",&muon_eta);
  trackTree->Branch("muon_phi",&muon_phi);
  trackTree->Branch("muon_pt",&muon_pt);
  trackTree->Branch("muon_e",&muon_e);
  trackTree->Branch("muon_Q2",&muon_Q2);
  trackTree->Branch("muon_x",&muon_x);

  trackTree->Branch("proton_eta",&proton_eta);
  trackTree->Branch("proton_phi",&proton_phi);
  trackTree->Branch("proton_pt",&proton_pt);
  trackTree->Branch("proton_e",&proton_e);
  trackTree->Branch("proton_Q2",&proton_Q2);
  trackTree->Branch("proton_x",&proton_x);
  trackTree->Branch("proton_num",&proton_num);
  trackTree->Branch("antiproton_eta",&antiproton_eta);
  trackTree->Branch("antiproton_phi",&antiproton_phi);
  trackTree->Branch("antiproton_pt",&antiproton_pt);
  trackTree->Branch("antiproton_e",&antiproton_e);
  trackTree->Branch("antiproton_Q2",&antiproton_Q2);
  trackTree->Branch("antiproton_x",&antiproton_x);
  trackTree->Branch("antiproton_num",&antiproton_num);
  
  trackTree->Branch("n0_eta",&n0_eta);
  trackTree->Branch("n0_phi",&n0_phi);
  trackTree->Branch("n0_pt",&n0_pt);
  trackTree->Branch("n0_e",&n0_e);
  trackTree->Branch("n0_Q2",&n0_Q2);
  trackTree->Branch("n0_x",&n0_x);
  trackTree->Branch("n0_num",&n0_num);
  trackTree->Branch("antin0_eta",&antin0_eta);
  trackTree->Branch("antin0_phi",&antin0_phi);
  trackTree->Branch("antin0_pt",&antin0_pt);
  trackTree->Branch("antin0_e",&antin0_e);
  trackTree->Branch("antin0_Q2",&antin0_Q2);
  trackTree->Branch("antin0_x",&antin0_x);
  trackTree->Branch("antin0_num",&antin0_num);
  
  trackTree->Branch("pionp_eta",&pionp_eta);
  trackTree->Branch("pionp_phi",&pionp_phi);
  trackTree->Branch("pionp_pt",&pionp_pt);
  trackTree->Branch("pionp_e",&pionp_e);
  trackTree->Branch("pionp_Q2",&pionp_Q2);
  trackTree->Branch("pionp_x",&pionp_x);
  trackTree->Branch("pionp_num",&pionp_num);
  trackTree->Branch("pionn_eta",&pionn_eta);
  trackTree->Branch("pionn_phi",&pionn_phi);
  trackTree->Branch("pionn_pt",&pionn_pt);
  trackTree->Branch("pionn_e",&pionn_e);
  trackTree->Branch("pionn_Q2",&pionn_Q2);
  trackTree->Branch("pionn_x",&pionn_x);
  trackTree->Branch("pionn_num",&pionn_num);
  
  trackTree->Branch("Kp_eta",&Kp_eta);
  trackTree->Branch("Kp_phi",&Kp_phi);
  trackTree->Branch("Kp_pt",&Kp_pt);
  trackTree->Branch("Kp_e",&Kp_e);
  trackTree->Branch("Kp_Q2",&Kp_Q2);
  trackTree->Branch("Kp_x",&Kp_x);
  trackTree->Branch("Kp_num",&Kp_num);
  trackTree->Branch("Kn_eta",&Kn_eta);
  trackTree->Branch("Kn_phi",&Kn_phi);
  trackTree->Branch("Kn_pt",&Kn_pt);
  trackTree->Branch("Kn_e",&Kn_e);
  trackTree->Branch("Kn_Q2",&Kn_Q2);
  trackTree->Branch("Kn_x",&Kn_x);
  trackTree->Branch("Kn_num",&Kn_num);
  trackTree->Branch("KL_eta",&KL_eta);
  trackTree->Branch("KL_phi",&KL_phi);
  trackTree->Branch("KL_pt",&KL_pt);
  trackTree->Branch("KL_e",&KL_e);
  trackTree->Branch("KL_Q2",&KL_Q2);
  trackTree->Branch("KL_x",&KL_x);
  trackTree->Branch("KL_num",&KL_num);
  
  trackTree->Branch("gamma_eta",&gamma_eta);
  trackTree->Branch("gamma_phi",&gamma_phi);
  trackTree->Branch("gamma_pt",&gamma_pt);
  trackTree->Branch("gamma_e",&gamma_e);
  trackTree->Branch("gamma_Q2",&gamma_Q2);
  trackTree->Branch("gamma_x",&gamma_x);
  trackTree->Branch("gamma_num",&gamma_num);

  trackTree->Branch("ep_eta",&ep_eta);
  trackTree->Branch("ep_phi",&ep_phi);
  trackTree->Branch("ep_pt",&ep_pt);
  trackTree->Branch("ep_e",&ep_e);
  trackTree->Branch("ep_Q2",&ep_Q2);
  trackTree->Branch("ep_x",&ep_x);
  trackTree->Branch("ep_num",&ep_num);
  trackTree->Branch("en_eta",&en_eta);
  trackTree->Branch("en_phi",&en_phi);
  trackTree->Branch("en_pt",&en_pt);
  trackTree->Branch("en_e",&en_e);
  trackTree->Branch("en_Q2",&en_Q2);
  trackTree->Branch("en_x",&en_x);
  trackTree->Branch("en_num",&en_num);
  
  trackTree->Branch("Q2",&Q2);
  trackTree->Branch("W2",&W2);
  trackTree->Branch("x",&x);
  trackTree->Branch("y",&y);
  
  trackTree->Branch("F2Q2",&F2Q2);
  trackTree->Branch("F2x",&F2x);
  trackTree->Branch("F2y",&F2y);
  trackTree->Branch("F2W2",&F2W2);
  
  trackTree->Branch("jF2Q2",&jF2Q2);
  trackTree->Branch("jF2x",&jF2x);
  trackTree->Branch("jF2y",&jF2y);
  
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
  trackTree->Branch("jh2dformuouteta2e_e",&jh2dformuouteta2e_e);
  
  trackTree->Branch("ph2dformuouteta2e_eta",&ph2dformuouteta2e_eta);
  trackTree->Branch("ph2dformuouteta2e_pt",&ph2dformuouteta2e_pt);
  trackTree->Branch("ph2dformuouteta2e_e",&ph2dformuouteta2e_e);

  
  //******************************************** ANALYZER ********************************************************

  //jet clustering
  // choose a jet definition
  double R = 1;
  JetDefinition jet_def(antikt_algorithm, R);

  // Begin event loop.
  int zoombeta=1;
  int nEvent = 5000*zoombeta;
  
for (int iQ=0;iQ<=5;iQ++) {
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
    for (int i = 7; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal()&&((pythia.event[i].id()!=13)||(pythia.event[i].id()==13&&pythia.event[pythia.event[i].mother1()].id()!=13))){
        //basic kinematic cuts
        if(  pythia.event[i].eta() < -etacut1 || pythia.event[i].eta() > etacut2 || pythia.event[i].pT() < Parptcut){
          continue;
        }
        
        PseudoJet pj = PseudoJet(   pythia.event[i].px(),  pythia.event[i].py(),  pythia.event[i].pz(), pythia.event[i].e() );
        pj.set_user_index(i);
        particles.push_back( pj );
        
        
        
        
        //proton
        if (pythia.event[i].id()==2212) {
          proton_eta.push_back(pythia.event[i].eta());
          proton_e.push_back(pythia.event[i].e());
          proton_phi.push_back(pythia.event[i].phi());
          proton_pt.push_back(pythia.event[i].pT());
          proton_Q2.push_back(Q2);
          proton_x.push_back(x);
          proton_numt++;
        }
        //antiproton
        else if (pythia.event[i].id()==-2212) {
          antiproton_eta.push_back(pythia.event[i].eta());
          antiproton_e.push_back(pythia.event[i].e());
          antiproton_phi.push_back(pythia.event[i].phi());
          antiproton_pt.push_back(pythia.event[i].pT());
          antiproton_Q2.push_back(Q2);
          antiproton_x.push_back(x);
          antiproton_numt++;
        }
        //n0
        else if (pythia.event[i].id()==2112) {
          n0_eta.push_back(pythia.event[i].eta());
          n0_e.push_back(pythia.event[i].e());
          n0_phi.push_back(pythia.event[i].phi());
          n0_pt.push_back(pythia.event[i].pT());
          n0_Q2.push_back(Q2);
          n0_x.push_back(x);
          n0_numt++;
        }
        //antin0
        else if (pythia.event[i].id()==-2112) {
          antin0_eta.push_back(pythia.event[i].eta());
          antin0_e.push_back(pythia.event[i].e());
          antin0_phi.push_back(pythia.event[i].phi());
          antin0_pt.push_back(pythia.event[i].pT());
          antin0_Q2.push_back(Q2);
          antin0_x.push_back(x);
          antin0_numt++;
        }
        //pionp
        else if (pythia.event[i].id()==211) {
          pionp_eta.push_back(pythia.event[i].eta());
          pionp_e.push_back(pythia.event[i].e());
          pionp_phi.push_back(pythia.event[i].phi());
          pionp_pt.push_back(pythia.event[i].pT());
          pionp_Q2.push_back(Q2);
          pionp_x.push_back(x);
          pionp_numt++;
        }
        //pionn
        else if (pythia.event[i].id()==-211) {
          pionn_eta.push_back(pythia.event[i].eta());
          pionn_e.push_back(pythia.event[i].e());
          pionn_phi.push_back(pythia.event[i].phi());
          pionn_pt.push_back(pythia.event[i].pT());
          pionn_Q2.push_back(Q2);
          pionn_x.push_back(x);
          pionn_numt++;
        }
        //Kp
        else if (pythia.event[i].id()==321) {
          Kp_eta.push_back(pythia.event[i].eta());
          Kp_e.push_back(pythia.event[i].e());
          Kp_phi.push_back(pythia.event[i].phi());
          Kp_pt.push_back(pythia.event[i].pT());
          Kp_Q2.push_back(Q2);
          Kp_x.push_back(x);
          Kp_numt++;
        }
        //Kn
        else if (pythia.event[i].id()==-321) {
          Kn_eta.push_back(pythia.event[i].eta());
          Kn_e.push_back(pythia.event[i].e());
          Kn_phi.push_back(pythia.event[i].phi());
          Kn_pt.push_back(pythia.event[i].pT());
          Kn_Q2.push_back(Q2);
          Kn_x.push_back(x);
          Kn_numt++;
        }
        //KL
        else if (pythia.event[i].id()==130) {
          KL_eta.push_back(pythia.event[i].eta());
          KL_e.push_back(pythia.event[i].e());
          KL_phi.push_back(pythia.event[i].phi());
          KL_pt.push_back(pythia.event[i].pT());
          KL_Q2.push_back(Q2);
          KL_x.push_back(x);
          KL_numt++;
        }
        //gamma
        else if (pythia.event[i].id()==22) {
          gamma_eta.push_back(pythia.event[i].eta());
          gamma_e.push_back(pythia.event[i].e());
          gamma_phi.push_back(pythia.event[i].phi());
          gamma_pt.push_back(pythia.event[i].pT());
          gamma_Q2.push_back(Q2);
          gamma_x.push_back(x);
          gamma_numt++;
        }
        else if (pythia.event[i].id()==11) {
          ep_eta.push_back(pythia.event[i].eta());
          ep_e.push_back(pythia.event[i].e());
          ep_phi.push_back(pythia.event[i].phi());
          ep_pt.push_back(pythia.event[i].pT());
          ep_Q2.push_back(Q2);
          ep_x.push_back(x);
          ep_numt++;
        }
        else if (pythia.event[i].id()==-11) {
          en_eta.push_back(pythia.event[i].eta());
          en_e.push_back(pythia.event[i].e());
          en_phi.push_back(pythia.event[i].phi());
          en_pt.push_back(pythia.event[i].pT());
          en_Q2.push_back(Q2);
          en_x.push_back(x);
          en_numt++;
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
    
    muon_eta.push_back(pythia.event[6].eta());
    muon_phi.push_back(pythia.event[6].phi());
    muon_pt.push_back(pythia.event[6].pT());
    muon_e.push_back(pythia.event[6].e());
    muon_Q2.push_back(Q2);
    muon_x.push_back(x);
      
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
        
        jF2x.push_back(x);
        jF2y.push_back(y);
        jF2Q2.push_back(Q2);
        
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
        
        gendPhiej.push_back(TMath::Abs(Absdphi(jets[0].phi(),pythia.event[6].phi())-pi));
        h1dPhiej->Fill(TMath::Abs(Absdphi(jets[0].phi(),pythia.event[6].phi())-pi));
        
        
        //$$$$$$$$$$$$$$$$$$$$$$QXcut2squre$$$$$$$$$$$$$$$$$$$$$$$
        
        float jeteta=jets[0].eta();
        float jetnum2=jets[0].pt();
        float jetnum3=jets[0].e();
        
        jh2dformuouteta2e->Fill(jeteta,jetnum2);
        jh2dformuouteta2e_eta.push_back(jeteta);
        jh2dformuouteta2e_pt.push_back(jetnum2);
        jh2dformuouteta2e_e.push_back(jetnum3);
        
        float peta=pythia.event[5].eta();
        float pnum2=pythia.event[5].pT();
        float pnum3=pythia.event[5].e();
        
        ph2dformuouteta2e_eta.push_back(peta);
        ph2dformuouteta2e_pt.push_back(pnum2);
        ph2dformuouteta2e_e.push_back(pnum3);
        
        
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
    proton_num.push_back(proton_numt);
    antiproton_num.push_back(antiproton_numt);
    n0_num.push_back(n0_numt);
    antin0_num.push_back(antin0_numt);
    pionp_num.push_back(pionp_numt);
    pionn_num.push_back(pionn_numt);
    Kp_num.push_back(Kp_numt);
    Kn_num.push_back(Kn_numt);
    ep_num.push_back(ep_numt);
    en_num.push_back(en_numt);
    KL_num.push_back(KL_numt);
    gamma_num.push_back(gamma_numt);
    trackTree->Fill();
    x=0;
    y=0;
    Q2=0;
    W2=0;
    Qre=0;
    Qp=0;
    
    proton_numt=0;
    antiproton_numt=0;
    n0_numt=0;
    antin0_numt=0;
    pionp_numt=0;
    pionn_numt=0;
    Kp_numt=0;
    Kn_numt=0;
    KL_numt=0;
    gamma_numt=0;
    ep_numt=0;
    en_numt=0;
    
    genJetPt.clear();
    genJetEta.clear();
    genJetPhi.clear();
    genJetChargedMultiplicity.clear();
    genoneJetPt.clear();
    
    F2Q2.clear();
    F2W2.clear();
    F2x.clear();
    F2y.clear();
    
    jF2Q2.clear();
    jF2x.clear();
    jF2y.clear();
    
    genJetNum.clear();  //**************
    gendiJetdR.clear();
    genParticleNum.clear();
    gendijetdPhi.clear();
    gendPhiej.clear();
    
    proton_eta.clear();
    proton_e.clear();
    proton_phi.clear();
    proton_pt.clear();
    proton_Q2.clear();
    proton_x.clear();
    proton_num.clear();
    antiproton_eta.clear();
    antiproton_e.clear();
    antiproton_phi.clear();
    antiproton_pt.clear();
    antiproton_Q2.clear();
    antiproton_x.clear();
    antiproton_num.clear();
    
    n0_eta.clear();
    n0_e.clear();
    n0_phi.clear();
    n0_pt.clear();
    n0_Q2.clear();
    n0_x.clear();
    n0_num.clear();
    antin0_eta.clear();
    antin0_e.clear();
    antin0_phi.clear();
    antin0_pt.clear();
    antin0_Q2.clear();
    antin0_x.clear();
    antin0_num.clear();
    
    pionp_eta.clear();
    pionp_e.clear();
    pionp_phi.clear();
    pionp_pt.clear();
    pionp_Q2.clear();
    pionp_x.clear();
    pionp_num.clear();
    pionn_eta.clear();
    pionn_e.clear();
    pionn_phi.clear();
    pionn_pt.clear();
    pionn_Q2.clear();
    pionn_x.clear();
    pionn_num.clear();
    
    Kp_eta.clear();
    Kp_e.clear();
    Kp_phi.clear();
    Kp_pt.clear();
    Kp_Q2.clear();
    Kp_x.clear();
    Kp_num.clear();
    Kn_eta.clear();
    Kn_e.clear();
    Kn_phi.clear();
    Kn_pt.clear();
    Kn_Q2.clear();
    Kn_x.clear();
    Kn_num.clear();
    KL_eta.clear();
    KL_e.clear();
    KL_phi.clear();
    KL_pt.clear();
    KL_Q2.clear();
    KL_x.clear();
    KL_num.clear();
    
    gamma_eta.clear();
    gamma_e.clear();
    gamma_phi.clear();
    gamma_pt.clear();
    gamma_Q2.clear();
    gamma_x.clear();
    gamma_num.clear();
    
    ep_eta.clear();
    ep_e.clear();
    ep_phi.clear();
    ep_pt.clear();
    ep_Q2.clear();
    ep_x.clear();
    ep_num.clear();
    en_eta.clear();
    en_e.clear();
    en_phi.clear();
    en_pt.clear();
    en_Q2.clear();
    en_x.clear();
    en_num.clear();
    
    muon_eta.clear();
    muon_phi.clear();
    muon_pt.clear();
    muon_e.clear();
    muon_Q2.clear();
    muon_x.clear();
    
    h2dformuouteta2e_eta.clear();
    h2dformuouteta2e_pt.clear();
    
    jh2dformuouteta2e_eta.clear();
    jh2dformuouteta2e_pt.clear();
    jh2dformuouteta2e_e.clear();
    
    ph2dformuouteta2e_eta.clear();
    ph2dformuouteta2e_pt.clear();
    ph2dformuouteta2e_e.clear();

  }
}


    
  trackTree->Write();
  
  f->Close();
  return 0;
}
    