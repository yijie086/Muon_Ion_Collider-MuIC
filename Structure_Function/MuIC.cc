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
  const float ycutmin=0.001;
  const float ycutmax=0.95;

  float pi;
  
  pi=3.14159265;

  float Q2,W2,x,y,tempjetnum;
  
  std::vector< int > genJetNum; //**************

  
  std::vector< float > F2x;
  std::vector< float > F2y;
  std::vector< float > F2Q2;
  


  TFile * f = TFile::Open("MuIC.root","recreate");
  TTree * trackTree = new TTree("trackTree","v1");


  
  //trackTree->Branch("Q2",&Q2);
  //trackTree->Branch("W2",&W2);
  //trackTree->Branch("x",&x);
  //trackTree->Branch("y",&y);
  
  trackTree->Branch("F2Q2",&F2Q2);
  trackTree->Branch("F2x",&F2x);
  trackTree->Branch("F2y",&F2y);
  
  //trackTree->Branch("genJetNum",&genJetNum);  //**************

  
  
  

  
  //******************************************** ANALYZER ********************************************************

  //jet clustering
  // choose a jet definition
  double R = 1;
  JetDefinition jet_def(antikt_algorithm, R);

  // Begin event loop.
  int zoombeta=20;
  int nEvent = 50000000*zoombeta;
  
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
 /*
    //getting pythia particle info
    vector<PseudoJet> particles;
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
        
      }
    }
    */
    F2x.push_back(x);
    F2y.push_back(y);
    F2Q2.push_back(Q2);
  /*    
    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
*/
    //loop over jets
    tempjetnum=0;
    /*
    for (unsigned i = 0; i < jets.size(); i++) {
      if(jets[i].pt() < jetptcut ) continue;
      tempjetnum++;
    }
    */
    genJetNum.push_back(tempjetnum);

    //fill and clear branches
    
    trackTree->Fill();
    x=0;
    y=0;
    Q2=0;
    W2=0;
    
    F2Q2.clear();
    F2x.clear();
    F2y.clear();
    
    genJetNum.clear();  //**************
    
  }
}
  
  gStyle->SetOptStat(0);
  
  
    
  trackTree->Write();
  
  f->Close();
  return 0;
}
    
