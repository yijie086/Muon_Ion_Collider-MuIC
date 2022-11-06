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

#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif

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
  Pythia8::Pythia8ToHepMC MuICHepMC("MuIC-ncbkg.hepmc");
  
  Pythia pythia;
  Event& event = pythia.event;



  pythia.readString("Beams:frameType = 4");
  
  pythia.readString("Beams:LHEF = ncbkg_events.lhe" );
  // Specify one must read inputs from the MadGraph banner.
  pythia.readString("JetMatching:setMad=off");
  // Disable match in and out
  pythia.readString("LesHouches:matchInOut = off");

  

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
  
  // For consistency with MG5 Studies
  // Can be commented out
  pythia.readString("PDF:pSet=LHAPDF6:PDF4LHC15_nlo_mc_pdfas");
  
  // These options are sometimes required when reading from LHE File.
  // Seems to not be needed when simulating directly from pythia,
  // but sometimes needed when reading from LHE File
  // Can be commented out
  pythia.readString("BeamRemnants:remnantMode=1");
  pythia.readString("ColourReconnection:mode=1");
  
  
  //cout<<"sadkljkfghodsalkfhdlk"<<endl;
  
  pythia.init();

  //cout<<"1132145678985764357687980"<<endl;
  //*********************************************** ROOT SETUP  **************************************************

  int nEvent = 10000;
  int events_processed = 0;
  
  while (pythia.next() && events_processed < nEvent) {
    //pythia.event.list();
      
    auto &proton_in = event[1];
    auto &muon_in = event[4];
    auto &muon_out = event[6];
    
    auto p_proton_in = proton_in.p();
    auto p_muon_in = muon_in.p();
    auto p_muon_out = muon_out.p();
    auto p_virtual = p_muon_in - p_muon_out;
    
    // Calculate Q2, W2, Bjorken x, y.
    double Q2 = -p_virtual.m2Calc();
    double W2 = (p_proton_in + p_virtual).m2Calc();
    double x = Q2 / (2. * p_proton_in * p_virtual);
    double y = (p_proton_in * p_virtual) / (p_proton_in * p_muon_in);
    
    if (Q2>400) {
      MuICHepMC.writeNextEvent( pythia );
    }
    events_processed++;
  }
  return 0;
}
    
