/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <string>
//#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCutG.h"
#include "TCanvas.h"
#include <numeric>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>

int findmax (double num0,double num1,double num2) {
  if (num0>=num1){
    if (num0>=num2) {
      return 0;
    }
    else {
      return 2;
    }
  }
  else {
    if (num1>=num2) {
      return 1;
    }
    else {
      return 2;
    }
  }
}
int findmin (double num0,double num1,double num2) {
  if (num0<=num1){
    if (num0<=num2) {
      return 0;
    }
    else {
      return 2;
    }
  }
  else {
    if (num1<=num2) {
      return 1;
    }
    else {
      return 2;
    }
  }
}
double eta2theta(double eta) {
  return 2.0*TMath::ATan(TMath::Exp(-eta));
}

double P2E(double mass, double P) {
  return sqrt(mass*mass + P*P);
}

void Mjjana(const char *inputFile)
{
  double eta[3];
  int max_num;
  double position[2];
  
  double jet1mass;
  double jet1theta;
  double jet1phi;
  double jet1pt;
  
  double jet2mass;
  double jet2theta;
  double jet2phi;
  double jet2pt;
  
  double jet1PX;
  double jet1PY;
  double jet1PZ;
  double jet1E;
  double jet1P;
  
  double jet2PX;
  double jet2PY;
  double jet2PZ;
  double jet2E;
  double jet2P;
  
  double Mjj;
  double jjPX;
  double jjPY;
  double jjPZ;
  double jjE;
  
  double MET;
  
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  // Book histograms
  TH1 *histMjj = new TH1F("Mjj", "M_{jj}", 30, 0.0, 300.0);
  TH1 *histfwdJeteta = new TH1F("eta", "forward #eta", 30, -6.0, 3.0);
  TH1 *histmissingET = new TH1F("missingET", "missing E_{T}", 25, 0.0, 100.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    //cout << "Num:"<<entry<< endl;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    for (int i=0;i<branchJet->GetEntries();i++) {
      Jet *jet = (Jet*) branchJet->At(i);
      //printf("%f\t",jet->Eta);
    }
    //printf("\n");
    MissingET *missinget = (MissingET*) branchMET->At(0);
    histmissingET->Fill(missinget->MET);
    MET=missinget->MET;
    if (branchJet->GetEntries() == 3 && MET >= 20) {
      Jet *jet = (Jet*) branchJet->At(0);
      eta[0] = jet->Eta;
      jet = (Jet*) branchJet->At(1);
      eta[1] = jet->Eta;
      jet = (Jet*) branchJet->At(2);
      eta[2] = jet->Eta;
      max_num = findmin(eta[0],eta[1],eta[2]);
      if (max_num==0) {
        position[0]=1;
        position[1]=2;
      }
      else if (max_num==1) {
        position[0]=0;
        position[1]=2;
      }
      else {
        position[0]=0;
        position[1]=1;
      }
      
      histfwdJeteta->Fill(eta[max_num]);
      //printf("%f\t%f\t%f\t%d\n",eta[0],eta[1],eta[2],max_num);
      
      jet = (Jet*) branchJet->At(position[0]);
      jet1mass=jet->Mass;
      jet1theta=eta2theta(jet->Eta);
      jet1phi=jet->Phi;
      jet1pt=jet->PT;
      
      jet = (Jet*) branchJet->At(position[1]);
      jet2mass=jet->Mass;
      jet2theta=eta2theta(jet->Eta);
      jet2phi=jet->Phi;
      jet2pt=jet->PT;
      
      jet1PX=jet1pt*cos(jet1phi);
      jet1PY=jet1pt*sin(jet1phi);
      jet1PZ=jet1pt/tan(jet1theta);
      jet1P=jet1pt/sin(jet1theta);
      jet1E=P2E(jet1mass,jet1P);
      
      jet2PX=jet2pt*cos(jet2phi);
      jet2PY=jet2pt*sin(jet2phi);
      jet2PZ=jet2pt/tan(jet2theta);
      jet2P=jet2pt/sin(jet2theta);
      jet2E=P2E(jet2mass,jet2P);
      
      jjPX=jet1PX+jet2PX;
      jjPY=jet1PY+jet2PY;
      jjPZ=jet1PZ+jet2PZ;
      jjE=jet1E+jet1E;
      
      Mjj=sqrt(jjE*jjE-jjPX*jjPX-jjPY*jjPY-jjPZ*jjPZ);
      //printf("%f\n",Mjj);
    histMjj->Fill(Mjj);
    }
    
  }
  
  auto hist0=new TCanvas();
  hist0->SetLogz();
  hist0->Divide(3,1);

  // Show resulting histograms
  hist0->cd(1);
  histMjj->SetXTitle("M_{jj}/GeV");
  histMjj->SetYTitle("Events");
  histMjj->Draw();
  
  hist0->cd(2);
  histfwdJeteta->SetXTitle("#eta");
  histfwdJeteta->SetYTitle("Events");
  histfwdJeteta->Draw();
  
  hist0->cd(3);
  histmissingET->SetXTitle("E_T");
  histmissingET->SetYTitle("Events");
  histmissingET->Draw();
  //histMass->Draw();
  //TFile *sf=TFile::Open("cards/MuIC/hist-higgs-signal.root","recreate");
  //histMjj->Write();
}
