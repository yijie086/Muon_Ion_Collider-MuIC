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
#include "TCutG.h"
#include "TCanvas.h"
#include <numeric>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>
//#include "include/coordinateTools.h"

double eta2theta(double eta) {
	return 2.0*TMath::ATan(TMath::Exp(-eta));
}

void readout_K() {
	
	TFile * f = TFile::Open("MuIC111.root","read");
	TTree * t = (TTree*) f->Get("trackTree");
	
	std::vector< float > * do_eta=0;
	std::vector< float > * do_pt=0;
	std::vector< float > * do_Q2=0;
	std::vector< float > * do_x=0;
	
	ofstream phifile;
	phifile.open("KKfull.csv",ios::out|ios::trunc);
	t->SetBranchAddress("Kn_eta",&do_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_eta->size();j++) {
			phifile<<eta2theta(do_eta->at(j))<<"\t";
		}
	}
	t->SetBranchAddress("Kp_eta",&do_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_eta->size();j++) {
			phifile<<eta2theta(do_eta->at(j))<<"\t";
		}
	}
	t->SetBranchAddress("Kn_pt",&do_pt);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_pt->size();j++) {
			phifile<<do_pt->at(j)<<"\t";
		}
	}
	t->SetBranchAddress("Kp_pt",&do_pt);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_pt->size();j++) {
			phifile<<do_pt->at(j)<<"\t";
		}
	}
	t->SetBranchAddress("Kn_Q2",&do_Q2);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_Q2->size();j++) {
			phifile<<do_Q2->at(j)<<"\t";
		}
	}
	t->SetBranchAddress("Kp_Q2",&do_Q2);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_Q2->size();j++) {
			phifile<<do_Q2->at(j)<<"\t";
		}
	}
	t->SetBranchAddress("Kn_x",&do_x);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_x->size();j++) {
			phifile<<do_x->at(j)<<"\t";
		}
	}
	t->SetBranchAddress("Kp_x",&do_x);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_x->size();j++) {
			phifile<<do_x->at(j)<<"\t";
		}
	}
	phifile<<"\n";
	phifile.close();
	
}

