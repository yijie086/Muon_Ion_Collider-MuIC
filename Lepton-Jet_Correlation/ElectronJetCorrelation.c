#include <string>
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
#include "TLegend.h"
#include <numeric>
#include <iomanip>
#include <vector>
#include <iostream>
#include <fstream>

//Number of bins
const int num=50;

void Binnormal(TH1 *h1, double c) {
	// ***** This function return hist with the total=1 *****
	for (int i=1; i<=num; i++) {
		if (h1->GetBinContent(i)!=0) {
			h1->SetBinContent(i,h1->GetBinContent(i)/c);
		}
	}
	return;
}

void ElectronJetCorrelation() {
	// ***** This is the main function *****
	
	//Draw four pt level hist pt \in (9,11) (18,22) (45,55) (80,120)
	auto hist9t11=new TH1D("","",num,0.0,0.5);
	auto hist18t22=new TH1D("hist18t22","hist18t22",num,0.0,0.5);
	auto hist45t55=new TH1D("hist45t55","hist45t55",num,0.0,0.5);
	auto hist80t120=new TH1D("hist80t120","hist80t120",num,0.0,0.5);
	
	//Load the Pythia MC data file
	TFile * f = TFile::Open("EIC.root","read");
	TTree * t = (TTree*) f->Get("trackTree");
	
	//Define some varables
	std::vector< float > * do_pt=0;
	std::vector< float > * do_pt2=0;
	std::vector< float > do_pttemp;
	std::vector< float > do_pt2temp;
	
	//Load the varables
	t->SetBranchAddress("gendPhiej",&do_pt);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_pt->size();j++) {
			do_pttemp.push_back(do_pt->at(j));
		}
	}
	t->SetBranchAddress("genoneJetPt",&do_pt2);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_pt2->size();j++) {
			do_pt2temp.push_back(do_pt2->at(j));
		}
	}
	
	// Define the total num and record for it to normalize
	double N9t11=0;
	double N18t22=0;
	double N45t55=0;
	double N80t120=0;
	for(int i=0; i<do_pttemp.size(); i++) {
		if (do_pt2temp[i]>9&&do_pt2temp[i]<11) {
			hist9t11->Fill(do_pttemp[i]);
			N9t11++;
		}
		if (do_pt2temp[i]>18&&do_pt2temp[i]<22) {
			hist18t22->Fill(do_pttemp[i]);
			N18t22++;
		}
		if (do_pt2temp[i]>45&&do_pt2temp[i]<55) {
			hist45t55->Fill(do_pttemp[i]);
			N45t55++;
		}
		if (do_pt2temp[i]>80&&do_pt2temp[i]<120) {
			hist80t120->Fill(do_pttemp[i]);
			N80t120++;
		}
	}
	
	//Normalize the hist
	Binnormal(hist9t11,N9t11);
	Binnormal(hist18t22,N18t22);
	Binnormal(hist45t55,N45t55);
	Binnormal(hist80t120,N80t120);
	
	//Draw plots
	auto jhist=new TCanvas();
	jhist->SetLogy();
	jhist->Divide(1,1);
	jhist->cd(1);
	
	//Set the Axis
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	hist9t11->SetAxisRange(0.001, 1.0,"Y");
	hist9t11->GetYaxis()->SetLabelSize(0.04);
	hist9t11->GetXaxis()->SetLabelSize(0.04);
	hist9t11->GetYaxis()->SetTitleSize(0.04);
	hist9t11->GetXaxis()->SetTitleSize(0.04);
	
	//Set the Lines
	hist9t11->SetLineWidth(4);
	hist18t22->SetLineWidth(4);
	hist45t55->SetLineWidth(4);
	hist80t120->SetLineWidth(4);
	hist9t11->SetLineColor(kRed);
	hist18t22->SetLineColor(kBlue);
	hist45t55->SetLineColor(kGreen);
	hist80t120->SetLineColor(kBlack);
	
	//Set the title
	hist9t11->Draw("HIST");
	hist9t11->SetYTitle("Normalized counts");
	hist9t11->SetXTitle("|#phi^{jet}-#phi^{electron}-#pi|");
	gStyle->SetLabelSize(.09, "XY");
	hist9t11->SetTitleSize(0.06,"xyz");
	
	//Draw the remaining hist
	hist18t22->Draw("HIST SAME");
	hist45t55->Draw("HIST SAME");
	hist80t120->Draw("HIST SAME");
	
	//Draw the legend
	auto legend = new TLegend(0.47,0.6,0.88,0.88);
	legend->SetTextSize(0.038);
	legend->AddEntry(hist9t11,"9GeV/c<p_{T}^{electron}<11GeV/c","l");
	legend->AddEntry(hist18t22,"18GeV/c<p_{T}^{electron}<22GeV/c","l");
	legend->AddEntry(hist45t55,"45GeV/c<p_{T}^{electron}<55GeV/c","l");
	legend->AddEntry(hist80t120,"80GeV/c<p_{T}^{electron}<120GeV/c","l");
	legend->Draw();
}

