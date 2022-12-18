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
#include "TLegend.h"
//#include "include/coordinateTools.h"

void Errorreset(TH1 *h1, double c) {
	for (int i=1; i<=25; i++) {
		if (h1->GetBinContent(i)!=0) {
			h1->SetBinError(i,h1->GetBinError(i)/c);
		}
		//printf("%f\n",E[i]);
	}
	return;
}

void Clear(TH1 *h1) {
	int flag=0;
	for (int i=25; h1->GetBinContent(i)==0||flag==0; i--) {
		if (h1->GetBinContent(i)!=0) {
			flag=1;
		}
		h1->SetBinContent(i,0);
		h1->SetBinError(i,0);
		//printf("%d\n",i);
	}
	return;
}

void Clear0(TH1 *h1) {
	int flag=0;
	for (int i=1; h1->GetBinContent(i)==0||flag==0; i++) {
		if (h1->GetBinContent(i)!=0) {
			flag=1;
		}
		h1->SetBinContent(i,0);
		h1->SetBinError(i,0);
		//printf("%d\n",i);
	}
	return;
}

void Jetcrosssectionanaunited() {
	
	double xbins[26];
	
	for (int i=0;i<=25;i++) {
		xbins[i]=TMath::Power(10,i*0.12);
	}
	
	double Nx=3371892;
	double Nw=94;
	double Px=Nx/(Nx+Nw);
	double Pw=Nw/(Nx+Nw);
	
	TH1 *Jetpthist=new TH1D("Jetpt","",25,xbins);
	TH1 *EJetpthist=new TH1D("Jetpt","",25,xbins);
	
	TFile * f = TFile::Open("MuICx4.root","read");
	TTree * t = (TTree*) f->Get("trackTree");
	
	std::vector< float > * genJetPt=0;
	std::vector< float > * genJetNum=0;
	
	std::vector< float > jetpttemp;
	
	double N;
	
	double binlen;
	
	t->SetBranchAddress("genJetPt",&genJetPt);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<genJetPt->size();j++) {
			jetpttemp.push_back(genJetPt->at(j));
		}
	}
	
	t->SetBranchAddress("genJetNum",&genJetNum);
	
	N=1.48*t->GetEntries();
	
	//h2dpionp->Fill(pionp_etatemp[i],pionp_pttemp[i]);
	
	for(int i=0; i<jetpttemp.size(); i++) {
		binlen=(TMath::Power(10,0.12)-1)*TMath::Power(10,0.12*ceil(TMath::Log10(jetpttemp[i])/0.12));
		Jetpthist->Fill(jetpttemp[i],Px*5.109*TMath::Power(10,6)/(N*binlen));
	}
	Clear(Jetpthist);
	Clear(Jetpthist);
	Errorreset(Jetpthist, 100);
	
	TFile * f1 = TFile::Open("MuICw4.root","read");
	TTree * t1 = (TTree*) f1->Get("trackTree");
	jetpttemp.clear();
	t1->SetBranchAddress("genJetPt",&genJetPt);
	for(int i=0;i<t1->GetEntries();i++){
		t1->GetEntry(i);
		for (int j=0;j<genJetPt->size();j++) {
			jetpttemp.push_back(genJetPt->at(j));
		}
	}
	t1->SetBranchAddress("genJetNum",&genJetNum);
	
	N=1.48*t1->GetEntries();
	
	//h2dpionp->Fill(pionp_etatemp[i],pionp_pttemp[i]);
	
	for(int i=0; i<jetpttemp.size(); i++) {
		binlen=(TMath::Power(10,0.12)-1)*TMath::Power(10,0.12*ceil(TMath::Log10(jetpttemp[i])/0.12));
		Jetpthist->Fill(jetpttemp[i],Pw*5.109*TMath::Power(10,6)/(N*binlen));
	}
	
	
	TFile * f2 = TFile::Open("EICx4.root","read");
	TTree * t2 = (TTree*) f2->Get("trackTree");
	jetpttemp.clear();
	t2->SetBranchAddress("genJetPt",&genJetPt);
	for(int i=0;i<t2->GetEntries();i++){
		t2->GetEntry(i);
		for (int j=0;j<genJetPt->size();j++) {
			jetpttemp.push_back(genJetPt->at(j));
		}
	}
	t2->SetBranchAddress("genJetNum",&genJetNum);
	N=1.23*t2->GetEntries();
	
	for(int i=0; i<jetpttemp.size(); i++) {
		binlen=(TMath::Power(10,0.12)-1)*TMath::Power(10,0.12*ceil(TMath::Log10(jetpttemp[i])/0.12));
		EJetpthist->Fill(jetpttemp[i],1.628*TMath::Power(10,6)/(N*binlen));
	}
	Errorreset(EJetpthist, 31.62);
	
	Clear0(Jetpthist);
	Clear0(EJetpthist);
	
	gStyle->SetOptStat(0);
	
	auto hist2d0=new TCanvas();
	hist2d0->SetLogz();
	hist2d0->Divide(1,1);
	hist2d0->cd(1);
	gPad->SetLogy();
	gPad->SetLogx();
	Jetpthist->Draw();
	EJetpthist->Draw("SAME");
	Jetpthist->SetAxisRange(5,800.,"X");
	Jetpthist->SetXTitle("p_{T} [GeV]");
	Jetpthist->SetYTitle("d #sigma/d p_{T} [pb/GeV]");
	Jetpthist->GetYaxis()->SetLabelSize(0.04);
	Jetpthist->GetXaxis()->SetLabelSize(0.04);
	gStyle->SetLabelSize(.08, "XY");
	Jetpthist->GetYaxis()->SetLabelSize(0.04);
	Jetpthist->GetXaxis()->SetLabelSize(0.04);
	Jetpthist->SetLineWidth(3);
	EJetpthist->SetLineWidth(3);
	EJetpthist->SetLineColor(kRed);
	Jetpthist->SetLineColor(kBlue);
	
	auto legend = new TLegend(0.6,0.6,0.88,0.88);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->SetTextSize(0.04);
	legend->AddEntry(Jetpthist,"MuIC","l");
	legend->AddEntry(EJetpthist,"EIC","l");
	legend->Draw();
	
	hist2d0->SaveAs("ptofjet.pdf");
}

