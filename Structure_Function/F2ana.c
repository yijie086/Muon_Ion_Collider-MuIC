#include <string>
//#include "include/Timer.h"
#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLatex.h"
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
#include "TAxis.h"
#include "TH1.h"
#include "TArrayD.h"

const int sz=20;

void Copy(double * E, const TH1 *h1) {
	for (int i=1; i<=sz; i++) {
		E[i]=h1->GetBinContent(i);
		//printf("%f\n",E[i]);
	}
	return;
}

void Clear(TH1 *h1) {
	int flag=0;
	for (int i=sz; h1->GetBinContent(i)==0||flag==0; i--) {
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
	for (int i=0; h1->GetBinContent(i)==0||flag==0; i++) {
		if (h1->GetBinContent(i)!=0) {
			flag=1;
		}
		h1->SetBinContent(i,0);
		h1->SetBinError(i,0);
		//printf("%d\n",i);
	}
	return;
}

void Set(double * E, TH1 *h1,double c1) {
	for (int i=1; i<=sz; i++) {
		if (h1->GetBinContent(i)!=0) {
			h1->SetBinContent(i,E[i]+c1);
		}
		//printf("%f\n",E[i]);
	}
	return;
}

void Print(TH1 *h1, double c) {
	for (int i=1; i<=sz; i++) {
		if (h1->GetBinContent(i)!=0) {
			printf("%f\t%f\t",TMath::Power(10,(6.0/sz)*(i-1)),TMath::Power(10,(6.0/sz)*i));
			printf("%f\t%f\n",h1->GetBinContent(i)-c,h1->GetBinError(i));
		}
		else {
			printf("%f\t%f\t",TMath::Power(10,(6.0/sz)*(i-1)),TMath::Power(10,(6.0/sz)*i));
			printf("%f\t%f\n",h1->GetBinContent(i),h1->GetBinError(i));
		}
	}
}

void Errorreset(TH1 *h1, double c) {
	for (int i=1; i<=sz; i++) {
		if (h1->GetBinContent(i)!=0) {
			h1->SetBinError(i,h1->GetBinError(i)/c);
		}
		//printf("%f\n",E[i]);
	}
	return;
}

void F2ana() {
	
	double xbins[sz+1];
	double error[sz+1];
	
	for (int i=0;i<=sz;i++) {
		xbins[i]=TMath::Power(10,i*(6.0/sz));
	}
	
	//TH2 *Q2tx=new TH2D("F2xtQ2","F2xtQ2",100,0,1,100,0,10000);
	TH1 *F2tx=new TH1D("","x=0.3",sz,xbins);
	TH1 *F2tx0=new TH1D("","x=0.1",sz,xbins);
	TH1 *F2tx05=new TH1D("","x=0.03",sz,xbins);
	TH1 *F2tx1=new TH1D("","x=0.01",sz,xbins);
	TH1 *F2tx15=new TH1D("","x=0.003",sz,xbins);
	TH1 *F2tx2=new TH1D("","x=0.001",sz,xbins);
	TH1 *F2tx25=new TH1D("","x=0.0003",sz,xbins);
	TH1 *F2tx3=new TH1D("","x=0.0001",sz,xbins);
	TH1 *F2tx4=new TH1D("","x=0.00003",sz,xbins);
	
	TH1 *mF2tx=new TH1D("","x=0.3",sz,xbins);
	TH1 *mF2tx0=new TH1D("","x=0.1",sz,xbins);
	TH1 *mF2tx05=new TH1D("","x=0.03",sz,xbins);
	TH1 *mF2tx1=new TH1D("","x=0.01",sz,xbins);
	TH1 *mF2tx15=new TH1D("","x=0.003",sz,xbins);
	TH1 *mF2tx2=new TH1D("","x=0.001",sz,xbins);
	TH1 *mF2tx25=new TH1D("","x=0.0003",sz,xbins);
	TH1 *mF2tx3=new TH1D("","x=0.0001",sz,xbins);
	TH1 *mF2tx4=new TH1D("","x=0.00003",sz,xbins);
	
	TFile * f = TFile::Open("EIC.root","read");
	TTree * t = (TTree*) f->Get("trackTree");
	
	TFile * f1 = TFile::Open("MuIC01.root","read");
	TTree * t1 = (TTree*) f1->Get("trackTree");
	
	float weight;
	
	std::vector< float > xtemp;
	std::vector< float > ytemp;
	std::vector< float > Q2temp;
	
	std::vector< float > * F2x=0;
	std::vector< float > * F2y=0;
	std::vector< float > * F2Q2=0;
	std::vector< float > * genJetNum=0;
	
	std::vector< float > mxtemp;
	std::vector< float > mytemp;
	std::vector< float > mQ2temp;
	
	std::vector< float > * mF2x=0;
	std::vector< float > * mF2y=0;
	std::vector< float > * mF2Q2=0;
	std::vector< float > * mgenJetNum=0;
	
	double Q2binlen;
	//double xbinset=0.01;
	double xbe;
	double xen;
	double xbinlen;
	double simgatotal=0.0131222744016*2/5.109;
	double msimgatotal=0.0131222744016;
	
	t->SetBranchAddress("F2x",&F2x);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<F2x->size();j++) {
			xtemp.push_back(F2x->at(j));
		}
	}
	
	t->SetBranchAddress("F2y",&F2y);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<F2y->size();j++) {
			ytemp.push_back(F2y->at(j));
		}
	}
	
	t->SetBranchAddress("F2Q2",&F2Q2);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<F2Q2->size();j++) {
			Q2temp.push_back(F2Q2->at(j));
		}
	}
	
	t1->SetBranchAddress("F2x",&mF2x);
	for(int i=0;i<t1->GetEntries();i++){
		t1->GetEntry(i);
		for (int j=0;j<mF2x->size();j++) {
			mxtemp.push_back(mF2x->at(j));
		}
	}
	
	t1->SetBranchAddress("F2y",&mF2y);
	for(int i=0;i<t1->GetEntries();i++){
		t1->GetEntry(i);
		for (int j=0;j<mF2y->size();j++) {
			mytemp.push_back(mF2y->at(j));
		}
	}
	
	t1->SetBranchAddress("F2Q2",&mF2Q2);
	for(int i=0;i<t1->GetEntries();i++){
		t1->GetEntry(i);
		for (int j=0;j<mF2Q2->size();j++) {
			mQ2temp.push_back(mF2Q2->at(j));
		}
	}
	
	double n;
	double mn;
	
	n=100000000;
	mn=50000000;
	//h2dpionp->Fill(pionp_etatemp[i],pionp_pttemp[i]);
	
	double scale=1.0/3.0;
	
	for(int i=0; i<xtemp.size(); i++) {
		
		//weight=((xtemp[i]*Q2temp[i]*Q2temp[i])/(4.0*3.14159265/(137.0*137.0)))/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
		
		//weight=1.0;
		
		//Q2tx->Fill(xtemp[i],Q2temp[i],weight);
		xbe=0.08;
		xen=0.12;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx0->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.028;
		xen=0.032;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx05->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.008;
		xen=0.012;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx1->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.0028;
		xen=0.0032;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx15->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.0008;
		xen=0.0012;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx2->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.00028;
		xen=0.00032;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx25->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.00008;
		xen=0.00012;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx3->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.000028;
		xen=0.000032;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx4->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
		xbe=0.28;
		xen=0.32;
		xbinlen=xen-xbe;
		if (xbe<xtemp[i]&&xtemp[i]<xen) {
			weight=xtemp[i]*Q2temp[i]*Q2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(Q2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			F2tx->Fill(Q2temp[i],scale*weight*simgatotal/(n*Q2binlen*xbinlen));
		}
	}
	for(int i=0; i<mxtemp.size(); i++) {
		
		//weight=((xtemp[i]*Q2temp[i]*Q2temp[i])/(4.0*3.14159265/(137.0*137.0)))/(1-ytemp[i]+ytemp[i]*ytemp[i]/2.0-(ytemp[i]*ytemp[i]*0.938*0.938)/(Q2temp[i]));
		
		//weight=1.0;
		
		//Q2tx->Fill(xtemp[i],Q2temp[i],weight);
		xbe=0.08;
		xen=0.12;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx0->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.028;
		xen=0.032;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx05->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.008;
		xen=0.012;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx1->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.0028;
		xen=0.0032;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx15->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.0008;
		xen=0.0012;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx2->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.00028;
		xen=0.00032;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx25->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.00008;
		xen=0.00012;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx3->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.000028;
		xen=0.000032;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx4->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
		xbe=0.28;
		xen=0.32;
		xbinlen=xen-xbe;
		if (xbe<mxtemp[i]&&mxtemp[i]<xen) {
			weight=mxtemp[i]*mQ2temp[i]*mQ2temp[i]/(4.0*3.14159265/(137.0*137.0));
			weight=weight/(1-mytemp[i]+mytemp[i]*mytemp[i]/2.0-(mytemp[i]*mytemp[i]*0.938*0.938)/(mQ2temp[i]));
			Q2binlen=(TMath::Power(10,(6.0/sz))-1)*TMath::Power(10,(6.0/sz)*floor(TMath::Log10(mQ2temp[i])/(6.0/sz)));
			//printf("%f\n",weight*simgatotal/(n*Q2binlen*xbinlen));
			mF2tx->Fill(mQ2temp[i],scale*weight*msimgatotal/(mn*Q2binlen*xbinlen));
		}
	}
	//gStyle->SetOptStat(0);
	
	//auto hist2d0=new TCanvas();
	//hist2d0->SetLogz();
	//hist2d0->Divide(1,1);
	//hist2d0->cd(1);
	//gPad->SetLogz();
	//gStyle->SetPalette(kBird);
	//Q2tx->SetMaximum(100);
	//Q2tx->SetMinimum(1000);
	//Q2tx->Draw("COLZ");
	//Q2tx->SetTitle("");
	//hist2d0->Draw();
	
	//hist2d0->SaveAs("F2.pdf");
	
	
	double resetnum=20;
	
	auto hist2d1=new TCanvas();
	hist2d1->Divide(1,1);
	hist2d1->cd(1);
	gPad->SetLogx();
	gStyle->SetErrorX(0.);
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	Copy(error,F2tx0);
	F2tx0->SetAxisRange(0.5, 7.,"Y");
	F2tx0->SetXTitle("Q^{2} [GeV^{2}]");
	F2tx0->SetYTitle("F_{2}(x,Q^{2})-log_{10}(x)");
	//F2tx->Scale(1,"");
	Set(error,F2tx0,1);
	Errorreset(F2tx0,3.16);
	
	Copy(error,F2tx05);
	Set(error,F2tx05,1.5);
	Errorreset(F2tx05,3.16);
	
	
	Copy(error,F2tx1);
	Set(error,F2tx1,2);
	Errorreset(F2tx1,3.16);
	
	
	Copy(error,F2tx15);
	Set(error,F2tx15,2.5);
	Errorreset(F2tx15,3.16);
	
	
	Copy(error,F2tx2);
	Set(error,F2tx2,3);
	Errorreset(F2tx2,3.16);
	
	
	Copy(error,F2tx25);
	Set(error,F2tx25,3.5);
	Errorreset(F2tx25,3.16);
	
	
	Copy(error,F2tx3);
	Set(error,F2tx3,4);
	Errorreset(F2tx3,3.16);
	
	
	Copy(error,F2tx4);
	Set(error,F2tx4,4.5);
	Errorreset(F2tx4,3.16);
	
	
	Copy(error,F2tx);
	Set(error,F2tx,0.3);
	Errorreset(F2tx,3.16);
	
	F2tx->SetBinError(2,0.1);
	
	Copy(error,mF2tx0);
	Set(error,mF2tx0,1);
	Errorreset(F2tx0,14);
	
	Copy(error,mF2tx05);
	Set(error,mF2tx05,1.5);
	Errorreset(F2tx05,14);
	
	
	Copy(error,mF2tx1);
	Set(error,mF2tx1,2);
	Errorreset(F2tx1,14);
	
	
	Copy(error,mF2tx15);
	Set(error,mF2tx15,2.5);
	Errorreset(F2tx15,14);
	
	
	Copy(error,mF2tx2);
	Set(error,mF2tx2,3);
	Errorreset(F2tx2,14);
	
	
	Copy(error,mF2tx25);
	Set(error,mF2tx25,3.5);
	Errorreset(F2tx25,14);
	
	
	Copy(error,mF2tx3);
	Set(error,mF2tx3,4);
	Errorreset(F2tx3,14);
	
	
	Copy(error,mF2tx4);
	Set(error,mF2tx4,4.5);
	Errorreset(F2tx4,14);
	
	
	Copy(error,mF2tx);
	Set(error,mF2tx,0.3);
	Errorreset(F2tx,14);
	
	
	
	Clear0(F2tx0);
	Clear(F2tx0);
	Clear0(mF2tx0);
	Clear0(mF2tx0);
	
	Clear(F2tx05);
	
	Clear(F2tx1);
	Clear(mF2tx1);
	//Clear0(mF2tx1);
	
	Clear(F2tx15);
	Clear0(mF2tx15);
	Clear(mF2tx15);
	
	Clear(F2tx2);
	Clear(mF2tx2);
	Clear(mF2tx2);
	
	Clear(F2tx25);
	Clear(mF2tx25);
	Clear(mF2tx25);
	
	Clear(mF2tx3);
	
	
	//F2tx->SetMarkerStyle(kFullCircle);
	F2tx0->SetMarkerStyle(kFullCircle);
	F2tx05->SetMarkerStyle(kFullCircle);
	F2tx1->SetMarkerStyle(kFullCircle);
	F2tx15->SetMarkerStyle(kFullCircle);
	F2tx2->SetMarkerStyle(kFullCircle);
	F2tx25->SetMarkerStyle(kFullCircle);
	F2tx3->SetMarkerStyle(kFullCircle);
	F2tx4->SetMarkerStyle(kFullCircle);
	
	F2tx0->SetMarkerColor(kRed);
	F2tx05->SetMarkerColor(kRed);
	F2tx1->SetMarkerColor(kRed);
	F2tx15->SetMarkerColor(kRed);
	F2tx2->SetMarkerColor(kRed);
	F2tx25->SetMarkerColor(kRed);
	F2tx3->SetMarkerColor(kRed);
	F2tx4->SetMarkerColor(kRed);
	
	mF2tx0->SetMarkerStyle(kFullCircle);
	mF2tx05->SetMarkerStyle(kFullCircle);
	mF2tx1->SetMarkerStyle(kFullCircle);
	mF2tx15->SetMarkerStyle(kFullCircle);
	mF2tx2->SetMarkerStyle(kFullCircle);
	mF2tx25->SetMarkerStyle(kFullCircle);
	mF2tx3->SetMarkerStyle(kFullCircle);
	mF2tx4->SetMarkerStyle(kFullCircle);
	
	mF2tx0->SetMarkerColor(kBlue);
	mF2tx05->SetMarkerColor(kBlue);
	mF2tx1->SetMarkerColor(kBlue);
	mF2tx15->SetMarkerColor(kBlue);
	mF2tx2->SetMarkerColor(kBlue);
	mF2tx25->SetMarkerColor(kBlue);
	mF2tx3->SetMarkerColor(kBlue);
	mF2tx4->SetMarkerColor(kBlue);
	
	F2tx0->Draw();
	//F2tx->Draw("SAME PLC PMC");
	F2tx05->Draw("SAME");
	F2tx1->Draw("SAME");
	F2tx15->Draw("SAME");
	F2tx2->Draw("SAME");
	F2tx25->Draw("SAME");
	F2tx3->Draw("SAME");
	F2tx4->Draw("SAME");
	
	mF2tx0->Draw("SAME");
	mF2tx05->Draw("SAME");
	mF2tx1->Draw("SAME");
	mF2tx15->Draw("SAME");
	mF2tx2->Draw("SAME");
	mF2tx25->Draw("SAME");
	mF2tx3->Draw("SAME");
	mF2tx4->Draw("SAME");
	
	printf("F2tx:x=0.3\n");
	Print(F2tx,0.5);
	
	printf("F2tx0:x=0.1\n");
	Print(F2tx0,1);
	
	printf("F2tx05:x=0.03\n");
	Print(F2tx05,1.5);
	
	printf("F2tx1:x=0.01\n");
	Print(F2tx1,2);
	
	printf("F2tx15:x=0.003\n");
	Print(F2tx15,2.5);
	
	printf("F2tx2:x=0.001\n");
	Print(F2tx2,3);
	
	printf("F2tx25:x=0.0003\n");
	Print(F2tx25,3.5);
	
	printf("F2tx3:x=0.0001\n");
	Print(F2tx3,4);
	
	printf("F2tx4:x=0.00003\n");
	Print(F2tx4,4.5);
	
	printf("mF2tx:x=0.3\n");
	Print(mF2tx,0.5);
	
	printf("mF2tx0:x=0.1\n");
	Print(mF2tx0,1);
	
	printf("mF2tx05:x=0.03\n");
	Print(mF2tx05,1.5);
	
	printf("mF2tx1:x=0.01\n");
	Print(mF2tx1,2);
	
	printf("mF2tx15:x=0.003\n");
	Print(mF2tx15,2.5);
	
	printf("mF2tx2:x=0.001\n");
	Print(mF2tx2,3);
	
	printf("mF2tx25:x=0.0003\n");
	Print(mF2tx25,3.5);
	
	printf("mF2tx3:x=0.0001\n");
	Print(mF2tx3,4);
	
	printf("mF2tx4:x=0.00003\n");
	Print(mF2tx4,4.5);
	
	TLatex tt; tt.SetNDC();
	tt.SetTextSize(0.035);
	tt.DrawLatex(0.28,  0.86, "x=0.00003");
	tt.DrawLatex(0.34,  0.82, "x=0.0001");
	tt.DrawLatex(0.42,  0.78, "x=0.0003");
	tt.DrawLatex(0.5,  0.70, "x=0.001");
	tt.DrawLatex(0.56,  0.56, "x=0.003");
	tt.DrawLatex(0.64,  0.47, "x=0.01");
	tt.DrawLatex(0.70,  0.38, "x=0.03");
	tt.DrawLatex(0.79,  0.3, "x=0.1");
	//tt.DrawLatex(0.82,  0.17, "x=0.3");
}

