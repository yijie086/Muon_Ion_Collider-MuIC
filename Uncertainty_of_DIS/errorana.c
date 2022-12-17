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
//#include "include/coordinateTools.h"
int readindec(std::vector< float > reading, int postion) {
	double flag=0;
	int i;
	for (i=postion; i<reading.size()&&flag==0; i++) {
		printf("%d\t %lf\n",i,reading[i]);
		if (reading[i+1]!=reading[postion]) {
			flag=1;
		}
	}
	printf("!!%d!!\n",i);
	return(i-1);
}

double eta2theta(double eta) {
	return 2.0*TMath::ATan(TMath::Exp(-eta));
}

double E2p(double mass, double energy) {
	return TMath::Power(energy*energy-mass*mass,0.5);
}

double error_combine(double error_1, double error_2) {
	return TMath::Power(error_1*error_1+error_2*error_2,0.5);
}
void MEAN(TH2 *num, TH2 *Data,double binx,double biny) {
	for (int i=1; i<=binx; i++) {
		for (int j=1; j<=biny; j++) {
			if (Data->GetBinContent(i,j)!=0&&num->GetBinContent(i,j)!=0) {
				printf("(%d,%d),%f,%f,%f\t\t!!!",i,j,Data->GetBinContent(i,j),num->GetBinContent(i,j),Data->GetBinContent(i,j)/num->GetBinContent(i,j));
				Data->SetBinContent(i,j,Data->GetBinContent(i,j)/num->GetBinContent(i,j));
				printf("%f\n",Data->GetBinContent(i,j));
			}
		}
	}
}
void Errorreset(TH1 *h1, double c, int num) {
	for (int i=1; i<=num; i++) {
		if (h1->GetBinContent(i)!=0) {
			h1->SetBinError(i,h1->GetBinError(i)*c);
		}
		//printf("%f\n",E[i]);
	}
	return;
}
void Binnormal(TH1 *h1, double c, int num) {
	for (int i=1; i<=num; i++) {
		if (h1->GetBinContent(i)!=0) {
			h1->SetBinContent(i,h1->GetBinContent(i)/c);
		}
		//printf("%f\n",E[i]);
	}
	return;
}


void errorana() {

	double tempi;
	
	int Q2xxnum=50;
	double xbe=-5;
	double ybe=0;
	double xen=0;
	double yen=6;
	int Q2xynum=50;
	double xbin[Q2xxnum+1];
	double ybin[Q2xynum+1];
	int D1num=50;
	
	for (int i=0; i<=Q2xxnum; i++) {
		xbin[i]=TMath::Power(10.0,xbe+i*((xen-xbe)/Q2xxnum));
	}
	for (int i=0; i<=Q2xynum; i++) {
		ybin[i]=TMath::Power(10.0,ybe+i*((yen-ybe)/Q2xynum));
	}
	
	TH2 *Q2txQ2j=new TH2D("Q2txQ2j","Q2 JB",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txyj=new TH2D("Q2txyj","y JB",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxj=new TH2D("Q2txxj","x JB",Q2xxnum,xbin,Q2xynum,ybin);
	
	TH2 *Q2txQ2l=new TH2D("Q2txQ2l","Q2 LO",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txyl=new TH2D("Q2txyl","y LO",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxl=new TH2D("Q2txxl","x LO",Q2xxnum,xbin,Q2xynum,ybin);
	
	TH2 *Q2txQ2d=new TH2D("Q2txQ2d","Q2 DA",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txyd=new TH2D("Q2txyd","y DA",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxd=new TH2D("Q2txxd","x DA",Q2xxnum,xbin,Q2xynum,ybin);
	
	TH2 *Q2txQ2s=new TH2D("Q2txQ2s","Q2 Isigma",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txys=new TH2D("Q2txys","y Isigma",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxs=new TH2D("Q2txxs","x Isigma",Q2xxnum,xbin,Q2xynum,ybin);
	
	TH2 *Q2txnum=new TH2D("Q2tx","Q2tx",Q2xxnum,xbin,Q2xynum,ybin);
	
	TH2 *Q2txxsj=new TH2D("Q2txxsj","x=Q2_Isigma/(s*y_JB)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxsl=new TH2D("Q2txxsl","x=Q2_Isigma/(s*y_LO)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxsd=new TH2D("Q2txxsd","x=Q2_Isigma/(s*y_DA)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxjs=new TH2D("Q2txxjs","x=Q2_JB/(s*y_Isigma)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxjl=new TH2D("Q2txxjl","x=Q2_JB/(s*y_LO)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxjd=new TH2D("Q2txxjd","x=Q2_JB/(s*y_DA)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxls=new TH2D("Q2txxls","x=Q2_LO/(s*y_Isigma)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxlj=new TH2D("Q2txxlj","x=Q2_LO/(s*y_JB)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxld=new TH2D("Q2txxld","x=Q2_LO/(s*y_DA)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxds=new TH2D("Q2txxds","x=Q2_DA/(s*y_Isigma)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxdj=new TH2D("Q2txxdj","x=Q2_DA/(s*y_JB)",Q2xxnum,xbin,Q2xynum,ybin);
	TH2 *Q2txxdl=new TH2D("Q2txxdl","x=Q2_DA/(s*y_LO)",Q2xxnum,xbin,Q2xynum,ybin);
	
	TFile * f = TFile::Open("MUIC.root","read");
	TTree * t = (TTree*) f->Get("trackTree");
	
	std::vector< float > true_eta;
	
	std::vector< float > * do_eta=0;
	std::vector< float > * do_phi=0;
	std::vector< float > * do_e=0;
	std::vector< float > * do_Q2=0;
	std::vector< float > * do_x=0;
	std::vector< float > * do_y=0;
	
	std::vector< float > do_etatemp;
	std::vector< float > do_phitemp;
	std::vector< float > do_etemp;
	std::vector< float > do_Q2temp;
	std::vector< float > do_xtemp;
	std::vector< float > do_ytemp;
	
	
	std::vector< float > * pp_eta=0;
	std::vector< float > * pp_phi=0;
	std::vector< float > * pp_e=0;
	std::vector< float > pp_etatemp;
	std::vector< float > pp_phitemp;
	std::vector< float > pp_etemp;
	std::vector< float > * pp_num=0;
	std::vector< int > pp_numtemp;
	
	std::vector< float > * pn_eta=0;
	std::vector< float > * pn_phi=0;
	std::vector< float > * pn_e=0;
	std::vector< float > pn_etatemp;
	std::vector< float > pn_phitemp;
	std::vector< float > pn_etemp;
	std::vector< float > * pn_num=0;
	std::vector< int > pn_numtemp;
	
	std::vector< float > * np_eta=0;
	std::vector< float > * np_phi=0;
	std::vector< float > * np_e=0;
	std::vector< float > np_etatemp;
	std::vector< float > np_phitemp;
	std::vector< float > np_etemp;
	std::vector< float > * np_num=0;
	std::vector< int > np_numtemp;
	
	std::vector< float > * nn_eta=0;
	std::vector< float > * nn_phi=0;
	std::vector< float > * nn_e=0;
	std::vector< float > nn_etatemp;
	std::vector< float > nn_phitemp;
	std::vector< float > nn_etemp;
	std::vector< float > * nn_num=0;
	std::vector< int > nn_numtemp;
	
	std::vector< float > * pip_eta=0;
	std::vector< float > * pip_phi=0;
	std::vector< float > * pip_e=0;
	std::vector< float > pip_etatemp;
	std::vector< float > pip_phitemp;
	std::vector< float > pip_etemp;
	std::vector< float > * pip_num=0;
	std::vector< int > pip_numtemp;
	
	std::vector< float > * pin_eta=0;
	std::vector< float > * pin_phi=0;
	std::vector< float > * pin_e=0;
	std::vector< float > pin_etatemp;
	std::vector< float > pin_phitemp;
	std::vector< float > pin_etemp;
	std::vector< float > * pin_num=0;
	std::vector< int > pin_numtemp;
	
	std::vector< float > * Kp_eta=0;
	std::vector< float > * Kp_phi=0;
	std::vector< float > * Kp_e=0;
	std::vector< float > Kp_etatemp;
	std::vector< float > Kp_phitemp;
	std::vector< float > Kp_etemp;
	std::vector< float > * Kp_num=0;
	std::vector< int > Kp_numtemp;
	
	std::vector< float > * Kn_eta=0;
	std::vector< float > * Kn_phi=0;
	std::vector< float > * Kn_e=0;
	std::vector< float > Kn_etatemp;
	std::vector< float > Kn_phitemp;
	std::vector< float > Kn_etemp;
	std::vector< float > * Kn_num=0;
	std::vector< int > Kn_numtemp;
	
	std::vector< float > * KL_eta=0;
	std::vector< float > * KL_phi=0;
	std::vector< float > * KL_e=0;
	std::vector< float > KL_etatemp;
	std::vector< float > KL_phitemp;
	std::vector< float > KL_etemp;
	std::vector< float > * KL_num=0;
	std::vector< int > KL_numtemp;
	
	std::vector< float > * g_eta=0;
	std::vector< float > * g_phi=0;
	std::vector< float > * g_e=0;
	std::vector< float > g_etatemp;
	std::vector< float > g_phitemp;
	std::vector< float > g_etemp;
	std::vector< float > * g_num=0;
	std::vector< int > g_numtemp;
	
	std::vector< float > * ep_eta=0;
	std::vector< float > * ep_phi=0;
	std::vector< float > * ep_e=0;
	std::vector< float > ep_etatemp;
	std::vector< float > ep_phitemp;
	std::vector< float > ep_etemp;
	std::vector< float > * ep_num=0;
	std::vector< int > ep_numtemp;
	
	std::vector< float > * en_eta=0;
	std::vector< float > * en_phi=0;
	std::vector< float > * en_e=0;
	std::vector< float > en_etatemp;
	std::vector< float > en_phitemp;
	std::vector< float > en_etemp;
	std::vector< float > * en_num=0;
	std::vector< int > en_numtemp;
	
	double Q2;
	double y;
	double x;
	double theta;
	double yevalu;
	double yevalu2;
	double ytemp2;
	double wy;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	t->SetBranchAddress("muon_Q2",&do_Q2);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_Q2->size();j++) {
			do_Q2temp.push_back(do_Q2->at(j));
		}
	}
	t->SetBranchAddress("muon_x",&do_x);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_x->size();j++) {
			do_xtemp.push_back(do_x->at(j));
		}
	}
	t->SetBranchAddress("F2y",&do_y);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_x->size();j++) {
			do_ytemp.push_back(do_y->at(j));
		}
	}
	t->SetBranchAddress("muon_eta",&do_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_eta->size();j++) {
			true_eta.push_back(do_eta->at(j));
		}
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	double mrandomsigma_E=0.01;
	double mrandomsigma2_E=0.0001;
	double mrandomsigma_angle=0.0002;
	double crandomsigma_E=0.01;
	double crandomsigma2_E=0.001;
	double crandomsigma_angle=0.0002;
	double crandomsigma2_angle=0.002;
	double nrandomsigma_E=0.1;
	double nrandomsigma2_E=0.5;
	double nrandomsigma_angle=0.0251;
	double grandomsigma_E=0.02;
	double grandomsigma2_E=0.1;
	double grandomsigma_angle=0.0251;
	
	t->SetBranchAddress("muon_eta",&do_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_eta->size();j++) {
			do_etatemp.push_back(eta2theta(do_eta->at(j))+gRandom->Gaus(0,mrandomsigma_angle));
		}
	}
	t->SetBranchAddress("muon_phi",&do_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_phi->size();j++) {
			do_phitemp.push_back(do_phi->at(j)+gRandom->Gaus(0,mrandomsigma_angle));
		}
	}
	t->SetBranchAddress("muon_e",&do_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<do_e->size();j++) {
			do_etemp.push_back(do_e->at(j)*(1+gRandom->Gaus(0,error_combine(mrandomsigma_E, mrandomsigma2_E*E2p(0.10566,do_e->at(j))))));
		}
	}
	
	t->SetBranchAddress("proton_e",&pp_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pp_e->size();j++) {
			pp_etemp.push_back(pp_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.93827,pp_e->at(j))))));
		}
	}
	t->SetBranchAddress("proton_phi",&pp_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pp_phi->size();j++) {
			pp_phitemp.push_back(pp_phi->at(j));
		}
	}
	t->SetBranchAddress("proton_eta",&pp_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pp_eta->size();j++) {
			pp_etatemp.push_back(pp_eta->at(j));
		}
	}
	t->SetBranchAddress("proton_num",&pp_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pp_num->size();j++) {
			pp_numtemp.push_back(pp_num->at(j));
		}
	}
	for (int i=0; i<pp_etemp.size(); i++) {
		pp_etatemp[i]=pp_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.93827,pp_etemp[i])));
		pp_phitemp[i]=pp_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.93827,pp_etemp[i])));
	}
	
	t->SetBranchAddress("antiproton_e",&pn_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pn_e->size();j++) {
			pn_etemp.push_back(pn_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.93827,pn_e->at(j))))));
		}
	}
	t->SetBranchAddress("antiproton_phi",&pn_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pn_phi->size();j++) {
			pn_phitemp.push_back(pn_phi->at(j));
		}
	}
	t->SetBranchAddress("antiproton_eta",&pn_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pn_eta->size();j++) {
			pn_etatemp.push_back(pn_eta->at(j));
		}
	}
	t->SetBranchAddress("antiproton_num",&pn_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pn_num->size();j++) {
			pn_numtemp.push_back(pn_num->at(j));
		}
	}
	for (int i=0; i<pn_etemp.size(); i++) {
		pn_etatemp[i]=pn_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.93827,pn_etemp[i])));
		pn_phitemp[i]=pn_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.93827,pn_etemp[i])));
	}
	
	t->SetBranchAddress("n0_e",&np_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<np_e->size();j++) {
			np_etemp.push_back(np_e->at(j)*(1+gRandom->Gaus(0,error_combine(nrandomsigma_E, nrandomsigma2_E/sqrt(np_e->at(j))))));
		}
	}
	t->SetBranchAddress("n0_phi",&np_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<np_phi->size();j++) {
			np_phitemp.push_back(np_phi->at(j)+gRandom->Gaus(0,nrandomsigma_angle));
		}
	}
	t->SetBranchAddress("n0_eta",&np_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<np_eta->size();j++) {
			np_etatemp.push_back(np_eta->at(j)+gRandom->Gaus(0,nrandomsigma_angle));
		}
	}
	t->SetBranchAddress("n0_num",&np_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<np_num->size();j++) {
			np_numtemp.push_back(np_num->at(j));
		}
	}
	
	t->SetBranchAddress("antin0_e",&nn_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<nn_e->size();j++) {
			nn_etemp.push_back(nn_e->at(j)*(1+gRandom->Gaus(0,error_combine(nrandomsigma_E, nrandomsigma2_E/sqrt(nn_e->at(j))))));
		}
	}
	t->SetBranchAddress("antin0_phi",&nn_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<nn_phi->size();j++) {
			nn_phitemp.push_back(nn_phi->at(j)+gRandom->Gaus(0,nrandomsigma_angle));
		}
	}
	t->SetBranchAddress("antin0_eta",&nn_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<nn_eta->size();j++) {
			nn_etatemp.push_back(nn_eta->at(j)+gRandom->Gaus(0,nrandomsigma_angle));
		}
	}
	t->SetBranchAddress("antin0_num",&nn_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<nn_num->size();j++) {
			nn_numtemp.push_back(nn_num->at(j));
		}
	}
	
	t->SetBranchAddress("pionp_e",&pip_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pip_e->size();j++) {
			pip_etemp.push_back(pip_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.13957,pip_e->at(j))))));
		}
	}
	t->SetBranchAddress("pionp_phi",&pip_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pip_phi->size();j++) {
			pip_phitemp.push_back(pip_phi->at(j));
		}
	}
	t->SetBranchAddress("pionp_eta",&pip_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pip_eta->size();j++) {
			pip_etatemp.push_back(pip_eta->at(j));
		}
	}
	t->SetBranchAddress("pionp_num",&pip_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pip_num->size();j++) {
			pip_numtemp.push_back(pip_num->at(j));
		}
	}
	for (int i=0; i<pip_etemp.size(); i++) {
		pip_etatemp[i]=pip_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.13957,pip_etemp[i])));
		pip_phitemp[i]=pip_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.13957,pip_etemp[i])));
	}
	
	t->SetBranchAddress("pionn_e",&pin_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pin_e->size();j++) {
			pin_etemp.push_back(pin_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.13957,pin_e->at(j))))));
		}
	}
	t->SetBranchAddress("pionn_phi",&pin_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pin_phi->size();j++) {
			pin_phitemp.push_back(pin_phi->at(j));
		}
	}
	t->SetBranchAddress("pionn_eta",&pin_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pin_eta->size();j++) {
			pin_etatemp.push_back(pin_eta->at(j));
		}
	}
	t->SetBranchAddress("pionn_num",&pin_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<pin_num->size();j++) {
			pin_numtemp.push_back(pin_num->at(j));
		}
	}
	for (int i=0; i<pin_etemp.size(); i++) {
		pin_etatemp[i]=pin_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.13957,pin_etemp[i])));
		pin_phitemp[i]=pin_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.13957,pin_etemp[i])));
	}
	
	t->SetBranchAddress("Kp_e",&Kp_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kp_e->size();j++) {
			Kp_etemp.push_back(Kp_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.49368,Kp_e->at(j))))));
		}
	}
	t->SetBranchAddress("Kp_phi",&Kp_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kp_phi->size();j++) {
			Kp_phitemp.push_back(Kp_phi->at(j));
		}
	}
	t->SetBranchAddress("Kp_eta",&Kp_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kp_eta->size();j++) {
			Kp_etatemp.push_back(Kp_eta->at(j));
		}
	}
	t->SetBranchAddress("Kp_num",&Kp_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kp_num->size();j++) {
			Kp_numtemp.push_back(Kp_num->at(j));
		}
	}
	for (int i=0; i<Kp_etemp.size(); i++) {
		Kp_etatemp[i]=Kp_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.49368,Kp_etemp[i])));
		Kp_phitemp[i]=Kp_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.49368,Kp_etemp[i])));
	}
	
	t->SetBranchAddress("Kn_e",&Kn_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kn_e->size();j++) {
			Kn_etemp.push_back(Kn_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.49368,Kn_e->at(j))))));
		}
	}
	t->SetBranchAddress("Kn_phi",&Kn_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kn_phi->size();j++) {
			Kn_phitemp.push_back(Kn_phi->at(j));
		}
	}
	t->SetBranchAddress("Kn_eta",&Kn_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kn_eta->size();j++) {
			Kn_etatemp.push_back(Kn_eta->at(j));
		}
	}
	t->SetBranchAddress("Kn_num",&Kn_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<Kn_num->size();j++) {
			Kn_numtemp.push_back(Kn_num->at(j));
		}
	}
	for (int i=0; i<Kn_etemp.size(); i++) {
		Kn_etatemp[i]=Kn_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.49368,Kn_etemp[i])));
		Kn_phitemp[i]=Kn_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.49368,Kn_etemp[i])));
	}
	
	t->SetBranchAddress("KL_e",&KL_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<KL_e->size();j++) {
			KL_etemp.push_back(KL_e->at(j)*(1+gRandom->Gaus(0,error_combine(nrandomsigma_E, nrandomsigma2_E/sqrt(KL_e->at(j))))));
		}
	}
	t->SetBranchAddress("KL_phi",&KL_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<KL_phi->size();j++) {
			KL_phitemp.push_back(KL_phi->at(j)+gRandom->Gaus(0,nrandomsigma_angle));
		}
	}
	t->SetBranchAddress("KL_eta",&KL_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<KL_eta->size();j++) {
			KL_etatemp.push_back(KL_eta->at(j)+gRandom->Gaus(0,nrandomsigma_angle));
		}
	}
	t->SetBranchAddress("KL_num",&KL_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<KL_num->size();j++) {
			KL_numtemp.push_back(KL_num->at(j));
		}
	}
	
	t->SetBranchAddress("gamma_e",&g_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<g_e->size();j++) {
			g_etemp.push_back(g_e->at(j)*(1+gRandom->Gaus(0,error_combine(grandomsigma_E, grandomsigma2_E/sqrt(g_e->at(j))))));
		}
	}
	t->SetBranchAddress("gamma_phi",&g_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<g_phi->size();j++) {
			g_phitemp.push_back(g_phi->at(j)+gRandom->Gaus(0,grandomsigma_angle));
		}
	}
	t->SetBranchAddress("gamma_eta",&g_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<g_eta->size();j++) {
			g_etatemp.push_back(g_eta->at(j)+gRandom->Gaus(0,grandomsigma_angle));
		}
	}
	t->SetBranchAddress("gamma_num",&g_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<g_num->size();j++) {
			g_numtemp.push_back(g_num->at(j));
		}
	}
	
	t->SetBranchAddress("ep_e",&ep_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<ep_e->size();j++) {
			ep_etemp.push_back(ep_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.000511,ep_e->at(j))))));
		}
	}
	t->SetBranchAddress("ep_phi",&ep_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<ep_phi->size();j++) {
			ep_phitemp.push_back(ep_phi->at(j));
		}
	}
	t->SetBranchAddress("ep_eta",&ep_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<ep_eta->size();j++) {
			ep_etatemp.push_back(ep_eta->at(j));
		}
	}
	t->SetBranchAddress("ep_num",&ep_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<ep_num->size();j++) {
			ep_numtemp.push_back(ep_num->at(j));
		}
	}
	for (int i=0; i<ep_etemp.size(); i++) {
		ep_etatemp[i]=ep_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.000511,ep_etemp[i])));
		ep_phitemp[i]=ep_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.000511,ep_etemp[i])));
	}
	
	t->SetBranchAddress("en_e",&en_e);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<en_e->size();j++) {
			en_etemp.push_back(en_e->at(j)*(1+gRandom->Gaus(0,error_combine(crandomsigma_E, crandomsigma2_E*E2p(0.000511,en_e->at(j))))));
		}
	}
	t->SetBranchAddress("en_phi",&en_phi);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<en_phi->size();j++) {
			en_phitemp.push_back(en_phi->at(j));
		}
	}
	t->SetBranchAddress("en_eta",&en_eta);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<en_eta->size();j++) {
			en_etatemp.push_back(en_eta->at(j));
		}
	}
	t->SetBranchAddress("en_num",&en_num);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<en_num->size();j++) {
			en_numtemp.push_back(en_num->at(j));
		}
	}
	for (int i=0; i<en_etemp.size(); i++) {
		en_etatemp[i]=en_etatemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.000511,en_etemp[i])));
		en_phitemp[i]=en_phitemp[i]+gRandom->Gaus(0,error_combine(crandomsigma_angle, crandomsigma2_angle/E2p(0.000511,en_etemp[i])));
	}
	
	int ppp=0;
	int pnp=0;
	int npp=0;
	int nnp=0;
	int pip=0;
	int pin=0;
	int Kpp=0;
	int Knp=0;
	int KLp=0;
	int gp=0;
	int epp=0;
	int enp=0;
	
	int ppp2=0;
	int pnp2=0;
	int npp2=0;
	int nnp2=0;
	int pip2=0;
	int pin2=0;
	int Kpp2=0;
	int Knp2=0;
	int KLp2=0;
	int gp2=0;
	int epp2=0;
	int enp2=0;
	
	double ytemp;
	double ytemp1;
	double Q2JBtemp;
	double Q2JBtempABS;
	double Q2JBtempX;
	double Q2JBtempY;
	double Q2JBtempZ;
	double Q2s=0;
	double xs=0;
	double ys=0;
	double Q2D=0;
	double xD=0;
	double yD=0;
	double Q2l=0;
	double xl=0;
	double yl=0;
	double Q2j=0;
	double xj=0;
	double yj=0;
	
	double xsj=0;
	double xsl=0;
	double xsd=0;
	double xjs=0;
	double xjl=0;
	double xjd=0;
	double xls=0;
	double xlj=0;
	double xld=0;
	double xds=0;
	double xdj=0;
	double xdl=0;
	
	double tanthetahalf=0;
	double ei=960;
	double eii=275.0;
	double s=(ei+eii)*(ei+eii)-(ei-eii)*(ei-eii);

	for (int i=0; i<do_Q2temp.size();i++) {
		ppp2=ppp+pp_numtemp[i];
		pnp2=pnp+pn_numtemp[i];
		npp2=npp+np_numtemp[i];
		nnp2=nnp+nn_numtemp[i];
		pip2=pip+pip_numtemp[i];
		pin2=pin+pin_numtemp[i];
		Kpp2=Kpp+Kp_numtemp[i];
		Knp2=Knp+Kn_numtemp[i];
		KLp2=KLp+KL_numtemp[i];
		gp2=gp+g_numtemp[i];
		epp2=epp+ep_numtemp[i];
		enp2=enp+en_numtemp[i];
		ytemp=0;
		Q2JBtemp=0;
		Q2JBtempX=0;
		Q2JBtempY=0;
		Q2JBtempZ=0;
		Q2JBtempABS=0;
		for (int j=ppp; j<ppp2; j++) {
			ytemp=ytemp+pp_etemp[j]-E2p(0.93827,pp_etemp[j])*TMath::Cos(eta2theta(pp_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.93827,pp_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.93827,pp_etemp[j])*TMath::Sin(eta2theta(pp_etatemp[j]))*TMath::Cos(pp_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.93827,pp_etemp[j])*TMath::Sin(eta2theta(pp_etatemp[j]))*TMath::Sin(pp_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.93827,pp_etemp[j])*TMath::Sin(eta2theta(pp_etatemp[j]));
		}
		for (int j=pnp; j<pnp2; j++) {
			ytemp=ytemp+pn_etemp[j]-E2p(0.93827,pn_etemp[j])*TMath::Cos(eta2theta(pn_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.93827,pn_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.93827,pn_etemp[j])*TMath::Sin(eta2theta(pn_etatemp[j]))*TMath::Cos(pn_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.93827,pn_etemp[j])*TMath::Sin(eta2theta(pn_etatemp[j]))*TMath::Sin(pn_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.93827,pn_etemp[j])*TMath::Sin(eta2theta(pn_etatemp[j]));
		}
		for (int j=npp; j<npp2; j++) {
			ytemp=ytemp+np_etemp[j]-E2p(0.93957,np_etemp[j])*TMath::Cos(eta2theta(np_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.93957,np_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.93957,np_etemp[j])*TMath::Sin(eta2theta(np_etatemp[j]))*TMath::Cos(np_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.93957,np_etemp[j])*TMath::Sin(eta2theta(np_etatemp[j]))*TMath::Sin(np_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.93957,np_etemp[j])*TMath::Sin(eta2theta(np_etatemp[j]));
		}
		for (int j=nnp; j<nnp2; j++) {
			ytemp=ytemp+nn_etemp[j]-E2p(0.93957,nn_etemp[j])*TMath::Cos(eta2theta(nn_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.93957,nn_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.93957,nn_etemp[j])*TMath::Sin(eta2theta(nn_etatemp[j]))*TMath::Cos(nn_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.93957,nn_etemp[j])*TMath::Sin(eta2theta(nn_etatemp[j]))*TMath::Sin(nn_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.93957,nn_etemp[j])*TMath::Sin(eta2theta(nn_etatemp[j]));
		}
		for (int j=pip; j<pip2; j++) {
			ytemp=ytemp+pip_etemp[j]-E2p(0.13957,pip_etemp[j])*TMath::Cos(eta2theta(pip_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.13957,pip_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.13957,pip_etemp[j])*TMath::Sin(eta2theta(pip_etatemp[j]))*TMath::Cos(pip_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.13957,pip_etemp[j])*TMath::Sin(eta2theta(pip_etatemp[j]))*TMath::Sin(pip_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.13957,pip_etemp[j])*TMath::Sin(eta2theta(pip_etatemp[j]));
		}
		for (int j=pin; j<pin2; j++) {
			ytemp=ytemp+pin_etemp[j]-E2p(0.13957,pin_etemp[j])*TMath::Cos(eta2theta(pin_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.13957,pin_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.13957,pin_etemp[j])*TMath::Sin(eta2theta(pin_etatemp[j]))*TMath::Cos(pin_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.13957,pin_etemp[j])*TMath::Sin(eta2theta(pin_etatemp[j]))*TMath::Sin(pin_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.13957,pin_etemp[j])*TMath::Sin(eta2theta(pin_etatemp[j]));
		}
		for (int j=Kpp; j<Kpp2; j++) {
			ytemp=ytemp+Kp_etemp[j]-E2p(0.49368,Kp_etemp[j])*TMath::Cos(eta2theta(Kp_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.49368,Kp_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.49368,Kp_etemp[j])*TMath::Sin(eta2theta(Kp_etatemp[j]))*TMath::Cos(Kp_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.49368,Kp_etemp[j])*TMath::Sin(eta2theta(Kp_etatemp[j]))*TMath::Sin(Kp_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.49368,Kp_etemp[j])*TMath::Sin(eta2theta(Kp_etatemp[j]));
		}
		for (int j=Knp; j<Knp2; j++) {
			ytemp=ytemp+Kn_etemp[j]-E2p(0.49368,Kn_etemp[j])*TMath::Cos(eta2theta(Kn_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.49368,Kn_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.49368,Kn_etemp[j])*TMath::Sin(eta2theta(Kn_etatemp[j]))*TMath::Cos(Kn_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.49368,Kn_etemp[j])*TMath::Sin(eta2theta(Kn_etatemp[j]))*TMath::Sin(Kn_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.49368,Kn_etemp[j])*TMath::Sin(eta2theta(Kn_etatemp[j]));
		}
		for (int j=KLp; j<KLp2; j++) {
			ytemp=ytemp+KL_etemp[j]-E2p(0.49761,KL_etemp[j])*TMath::Cos(eta2theta(KL_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.49761,KL_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.49761,KL_etemp[j])*TMath::Sin(eta2theta(KL_etatemp[j]))*TMath::Cos(KL_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.49761,KL_etemp[j])*TMath::Sin(eta2theta(KL_etatemp[j]))*TMath::Sin(KL_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.49761,KL_etemp[j])*TMath::Sin(eta2theta(KL_etatemp[j]));
		}
		for (int j=gp; j<gp2; j++) {
			ytemp=ytemp+g_etemp[j]-E2p(0.0,g_etemp[j])*TMath::Cos(eta2theta(g_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.0,g_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.0,g_etemp[j])*TMath::Sin(eta2theta(g_etatemp[j]))*TMath::Cos(g_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.0,g_etemp[j])*TMath::Sin(eta2theta(g_etatemp[j]))*TMath::Sin(g_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.0,g_etemp[j])*TMath::Sin(eta2theta(g_etatemp[j]));
		}
		for (int j=epp; j<epp2; j++) {
			ytemp=ytemp+ep_etemp[j]-E2p(0.000511,ep_etemp[j])*TMath::Cos(eta2theta(ep_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.000511,ep_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.000511,ep_etemp[j])*TMath::Sin(eta2theta(ep_etatemp[j]))*TMath::Cos(ep_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.000511,ep_etemp[j])*TMath::Sin(eta2theta(ep_etatemp[j]))*TMath::Sin(ep_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.000511,ep_etemp[j])*TMath::Sin(eta2theta(ep_etatemp[j]));
		}
		for (int j=enp; j<enp2; j++) {
			ytemp=ytemp+en_etemp[j]-E2p(0.000511,en_etemp[j])*TMath::Cos(eta2theta(en_etatemp[j]));
			Q2JBtempZ=Q2JBtempZ+E2p(0.000511,en_etemp[j]);
			Q2JBtempX=Q2JBtempX+E2p(0.000511,en_etemp[j])*TMath::Sin(eta2theta(en_etatemp[j]))*TMath::Cos(en_phitemp[j]);
			Q2JBtempY=Q2JBtempY+E2p(0.000511,en_etemp[j])*TMath::Sin(eta2theta(en_etatemp[j]))*TMath::Sin(en_phitemp[j]);
			Q2JBtempABS=Q2JBtempABS+E2p(0.000511,en_etemp[j])*TMath::Sin(eta2theta(en_etatemp[j]));
		}

		if (!isnan(ytemp)&&!isnan(Q2JBtempX)&&!isnan(Q2JBtempY)) {
			
			tanthetahalf=ytemp/sqrt(Q2JBtempX*Q2JBtempX+Q2JBtempY*Q2JBtempY);
			
			Q2D=4*ei*ei*(1/tan(do_etatemp[i]/2))/(tanthetahalf+tan(do_etatemp[i]/2));
			yD=(tanthetahalf)/(tanthetahalf+tan(do_etatemp[i]/2));
			xD=Q2D/(s*yD);
			
			ys=ytemp/(ytemp+do_etemp[i]*(1-cos(do_etatemp[i])));
			Q2s=do_etemp[i]*do_etemp[i]*sin(do_etatemp[i])*sin(do_etatemp[i])/(1-ys);
			xs=Q2s/(s*ys);
			
			ytemp2=(2*ei-do_etemp[i]*(1-TMath::Cos(do_etatemp[i])));
			ytemp1=ytemp;
			
			yl=ytemp2/(2.0*ei);
			Q2l=2.0*ei*do_etemp[i]*(1+TMath::Cos(do_etatemp[i]));
			xl=Q2l/(s*yl);
			
			yj=ytemp/(2.0*ei);
			Q2j=(Q2JBtempX*Q2JBtempX+Q2JBtempY*Q2JBtempY)/(1-ytemp/(2.0*ei));
			xj=Q2j/(s*yj);
			
			xjl=Q2j/(s*yl);
			xjs=Q2j/(s*ys);
			xjd=Q2j/(s*yD);
			
			xlj=Q2l/(s*yj);
			xls=Q2l/(s*ys);
			xld=Q2l/(s*yD);
			
			xsj=Q2s/(s*yj);
			xsl=Q2s/(s*yl);
			xsd=Q2s/(s*yD);
			
			xdj=Q2D/(s*yj);
			xdl=Q2D/(s*yl);
			xds=Q2D/(s*ys);
			
			Q2txnum->Fill(do_xtemp[i],do_Q2temp[i]);
			
			Q2txyj  ->Fill(do_xtemp[i],do_Q2temp[i],fabs(yj-do_ytemp[i])/do_ytemp[i]);
			Q2txQ2j ->Fill(do_xtemp[i],do_Q2temp[i],fabs(Q2j-do_Q2temp[i])/do_Q2temp[i]);
			if (!isinf((xj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxj->Fill(do_xtemp[i],do_Q2temp[i],fabs(xj-do_xtemp[i])/do_xtemp[i]);
			}
			
			Q2txys  ->Fill(do_xtemp[i],do_Q2temp[i],fabs(ys-do_ytemp[i])/do_ytemp[i]);
			Q2txQ2s ->Fill(do_xtemp[i],do_Q2temp[i],fabs(Q2s-do_Q2temp[i])/do_Q2temp[i]);
			if (!isinf((xj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxs->Fill(do_xtemp[i],do_Q2temp[i],fabs(xs-do_xtemp[i])/do_xtemp[i]);
			}
			
			Q2txyl  ->Fill(do_xtemp[i],do_Q2temp[i],fabs(yl-do_ytemp[i])/do_ytemp[i]);
			Q2txQ2l ->Fill(do_xtemp[i],do_Q2temp[i],fabs(Q2l-do_Q2temp[i])/do_Q2temp[i]);
			if (!isinf((xj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxl->Fill(do_xtemp[i],do_Q2temp[i],fabs(xl-do_xtemp[i])/do_xtemp[i]);
			}
			
			Q2txyd  ->Fill(do_xtemp[i],do_Q2temp[i],fabs(yD-do_ytemp[i])/do_ytemp[i]);
			Q2txQ2d ->Fill(do_xtemp[i],do_Q2temp[i],fabs(Q2D-do_Q2temp[i])/do_Q2temp[i]);
			if (!isinf((xj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxd->Fill(do_xtemp[i],do_Q2temp[i],fabs(xD-do_xtemp[i])/do_xtemp[i]);
			}
			
			if (!isinf((xdl-do_xtemp[i])/do_xtemp[i])) {
				Q2txxdl->Fill(do_xtemp[i],do_Q2temp[i],fabs(xdl-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xds-do_xtemp[i])/do_xtemp[i])) {
				Q2txxds->Fill(do_xtemp[i],do_Q2temp[i],fabs(xds-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xdj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxdj->Fill(do_xtemp[i],do_Q2temp[i],fabs(xdj-do_xtemp[i])/do_xtemp[i]);
			}
			
			if (!isinf((xld-do_xtemp[i])/do_xtemp[i])) {
				Q2txxld->Fill(do_xtemp[i],do_Q2temp[i],fabs(xld-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xls-do_xtemp[i])/do_xtemp[i])) {
				Q2txxls->Fill(do_xtemp[i],do_Q2temp[i],fabs(xls-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xlj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxlj->Fill(do_xtemp[i],do_Q2temp[i],fabs(xlj-do_xtemp[i])/do_xtemp[i]);
			}
			
			if (!isinf((xsd-do_xtemp[i])/do_xtemp[i])) {
				Q2txxsd->Fill(do_xtemp[i],do_Q2temp[i],fabs(xsd-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xsl-do_xtemp[i])/do_xtemp[i])) {
				Q2txxsl->Fill(do_xtemp[i],do_Q2temp[i],fabs(xsl-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xsj-do_xtemp[i])/do_xtemp[i])) {
				Q2txxsj->Fill(do_xtemp[i],do_Q2temp[i],fabs(xsj-do_xtemp[i])/do_xtemp[i]);
			}
			
			if (!isinf((xjd-do_xtemp[i])/do_xtemp[i])) {
				Q2txxjd->Fill(do_xtemp[i],do_Q2temp[i],fabs(xjd-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xjl-do_xtemp[i])/do_xtemp[i])) {
				Q2txxjl->Fill(do_xtemp[i],do_Q2temp[i],fabs(xjl-do_xtemp[i])/do_xtemp[i]);
			}
			if (!isinf((xjs-do_xtemp[i])/do_xtemp[i])) {
				Q2txxjs->Fill(do_xtemp[i],do_Q2temp[i],fabs(xjs-do_xtemp[i])/do_xtemp[i]);
			}
			
		}
		ppp=ppp2;
		pnp=pnp2;
		npp=npp2;
		nnp=nnp2;
		pip=pip2;
		pin=pin2;
		Kpp=Kpp2;
		Knp=Knp2;
		KLp=KLp2;
		gp=gp2;
		epp=epp2;
		enp=enp2;
	}
	gStyle->SetOptStat(0);
	
	auto hist2d0=new TCanvas();
	hist2d0->SetLogz();
	gStyle->SetNumberContours(55);
	hist2d0->Divide(3,4);
	
	Int_t palette[20];
	palette[0] = 887;
	palette[1] = 890;
	palette[2] = 600;
	palette[3] = 860;
	palette[4] = 867;
	palette[5] = 870;
	palette[6] = 432;
	palette[7] = 836;
	palette[8] = 847;
	palette[9] = 837;
	palette[10] = 820;
	palette[11] = 412;
	palette[12] = 817;
	palette[13] = 811;
	palette[14] = 828;
	palette[15] = 396;
	palette[16] = 798;
	palette[17] = 807;
	palette[18] = 810;
	palette[19] = 900;
	
	hist2d0->cd(1);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txQ2j->Divide(Q2txnum);
	Q2txQ2j->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txQ2j->Draw("COLZ");
	
	hist2d0->cd(2);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txyj->Divide(Q2txnum);
	Q2txyj->SetAxisRange(0.0,1.0,"Z");
	Q2txyj->Draw("COLZ");
	
	hist2d0->cd(3);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxj->Divide(Q2txnum);
	Q2txxj->SetAxisRange(0.0,1.0,"Z");
	Q2txxj->Draw("COLZ");
	
	hist2d0->cd(4);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txQ2l->Divide(Q2txnum);
	Q2txQ2l->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txQ2l->Draw("COLZ");
	
	hist2d0->cd(5);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txyl->Divide(Q2txnum);
	Q2txyl->SetAxisRange(0.0,1.0,"Z");
	Q2txyl->Draw("COLZ");
	
	hist2d0->cd(6);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxl->Divide(Q2txnum);
	Q2txxl->SetAxisRange(0.0,1.0,"Z");
	Q2txxl->Draw("COLZ");
	
	hist2d0->cd(7);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txQ2d->Divide(Q2txnum);
	Q2txQ2d->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txQ2d->Draw("COLZ");
	
	hist2d0->cd(8);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txyd->Divide(Q2txnum);
	Q2txyd->SetAxisRange(0.0,1.0,"Z");
	Q2txyd->Draw("COLZ");
	
	hist2d0->cd(9);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxd->Divide(Q2txnum);
	Q2txxd->SetAxisRange(0.0,1.0,"Z");
	Q2txxd->Draw("COLZ");
	
	hist2d0->cd(10);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txQ2s->Divide(Q2txnum);
	Q2txQ2s->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txQ2s->Draw("COLZ");
	
	hist2d0->cd(11);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txys->Divide(Q2txnum);
	Q2txys->SetAxisRange(0.0,1.0,"Z");
	Q2txys->Draw("COLZ");
	
	hist2d0->cd(12);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxs->Divide(Q2txnum);
	Q2txxs->SetAxisRange(0.0,1.0,"Z");
	Q2txxs->Draw("COLZ");
	
	auto hist2d1=new TCanvas();
	hist2d1->SetLogz();
	gStyle->SetNumberContours(55);
	hist2d1->Divide(3,4);
	
	hist2d1->cd(1);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxld->Divide(Q2txnum);
	Q2txxld->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txxld->Draw("COLZ");
	
	hist2d1->cd(2);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxls->Divide(Q2txnum);
	Q2txxls->SetAxisRange(0.0,1.0,"Z");
	Q2txxls->Draw("COLZ");
	
	hist2d1->cd(3);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxlj->Divide(Q2txnum);
	Q2txxlj->SetAxisRange(0.0,1.0,"Z");
	Q2txxlj->Draw("COLZ");
	
	hist2d1->cd(4);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxdl->Divide(Q2txnum);
	Q2txxdl->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txxdl->Draw("COLZ");
	
	hist2d1->cd(5);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxds->Divide(Q2txnum);
	Q2txxds->SetAxisRange(0.0,1.0,"Z");
	Q2txxds->Draw("COLZ");
	
	hist2d1->cd(6);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxdj->Divide(Q2txnum);
	Q2txxdj->SetAxisRange(0.0,1.0,"Z");
	Q2txxdj->Draw("COLZ");
	
	hist2d1->cd(7);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxsl->Divide(Q2txnum);
	Q2txxsl->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txxsl->Draw("COLZ");
	
	hist2d1->cd(8);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxsd->Divide(Q2txnum);
	Q2txxsd->SetAxisRange(0.0,1.0,"Z");
	Q2txxsd->Draw("COLZ");
	
	hist2d1->cd(9);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxsj->Divide(Q2txnum);
	Q2txxsj->SetAxisRange(0.0,1.0,"Z");
	Q2txxsj->Draw("COLZ");
	
	hist2d1->cd(10);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxjl->Divide(Q2txnum);
	Q2txxjl->SetAxisRange(0.0,1.0,"Z");
	gStyle->SetPalette(20,palette);
	Q2txxjl->Draw("COLZ");
	
	hist2d1->cd(11);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxjd->Divide(Q2txnum);
	Q2txxjd->SetAxisRange(0.0,1.0,"Z");
	Q2txxjd->Draw("COLZ");
	
	hist2d1->cd(12);
	gPad->SetLogy();
	gPad->SetLogx();
	Q2txxjs->Divide(Q2txnum);
	Q2txxjs->SetAxisRange(0.0,1.0,"Z");
	Q2txxjs->Draw("COLZ");
}

