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

void readout_W2() {
	
	TFile * f = TFile::Open("EIC111.root","read");
	TTree * t = (TTree*) f->Get("trackTree");
	
	std::vector< float > * F2W2=0;
	
	ofstream phifile;
	phifile.open("EIC_W2.csv",ios::out|ios::trunc);
	t->SetBranchAddress("F2W2",&F2W2);
	for(int i=0;i<t->GetEntries();i++){
		t->GetEntry(i);
		for (int j=0;j<F2W2->size();j++) {
			//printf("%f\t",F2W2->at(j));
			phifile<<F2W2->at(j)<<"\t";
		}
	}
	phifile<<"\n";
	phifile.close();
	
}

