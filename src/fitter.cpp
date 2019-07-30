#include "fitter.h"
using namespace std;

fitter::fitter(){//constructor
	RooAbsData::setDefaultStorageType(RooAbsData::Tree);//this allows to use the tree() method to get a tre from the roodataset at the end
}

void fitter::makeDataSet(TChain* chain){

	//Declare TTreeReader and the necessary variables
	TTreeReader r(chain);
	TTreeReaderValue<UInt_t> mus(r, "nMuon");
	TTreeReaderArray<Float_t> pt(r, "Muon_pt");
	TTreeReaderArray<Float_t> eta(r, "Muon_eta");
	TTreeReaderArray<Float_t> phi(r, "Muon_phi");
	TTreeReaderArray<Float_t> mass(r, "Muon_mass");
	TTreeReaderArray<Int_t> c(r, "Muon_charge");
	TTreeReaderValue<bool> HLT_sig(r, "HLT_Dimuon0_LowMass");

	///RooFit stuff
	const float M_min =2.91;
	const float M_max =3.29;
	RooRealVar *M = new RooRealVar("M","m(#mu#mu)",M_min, M_max);
	RooRealVar *pt1 = new RooRealVar("mu1_pt","mu1_pt", -1);
	RooRealVar *pt2 = new RooRealVar("mu2_pt","mu2_pt", -1);
	RooRealVar *HLT_sig_roo = new RooRealVar("HLT_Dimuon0_LowMass","HLT_Dimuon0_LowMass", -1);
	data = new RooDataSet("data", "data", RooArgSet(*M,*pt1,*pt2,*HLT_sig_roo));


	//Loop on events...
	unsigned int max_evts=0;
	unsigned int i=0;

	while (r.Next()) {
		int c_prod=c[0]; //charge of the first muon
		int index_second_muon=-1;
		for (int i=1; i<(int)*mus; i++){
			if (c_prod*c[i]==-1){//opposite charge requirement. That"s going to be our second muon
				index_second_muon=i;
				break;
			}
		}

		if (index_second_muon>0){

			TLorentzVector P1;
			TLorentzVector P2;
			P1.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
			P2.SetPtEtaPhiM(pt[index_second_muon], eta[index_second_muon], phi[index_second_muon], mass[index_second_muon]);
			TLorentzVector P = P1+P2;

			///add RooDataSet
			float DiMuon_mass=P.M();
			if (DiMuon_mass>=M_min && DiMuon_mass<=M_max){
				*M=DiMuon_mass;
				*pt1=pt[0];
				*pt2=pt[index_second_muon];
				*HLT_sig_roo=*HLT_sig;
				data->add(RooArgSet(*M, *pt1, *pt2, *HLT_sig_roo));
			}
		}
		if (++i>max_evts and max_evts>0) break;
	}

	//Now we can create our pdf and fit it
	w.import(*M);
	w.import(*data);
}

void fitter::reduceDataSet(string cut){
	data = (RooDataSet*)w.data("data")->reduce(cut.c_str());
}

void fitter::fit(){
	if (data->sumEntries() == 0) return;
	w.factory("RooCBShape::cb(M,mu[3.1,3,3.2],sigma0[0.01,0.005,0.05], alpha[0.1,3],n[3])");
	w.factory("Gaussian::g1(M,mu,sigma1[0.04,0.01,0.15])");
	w.factory("Gaussian::g2(M,mu,sigma2[0.04,0.01,0.15])");
	//w.factory("Gaussian::g3(M,mu,sigma3[0.01,0.005,0.1])");
	w.factory("SUM::2gau(gf1[0.05,1.0]*g1, g2)");
	//w.factory("SUM::3gau(gf2[0.05,1.0]*2gau, g3)");
	w.factory("SUM::sig(gf3[0.05,1.0]*2gau, cb)");
	w.factory("Exponential::e(M,tau1[-0.05,-0.1,-0.005])");
	float nentries = data->sumEntries();
	RooRealVar s("s", "signal yield", 1,0,2); //signal yield
	w.import(s);
	w.var("s")->setRange(0.3*nentries, 1.01*nentries);
	w.var("s")->setVal(0.9*nentries);
	RooRealVar b("b", "background yield", 1,0,2); //background yield
	w.import(b);
	w.var("b")->setRange(-0.01*nentries, 0.3*nentries);
	w.var("b")->setVal(0.05*nentries);
	w.factory("SUM::model(s*sig,b*e)");

	w.pdf("model")->fitTo(*data, RooFit::Save());
}

void fitter::saveFitPdf(string pdffname){
	if (data->sumEntries() == 0) return;
	auto M = w.var("M");
	auto frame = M->frame();
	frame->SetTitle("");
	data->plotOn(frame);
	w.pdf("model")->plotOn(frame);
	auto hpull = frame->pullHist();
	w.pdf("model")->plotOn(frame, RooFit::Components("e"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
	w.pdf("model")->plotOn(frame, RooFit::Components("sig"), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

	frame->GetXaxis()->SetTitleSize(0);
	frame->GetXaxis()->SetLabelSize(0);
	auto pframe = M->frame();
	pframe->SetTitle("");
	pframe->GetXaxis()->SetLabelSize(0.1);
	pframe->GetXaxis()->SetTitle("m(#mu#mu)");
	pframe->GetXaxis()->SetTitleSize(0.15);
	pframe->GetYaxis()->SetLabelSize(0.1);
	pframe->GetYaxis()->SetTitle("Pool  ");
	pframe->GetYaxis()->SetTitleSize(0.15);
	pframe->GetYaxis()->SetTitleOffset(0.3);
	pframe->addPlotable(hpull,"P");

	TCanvas *canvas = new TCanvas("fit", "fit", 800, 600);;

	TPad* pad1 = new TPad("pad1", "pad1",0,0.25,1,1);
	pad1->SetBottomMargin(0.02);
	pad1->cd();
	frame->Draw();

	TPad* pad2 = new TPad("pad2", "pad2",0,0,1,0.25);
	pad2->SetTopMargin(0.05);
	pad2->SetBottomMargin(0.4);
	pad2->cd();
	pframe->Draw();

	canvas->cd();
	pad1->Draw();
	pad2->Draw();
	canvas->Draw();
	canvas->SaveAs(pdffname.c_str());
}

float fitter::getSignalYield() {
	if (data->sumEntries() == 0) return 0;
	return w.var("s")->getVal();
}

float fitter::getSignalYieldError() {
	if (data->sumEntries() == 0) return 0;
	return w.var("s")->getError();
}