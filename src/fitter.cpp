#include "fitter.h"
using namespace std;

fitter::fitter(){//constructor
	RooAbsData::setDefaultStorageType(RooAbsData::Tree);//this allows to use the tree() method to get a tre from the roodataset at the end
}

void fitter::makeDataSet(TChain* chain, float M_min, float M_max){

	cout<<"makeDataSet called"<<endl;
	//Declare TTreeReader and the necessary variables
	TTreeReader r(chain);
	TTreeReaderArray<Float_t> vtx_prob(r, "mm_kin_vtx_prob");
	TTreeReaderArray<Float_t> DiMuon_mass(r, "mm_mass");
	TTreeReaderValue<Float_t> mm_kin_cosAlpha(r, "mm_kin_cosAlpha");
	TTreeReaderValue<bool> HLT_ref(r, "HLT_Dimuon0_LowMass");
	TTreeReaderValue<bool> HLT_sig(r, "HLT_Dimuon0_Jpsi");
	TTreeReaderValue<bool> HLT_norm(r, "HLT_Dimuon0_Jpsi_NoVertexing");

	RooRealVar *M = new RooRealVar("M","m(#mu#mu)",M_min, M_max);
	RooRealVar *HLT_ref_roo = new RooRealVar("HLT_Dimuon0_LowMass","HLT_Dimuon0_LowMass", -1);
	RooRealVar *HLT_sig_roo = new RooRealVar("HLT_Dimuon0_Jpsi","HLT_Dimuon0_Jpsi", -1);
	RooRealVar *HLT_norm_roo = new RooRealVar("HLT_Dimuon0_Jpsi_NoVertexing","HLT_Dimuon0_Jpsi_NoVertexing", -1);
	RooRealVar *vtx_prob_roo = new RooRealVar("mm_kin_vtx_prob","mm_kin_vtx_prob", -1);
	data = new RooDataSet("data", "data", RooArgSet(*M,*HLT_ref_roo,*HLT_sig_roo,*HLT_norm_roo,*vtx_prob_roo));


	//Loop on events...
	unsigned int max_evts=0;
	unsigned int i=0;
	chain->LoadTree(-1); //load the first tree in the chain
	unsigned long long int entries = chain->GetEntries();

	string str="Analyzed "+std::to_string(100*float(i)/entries)+"\% of the events.";
	while (r.Next()) {
		if (i%100000==0) {
			str="Analyzed "+std::to_string(100*float(i)/entries)+"\% of the events.";
			cout << string(str.length(),'\b');
			cout << str;
		}

		///add point int RooDataSet
		if (DiMuon_mass[0]>=M_min && DiMuon_mass[0]<=M_max && *mm_kin_cosAlpha>-2){
			*M=DiMuon_mass[0];
			*HLT_sig_roo=*HLT_sig;
			*HLT_norm_roo=*HLT_norm;
			*HLT_ref_roo=*HLT_ref;
			*vtx_prob_roo=vtx_prob[0];
			data->add(RooArgSet(*M,*HLT_ref_roo,*HLT_sig_roo,*HLT_norm_roo,*vtx_prob_roo));
		}
		if (++i>max_evts and max_evts>0) break;
	}

	if (w.var("M") != 0){
		w.var("M")->setRange(M_min, M_max);
	} else {
		w.import(*M);
	}
	data_full = (RooDataSet*)data->Clone();
	M_min_def = M_min;
	M_max_def = M_max;
}

void fitter::reduceDataSet(string cut, float M_min, float M_max){
	cout<<"reduceDataSet called"<<endl;
	w.var("M")->setMin(M_min);
	w.var("M")->setMax(M_max);
	data = (RooDataSet*)data_full->reduce(cut.c_str());
}

void fitter::resetDataSet(){
	cout<<"resetDataSet called"<<endl;
	w.var("M")->setMin(M_min_def);
	w.var("M")->setMax(M_max_def);
	delete data;
	data = data_full;
}

void fitter::preparePDF(bool do_preliminary_fit = false){

	binned_data = new RooDataHist("dh", "binned data", RooArgSet(*w.var("M")), *data);

	w.factory("RooCBShape::cb(M,mu[3.1,3,3.2],sigma0[0.01,0.005,0.05], alpha[0.1,3],n[1,5])");
	w.factory("Gaussian::g1(M,mu,sigma1[0.04,0.01,0.15])");
	w.factory("Gaussian::g2(M,mu,sigma2[0.04,0.01,0.15])");
	//w.factory("Gaussian::g3(M,mu,sigma3[0.01,0.005,0.1])");
	w.factory("SUM::2gau(gf1[0.2,0.7]*g1, g2)");
	//w.factory("SUM::3gau(gf2[0.05,1.0]*2gau, g3)");
	w.factory("SUM::sig(gf3[0.1,0.5]*2gau, cb)");
	w.factory("Exponential::e(M,tau1[-0.5,-3,-0.05])");
	float nentries = binned_data->sumEntries();
	RooRealVar s("s", "signal yield", 1,0,2); //signal yield
	w.import(s);
	w.var("s")->setRange(0.3*nentries, 1.01*nentries);
	w.var("s")->setVal(0.9*nentries);
	RooRealVar b("b", "background yield", 1,0,2); //background yield
	w.import(b);
	w.var("b")->setRange(-0.01*nentries, 0.3*nentries);
	w.var("b")->setVal(0.05*nentries);
	w.factory("SUM::model(s*sig,b*e)");

	if (do_preliminary_fit) {
		w.pdf("model")->fitTo(*binned_data, RooFit::Save()); //binned fit
		w.var("tau1")->setConstant(kTRUE);
		w.var("mu")->setConstant(kTRUE);
		w.var("alpha")->setConstant(kTRUE);
		w.var("n")->setConstant(kTRUE);
		w.var("gf1")->setConstant(kTRUE);
		w.var("gf3")->setConstant(kTRUE);

		//w.var("sigma0")->setConstant(kTRUE);
		//w.var("sigma1")->setConstant(kTRUE);
		//w.var("sigma2")->setConstant(kTRUE);
		w.saveSnapshot("default", w.allVars());
	}
}

void fitter::fit(){

	binned_data = new RooDataHist("dh", "binned data", RooArgSet(*w.var("M")), *data);

	float nentries = data->sumEntries();
	if (nentries == 0) return;
	if (w.getSnapshot("default") != 0) w.loadSnapshot("default");
	w.var("s")->setConstant(kFALSE);
	w.var("s")->setRange(0.3*nentries, 1.01*nentries);
	w.var("s")->setVal(0.9*nentries);
	w.var("b")->setConstant(kFALSE);
	w.var("b")->setRange(-0.01*nentries, 0.3*nentries);
	w.var("b")->setVal(0.05*nentries);
	w.factory("SUM::model(s*sig,b*e)");

	RooFitResult *result = 0;
	int n_attempts = 0;
	while ((result==0 || result->covQual()<3 || result->status()!=0) && n_attempts<10) { // loop until convergence
		result = w.pdf("model")->fitTo(*binned_data, RooFit::Save());
		n_attempts++;
	}
	if (n_attempts == 10) {
		cout<<"Fit did not converge after many attempts. Giving up."<<endl;
		exit(0);
	}
}

void fitter::saveFitPdf(string pdffname){
	if (data->sumEntries() == 0) {
		cout<<"Empty dataset. Nothing to plot"<<endl;
		return;
	}
	auto M = w.var("M");
	auto frame = M->frame();
	frame->SetTitle("");
	binned_data->plotOn(frame);
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