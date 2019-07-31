#include "fitter.cpp"

//data in "/mnt/hadoop/scratch/gandreas/Charmonium_Dimuon0_Jpsi_NoVertexing/*.root"

int main (int argc, char *argv[]) {

	TChain* chain0 = new TChain("Events");
	chain0->Add("/mnt/hadoop/scratch/gandreas/Charmonium_Dimuon0_LowMass_skimmed/*.root");
	fitter f;

	///prepare pdf by fitting it to the whole DiMuon0_LowMass dataset
	f.makeDataSet(chain0, 2.5, 3.5);
	f.preparePDF();
	f.saveFitPdf("plots/fit_prep.pdf");
  
	TChain* chain = new TChain("Events");
	chain->Add("/mnt/hadoop/scratch/gandreas/Charmonium_Dimuon0_Jpsi_NoVertexing/*.root");

	f.makeDataSet(chain, 2.95, 3.25);

	float pt_bin_boundaries[] = {3, 3.5, 4, 5, 6, 7, 8, 10, 15, 20, 30, 100};
	//float pt_bin_boundaries[] = {3, 3.5, 30, 100};
	int nbins1 = (float)sizeof(pt_bin_boundaries)/sizeof(pt_bin_boundaries[0])-1;

	TH1F hpass = TH1F("hpass", "hpass", nbins1, pt_bin_boundaries);
	TH1F htot = TH1F("htot", "htot", nbins1, pt_bin_boundaries);


	for (int i=0; i<nbins1; i++){
		string bincut = Form("mu1_pt>%f && mu1_pt<%f", pt_bin_boundaries[i], pt_bin_boundaries[i+1]);
		f.reduceDataSet(bincut);
	  	f.fit();
	  	f.saveFitPdf(Form("plots/fit_%d.pdf", i));
	  	htot.SetBinContent(i+1, f.getSignalYield());
	  	htot.SetBinError(i+1, f.getSignalYieldError());

	  	double tot =  f.getSignalYield();

		f.reduceDataSet(bincut+" && HLT_Dimuon0_LowMass");
	  	f.fit();
	  	f.saveFitPdf(Form("plots/fit_pass_%d.pdf", i));
	  	hpass.SetBinContent(i+1, f.getSignalYield());
	  	hpass.SetBinError(i+1, f.getSignalYieldError());

	  	std::cout<<"tot: "<<tot<<" pass: "<<f.getSignalYield()<<endl;
	}

	TCanvas* cc = new TCanvas("cc","",600,400);
	htot.SetMinimum(0);
	htot.Draw();
	hpass.SetLineColor(kRed);
	hpass.Draw("same");
	cc->SaveAs("plots/tot-pass.pdf");

	TEfficiency eff = TEfficiency(hpass, htot);
	TCanvas* c = new TCanvas("c","",600,400);
	eff.Draw("AP");
	c->SaveAs("plots/efficiency.pdf");

  return 0;
}
