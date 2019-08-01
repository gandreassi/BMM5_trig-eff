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

	float pt_bin_boundaries[] = {3, 4, 5, 6, 8, 12, 100};
	//float pt_bin_boundaries[] = {3, 3.5, 30, 100};
	int nbins1 = (float)sizeof(pt_bin_boundaries)/sizeof(pt_bin_boundaries[0])-1;

	TH2F hpass = TH2F("hpass", "hpass", nbins1, pt_bin_boundaries, nbins1, pt_bin_boundaries);
	TH2F htot = TH2F("htot", "htot", nbins1, pt_bin_boundaries, nbins1, pt_bin_boundaries);


	for (int i=0; i<nbins1; i++){
		for (int j=0; j<nbins1; j++){
			string bincut = Form("mu1_pt>%f && mu1_pt<%f && mu2_pt>%f && mu2_pt<%f", pt_bin_boundaries[i], pt_bin_boundaries[i+1], pt_bin_boundaries[j], pt_bin_boundaries[j+1]);
			f.reduceDataSet(bincut);
		  	f.fit();
		  	f.saveFitPdf(Form("plots/fit_%d_%d.pdf", i, j));
		  	htot.SetBinContent(i+1, j+1, f.getSignalYield());
		  	htot.SetBinError(i+1, j+1, f.getSignalYieldError());

		  	double tot =  f.getSignalYield();

			f.reduceDataSet(bincut+" && HLT_Dimuon0_LowMass");
		  	f.fit();
		  	f.saveFitPdf(Form("plots/fit_pass_%d_%d.pdf", i, j));
		  	hpass.SetBinContent(i+1, j+1, f.getSignalYield());
		  	hpass.SetBinError(i+1, j+1, f.getSignalYieldError());

		  	std::cout<<"tot: "<<tot<<" pass: "<<f.getSignalYield()<<endl;
		}
	}

	htot.GetXaxis()->SetTitle("p_{T} (#mu_{1})");
	htot.GetYaxis()->SetTitle("p_{T} (#mu_{2})");
	hpass.GetXaxis()->SetTitle("p_{T} (#mu_{1})");
	hpass.GetYaxis()->SetTitle("p_{T} (#mu_{2})");
	TCanvas* cc = new TCanvas("cc","",600,400);
	htot.Draw("colz");
	cc->SaveAs("plots/tot.pdf");
	TCanvas* ccc = new TCanvas("ccc","",600,400);
	hpass.Draw("colz");
	ccc->SaveAs("plots/pass.pdf");


	TEfficiency* eff = new TEfficiency(hpass, htot);
	TCanvas* c = new TCanvas("c","",600,400);
	eff->Draw("colz");
	gPad->Update();

	c->SaveAs("plots/efficiency.pdf");

  return 0;
}
