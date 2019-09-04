#include "fitter.cpp"

//data in "/mnt/hadoop/scratch/gandreas/Charmonium_Dimuon0_Jpsi_NoVertexing/*.root"

int main (int argc, char *argv[]) {

	std::map<int, std::map<std::string, std::string>> files;

	//files[2017]["B"] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017B-31Mar2018-v1+MINIAOD/*.root";
	files[2017]["C"] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017C-31Mar2018-v1+MINIAOD/*.root";
	files[2017]["D"] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017D-31Mar2018-v1+MINIAOD/*.root";
	files[2017]["E"] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017E-31Mar2018-v1+MINIAOD/*.root";
	files[2017]["F"] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017F-31Mar2018-v1+MINIAOD/*.root";

	fitter f;

	TChain* chain0 = new TChain("Events");
	int year = std::atoi(argv[1]);
	string section = argv[2];
	chain0->Add(files[year][section].c_str());

	///prepare pdf by fitting it to the whole DiMuon0_LowMass dataset
	f.makeDataSet(chain0, 2.95, 3.35);
	f.reduceDataSet("HLT_Dimuon0_LowMass", 2.95, 3.35);
	delete chain0;
	f.preparePDF(true);
	string plots_folder = "plots/"+to_string(year)+section;
	gSystem->Exec(("mkdir "+plots_folder).c_str());
	f.saveFitPdf(plots_folder+"/fit_prep.pdf");
	f.resetDataSet();
	//f.preparePDF(false);

	float vtxprob_bin_boundaries[] = {0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0};
	//float pt_bin_boundaries[] = {3, 3.5, 30, 100};
	int nbins1 = (float)sizeof(vtxprob_bin_boundaries)/sizeof(vtxprob_bin_boundaries[0])-1;

	TH1F hpass = TH1F("hpass", "hpass", nbins1, vtxprob_bin_boundaries);
	TH1F htot = TH1F("htot", "htot", nbins1, vtxprob_bin_boundaries);
	hpass.Sumw2();
	htot.Sumw2();


	int i =0;

	string mass_cut = " && (M>2.95 && M<3.25)";

	for (int j=0; j<nbins1; j++){
		string bincut = Form("mm_kin_vtx_prob>%f && mm_kin_vtx_prob<%f", vtxprob_bin_boundaries[j], vtxprob_bin_boundaries[j+1]);
		f.reduceDataSet(bincut+mass_cut+" && HLT_Dimuon0_Jpsi_NoVertexing", 2.95, 3.25);
		f.fit();
		f.saveFitPdf(plots_folder+Form("/fit_%d_%d.pdf", i, j));
		htot.SetBinContent(j+1, f.getSignalYield());
		htot.SetBinError(j+1, f.getSignalYieldError());

		double tot =  f.getSignalYield();

		f.reduceDataSet(bincut+mass_cut+" && HLT_Dimuon0_Jpsi", 2.95, 3.25); //It's right as it is, trust me. You don't have to apply also HLT_Dimuon0_Jpsi_NoVertexing, because of mis-aligned prescales!
		f.fit();
		f.saveFitPdf(plots_folder+Form("/fit_pass_%d_%d.pdf", i, j));
		hpass.SetBinContent(j+1, f.getSignalYield());
		hpass.SetBinError(j+1, f.getSignalYieldError());

		std::cout<<"tot: "<<tot<<" pass: "<<f.getSignalYield()<<endl;
	}

	htot.GetXaxis()->SetTitle("vertex probability");
	hpass.GetXaxis()->SetTitle("vertex probability");
	TCanvas* cc = new TCanvas("cc","",600,400);
	htot.Draw();
	cc->SaveAs("plots/tot.pdf");
	TCanvas* ccc = new TCanvas("ccc","",600,400);
	hpass.Draw();
	ccc->SaveAs("plots/pass.pdf");

	TFile *fout = new TFile(Form("hists%d%s.root",year,section.c_str()), "recreate");
	fout->cd();
	hpass.Write();
	htot.Write();
	//eff->Write();
	fout->Write();
	fout->Close();


  return 0;
}
