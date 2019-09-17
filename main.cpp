#include "fitter.cpp"

//data in "/mnt/hadoop/scratch/gandreas/Charmonium_Dimuon0_Jpsi_NoVertexing/*.root"

int main (int argc, char *argv[]) {

	std::map<int, std::map<char, std::string>> files;

	//files[2017]["B"] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017B-31Mar2018-v1+MINIAOD/*.root";
	//files[2017]['C'] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017C-31Mar2018-v1+MINIAOD/*.root";
	files[2017]['D'] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017D-31Mar2018-v1+MINIAOD/*.root";
	files[2017]['E'] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017E-31Mar2018-v1+MINIAOD/*.root";
	files[2017]['F'] = "/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD/501/Charmonium+Run2017F-31Mar2018-v1+MINIAOD/*.root";

	files[2018]['A'] = "/eos/user/g/gandreas/jpsikdata/2018/Charmonium+Run2018A-17Sep2018-v1+MINIAOD/*.root";
	files[2018]['B'] = "/eos/user/g/gandreas/jpsikdata/2018/Charmonium+Run2018B-17Sep2018-v1+MINIAOD/*.root";
	files[2018]['C'] = "/eos/user/g/gandreas/jpsikdata/2018/Charmonium+Run2018C-17Sep2018-v1+MINIAOD/*.root";
	files[2018]['D'] = "/eos/user/g/gandreas/jpsikdata/2018/Charmonium+Run2018D-17Sep2018-v1+MINIAOD/*.root";

	fitter f;

	TChain* chain0 = new TChain("Events");
	int year = std::atoi(argv[1]);
	string section = argv[2];
	for(char& c : section) {
		chain0->Add(files[year][c].c_str());
	}

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

	float vtxprob_bin_boundaries[] = {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.7, 1.0};
	float displ_bin_boundaries[] = {0, 3, 6, 9, 20, 40, 150};
	//float pt_bin_boundaries[] = {3, 3.5, 30, 100};
	int nbins_vtx = (float)sizeof(vtxprob_bin_boundaries)/sizeof(vtxprob_bin_boundaries[0])-1;
	int nbins_displ = (float)sizeof(displ_bin_boundaries)/sizeof(displ_bin_boundaries[0])-1;

	TH2F hpass = TH2F("hpass", "hpass", nbins_vtx, vtxprob_bin_boundaries, nbins_displ, displ_bin_boundaries);
	TH2F htot = TH2F("htot", "htot", nbins_vtx, vtxprob_bin_boundaries, nbins_displ, displ_bin_boundaries);
	hpass.Sumw2();
	htot.Sumw2();


	string mass_cut = " && (M>2.95 && M<3.25)";

	for (int j=0; j<nbins_vtx; j++){ 
		string vtx_bincut = Form("mm_kin_vtx_prob>%f && mm_kin_vtx_prob<%f", vtxprob_bin_boundaries[j], vtxprob_bin_boundaries[j+1]);
		for (int k=0; k<nbins_displ; k++){
			string displ_bincut = Form("&& mm_kin_slxy>%f && mm_kin_slxy<%f", displ_bin_boundaries[k], displ_bin_boundaries[k+1]);

			f.reduceDataSet(vtx_bincut+displ_bincut+mass_cut+" && HLT_DoubleMu4_Jpsi_NoVertexing", 2.95, 3.25);
			f.fit();
			f.saveFitPdf(plots_folder+Form("/fit_%d_%d.pdf", j, k));
			htot.SetBinContent(j+1, k+1, f.getSignalYield());
			htot.SetBinError(j+1, k+1, f.getSignalYieldError());

			double tot =  f.getSignalYield();

			f.reduceDataSet(vtx_bincut+displ_bincut+mass_cut+" && HLT_DoubleMu4_Jpsi_Displaced", 2.95, 3.25); //It's right as it is, trust me. You don't have to apply also HLT_Dimuon0_Jpsi_NoVertexing, because of mis-aligned prescales!
			f.fit();
			f.saveFitPdf(plots_folder+Form("/fit_pass_%d_%d.pdf", j, k));
			hpass.SetBinContent(j+1, k+1, f.getSignalYield());
			hpass.SetBinError(j+1, k+1, f.getSignalYieldError());

			std::cout<<"tot: "<<tot<<" pass: "<<f.getSignalYield()<<endl;
		}

	}

	htot.GetXaxis()->SetTitle("vertex probability");
	hpass.GetXaxis()->SetTitle("vertex probability");
	TCanvas* cc = new TCanvas("cc","",600,400);
	htot.Draw("colz");
	cc->SaveAs((plots_folder+"/tot.pdf").c_str());
	TCanvas* ccc = new TCanvas("ccc","",600,400);
	hpass.Draw("colz");
	ccc->SaveAs((plots_folder+"/pass.pdf").c_str());

	TFile *fout = new TFile(Form("hists%d%s.root",year,section.c_str()), "recreate");
	fout->cd();
	hpass.Write();
	htot.Write();
	//eff->Write();
	fout->Write();
	fout->Close();


  return 0;
}
